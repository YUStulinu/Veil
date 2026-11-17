// CUDA kernels for the prover: the G1 multi-scalar multiplication (Pippenger
// with buckets sorted by CUB) and the NTT over Fr. Field elements stay in the
// CPU's Montgomery form (R = 2^384 for Fq, 2^256 for Fr), read as 32-bit limbs.
//
// C interface (see src/gpu.rs):
//   veil_gpu_msm_g1(points, stride, scalars, n, out)  -> 0 on success
//   veil_gpu_ntt(data, n, omega)                      -> 0 on success

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cub/cub.cuh>
#include <chrono>
#include <cstdlib>

static bool g_prof = getenv("VEIL_GPU_PROFILE") != nullptr;
static auto g_t0 = std::chrono::steady_clock::now();
static void prof(const char* what) {
    if (!g_prof) return;
    cudaDeviceSynchronize();
    auto t = std::chrono::steady_clock::now();
    fprintf(stderr, "  %-10s %8.3f ms\n", what, std::chrono::duration<double, std::milli>(t - g_t0).count());
    g_t0 = t;
}

typedef uint32_t u32;
typedef uint64_t u64;

#define CHECK(x)                                                                     \
    do {                                                                             \
        cudaError_t e_ = (x);                                                        \
        if (e_ != cudaSuccess) {                                                     \
            fprintf(stderr, "veil cuda: %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); \
            return 1;                                                                \
        }                                                                            \
    } while (0)

// ---------------- fields ----------------

template <int N>
struct Fp {
    u32 v[N];
};
typedef Fp<12> Fq;
typedef Fp<8> Fr;

__constant__ u32 FQ_P[12] = {0xffffaaab, 0xb9feffff, 0xb153ffff, 0x1eabfffe, 0xf6b0f624, 0x6730d2a0, 0xf38512bf, 0x64774b84, 0x434bacd7, 0x4b1ba7b6, 0x397fe69a, 0x1a0111ea};
__constant__ u32 FR_P[8] = {0x00000001, 0xffffffff, 0xfffe5bfe, 0x53bda402, 0x09a1d805, 0x3339d808, 0x299d7d48, 0x73eda753};
#define FQ_INV 0xfffcfffdu
#define FR_INV 0xffffffffu

template <int N>
__device__ __forceinline__ bool geq(const u32* a, const u32* p) {
#pragma unroll
    for (int i = N - 1; i >= 0; i--) {
        if (a[i] != p[i]) return a[i] > p[i];
    }
    return true;
}

template <int N>
__device__ __forceinline__ void sub_p(u32* a, const u32* p) {
    u64 borrow = 0;
#pragma unroll
    for (int i = 0; i < N; i++) {
        u64 d = (u64)a[i] - p[i] - borrow;
        a[i] = (u32)d;
        borrow = (d >> 63) & 1;
    }
}

template <int N>
__device__ __forceinline__ Fp<N> add(const Fp<N>& a, const Fp<N>& b, const u32* p) {
    Fp<N> r;
    u64 c = 0;
#pragma unroll
    for (int i = 0; i < N; i++) {
        u64 s = (u64)a.v[i] + b.v[i] + c;
        r.v[i] = (u32)s;
        c = s >> 32;
    }
    if (geq<N>(r.v, p)) sub_p<N>(r.v, p);
    return r;
}

template <int N>
__device__ __forceinline__ Fp<N> sub(const Fp<N>& a, const Fp<N>& b, const u32* p) {
    Fp<N> r;
    u64 borrow = 0;
#pragma unroll
    for (int i = 0; i < N; i++) {
        u64 d = (u64)a.v[i] - b.v[i] - borrow;
        r.v[i] = (u32)d;
        borrow = (d >> 63) & 1;
    }
    if (borrow) {
        u64 c = 0;
#pragma unroll
        for (int i = 0; i < N; i++) {
            u64 s = (u64)r.v[i] + p[i] + c;
            r.v[i] = (u32)s;
            c = s >> 32;
        }
    }
    return r;
}

// Montgomery multiplication, CIOS
template <int N>
__device__ __forceinline__ Fp<N> mul(const Fp<N>& a, const Fp<N>& b, const u32* p, u32 inv) {
    u32 t[N + 2];
#pragma unroll
    for (int i = 0; i < N + 2; i++) t[i] = 0;
#pragma unroll
    for (int i = 0; i < N; i++) {
        u64 c = 0;
#pragma unroll
        for (int j = 0; j < N; j++) {
            u64 s = (u64)t[j] + (u64)a.v[j] * b.v[i] + c;
            t[j] = (u32)s;
            c = s >> 32;
        }
        u64 s = (u64)t[N] + c;
        t[N] = (u32)s;
        t[N + 1] = (u32)(s >> 32);
        u32 m = t[0] * inv;
        s = (u64)t[0] + (u64)m * p[0];
        c = s >> 32;
#pragma unroll
        for (int j = 1; j < N; j++) {
            s = (u64)t[j] + (u64)m * p[j] + c;
            t[j - 1] = (u32)s;
            c = s >> 32;
        }
        s = (u64)t[N] + c;
        t[N - 1] = (u32)s;
        t[N] = t[N + 1] + (u32)(s >> 32);
    }
    Fp<N> r;
#pragma unroll
    for (int i = 0; i < N; i++) r.v[i] = t[i];
    if (t[N] || geq<N>(r.v, p)) sub_p<N>(r.v, p);
    return r;
}

template <int N>
__device__ __forceinline__ bool is_zero(const Fp<N>& a) {
    u32 x = 0;
#pragma unroll
    for (int i = 0; i < N; i++) x |= a.v[i];
    return x == 0;
}

__device__ __forceinline__ Fq qadd(const Fq& a, const Fq& b) { return add<12>(a, b, FQ_P); }
__device__ __forceinline__ Fq qsub(const Fq& a, const Fq& b) { return sub<12>(a, b, FQ_P); }
__device__ __forceinline__ Fq qmul(const Fq& a, const Fq& b) { return mul<12>(a, b, FQ_P, FQ_INV); }
__device__ __forceinline__ Fq qsqr(const Fq& a) { return mul<12>(a, a, FQ_P, FQ_INV); }
__device__ __forceinline__ Fq qdbl(const Fq& a) { return add<12>(a, a, FQ_P); }

__device__ __forceinline__ Fr radd(const Fr& a, const Fr& b) { return add<8>(a, b, FR_P); }
__device__ __forceinline__ Fr rsub(const Fr& a, const Fr& b) { return sub<8>(a, b, FR_P); }
__device__ __forceinline__ Fr rmul(const Fr& a, const Fr& b) { return mul<8>(a, b, FR_P, FR_INV); }

// ---------------- G1 (Jacobian, a = 0) ----------------

struct Jac {
    Fq x, y, z;
};
struct Aff {
    Fq x, y;
};

__device__ __forceinline__ Jac jac_zero() {
    Jac r;
    for (int i = 0; i < 12; i++) r.x.v[i] = r.y.v[i] = r.z.v[i] = 0;
    return r;
}

__device__ Jac jac_double(const Jac& p) {
    if (is_zero(p.z)) return p;
    Fq a = qsqr(p.x), b = qsqr(p.y), c = qsqr(b);
    Fq d = qdbl(qsub(qsub(qsqr(qadd(p.x, b)), a), c));
    Fq e = qadd(qdbl(a), a);
    Fq f = qsqr(e);
    Jac r;
    r.x = qsub(f, qdbl(d));
    r.y = qsub(qmul(e, qsub(d, r.x)), qdbl(qdbl(qdbl(c))));
    r.z = qdbl(qmul(p.y, p.z));
    return r;
}

__device__ Jac jac_add(const Jac& p, const Jac& q) {
    if (is_zero(p.z)) return q;
    if (is_zero(q.z)) return p;
    Fq z1z1 = qsqr(p.z), z2z2 = qsqr(q.z);
    Fq u1 = qmul(p.x, z2z2), u2 = qmul(q.x, z1z1);
    Fq s1 = qmul(qmul(p.y, q.z), z2z2), s2 = qmul(qmul(q.y, p.z), z1z1);
    Fq h = qsub(u2, u1);
    Fq r = qdbl(qsub(s2, s1));
    if (is_zero(h)) {
        if (is_zero(r)) return jac_double(p);
        return jac_zero();
    }
    Fq i = qsqr(qdbl(h));
    Fq j = qmul(h, i);
    Fq v = qmul(u1, i);
    Jac o;
    o.x = qsub(qsub(qsqr(r), j), qdbl(v));
    o.y = qsub(qmul(r, qsub(v, o.x)), qdbl(qmul(s1, j)));
    o.z = qmul(qsub(qsub(qsqr(qadd(p.z, q.z)), z1z1), z2z2), h);
    return o;
}

// p + q with q affine (not the identity)
__device__ Jac jac_add_aff(const Jac& p, const Aff& q, const Fq& one) {
    if (is_zero(p.z)) {
        Jac r;
        r.x = q.x;
        r.y = q.y;
        r.z = one;
        return r;
    }
    Fq z1z1 = qsqr(p.z);
    Fq u2 = qmul(q.x, z1z1);
    Fq s2 = qmul(qmul(q.y, p.z), z1z1);
    Fq h = qsub(u2, p.x);
    Fq r = qdbl(qsub(s2, p.y));
    if (is_zero(h)) {
        if (is_zero(r)) return jac_double(p);
        return jac_zero();
    }
    Fq hh = qsqr(h);
    Fq i = qdbl(qdbl(hh));
    Fq j = qmul(h, i);
    Fq v = qmul(p.x, i);
    Jac o;
    o.x = qsub(qsub(qsqr(r), j), qdbl(v));
    o.y = qsub(qmul(r, qsub(v, o.x)), qdbl(qmul(p.y, j)));
    o.z = qsub(qsub(qsqr(qadd(p.z, h)), z1z1), hh);
    return o;
}

// ---------------- MSM ----------------

// the CPU's Affine<Fq>: x (6 u64), y (6 u64), then a bool "infinity"
__global__ void k_load_points(const unsigned char* src, size_t stride, size_t n, Aff* dst, u32* inf) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i >= n) return;
    const u32* p = (const u32*)(src + i * stride);
    for (int k = 0; k < 12; k++) {
        dst[i].x.v[k] = p[k];
        dst[i].y.v[k] = p[12 + k];
    }
    inf[i] = src[i * stride + 96];
}

// Fr from Montgomery to canonical, then (bucket key, point index) per window.
// Digit 0 and points at infinity go to bucket 0 of their window, which is ignored.
__global__ void k_keys(const Fr* scalars, const u32* inf, size_t n, int c, int windows, u32* keys, u32* vals) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i >= n) return;
    Fr one;
    for (int k = 0; k < 8; k++) one.v[k] = k == 0;
    Fr s = rmul(scalars[i], one);
    for (int w = 0; w < windows; w++) {
        int bit = w * c;
        int limb = bit / 32, off = bit % 32;
        u64 lo = limb < 8 ? s.v[limb] : 0;
        u64 hi = limb + 1 < 8 ? s.v[limb + 1] : 0;
        u32 d = (u32)(((lo | (hi << 32)) >> off) & ((1u << c) - 1));
        if (inf[i]) d = 0;
        keys[(size_t)w * n + i] = ((u32)w << c) | d;
        vals[(size_t)w * n + i] = (u32)i;
    }
}

__global__ void k_bounds(const u32* keys, size_t m, u32* start, u32* end) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i >= m) return;
    u32 k = keys[i];
    if (i == 0 || keys[i - 1] != k) start[k] = (u32)i;
    if (i == m - 1 || keys[i + 1] != k) end[k] = (u32)i + 1;
}

__global__ void k_buckets(const u32* start, const u32* end, const u32* vals, const Aff* pts, size_t nb, int c, Jac* buckets, Fq one) {
    size_t b = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (b >= nb) return;
    Jac acc = jac_zero();
    if ((b & ((1u << c) - 1)) != 0) {
        for (u32 k = start[b]; k < end[b]; k++) acc = jac_add_aff(acc, pts[vals[k]], one);
    }
    buckets[b] = acc;
}

// window w, segment t of `segs`: sum over its buckets d of d * B_d
__global__ void k_reduce(const Jac* buckets, int c, int windows, int segs, Jac* out) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= windows * segs) return;
    int w = id / segs, t = id % segs;
    u32 nb = 1u << c;
    u32 per = (nb + segs - 1) / segs;
    u32 lo = t * per, hi = min(nb, lo + per);
    if (lo == 0) lo = 1;
    Jac running = jac_zero(), acc = jac_zero();
    for (u32 d = hi; d > lo; d--) {
        running = jac_add(running, buckets[((size_t)w << c) + d - 1]);
        acc = jac_add(acc, running);
    }
    // acc = sum (d - lo + 1) B_d; add (lo - 1) * running
    u32 k = lo - 1;
    Jac kr = jac_zero();
    for (int bit = 31; bit >= 0; bit--) {
        kr = jac_double(kr);
        if ((k >> bit) & 1) kr = jac_add(kr, running);
    }
    out[id] = jac_add(acc, kr);
}

struct MsmCache {
    const void* host = nullptr;
    size_t n = 0;
    u64 tag[4] = {0, 0, 0, 0};
    Aff* pts = nullptr;
    u32* inf = nullptr;
};
static MsmCache g_cache;

// Windows that divide 255 bits exactly (255 = 15 * 17): a short last window
// would put all points in a few buckets and serialize them on one thread.
static int window_bits(size_t n) {
    return n >= ((size_t)1 << 21) ? 17 : 15;
}

extern "C" int veil_gpu_msm_g1(const unsigned char* points, size_t stride, const u64* scalars, size_t n, u64* out, size_t* out_count) {
    if (n == 0) return 1;
    // the SRS points stay on the GPU between calls
    u64 tag[4];
    memcpy(tag, points, sizeof(tag));
    if (g_cache.host != points || g_cache.n < n || memcmp(g_cache.tag, tag, sizeof(tag)) != 0) {
        cudaFree(g_cache.pts);
        cudaFree(g_cache.inf);
        unsigned char* raw;
        CHECK(cudaMalloc(&raw, n * stride));
        CHECK(cudaMemcpy(raw, points, n * stride, cudaMemcpyHostToDevice));
        CHECK(cudaMalloc(&g_cache.pts, n * sizeof(Aff)));
        CHECK(cudaMalloc(&g_cache.inf, n * sizeof(u32)));
        k_load_points<<<(n + 255) / 256, 256>>>(raw, stride, n, g_cache.pts, g_cache.inf);
        CHECK(cudaGetLastError());
        CHECK(cudaFree(raw));
        g_cache.host = points;
        g_cache.n = n;
        memcpy(g_cache.tag, tag, sizeof(tag));
    }
    prof("start");
    int c = window_bits(n);
    int windows = (255 + c - 1) / c;
    size_t m = n * windows;
    size_t nb = (size_t)windows << c;
    Fr* d_s;
    u32 *keys, *vals, *keys2, *vals2, *start, *end;
    Jac *buckets, *partial;
    int segs = 256;
    CHECK(cudaMalloc(&d_s, n * sizeof(Fr)));
    CHECK(cudaMemcpy(d_s, scalars, n * sizeof(Fr), cudaMemcpyHostToDevice));
    CHECK(cudaMalloc(&keys, m * 4));
    CHECK(cudaMalloc(&vals, m * 4));
    CHECK(cudaMalloc(&keys2, m * 4));
    CHECK(cudaMalloc(&vals2, m * 4));
    CHECK(cudaMalloc(&start, nb * 4));
    CHECK(cudaMalloc(&end, nb * 4));
    CHECK(cudaMemset(start, 0, nb * 4));
    CHECK(cudaMemset(end, 0, nb * 4));
    CHECK(cudaMalloc(&buckets, nb * sizeof(Jac)));
    CHECK(cudaMalloc(&partial, (size_t)windows * segs * sizeof(Jac)));
    prof("alloc");
    k_keys<<<(n + 255) / 256, 256>>>(d_s, g_cache.inf, n, c, windows, keys, vals);
    CHECK(cudaGetLastError());
    void* tmp = nullptr;
    size_t tmp_bytes = 0;
    int key_bits = c;
    while ((1 << (key_bits - c)) < windows) key_bits++;
    cub::DeviceRadixSort::SortPairs(tmp, tmp_bytes, keys, keys2, vals, vals2, (int)m, 0, key_bits);
    CHECK(cudaMalloc(&tmp, tmp_bytes));
    CHECK(cub::DeviceRadixSort::SortPairs(tmp, tmp_bytes, keys, keys2, vals, vals2, (int)m, 0, key_bits));
    prof("sort");
    k_bounds<<<(m + 255) / 256, 256>>>(keys2, m, start, end);
    CHECK(cudaGetLastError());
    // Montgomery one: R mod p
    Fq one = {{0x0002fffd, 0x76090000, 0xc40c0002, 0xebf4000b, 0x53c758ba, 0x5f489857, 0x70525745, 0x77ce5853, 0xa256ec6d, 0x5c071a97, 0xfa80e493, 0x15f65ec3}};
    prof("bounds");
    k_buckets<<<(nb + 127) / 128, 128>>>(start, end, vals2, g_cache.pts, nb, c, buckets, one);
    CHECK(cudaGetLastError());
    prof("buckets");
    k_reduce<<<(windows * segs + 63) / 64, 64>>>(buckets, c, windows, segs, partial);
    CHECK(cudaGetLastError());
    prof("reduce");
    // windows * segs Jacobian points back to the CPU, which adds them up
    CHECK(cudaMemcpy(out, partial, (size_t)windows * segs * sizeof(Jac), cudaMemcpyDeviceToHost));
    *out_count = (size_t)windows * segs;
    out_count[1] = (size_t)c;
    out_count[2] = (size_t)segs;
    cudaFree(tmp);
    cudaFree(d_s);
    cudaFree(keys);
    cudaFree(vals);
    cudaFree(keys2);
    cudaFree(vals2);
    cudaFree(start);
    cudaFree(end);
    cudaFree(buckets);
    cudaFree(partial);
    return 0;
}

extern "C" size_t veil_gpu_msm_out_points(size_t n) {
    int c = window_bits(n);
    return (size_t)((255 + c - 1) / c) * 256;
}

// ---------------- NTT ----------------

__global__ void k_twiddles(Fr* tw, size_t half, Fr w) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i >= half) return;
    // w^i by square-and-multiply; R mod r is the Montgomery one
    Fr r = {{0xfffffffe, 0x00000001, 0x00034802, 0x5884b7fa, 0xecbc4ff5, 0x998c4fef, 0xacc5056f, 0x1824b159}};
    Fr b = w;
    size_t e = i;
    while (e) {
        if (e & 1) r = rmul(r, b);
        b = rmul(b, b);
        e >>= 1;
    }
    tw[i] = r;
}

__global__ void k_bitrev(Fr* a, size_t n, int logn) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i >= n) return;
    size_t j = __brevll((unsigned long long)i) >> (64 - logn);
    if (i < j) {
        Fr t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
}

__global__ void k_stage(Fr* a, const Fr* tw, size_t n, size_t m) {
    size_t id = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (id >= n / 2) return;
    size_t j = id % m;
    size_t lo = (id / m) * 2 * m + j;
    Fr x = rmul(a[lo + m], tw[j * (n / (2 * m))]);
    Fr u = a[lo];
    a[lo] = radd(u, x);
    a[lo + m] = rsub(u, x);
}

struct NttCache {
    size_t n = 0;
    u64 w[4] = {0, 0, 0, 0};
    Fr* tw = nullptr;
};
static NttCache g_ntt[4];
static int g_ntt_next = 0;
static Fr* g_buf = nullptr;
static size_t g_buf_n = 0;

extern "C" int veil_gpu_ntt(u64* data, size_t n, const u64* omega) {
    int logn = 0;
    while (((size_t)1 << logn) < n) logn++;
    if (((size_t)1 << logn) != n || n < 2) return 1;
    Fr* tw = nullptr;
    for (auto& e : g_ntt) {
        if (e.n == n && memcmp(e.w, omega, 32) == 0) tw = e.tw;
    }
    if (!tw) {
        NttCache& e = g_ntt[g_ntt_next++ % 4];
        cudaFree(e.tw);
        CHECK(cudaMalloc(&e.tw, n / 2 * sizeof(Fr)));
        Fr w;
        memcpy(w.v, omega, 32);
        k_twiddles<<<(n / 2 + 255) / 256, 256>>>(e.tw, n / 2, w);
        CHECK(cudaGetLastError());
        e.n = n;
        memcpy(e.w, omega, 32);
        tw = e.tw;
    }
    if (g_buf_n < n) {
        cudaFree(g_buf);
        CHECK(cudaMalloc(&g_buf, n * sizeof(Fr)));
        g_buf_n = n;
    }
    CHECK(cudaMemcpy(g_buf, data, n * sizeof(Fr), cudaMemcpyHostToDevice));
    k_bitrev<<<(n + 255) / 256, 256>>>(g_buf, n, logn);
    for (size_t m = 1; m < n; m *= 2) k_stage<<<(n / 2 + 255) / 256, 256>>>(g_buf, tw, n, m);
    CHECK(cudaGetLastError());
    CHECK(cudaMemcpy(data, g_buf, n * sizeof(Fr), cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" int veil_gpu_device(char* name, size_t len) {
    cudaDeviceProp p;
    if (cudaGetDeviceProperties(&p, 0) != cudaSuccess) return 1;
    snprintf(name, len, "%s", p.name);
    return 0;
}
