// BLS12-381 in TypeScript with BigInt, written apart from the Rust code (the
// same tower and pairing algorithm, so the two check each other): Fq, Fr,
// Fq2 = Fq[u]/(u^2+1), Fq6 = Fq2[v]/(v^3-(u+1)), Fq12 = Fq6[w]/(w^2-v),
// G1/G2 in Jacobian coordinates, the optimal ate pairing.

export const P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaabn;
export const R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n;

export const mod = (a: bigint, m: bigint = P): bigint => {
  const r = a % m;
  return r < 0n ? r + m : r;
};

/** a^-1 mod m by the extended Euclidean algorithm */
export function inv(a: bigint, m: bigint = P): bigint {
  let [t, nt, r, nr] = [0n, 1n, m, mod(a, m)];
  while (nr !== 0n) {
    const q = r / nr;
    [t, nt] = [nt, t - q * nt];
    [r, nr] = [nr, r - q * nr];
  }
  if (r !== 1n) throw new Error("not invertible");
  return mod(t, m);
}

export function powMod(b: bigint, e: bigint, m: bigint): bigint {
  let r = 1n;
  b = mod(b, m);
  while (e > 0n) {
    if (e & 1n) r = (r * b) % m;
    b = (b * b) % m;
    e >>= 1n;
  }
  return r;
}

// ---- Fq2 ----
export type F2 = [bigint, bigint];
export const F2_ZERO: F2 = [0n, 0n];
export const F2_ONE: F2 = [1n, 0n];
export const f2add = (a: F2, b: F2): F2 => [mod(a[0] + b[0]), mod(a[1] + b[1])];
export const f2sub = (a: F2, b: F2): F2 => [mod(a[0] - b[0]), mod(a[1] - b[1])];
export const f2neg = (a: F2): F2 => [mod(-a[0]), mod(-a[1])];
export const f2mul = (a: F2, b: F2): F2 => {
  const aa = a[0] * b[0];
  const bb = a[1] * b[1];
  return [mod(aa - bb), mod((a[0] + a[1]) * (b[0] + b[1]) - aa - bb)];
};
export const f2sqr = (a: F2): F2 => [mod((a[0] + a[1]) * (a[0] - a[1])), mod(2n * a[0] * a[1])];
export const f2scale = (a: F2, k: bigint): F2 => [mod(a[0] * k), mod(a[1] * k)];
export const f2inv = (a: F2): F2 => {
  const t = inv(a[0] * a[0] + a[1] * a[1]);
  return [mod(a[0] * t), mod(-a[1] * t)];
};
/** times xi = u + 1 */
export const f2nr = (a: F2): F2 => [mod(a[0] - a[1]), mod(a[0] + a[1])];
export const f2eq = (a: F2, b: F2) => a[0] === b[0] && a[1] === b[1];
export function f2pow(a: F2, e: bigint): F2 {
  let r = F2_ONE;
  let b = a;
  while (e > 0n) {
    if (e & 1n) r = f2mul(r, b);
    b = f2sqr(b);
    e >>= 1n;
  }
  return r;
}

// ---- Fq6 ----
export type F6 = [F2, F2, F2];
const F6_ZERO: F6 = [F2_ZERO, F2_ZERO, F2_ZERO];
const F6_ONE: F6 = [F2_ONE, F2_ZERO, F2_ZERO];
const f6add = (a: F6, b: F6): F6 => [f2add(a[0], b[0]), f2add(a[1], b[1]), f2add(a[2], b[2])];
const f6sub = (a: F6, b: F6): F6 => [f2sub(a[0], b[0]), f2sub(a[1], b[1]), f2sub(a[2], b[2])];
const f6neg = (a: F6): F6 => [f2neg(a[0]), f2neg(a[1]), f2neg(a[2])];
/** times v */
const f6nr = (a: F6): F6 => [f2nr(a[2]), a[0], a[1]];
function f6mul(a: F6, b: F6): F6 {
  const aa = f2mul(a[0], b[0]);
  const bb = f2mul(a[1], b[1]);
  const cc = f2mul(a[2], b[2]);
  const t1 = f2add(f2nr(f2sub(f2sub(f2mul(f2add(a[1], a[2]), f2add(b[1], b[2])), bb), cc)), aa);
  const t2 = f2add(f2sub(f2sub(f2mul(f2add(a[0], a[1]), f2add(b[0], b[1])), aa), bb), f2nr(cc));
  const t3 = f2add(f2sub(f2sub(f2mul(f2add(a[0], a[2]), f2add(b[0], b[2])), aa), cc), bb);
  return [t1, t2, t3];
}
function f6inv(a: F6): F6 {
  const c0 = f2sub(f2sqr(a[0]), f2nr(f2mul(a[1], a[2])));
  const c1 = f2sub(f2nr(f2sqr(a[2])), f2mul(a[0], a[1]));
  const c2 = f2sub(f2sqr(a[1]), f2mul(a[0], a[2]));
  const t = f2inv(f2add(f2add(f2mul(a[0], c0), f2nr(f2mul(a[2], c1))), f2nr(f2mul(a[1], c2))));
  return [f2mul(c0, t), f2mul(c1, t), f2mul(c2, t)];
}

// ---- Fq12 ----
export type F12 = [F6, F6];
export const F12_ONE: F12 = [F6_ONE, F6_ZERO];
export function f12mul(a: F12, b: F12): F12 {
  const aa = f6mul(a[0], b[0]);
  const bb = f6mul(a[1], b[1]);
  return [f6add(aa, f6nr(bb)), f6sub(f6sub(f6mul(f6add(a[0], a[1]), f6add(b[0], b[1])), aa), bb)];
}
const f12inv = (a: F12): F12 => {
  const t = f6inv(f6sub(f6mul(a[0], a[0]), f6nr(f6mul(a[1], a[1]))));
  return [f6mul(a[0], t), f6neg(f6mul(a[1], t))];
};
const f12conj = (a: F12): F12 => [a[0], f6neg(a[1])];
export const f12eq = (a: F12, b: F12) => a.every((x, i) => x.every((y, j) => f2eq(y, b[i][j])));

let gamma2: F2[] | null = null;
/** f^(q^2): coefficient of w^i times xi^(i (q^2 - 1) / 6) */
function f12frob2(a: F12): F12 {
  if (!gamma2) {
    const g = f2pow([1n, 1n], (P * P - 1n) / 6n);
    gamma2 = [F2_ONE];
    for (let i = 1; i < 6; i++) gamma2.push(f2mul(gamma2[i - 1], g));
  }
  const g = gamma2;
  return [
    [a[0][0], f2mul(a[0][1], g[2]), f2mul(a[0][2], g[4])],
    [f2mul(a[1][0], g[1]), f2mul(a[1][1], g[3]), f2mul(a[1][2], g[5])],
  ];
}
export function f12pow(a: F12, e: bigint): F12 {
  let r = F12_ONE;
  for (let i = BigInt(e.toString(2).length) - 1n; i >= 0n; i--) {
    r = f12mul(r, r);
    if ((e >> i) & 1n) r = f12mul(r, a);
  }
  return r;
}

// ---- curves ----
export interface G1 {
  x: bigint;
  y: bigint;
  z: bigint;
} // Jacobian, z = 0: infinity
export interface G1Aff {
  x: bigint;
  y: bigint;
  inf: boolean;
}
export interface G2Aff {
  x: F2;
  y: F2;
  inf: boolean;
}

export const G1_ZERO: G1 = { x: 1n, y: 1n, z: 0n };
export const G1_GEN: G1Aff = {
  x: 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bbn,
  y: 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1n,
  inf: false,
};

export const toJac = (p: G1Aff): G1 => (p.inf ? G1_ZERO : { x: p.x, y: p.y, z: 1n });
export function toAff(p: G1): G1Aff {
  if (p.z === 0n) return { x: 0n, y: 0n, inf: true };
  const zi = inv(p.z);
  const zi2 = mod(zi * zi);
  return { x: mod(p.x * zi2), y: mod(p.y * zi2 * zi), inf: false };
}
export function g1dbl(p: G1): G1 {
  if (p.z === 0n) return p;
  const a = mod(p.x * p.x), b = mod(p.y * p.y), c = mod(b * b);
  const d = mod(2n * ((p.x + b) ** 2n - a - c));
  const e = 3n * a;
  const f = mod(e * e);
  const x3 = mod(f - 2n * d);
  return { x: x3, y: mod(e * (d - x3) - 8n * c), z: mod(2n * p.y * p.z) };
}
export function g1add(p: G1, q: G1): G1 {
  if (p.z === 0n) return q;
  if (q.z === 0n) return p;
  const z1z1 = mod(p.z * p.z), z2z2 = mod(q.z * q.z);
  const u1 = mod(p.x * z2z2), u2 = mod(q.x * z1z1);
  const s1 = mod(p.y * q.z * z2z2), s2 = mod(q.y * p.z * z1z1);
  const h = mod(u2 - u1), r = mod(2n * (s2 - s1));
  if (h === 0n) return r === 0n ? g1dbl(p) : G1_ZERO;
  const i = mod(4n * h * h), j = mod(h * i), v = mod(u1 * i);
  const x3 = mod(r * r - j - 2n * v);
  return { x: x3, y: mod(r * (v - x3) - 2n * s1 * j), z: mod(((p.z + q.z) ** 2n - z1z1 - z2z2) * h) };
}
export function g1mul(p: G1, k: bigint): G1 {
  let acc = G1_ZERO;
  for (let i = BigInt(k.toString(2).length) - 1n; i >= 0n; i--) {
    acc = g1dbl(acc);
    if ((k >> i) & 1n) acc = g1add(acc, p);
  }
  return k === 0n ? G1_ZERO : acc;
}
export const g1neg = (p: G1Aff): G1Aff => ({ x: p.x, y: p.inf ? 0n : mod(-p.y), inf: p.inf });
export const g1onCurve = (p: G1Aff) => p.inf || mod(p.y * p.y - p.x ** 3n - 4n) === 0n;
export const g1inSubgroup = (p: G1Aff) => g1mul(toJac(p), R).z === 0n;

const g2onCurve = (p: G2Aff) => p.inf || f2eq(f2sqr(p.y), f2add(f2mul(f2sqr(p.x), p.x), [4n, 4n]));

// ---- encodings ----
export const bytesToBig = (b: Uint8Array): bigint => b.reduce((a, x) => (a << 8n) | BigInt(x), 0n);
export function bigToBytes(x: bigint, len: number): Uint8Array {
  const out = new Uint8Array(len);
  for (let i = len - 1; i >= 0; i--) {
    out[i] = Number(x & 0xffn);
    x >>= 8n;
  }
  return out;
}

/** the ZCash compressed format; checks the curve and the subgroup */
export function g1FromCompressed(b: Uint8Array): G1Aff {
  if (b.length !== 48) throw new Error("length");
  const [c, inf, sign] = [(b[0] & 0x80) !== 0, (b[0] & 0x40) !== 0, (b[0] & 0x20) !== 0];
  if (!c) throw new Error("not compressed");
  const xb = b.slice();
  xb[0] &= 0x1f;
  const x = bytesToBig(xb);
  if (inf) {
    if (sign || x !== 0n) throw new Error("bad infinity");
    return { x: 0n, y: 0n, inf: true };
  }
  if (x >= P) throw new Error("x not canonical");
  const y2 = mod(x ** 3n + 4n);
  let y = powMod(y2, (P + 1n) / 4n, P);
  if (mod(y * y) !== y2) throw new Error("not on the curve");
  if (y > (P - 1n) / 2n !== sign) y = mod(-y);
  const p = { x, y, inf: false };
  if (!g1inSubgroup(p)) throw new Error("not in the subgroup");
  return p;
}

/** EIP-2537 (64-byte padded coordinates), as in the verifying key and the transcript */
export function g1FromEip2537(b: Uint8Array): G1Aff {
  if (b.every((v) => v === 0)) return { x: 0n, y: 0n, inf: true };
  const p = { x: bytesToBig(b.slice(16, 64)), y: bytesToBig(b.slice(80, 128)), inf: false };
  if (p.x >= P || p.y >= P || !g1onCurve(p) || !g1inSubgroup(p)) throw new Error("bad G1 point");
  return p;
}
export function g1ToEip2537(p: G1Aff): Uint8Array {
  const out = new Uint8Array(128);
  if (!p.inf) {
    out.set(bigToBytes(p.x, 48), 16);
    out.set(bigToBytes(p.y, 48), 80);
  }
  return out;
}
export function g2FromEip2537(b: Uint8Array): G2Aff {
  const f = (o: number) => bytesToBig(b.slice(o + 16, o + 64));
  const p: G2Aff = { x: [f(0), f(64)], y: [f(128), f(192)], inf: false };
  if (!g2onCurve(p)) throw new Error("bad G2 point");
  return p;
}

// ---- pairing ----
const X = 0xd201000000010000n; // |x|, x < 0
const HARD = (P ** 4n - P ** 2n + 1n) / R;

/** the line with slope l through t, at p, times w^3 */
function line(l: F2, t: G2Aff, p: G1Aff): F12 {
  return [
    [f2sub(f2mul(l, t.x), t.y), f2scale(l, mod(-p.x)), F2_ZERO],
    [F2_ZERO, [p.y, 0n], F2_ZERO],
  ];
}

function millerLoop(pairs: [G1Aff, G2Aff][]): F12 {
  const ps = pairs.filter(([p, q]) => !p.inf && !q.inf);
  const ts = ps.map(([, q]) => ({ ...q }));
  let f = F12_ONE;
  const bits = X.toString(2);
  for (let i = 1; i < bits.length; i++) {
    f = f12mul(f, f);
    ps.forEach(([p], k) => {
      const t = ts[k];
      const l = f2mul(f2scale(f2sqr(t.x), 3n), f2inv(f2scale(t.y, 2n)));
      f = f12mul(f, line(l, t, p));
      const x3 = f2sub(f2sqr(l), f2scale(t.x, 2n));
      ts[k] = { x: x3, y: f2sub(f2mul(l, f2sub(t.x, x3)), t.y), inf: false };
    });
    if (bits[i] === "1") {
      ps.forEach(([p, q], k) => {
        const t = ts[k];
        const l = f2mul(f2sub(q.y, t.y), f2inv(f2sub(q.x, t.x)));
        f = f12mul(f, line(l, t, p));
        const x3 = f2sub(f2sub(f2sqr(l), t.x), q.x);
        ts[k] = { x: x3, y: f2sub(f2mul(l, f2sub(t.x, x3)), t.y), inf: false };
      });
    }
  }
  return f12conj(f);
}

function finalExp(f: F12): F12 {
  const f1 = f12mul(f12conj(f), f12inv(f));
  const f2 = f12mul(f12frob2(f1), f1);
  return f12pow(f2, HARD);
}

/** prod e(p_i, q_i) == 1 */
export const pairingCheck = (pairs: [G1Aff, G2Aff][]) => f12eq(finalExp(millerLoop(pairs)), F12_ONE);
export const pairing = (p: G1Aff, q: G2Aff) => finalExp(millerLoop([[p, q]]));
