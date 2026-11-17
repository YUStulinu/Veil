//! Polynomials over Fr: the radix-2 NTT on a multiplicative subgroup of size
//! 2^k (and on its cosets), plus the few coefficient-form helpers PLONK needs.

use crate::ff::{Field, Fr};

#[derive(Clone, Debug)]
pub struct Domain {
    pub n: usize,
    pub log_n: u32,
    pub omega: Fr,
    pub omega_inv: Fr,
    pub n_inv: Fr,
}

impl Domain {
    pub fn new(n: usize) -> Domain {
        assert!(n.is_power_of_two(), "domain size must be a power of two");
        let log_n = n.trailing_zeros();
        let omega = Fr::root_of_unity(log_n);
        Domain { n, log_n, omega, omega_inv: omega.inverse().unwrap(), n_inv: Fr::from_u64(n as u64).inverse().unwrap() }
    }

    /// omega^0 .. omega^(n-1)
    pub fn elements(&self) -> Vec<Fr> {
        powers(self.omega, self.n)
    }

    /// coefficients -> evaluations at omega^i
    pub fn fft(&self, a: &mut [Fr]) {
        assert_eq!(a.len(), self.n);
        ntt(a, self.omega);
    }

    /// evaluations -> coefficients
    pub fn ifft(&self, a: &mut [Fr]) {
        assert_eq!(a.len(), self.n);
        ntt(a, self.omega_inv);
        scale(a, self.n_inv);
    }

    /// coefficients -> evaluations at g omega^i
    pub fn coset_fft(&self, a: &mut [Fr], g: Fr) {
        distribute_powers(a, g);
        self.fft(a);
    }

    /// evaluations at g omega^i -> coefficients
    pub fn coset_ifft(&self, a: &mut [Fr], g: Fr) {
        self.ifft(a);
        distribute_powers(a, g.inverse().unwrap());
    }

    /// Z_H(z) = z^n - 1
    pub fn vanishing_at(&self, z: Fr) -> Fr {
        pow_u64(z, self.n as u64) - Fr::ONE
    }

    /// L_i(z) = omega^i (z^n - 1) / (n (z - omega^i)), for z outside the domain
    pub fn lagrange_at(&self, i: usize, z: Fr) -> Fr {
        let wi = pow_u64(self.omega, i as u64);
        wi * self.vanishing_at(z) * self.n_inv * (z - wi).inverse().expect("z is in the domain")
    }
}

pub fn pow_u64(x: Fr, e: u64) -> Fr {
    x.pow(&[e])
}

pub fn powers(x: Fr, n: usize) -> Vec<Fr> {
    let mut v = Vec::with_capacity(n);
    let mut acc = Fr::ONE;
    for _ in 0..n {
        v.push(acc);
        acc *= x;
    }
    v
}

pub fn scale(a: &mut [Fr], k: Fr) {
    for x in a {
        *x *= k;
    }
}

/// a_i *= g^i
pub fn distribute_powers(a: &mut [Fr], g: Fr) {
    let mut acc = Fr::ONE;
    for x in a {
        *x *= acc;
        acc *= g;
    }
}

fn threads() -> usize {
    if cfg!(target_arch = "wasm32") {
        1
    } else {
        std::thread::available_parallelism().map_or(1, |t| t.get())
    }
}

/// in-place NTT with root `w` (natural order in and out)
pub fn ntt(a: &mut [Fr], w: Fr) {
    #[cfg(feature = "cuda")]
    if crate::gpu::enabled() && a.len() >= 1 << 14 {
        return crate::gpu::ntt(a, w);
    }
    ntt_cpu(a, w)
}

pub fn ntt_cpu(a: &mut [Fr], w: Fr) {
    let n = a.len();
    let log_n = n.trailing_zeros();
    if n <= 1 {
        return;
    }
    for i in 0..n {
        let j = i.reverse_bits() >> (usize::BITS - log_n);
        if i < j {
            a.swap(i, j);
        }
    }
    // twiddles of the last stage; the stage with half-size m uses every (n/2m)-th
    let tw = powers(w, n / 2);
    // stages with half-size m in [m0, m1) on a slice made of whole groups
    let stages = |chunk: &mut [Fr], m0: usize, m1: usize| {
        let mut m = m0;
        while m < m1 {
            let stride = n / (2 * m);
            for g in chunk.chunks_mut(2 * m) {
                let (lo, hi) = g.split_at_mut(m);
                for j in 0..m {
                    let x = hi[j] * tw[j * stride];
                    hi[j] = lo[j] - x;
                    lo[j] += x;
                }
            }
            m *= 2;
        }
    };
    let t = if n >= 1 << 14 { threads().next_power_of_two().min(n / 2) } else { 1 };
    if t == 1 {
        stages(a, 1, n);
        return;
    }
    // the early stages inside n/t-sized chunks: one thread per chunk, spawned once
    let chunk = n / t;
    std::thread::scope(|s| {
        for c in a.chunks_mut(chunk) {
            let stages = &stages;
            s.spawn(move || stages(c, 1, chunk));
        }
    });
    // the last log2(t) stages: split the j range of every group across threads
    let mut m = chunk;
    while m < n {
        let stride = n / (2 * m);
        let per = (n / (2 * t)).max(1); // t threads per stage in all
        std::thread::scope(|s| {
            for g in a.chunks_mut(2 * m) {
                let (lo, hi) = g.split_at_mut(m);
                for (k, (l, h)) in lo.chunks_mut(per).zip(hi.chunks_mut(per)).enumerate() {
                    let tw = &tw;
                    s.spawn(move || {
                        for j in 0..l.len() {
                            let x = h[j] * tw[(k * per + j) * stride];
                            h[j] = l[j] - x;
                            l[j] += x;
                        }
                    });
                }
            }
        });
        m *= 2;
    }
}

/// p(z) by Horner
pub fn evaluate(coeffs: &[Fr], z: Fr) -> Fr {
    coeffs.iter().rev().fold(Fr::ZERO, |acc, c| acc * z + *c)
}

/// p(X) / (X - z), dropping the remainder (exact when p(z) = 0)
pub fn divide_by_linear(coeffs: &[Fr], z: Fr) -> Vec<Fr> {
    let n = coeffs.len();
    if n <= 1 {
        return vec![];
    }
    let mut q = vec![Fr::ZERO; n - 1];
    let mut carry = Fr::ZERO;
    for i in (1..n).rev() {
        carry = coeffs[i] + carry * z;
        q[i - 1] = carry;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn ntt_matches_evaluation_and_inverts() {
        let mut rng = Rng::from_seed([1; 32]);
        for log_n in [1u32, 3, 8, 13, 15] {
            let d = Domain::new(1 << log_n);
            let p: Vec<Fr> = (0..d.n).map(|_| Fr::random(&mut rng)).collect();
            let mut e = p.clone();
            d.fft(&mut e);
            for i in [0usize, 1, d.n - 1] {
                assert_eq!(e[i], evaluate(&p, pow_u64(d.omega, i as u64)));
            }
            let g = Fr::from_u64(7);
            let mut c = p.clone();
            d.coset_fft(&mut c, g);
            assert_eq!(c[1], evaluate(&p, g * d.omega));
            d.coset_ifft(&mut c, g);
            assert_eq!(c, p);
            d.ifft(&mut e);
            assert_eq!(e, p);
        }
    }

    #[test]
    fn lagrange_and_division() {
        let d = Domain::new(8);
        let z = Fr::from_u64(12345);
        // sum_i L_i(z) = 1
        let s = (0..8).fold(Fr::ZERO, |a, i| a + d.lagrange_at(i, z));
        assert_eq!(s, Fr::ONE);
        let p = vec![Fr::from_u64(3), Fr::from_u64(1), Fr::from_u64(4), Fr::from_u64(1)];
        let mut shifted = p.clone();
        shifted[0] -= evaluate(&p, z);
        let q = divide_by_linear(&shifted, z);
        let w = Fr::from_u64(99);
        assert_eq!(evaluate(&q, w) * (w - z), evaluate(&shifted, w));
    }
}
