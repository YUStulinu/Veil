//! Multi-scalar multiplication sum k_i P_i: Pippenger's bucket method, one
//! thread per window, plus the CUDA version when it is built in.

use crate::curve::{Affine, CurveField, Projective};
use crate::ff::Fr;

fn window_bits(n: usize) -> usize {
    if n < 32 {
        3
    } else {
        ((n as f64).ln().ceil() as usize).clamp(4, 16)
    }
}

/// the c-bit digit of k starting at bit `start`
#[inline]
fn digit(k: &[u64; 4], start: usize, c: usize) -> usize {
    let limb = start / 64;
    let off = start % 64;
    if limb >= 4 {
        return 0;
    }
    let mut v = k[limb] >> off;
    if off + c > 64 && limb + 1 < 4 {
        v |= k[limb + 1] << (64 - off);
    }
    (v as usize) & ((1 << c) - 1)
}

pub fn msm_cpu<F: CurveField + Send + Sync>(points: &[Affine<F>], scalars: &[Fr]) -> Projective<F> {
    assert!(points.len() >= scalars.len(), "more scalars than points");
    let n = scalars.len();
    if n == 0 {
        return Projective::IDENTITY;
    }
    let ks: Vec<[u64; 4]> = scalars.iter().map(|s| s.to_limbs()).collect();
    let c = window_bits(n);
    let windows = 255usize.div_ceil(c);
    let window_sum = |w: usize| -> Projective<F> {
        let mut buckets = vec![Projective::<F>::IDENTITY; (1 << c) - 1];
        for i in 0..n {
            let d = digit(&ks[i], w * c, c);
            if d != 0 {
                buckets[d - 1] = buckets[d - 1].add_affine(&points[i]);
            }
        }
        // sum_j j * B_j as a running sum from the top
        let mut running = Projective::IDENTITY;
        let mut total = Projective::IDENTITY;
        for b in buckets.iter().rev() {
            running = running.add(b);
            total = total.add(&running);
        }
        total
    };
    let threads = std::thread::available_parallelism().map_or(1, |t| t.get()).min(windows);
    let sums: Vec<Projective<F>> = if threads <= 1 || n < 64 || cfg!(target_arch = "wasm32") {
        (0..windows).map(window_sum).collect()
    } else {
        let mut out = vec![Projective::IDENTITY; windows];
        std::thread::scope(|s| {
            for (t, chunk) in out.chunks_mut(windows.div_ceil(threads)).enumerate() {
                let ws = &window_sum;
                let base = t * windows.div_ceil(threads);
                s.spawn(move || {
                    for (j, o) in chunk.iter_mut().enumerate() {
                        *o = ws(base + j);
                    }
                });
            }
        });
        out
    };
    let mut acc = Projective::IDENTITY;
    for s in sums.iter().rev() {
        for _ in 0..c {
            acc = acc.double();
        }
        acc = acc.add(s);
    }
    acc
}

/// G1 MSM on the GPU when built with CUDA and asked for, else on the CPU
pub fn msm_g1(points: &[crate::curve::G1Affine], scalars: &[Fr]) -> crate::curve::G1 {
    #[cfg(feature = "cuda")]
    if crate::gpu::enabled() && scalars.len() >= 1 << 12 {
        return crate::gpu::msm_g1(points, scalars);
    }
    msm_cpu(points, scalars)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::g1_generator;
    use crate::ff::Field;
    use crate::rng::Rng;

    #[test]
    fn matches_naive() {
        let mut rng = Rng::from_seed([7; 32]);
        let g = g1_generator().to_projective();
        for n in [1usize, 5, 33, 200] {
            let pts: Vec<_> = (0..n).map(|_| g.mul(&Fr::random(&mut rng))).collect();
            let pts = Projective::batch_to_affine(&pts);
            let ks: Vec<Fr> = (0..n).map(|i| if i == 2 { Fr::ZERO } else { Fr::random(&mut rng) }).collect();
            let naive = pts.iter().zip(&ks).fold(Projective::IDENTITY, |a, (p, k)| a.add(&p.to_projective().mul(k)));
            assert_eq!(msm_cpu(&pts, &ks), naive);
        }
    }
}
