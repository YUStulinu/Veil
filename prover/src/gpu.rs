//! The CUDA kernels (cuda/veil.cu): G1 MSM and the Fr NTT. Field elements
//! cross the boundary in the CPU's own Montgomery layout; the GPU returns
//! per-window partial sums of the MSM and the CPU adds them up.

use crate::curve::{G1Affine, Projective, G1};
use crate::ff::{Fq, Fr};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

extern "C" {
    fn veil_gpu_msm_g1(points: *const u8, stride: usize, scalars: *const u64, n: usize, out: *mut u64, layout: *mut usize) -> i32;
    fn veil_gpu_msm_out_points(n: usize) -> usize;
    fn veil_gpu_ntt(data: *mut u64, n: usize, omega: *const u64) -> i32;
    fn veil_gpu_device(name: *mut u8, len: usize) -> i32;
}

static ENABLED: AtomicBool = AtomicBool::new(false);
/// the kernels keep device buffers between calls: one caller at a time
static LOCK: Mutex<()> = Mutex::new(());

pub fn enable() {
    ENABLED.store(true, Ordering::SeqCst);
}

pub fn enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

pub fn device_name() -> Option<String> {
    let mut b = [0u8; 256];
    (unsafe { veil_gpu_device(b.as_mut_ptr(), b.len()) } == 0).then(|| String::from_utf8_lossy(&b[..b.iter().position(|c| *c == 0).unwrap_or(0)]).into())
}

pub fn msm_g1(points: &[G1Affine], scalars: &[Fr]) -> G1 {
    let _g = LOCK.lock().unwrap();
    let n = scalars.len();
    let count = unsafe { veil_gpu_msm_out_points(n) };
    let mut out = vec![0u64; count * 18];
    let mut layout = [0usize; 3];
    let rc = unsafe {
        veil_gpu_msm_g1(points.as_ptr() as *const u8, std::mem::size_of::<G1Affine>(), scalars.as_ptr() as *const u64, n, out.as_mut_ptr(), layout.as_mut_ptr())
    };
    assert_eq!(rc, 0, "GPU MSM failed");
    let (c, segs) = (layout[1], layout[2]);
    let windows = layout[0] / segs;
    let point = |k: usize| {
        let o = &out[18 * k..18 * k + 18];
        let f = |i: usize| Fq(o[6 * i..6 * i + 6].try_into().unwrap());
        Projective { x: f(0), y: f(1), z: f(2) }
    };
    let mut acc = G1::IDENTITY;
    for w in (0..windows).rev() {
        for _ in 0..c {
            acc = acc.double();
        }
        for t in 0..segs {
            acc = acc.add(&point(w * segs + t));
        }
    }
    acc
}

pub fn ntt(a: &mut [Fr], w: Fr) {
    let _g = LOCK.lock().unwrap();
    let rc = unsafe { veil_gpu_ntt(a.as_mut_ptr() as *mut u64, a.len(), &w as *const Fr as *const u64) };
    assert_eq!(rc, 0, "GPU NTT failed");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::g1_generator;
    use crate::rng::Rng;

    #[test]
    fn gpu_matches_cpu() {
        let Some(name) = device_name() else {
            eprintln!("no CUDA device: skipped");
            return;
        };
        eprintln!("GPU: {name}");
        let mut rng = Rng::from_seed([11; 32]);
        let g = g1_generator().to_projective();
        for n in [1usize << 10, 3000, 1 << 13] {
            let pts: Vec<_> = (0..n).map(|i| g.mul(&Fr::from_u64(i as u64 * 7919 + 1))).collect();
            let mut pts = Projective::batch_to_affine(&pts);
            pts[5] = G1Affine::IDENTITY;
            let ks: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
            assert_eq!(msm_g1(&pts, &ks), crate::msm::msm_cpu(&pts, &ks), "MSM n = {n}");
            let w = Fr::root_of_unity(n.next_power_of_two().trailing_zeros());
            let mut a: Vec<Fr> = (0..n.next_power_of_two()).map(|_| Fr::random(&mut rng)).collect();
            let mut b = a.clone();
            ntt(&mut a, w);
            crate::poly::ntt_cpu(&mut b, w);
            assert_eq!(a, b, "NTT n = {n}");
        }
    }
}
