//! The PLONK prover (Gabizon, Williamson, Ciobotaru 2019, with the
//! linearization polynomial r(X) and r(zeta) = 0), round by round.

use super::proof::Proof;
use super::setup::{ProvingKey, COSET};
use super::transcript::Transcript;
use crate::curve::{G1Affine, Projective, G1};
use crate::ff::{Field, Fr};
use crate::msm::msm_g1;
use crate::poly::{divide_by_linear, evaluate, pow_u64};
use crate::rng::Rng;

/// a clock in seconds that also works in the browser
fn clock() -> f64 {
    #[cfg(target_arch = "wasm32")]
    return crate::wasm::now_ms() / 1e3;
    #[cfg(not(target_arch = "wasm32"))]
    {
        use std::sync::OnceLock;
        static START: OnceLock<std::time::Instant> = OnceLock::new();
        START.get_or_init(std::time::Instant::now).elapsed().as_secs_f64()
    }
}

#[derive(Clone, Debug, Default)]
pub struct Timings {
    /// seconds per round
    pub rounds: [f64; 5],
}

/// acc += k p
fn add_scaled(acc: &mut Vec<Fr>, p: &[Fr], k: Fr) {
    if acc.len() < p.len() {
        acc.resize(p.len(), Fr::ZERO);
    }
    for (a, x) in acc.iter_mut().zip(p) {
        *a += k * *x;
    }
}

fn commit(pk: &ProvingKey, p: &[Fr]) -> G1 {
    assert!(p.len() <= pk.g1.len(), "polynomial longer than the SRS");
    msm_g1(&pk.g1[..p.len()], p)
}

/// + (b_0 + b_1 X + ...) (X^n - 1)
fn blind(mut p: Vec<Fr>, n: usize, bs: &[Fr]) -> Vec<Fr> {
    p.resize(n + bs.len(), Fr::ZERO);
    for (i, b) in bs.iter().enumerate() {
        p[n + i] += *b;
        p[i] -= *b;
    }
    p
}

/// the values of `f` at 4n coset points, computed on all threads
fn par_map(len: usize, f: impl Fn(usize) -> Fr + Sync) -> Vec<Fr> {
    let threads = if cfg!(target_arch = "wasm32") { 1 } else { std::thread::available_parallelism().map_or(1, |t| t.get()) };
    let mut out = vec![Fr::ZERO; len];
    if threads == 1 || len < 1 << 12 {
        for (i, o) in out.iter_mut().enumerate() {
            *o = f(i);
        }
        return out;
    }
    let per = len.div_ceil(threads);
    std::thread::scope(|s| {
        for (k, chunk) in out.chunks_mut(per).enumerate() {
            let f = &f;
            s.spawn(move || {
                for (j, o) in chunk.iter_mut().enumerate() {
                    *o = f(k * per + j);
                }
            });
        }
    });
    out
}

/// Proves that `w` (all variable values) satisfies the circuit behind `pk`.
/// It does not check the witness itself: a wrong one gives a proof that does
/// not verify (the soundness tests rely on this).
pub fn prove(pk: &ProvingKey, public: &[Fr], w: &[Fr], rng: &mut Rng) -> (Proof, Timings) {
    assert_eq!(public.len(), pk.num_public);
    let n = pk.n;
    let n4 = 4 * n;
    let mut tm = Timings::default();
    let mut tr = Transcript::new(&pk.vk.to_bytes());
    for x in public {
        tr.absorb_fr(x);
    }
    let b: Vec<Fr> = (0..11).map(|_| Fr::random(rng)).collect();
    let g = Fr::from_u64(COSET);
    let coset = |p: &[Fr]| {
        let mut e = p.to_vec();
        e.resize(n4, Fr::ZERO);
        pk.domain4.coset_fft(&mut e, g);
        e
    };

    // round 1: the wire polynomials
    let t0 = clock();
    let wire_values: [Vec<Fr>; 3] = std::array::from_fn(|col| pk.wires[col].iter().map(|v| w[*v as usize]).collect());
    let wire_polys: [Vec<Fr>; 3] = std::array::from_fn(|col| {
        let mut p = wire_values[col].clone();
        pk.domain.ifft(&mut p);
        blind(p, n, &[b[2 * col + 1], b[2 * col]])
    });
    let wire_c: Vec<G1> = wire_polys.iter().map(|p| commit(pk, p)).collect();
    let wire_c = Projective::batch_to_affine(&wire_c);
    for c in &wire_c {
        tr.absorb_g1(c);
    }
    let beta = tr.challenge();
    let gamma = tr.challenge();
    tm.rounds[0] = clock() - t0;

    // round 2: the permutation accumulator z
    let t0 = clock();
    let omegas = pk.domain.elements();
    let ks = [Fr::ONE, pk.vk.k1, pk.vk.k2];
    let mut num = vec![Fr::ONE; n];
    let mut den = vec![Fr::ONE; n];
    for i in 0..n {
        for col in 0..3 {
            let v = wire_values[col][i] + gamma;
            num[i] *= v + beta * ks[col] * omegas[i];
            den[i] *= v + beta * pk.sigma_values[col][i];
        }
    }
    Fr::batch_inverse(&mut den);
    let mut z = vec![Fr::ONE; n];
    for i in 0..n - 1 {
        z[i + 1] = z[i] * num[i] * den[i];
    }
    pk.domain.ifft(&mut z);
    let z_poly = blind(z, n, &[b[8], b[7], b[6]]);
    let z_c = commit(pk, &z_poly).to_affine();
    tr.absorb_g1(&z_c);
    let alpha = tr.challenge();
    tm.rounds[1] = clock() - t0;

    // round 3: the quotient t on the 4n coset (t has degree 3n + 5 < 4n)
    let t0 = clock();
    let mut pi_poly = vec![Fr::ZERO; n];
    for (i, x) in public.iter().enumerate() {
        pi_poly[i] = -*x;
    }
    pk.domain.ifft(&mut pi_poly);
    let [ac, bc, cc] = [coset(&wire_polys[0]), coset(&wire_polys[1]), coset(&wire_polys[2])];
    let zc = coset(&z_poly);
    let pic = coset(&pi_poly);
    let [qm, ql, qr, qo, qcst, s1, s2, s3] = &pk.cosets;
    let (k1, k2) = (pk.vk.k1, pk.vk.k2);
    let alpha2 = alpha.square();
    let t_evals = par_map(n4, |i| {
        let x = pk.x_coset[i];
        let (a, bb, c, zz) = (ac[i], bc[i], cc[i], zc[i]);
        let zw = zc[(i + 4) % n4];
        let gate = qm[i] * a * bb + ql[i] * a + qr[i] * bb + qo[i] * c + qcst[i] + pic[i];
        let bx = beta * x;
        let p1 = (a + bx + gamma) * (bb + bx * k1 + gamma) * (c + bx * k2 + gamma) * zz;
        let p2 = (a + beta * s1[i] + gamma) * (bb + beta * s2[i] + gamma) * (c + beta * s3[i] + gamma) * zw;
        (gate + alpha * (p1 - p2) + alpha2 * (zz - Fr::ONE) * pk.l0_coset[i]) * pk.zh_inv_coset[i % 4]
    });
    let mut t = t_evals;
    pk.domain4.coset_ifft(&mut t, g);
    t.truncate(3 * n + 6); // nonzero beyond this only for a wrong witness
    let mut t_lo = t[..n].to_vec();
    let mut t_mid = t[n..2 * n].to_vec();
    let mut t_hi = t[2 * n..].to_vec();
    t_lo.push(b[9]);
    t_mid[0] -= b[9];
    t_mid.push(b[10]);
    t_hi[0] -= b[10];
    let t_c = Projective::batch_to_affine(&[commit(pk, &t_lo), commit(pk, &t_mid), commit(pk, &t_hi)]);
    for c in &t_c {
        tr.absorb_g1(c);
    }
    let zeta = tr.challenge();
    tm.rounds[2] = clock() - t0;

    // round 4: evaluations
    let t0 = clock();
    let zeta_omega = zeta * pk.domain.omega;
    let ev = [
        evaluate(&wire_polys[0], zeta),
        evaluate(&wire_polys[1], zeta),
        evaluate(&wire_polys[2], zeta),
        evaluate(&pk.polys[5], zeta),
        evaluate(&pk.polys[6], zeta),
        evaluate(&z_poly, zeta_omega),
    ];
    for e in &ev {
        tr.absorb_fr(e);
    }
    let v = tr.challenge();
    tm.rounds[3] = clock() - t0;

    // round 5: the linearization polynomial and the two opening proofs
    let t0 = clock();
    let [ea, eb, ec, es1, es2, ezw] = ev;
    let zh_zeta = pk.domain.vanishing_at(zeta);
    let l0_zeta = pk.domain.lagrange_at(0, zeta);
    let pi_zeta = evaluate(&pi_poly, zeta);
    let mut r = vec![];
    add_scaled(&mut r, &pk.polys[0], ea * eb);
    add_scaled(&mut r, &pk.polys[1], ea);
    add_scaled(&mut r, &pk.polys[2], eb);
    add_scaled(&mut r, &pk.polys[3], ec);
    add_scaled(&mut r, &pk.polys[4], Fr::ONE);
    r[0] += pi_zeta;
    let perm_z = alpha * (ea + beta * zeta + gamma) * (eb + beta * k1 * zeta + gamma) * (ec + beta * k2 * zeta + gamma) + alpha2 * l0_zeta;
    add_scaled(&mut r, &z_poly, perm_z);
    let perm_s = alpha * (ea + beta * es1 + gamma) * (eb + beta * es2 + gamma) * ezw;
    add_scaled(&mut r, &pk.polys[7], -perm_s * beta);
    r[0] -= perm_s * (ec + gamma) + alpha2 * l0_zeta;
    let zn = pow_u64(zeta, n as u64);
    add_scaled(&mut r, &t_lo, -zh_zeta);
    add_scaled(&mut r, &t_mid, -zh_zeta * zn);
    add_scaled(&mut r, &t_hi, -zh_zeta * zn * zn);
    // r(zeta) = 0 exactly when the witness is right
    let mut num = r;
    let mut vk_ = v;
    for (p, e) in [(&wire_polys[0], ea), (&wire_polys[1], eb), (&wire_polys[2], ec), (&pk.polys[5], es1), (&pk.polys[6], es2)] {
        add_scaled(&mut num, p, vk_);
        num[0] -= vk_ * e;
        vk_ *= v;
    }
    let w_zeta = divide_by_linear(&num, zeta);
    let mut zs = z_poly.clone();
    zs[0] -= ezw;
    let w_zeta_omega = divide_by_linear(&zs, zeta_omega);
    let w_c = Projective::batch_to_affine(&[commit(pk, &w_zeta), commit(pk, &w_zeta_omega)]);
    tm.rounds[4] = clock() - t0;

    let points: [G1Affine; 9] = [wire_c[0], wire_c[1], wire_c[2], z_c, t_c[0], t_c[1], t_c[2], w_c[0], w_c[1]];
    (Proof { points, evals: ev }, tm)
}
