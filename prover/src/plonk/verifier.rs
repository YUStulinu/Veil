//! The PLONK verifier: recompute the challenges, build [F] - [E] from the
//! commitments, and one pairing check (two Miller loops).

use super::proof::Proof;
use super::setup::VerifyingKey;
use super::transcript::Transcript;
use crate::curve::{g1_generator, G1Affine, G1};
use crate::ff::{Field, Fr};
use crate::msm::msm_cpu;
use crate::pairing::pairing_check;
use crate::poly::pow_u64;

#[derive(Debug, PartialEq)]
pub enum VerifyError {
    PublicInputs,
    /// a point of the proof is not in G1, or the identity where it may not be
    BadPoint,
    Pairing,
}

/// The challenges, also exported for the Solidity and TypeScript tests.
#[derive(Debug, Clone, PartialEq)]
pub struct Challenges {
    pub beta: Fr,
    pub gamma: Fr,
    pub alpha: Fr,
    pub zeta: Fr,
    pub v: Fr,
    pub u: Fr,
}

pub fn challenges(vk: &VerifyingKey, public: &[Fr], proof: &Proof) -> Challenges {
    let mut tr = Transcript::new(&vk.to_bytes());
    for x in public {
        tr.absorb_fr(x);
    }
    let p = &proof.points;
    for c in &p[0..3] {
        tr.absorb_g1(c);
    }
    let beta = tr.challenge();
    let gamma = tr.challenge();
    tr.absorb_g1(&p[3]);
    let alpha = tr.challenge();
    for c in &p[4..7] {
        tr.absorb_g1(c);
    }
    let zeta = tr.challenge();
    for e in &proof.evals {
        tr.absorb_fr(e);
    }
    let v = tr.challenge();
    tr.absorb_g1(&p[7]);
    tr.absorb_g1(&p[8]);
    let u = tr.challenge();
    Challenges { beta, gamma, alpha, zeta, v, u }
}

/// The proof's points must already be valid G1 elements (`Proof::from_bytes`
/// checks them); scalars are field elements by construction.
pub fn verify(vk: &VerifyingKey, public: &[Fr], proof: &Proof) -> Result<(), VerifyError> {
    if public.len() != vk.num_public {
        return Err(VerifyError::PublicInputs);
    }
    if proof.points.iter().any(|p| !p.is_on_curve() || !p.is_torsion_free()) {
        return Err(VerifyError::BadPoint);
    }
    let Challenges { beta, gamma, alpha, zeta, v, u } = challenges(vk, public, proof);
    let n = vk.n as u64;
    let [ea, eb, ec, es1, es2, ezw] = proof.evals;

    // Z_H(zeta), L_0(zeta) and PI(zeta) = -sum x_i L_i(zeta), with one inversion
    let zn = pow_u64(zeta, n);
    let zh = zn - Fr::ONE;
    let n_inv = Fr::from_u64(n).inverse().unwrap();
    let l = public.len().max(1);
    let mut w_i = Fr::ONE;
    let mut dens = Vec::with_capacity(l);
    let mut ws = Vec::with_capacity(l);
    for _ in 0..l {
        ws.push(w_i);
        dens.push(zeta - w_i);
        w_i *= vk.omega;
    }
    if dens.iter().any(|d| d.is_zero()) {
        return Err(VerifyError::Pairing); // zeta in the domain: negligible, reject
    }
    Fr::batch_inverse(&mut dens);
    let lag = |i: usize| ws[i] * zh * n_inv * dens[i];
    let l0 = lag(0);
    let pi = public.iter().enumerate().fold(Fr::ZERO, |acc, (i, x)| acc - *x * lag(i));

    let alpha2 = alpha.square();
    let perm_s = alpha * (ea + beta * es1 + gamma) * (eb + beta * es2 + gamma) * ezw;
    let r0 = pi - perm_s * (ec + gamma) - alpha2 * l0;
    let perm_z = alpha * (ea + beta * zeta + gamma) * (eb + beta * vk.k1 * zeta + gamma) * (ec + beta * vk.k2 * zeta + gamma) + alpha2 * l0;
    let [a_c, b_c, c_c, z_c, tlo, tmid, thi, w1, w2] = proof.points;
    let [qm, ql, qr, qo, qc, s1, s2, s3] = vk.commitments;
    let (v2, v3, v4, v5) = (v * v, v * v * v, v * v * v * v, v * v * v * v * v);

    // [F] - [E] + zeta [W] + u zeta omega [W'] in one MSM
    let e = -r0 + v * ea + v2 * eb + v3 * ec + v4 * es1 + v5 * es2 + u * ezw;
    let pts: [G1Affine; 17] = [qm, ql, qr, qo, qc, z_c, s3, tlo, tmid, thi, a_c, b_c, c_c, s1, s2, g1_generator(), w1];
    let ks: [Fr; 17] = [
        ea * eb,
        ea,
        eb,
        ec,
        Fr::ONE,
        perm_z + u,
        -perm_s * beta,
        -zh,
        -zh * zn,
        -zh * zn * zn,
        v,
        v2,
        v3,
        v4,
        v5,
        -e,
        zeta,
    ];
    let rhs: G1 = msm_cpu(&pts, &ks).add(&w2.to_projective().mul(&(u * zeta * vk.omega)));
    let lhs: G1 = w1.to_projective().add(&w2.to_projective().mul(&u));
    // e(W + u W', [tau]_2) = e(rhs, [1]_2)
    if pairing_check(&[(lhs.to_affine(), vk.g2[1]), (rhs.to_affine().neg(), vk.g2[0])]) {
        Ok(())
    } else {
        Err(VerifyError::Pairing)
    }
}
