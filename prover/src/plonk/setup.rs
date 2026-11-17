//! Preprocessing: the PLONK table from a circuit (public-input rows, gates,
//! padding), the selector and permutation polynomials, their commitments.

use super::circuit::{Circuit, Var};
use crate::curve::{g1_to_eip2537, g2_to_eip2537, G1Affine, G2Affine, Projective};
use crate::ff::{Field, Fr};
use crate::kzg::Srs;
use crate::poly::{powers, Domain};

/// the coset shifts of the second and third wire columns
pub const K1: u64 = 7;
pub const K2: u64 = 13;
/// the shift of the 4n evaluation coset
pub const COSET: u64 = 5;

pub struct ProvingKey {
    pub n: usize,
    pub domain: Domain,
    pub domain4: Domain,
    pub g1: Vec<G1Affine>,
    pub num_public: usize,
    /// variable of each wire, per row: [a, b, c]
    pub wires: [Vec<Var>; 3],
    /// q_M, q_L, q_R, q_O, q_C, S_sigma1, S_sigma2, S_sigma3 in coefficient form
    pub polys: [Vec<Fr>; 8],
    /// the same on the 4n coset
    pub cosets: [Vec<Fr>; 8],
    /// S_sigma_j(omega^i)
    pub sigma_values: [Vec<Fr>; 3],
    /// L_0 on the 4n coset
    pub l0_coset: Vec<Fr>,
    /// the points of the 4n coset
    pub x_coset: Vec<Fr>,
    /// 1 / Z_H on the 4n coset: it takes 4 values
    pub zh_inv_coset: [Fr; 4],
    pub vk: VerifyingKey,
}

#[derive(Clone, Debug, PartialEq)]
pub struct VerifyingKey {
    pub n: usize,
    pub num_public: usize,
    pub k1: Fr,
    pub k2: Fr,
    pub omega: Fr,
    /// [q_M], [q_L], [q_R], [q_O], [q_C], [S_sigma1], [S_sigma2], [S_sigma3]
    pub commitments: [G1Affine; 8],
    /// [1]_2, [tau]_2
    pub g2: [G2Affine; 2],
}

impl VerifyingKey {
    /// "VEILVK1", n, l (u32 big-endian), k1, k2, omega (32 bytes each),
    /// 8 G1 (128 bytes each), 2 G2 (256 bytes each)
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = b"VEILVK1".to_vec();
        v.extend_from_slice(&(self.n as u32).to_be_bytes());
        v.extend_from_slice(&(self.num_public as u32).to_be_bytes());
        for x in [self.k1, self.k2, self.omega] {
            v.extend_from_slice(&x.to_bytes_be());
        }
        for c in &self.commitments {
            v.extend_from_slice(&g1_to_eip2537(c));
        }
        for g in &self.g2 {
            v.extend_from_slice(&g2_to_eip2537(g));
        }
        v
    }

    pub fn from_bytes(b: &[u8]) -> Result<VerifyingKey, String> {
        use crate::ff::Fq;
        use crate::tower::Fq2;
        if b.len() != 7 + 8 + 96 + 8 * 128 + 2 * 256 || &b[..7] != b"VEILVK1" {
            return Err("not a verifying key".into());
        }
        let u32_at = |o: usize| u32::from_be_bytes(b[o..o + 4].try_into().unwrap()) as usize;
        let fr = |o: usize| Fr::from_bytes_be(&b[o..o + 32]).ok_or("scalar not canonical");
        let fq = |o: usize| -> Result<Fq, String> {
            if b[o..o + 16].iter().any(|x| *x != 0) {
                return Err("bad padding".into());
            }
            Fq::from_bytes_be(&b[o + 16..o + 64]).ok_or("coordinate not canonical".into())
        };
        let g1 = |o: usize| -> Result<G1Affine, String> {
            if b[o..o + 128].iter().all(|x| *x == 0) {
                return Ok(G1Affine::IDENTITY);
            }
            let p = G1Affine { x: fq(o)?, y: fq(o + 64)?, infinity: false };
            if p.is_on_curve() && p.is_torsion_free() {
                Ok(p)
            } else {
                Err("invalid G1 point".into())
            }
        };
        let mut commitments = [G1Affine::IDENTITY; 8];
        for (i, c) in commitments.iter_mut().enumerate() {
            *c = g1(7 + 8 + 96 + 128 * i)?;
        }
        let mut g2 = [G2Affine::IDENTITY; 2];
        for (i, g) in g2.iter_mut().enumerate() {
            let o = 7 + 8 + 96 + 1024 + 256 * i;
            *g = G2Affine { x: Fq2::new(fq(o)?, fq(o + 64)?), y: Fq2::new(fq(o + 128)?, fq(o + 192)?), infinity: false };
            if !g.is_on_curve() || !g.is_torsion_free() {
                return Err("invalid G2 point".into());
            }
        }
        let vk = VerifyingKey { n: u32_at(7), num_public: u32_at(11), k1: fr(15)?, k2: fr(47)?, omega: fr(79)?, commitments, g2 };
        if !vk.n.is_power_of_two() || vk.omega != Fr::root_of_unity(vk.n.trailing_zeros()) {
            return Err("inconsistent domain".into());
        }
        Ok(vk)
    }
}

/// the smallest power-of-two domain that holds the circuit's rows
pub fn domain_size(c: &Circuit) -> usize {
    c.rows().next_power_of_two().max(8)
}

pub fn preprocess(c: &Circuit, srs: &Srs) -> Result<ProvingKey, String> {
    let n = domain_size(c);
    if srs.g1.len() < n + 6 {
        return Err(format!("the circuit needs {} SRS powers (n = {n}), the SRS has {}", n + 6, srs.g1.len()));
    }
    let domain = Domain::new(n);
    let domain4 = Domain::new(4 * n);
    let (k1, k2) = (Fr::from_u64(K1), Fr::from_u64(K2));
    let nn = [n as u64];
    // H, k1 H, k2 H are disjoint
    assert!(k1.pow(&nn) != Fr::ONE && k2.pow(&nn) != Fr::ONE && (k2 * k1.inverse().unwrap()).pow(&nn) != Fr::ONE);

    // the table: public-input rows (a - x = 0 with q_L = 1), the gates, zero padding
    let l = c.public.len();
    let mut wires = [vec![0 as Var; n], vec![0; n], vec![0; n]];
    let mut sel = [vec![Fr::ZERO; n], vec![Fr::ZERO; n], vec![Fr::ZERO; n], vec![Fr::ZERO; n], vec![Fr::ZERO; n]];
    for (i, v) in c.public.iter().enumerate() {
        wires[0][i] = *v;
        sel[1][i] = Fr::ONE;
    }
    for (k, g) in c.gates.iter().enumerate() {
        let i = l + k;
        wires[0][i] = g.a;
        wires[1][i] = g.b;
        wires[2][i] = g.c;
        for (s, q) in sel.iter_mut().zip([g.q_m, g.q_l, g.q_r, g.q_o, g.q_c]) {
            s[i] = q;
        }
    }

    // the copy permutation: each variable's positions form one cycle
    let mut positions: Vec<Vec<usize>> = vec![vec![]; c.vars];
    for col in 0..3 {
        for row in 0..n {
            positions[wires[col][row] as usize].push(col * n + row);
        }
    }
    let mut sigma = vec![0usize; 3 * n];
    for ps in &positions {
        for i in 0..ps.len() {
            sigma[ps[i]] = ps[(i + 1) % ps.len()];
        }
    }
    let omegas = domain.elements();
    let shifts = [Fr::ONE, k1, k2];
    let id = |p: usize| shifts[p / n] * omegas[p % n];
    let sigma_values: [Vec<Fr>; 3] = std::array::from_fn(|col| (0..n).map(|row| id(sigma[col * n + row])).collect());

    let mut polys: [Vec<Fr>; 8] = std::array::from_fn(|i| if i < 5 { sel[i].clone() } else { sigma_values[i - 5].clone() });
    for p in polys.iter_mut() {
        domain.ifft(p);
    }
    let g = Fr::from_u64(COSET);
    let cosets: [Vec<Fr>; 8] = std::array::from_fn(|i| {
        let mut e = polys[i].clone();
        e.resize(4 * n, Fr::ZERO);
        domain4.coset_fft(&mut e, g);
        e
    });
    let commit = |p: &[Fr]| crate::msm::msm_g1(&srs.g1[..p.len()], p);
    let commitments_proj: Vec<_> = polys.iter().map(|p| commit(p)).collect();
    let commitments: [G1Affine; 8] = Projective::batch_to_affine(&commitments_proj).try_into().unwrap();

    // L_0 = (X^n - 1) / (n (X - 1)) and Z_H on the coset
    let x_coset: Vec<Fr> = powers(domain4.omega, 4 * n).into_iter().map(|w| w * g).collect();
    let gn = g.pow(&nn);
    let w4 = Fr::root_of_unity(2);
    let zh: [Fr; 4] = std::array::from_fn(|i| gn * w4.pow(&[i as u64]) - Fr::ONE);
    let mut zh_inv_coset = zh;
    Fr::batch_inverse(&mut zh_inv_coset);
    let n_fr = Fr::from_u64(n as u64);
    let mut den: Vec<Fr> = x_coset.iter().map(|x| n_fr * (*x - Fr::ONE)).collect();
    Fr::batch_inverse(&mut den);
    let l0_coset = (0..4 * n).map(|i| zh[i % 4] * den[i]).collect();

    let vk = VerifyingKey { n, num_public: l, k1, k2, omega: domain.omega, commitments, g2: srs.g2 };
    Ok(ProvingKey {
        n,
        domain,
        domain4,
        g1: srs.g1[..n + 6].to_vec(),
        num_public: l,
        wires,
        polys,
        cosets,
        sigma_values,
        l0_coset,
        x_coset,
        zh_inv_coset,
        vk,
    })
}
