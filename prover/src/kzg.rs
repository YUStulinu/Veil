//! KZG polynomial commitments: the structured reference string, commitments,
//! single-point openings and their pairing check.

use crate::curve::{g1_from_compressed, g1_generator, g2_from_compressed, g2_generator, G1Affine, G2Affine, Projective, G1};
use crate::ff::{Field, Fq, Fr};
use crate::hash::unhex;
use crate::msm::msm_g1;
use crate::pairing::pairing_check;
use crate::poly::{divide_by_linear, evaluate};
use crate::tower::Fq2;

#[derive(Clone)]
pub struct Srs {
    /// [tau^i]_1
    pub g1: Vec<G1Affine>,
    /// [1]_2 and [tau]_2
    pub g2: [G2Affine; 2],
    /// [L_i(tau)]_1 over the domain of size g1_lagrange.len(), bit-reversed (EIP-4844 order)
    pub g1_lagrange: Vec<G1Affine>,
    /// "ethereum-kzg-ceremony" or "insecure-local"
    pub origin: String,
}

impl Srs {
    /// From a known tau: anyone who knows tau can forge proofs. For benchmarks
    /// beyond the ceremony's 4096 powers and for tests only.
    pub fn insecure_from_tau(tau: Fr, size: usize) -> Srs {
        let g = g1_generator().to_projective();
        let mut pts = Vec::with_capacity(size);
        let mut acc = Fr::ONE;
        // fixed-base: a window table would be faster; this runs once per size
        let threads = std::thread::available_parallelism().map_or(1, |t| t.get());
        let scalars: Vec<Fr> = (0..size)
            .map(|_| {
                let s = acc;
                acc *= tau;
                s
            })
            .collect();
        if cfg!(target_arch = "wasm32") || threads == 1 {
            pts = scalars.iter().map(|s| g.mul(s)).collect();
        } else {
            pts.resize(size, Projective::IDENTITY);
            std::thread::scope(|sc| {
                for (out, ks) in pts.chunks_mut(size.div_ceil(threads)).zip(scalars.chunks(size.div_ceil(threads))) {
                    sc.spawn(move || {
                        for (o, k) in out.iter_mut().zip(ks) {
                            *o = g.mul(k);
                        }
                    });
                }
            });
        }
        let h = g2_generator();
        Srs {
            g1: Projective::batch_to_affine(&pts),
            g2: [h, h.to_projective().mul(&tau).to_affine()],
            g1_lagrange: vec![],
            origin: "insecure-local".into(),
        }
    }

    /// The Ethereum KZG ceremony output in c-kzg-4844's trusted_setup.txt:
    /// counts, the Lagrange G1 points, 65 G2 powers, the monomial G1 points.
    /// Points are checked to be on the curve; the subgroup check is skipped
    /// for this trusted file (it costs seconds).
    pub fn from_ceremony_file(text: &str) -> Result<Srs, String> {
        let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
        let n1: usize = lines.next().ok_or("empty")?.parse().map_err(|_| "count")?;
        let n2: usize = lines.next().ok_or("short")?.parse().map_err(|_| "count")?;
        let mut lag = Vec::with_capacity(n1);
        for _ in 0..n1 {
            lag.push(g1_from_compressed_trusted(&unhex(lines.next().ok_or("short")?).ok_or("hex")?)?);
        }
        let mut g2 = Vec::with_capacity(n2);
        for _ in 0..n2 {
            g2.push(g2_from_compressed(&unhex(lines.next().ok_or("short")?).ok_or("hex")?).map_err(|e| format!("{e:?}"))?);
        }
        let mut mono = Vec::with_capacity(n1);
        for l in lines.take(n1) {
            mono.push(g1_from_compressed_trusted(&unhex(l).ok_or("hex")?)?);
        }
        if mono.len() != n1 || g2.len() < 2 {
            return Err("this trusted setup file has no monomial G1 points".into());
        }
        // the file lists them in natural order; EIP-4844 blobs are bit-reversed
        let bits = n1.trailing_zeros();
        let lag = (0..n1).map(|i| lag[i.reverse_bits() >> (usize::BITS - bits)]).collect();
        Ok(Srs { g1: mono, g2: [g2[0], g2[1]], g1_lagrange: lag, origin: "ethereum-kzg-ceremony".into() })
    }

    pub fn max_degree(&self) -> usize {
        self.g1.len() - 1
    }

    pub fn commit(&self, coeffs: &[Fr]) -> G1 {
        assert!(coeffs.len() <= self.g1.len(), "polynomial of degree {} exceeds the SRS ({} powers)", coeffs.len() - 1, self.g1.len());
        msm_g1(&self.g1[..coeffs.len()], coeffs)
    }

    /// (p(z), [q(tau)]_1) with q = (p - p(z)) / (X - z)
    pub fn open(&self, coeffs: &[Fr], z: Fr) -> (Fr, G1) {
        let y = evaluate(coeffs, z);
        let mut shifted = coeffs.to_vec();
        if !shifted.is_empty() {
            shifted[0] -= y;
        }
        (y, self.commit(&divide_by_linear(&shifted, z)))
    }

    /// e(C - [y]_1, [1]_2) = e(pi, [tau - z]_2)
    pub fn verify(&self, c: &G1Affine, z: Fr, y: Fr, pi: &G1Affine) -> bool {
        let g = g1_generator().to_projective();
        let lhs = c.to_projective().add(&g.mul(&y).neg()).to_affine();
        let rhs_g2 = self.g2[1].to_projective().add(&self.g2[0].to_projective().mul(&z).neg()).to_affine();
        pairing_check(&[(lhs, self.g2[0].neg()), (*pi, rhs_g2)])
    }

    /// uncompressed little file for fast loading (the browser uses it):
    /// "VSRS" u32 count, G1 points as x||y (96 bytes), then the two G2 points
    pub fn to_bytes(&self, count: usize) -> Vec<u8> {
        let mut out = b"VSRS".to_vec();
        out.extend_from_slice(&(count as u32).to_le_bytes());
        for p in &self.g1[..count] {
            out.extend_from_slice(&p.x.to_bytes_be());
            out.extend_from_slice(&p.y.to_bytes_be());
        }
        for p in &self.g2 {
            for c in [p.x.c0, p.x.c1, p.y.c0, p.y.c1] {
                out.extend_from_slice(&c.to_bytes_be());
            }
        }
        out
    }

    pub fn from_bytes(b: &[u8], origin: &str) -> Result<Srs, String> {
        if b.len() < 8 || &b[..4] != b"VSRS" {
            return Err("not an SRS file".into());
        }
        let n = u32::from_le_bytes(b[4..8].try_into().unwrap()) as usize;
        if b.len() != 8 + 96 * n + 2 * 192 {
            return Err("SRS file has the wrong length".into());
        }
        let fq = |o: usize| Fq::from_bytes_be(&b[o..o + 48]).ok_or("coordinate not canonical");
        let mut g1 = Vec::with_capacity(n);
        for i in 0..n {
            let p = G1Affine { x: fq(8 + 96 * i)?, y: fq(8 + 96 * i + 48)?, infinity: false };
            if !p.is_on_curve() {
                return Err("SRS point not on the curve".into());
            }
            g1.push(p);
        }
        let mut g2 = [G2Affine::IDENTITY; 2];
        for (k, p) in g2.iter_mut().enumerate() {
            let o = 8 + 96 * n + 192 * k;
            *p = G2Affine { x: Fq2::new(fq(o)?, fq(o + 48)?), y: Fq2::new(fq(o + 96)?, fq(o + 144)?), infinity: false };
            if !p.is_on_curve() {
                return Err("SRS point not on the curve".into());
            }
        }
        Ok(Srs { g1, g2, g1_lagrange: vec![], origin: origin.into() })
    }
}

fn g1_from_compressed_trusted(b: &[u8]) -> Result<G1Affine, String> {
    // same as g1_from_compressed without the r P = 0 check
    if b.len() != 48 || b[0] & 0x80 == 0 || b[0] & 0x40 != 0 {
        return g1_from_compressed(b).map_err(|e| format!("{e:?}"));
    }
    let sign = b[0] & 0x20 != 0;
    let mut xb = [0u8; 48];
    xb.copy_from_slice(b);
    xb[0] &= 0x1f;
    let x = Fq::from_bytes_be(&xb).ok_or("x not canonical")?;
    let mut y = (x.square() * x + Fq::from_u64(4)).sqrt().ok_or("not on curve")?;
    if y.is_lexicographically_largest() != sign {
        y = -y;
    }
    Ok(G1Affine { x, y, infinity: false })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn open_and_verify() {
        let mut rng = Rng::from_seed([3; 32]);
        let srs = Srs::insecure_from_tau(Fr::random(&mut rng), 64);
        let p: Vec<Fr> = (0..50).map(|_| Fr::random(&mut rng)).collect();
        let c = srs.commit(&p).to_affine();
        let z = Fr::random(&mut rng);
        let (y, pi) = srs.open(&p, z);
        let pi = pi.to_affine();
        assert!(srs.verify(&c, z, y, &pi));
        assert!(!srs.verify(&c, z, y + Fr::ONE, &pi));
        assert!(!srs.verify(&c, z + Fr::ONE, y, &pi));
        let srs2 = Srs::from_bytes(&srs.to_bytes(64), "t").unwrap();
        assert!(srs2.verify(&c, z, y, &pi));
    }
}
