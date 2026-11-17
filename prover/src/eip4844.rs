//! The EIP-4844 KZG functions (polynomial-commitments.md of the consensus
//! specs), so the official test vectors can check the KZG layer and the
//! pairing against the Ethereum ceremony's setup.

use crate::curve::{g1_from_compressed, g1_to_compressed, G1Affine};
use crate::ff::{Field, Fr};
use crate::hash::sha256;
use crate::kzg::Srs;
use crate::msm::msm_g1;
use crate::poly::{pow_u64, Domain};

pub const N: usize = 4096;

#[derive(Debug, PartialEq)]
pub struct InvalidInput;

/// the domain in bit-reversed order
pub fn roots_brp() -> Vec<Fr> {
    let d = Domain::new(N);
    let e = d.elements();
    (0..N).map(|i| e[i.reverse_bits() >> (usize::BITS - 12)]).collect()
}

pub fn blob_to_poly(blob: &[u8]) -> Result<Vec<Fr>, InvalidInput> {
    if blob.len() != 32 * N {
        return Err(InvalidInput);
    }
    blob.chunks(32).map(|c| Fr::from_bytes_be(c).ok_or(InvalidInput)).collect()
}

fn field(b: &[u8]) -> Result<Fr, InvalidInput> {
    if b.len() != 32 {
        return Err(InvalidInput);
    }
    Fr::from_bytes_be(b).ok_or(InvalidInput)
}

fn point(b: &[u8]) -> Result<G1Affine, InvalidInput> {
    g1_from_compressed(b).map_err(|_| InvalidInput)
}

pub fn evaluate_in_evaluation_form(poly: &[Fr], z: Fr, roots: &[Fr]) -> Fr {
    if let Some(i) = roots.iter().position(|w| *w == z) {
        return poly[i];
    }
    let mut den: Vec<Fr> = roots.iter().map(|w| z - *w).collect();
    Fr::batch_inverse(&mut den);
    let mut acc = Fr::ZERO;
    for i in 0..N {
        acc += poly[i] * roots[i] * den[i];
    }
    acc * (pow_u64(z, N as u64) - Fr::ONE) * Fr::from_u64(N as u64).inverse().unwrap()
}

pub fn blob_to_kzg_commitment(srs: &Srs, blob: &[u8]) -> Result<[u8; 48], InvalidInput> {
    let p = blob_to_poly(blob)?;
    Ok(g1_to_compressed(&msm_g1(&srs.g1_lagrange, &p).to_affine()))
}

fn compute_kzg_proof_impl(srs: &Srs, poly: &[Fr], z: Fr) -> (G1Affine, Fr) {
    let roots = roots_brp();
    let y = evaluate_in_evaluation_form(poly, z, &roots);
    let mut den: Vec<Fr> = roots.iter().map(|w| *w - z).collect();
    Fr::batch_inverse(&mut den);
    let mut q = vec![Fr::ZERO; N];
    for i in 0..N {
        if roots[i] == z {
            // the quotient at z itself, from the other points
            let mut acc = Fr::ZERO;
            for j in 0..N {
                if j != i {
                    acc += (poly[j] - y) * roots[j] * (z * (z - roots[j])).inverse().unwrap();
                }
            }
            q[i] = acc;
        } else {
            q[i] = (poly[i] - y) * den[i];
        }
    }
    (msm_g1(&srs.g1_lagrange, &q).to_affine(), y)
}

pub fn compute_kzg_proof(srs: &Srs, blob: &[u8], z: &[u8]) -> Result<([u8; 48], [u8; 32]), InvalidInput> {
    let p = blob_to_poly(blob)?;
    let (pi, y) = compute_kzg_proof_impl(srs, &p, field(z)?);
    Ok((g1_to_compressed(&pi), y.to_bytes_be()))
}

pub fn verify_kzg_proof(srs: &Srs, c: &[u8], z: &[u8], y: &[u8], pi: &[u8]) -> Result<bool, InvalidInput> {
    let (c, z, y, pi) = (point(c)?, field(z)?, field(y)?, point(pi)?);
    Ok(srs.verify(&c, z, y, &pi))
}

pub fn compute_challenge(blob: &[u8], commitment: &[u8]) -> Fr {
    let mut data = b"FSBLOBVERIFY_V1_".to_vec();
    data.extend_from_slice(&(N as u128).to_be_bytes());
    data.extend_from_slice(blob);
    data.extend_from_slice(commitment);
    Fr::from_bytes_be_mod(&sha256(&data))
}

pub fn compute_blob_kzg_proof(srs: &Srs, blob: &[u8], commitment: &[u8]) -> Result<[u8; 48], InvalidInput> {
    let p = blob_to_poly(blob)?;
    point(commitment)?;
    let z = compute_challenge(blob, commitment);
    Ok(g1_to_compressed(&compute_kzg_proof_impl(srs, &p, z).0))
}

pub fn verify_blob_kzg_proof(srs: &Srs, blob: &[u8], commitment: &[u8], proof: &[u8]) -> Result<bool, InvalidInput> {
    let p = blob_to_poly(blob)?;
    let c = point(commitment)?;
    let pi = point(proof)?;
    let z = compute_challenge(blob, commitment);
    let y = evaluate_in_evaluation_form(&p, z, &roots_brp());
    Ok(srs.verify(&c, z, y, &pi))
}
