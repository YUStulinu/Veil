//! The proof: 9 G1 points and 6 scalars. Compressed: 9 x 48 + 6 x 32 = 624
//! bytes. For the EVM: 9 x 128 + 6 x 32 = 1344 bytes (EIP-2537 points, the
//! precompiles have no decompression).

use crate::curve::{g1_from_compressed, g1_to_compressed, g1_to_eip2537, G1Affine};
use crate::ff::Fr;

#[derive(Clone, Debug, PartialEq)]
pub struct Proof {
    /// [a], [b], [c], [z], [t_lo], [t_mid], [t_hi], [W_zeta], [W_zeta_omega]
    pub points: [G1Affine; 9],
    /// a(zeta), b(zeta), c(zeta), S_sigma1(zeta), S_sigma2(zeta), z(zeta omega)
    pub evals: [Fr; 6],
}

pub const PROOF_BYTES: usize = 9 * 48 + 6 * 32;
pub const PROOF_EVM_BYTES: usize = 9 * 128 + 6 * 32;

impl Proof {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(PROOF_BYTES);
        for p in &self.points {
            v.extend_from_slice(&g1_to_compressed(p));
        }
        for e in &self.evals {
            v.extend_from_slice(&e.to_bytes_be());
        }
        v
    }

    /// checks every point (on the curve, in the subgroup) and every scalar (< r)
    pub fn from_bytes(b: &[u8]) -> Result<Proof, String> {
        if b.len() != PROOF_BYTES {
            return Err(format!("a proof has {PROOF_BYTES} bytes, got {}", b.len()));
        }
        let mut points = [G1Affine::IDENTITY; 9];
        for (i, p) in points.iter_mut().enumerate() {
            *p = g1_from_compressed(&b[48 * i..48 * i + 48]).map_err(|e| format!("point {i}: {e:?}"))?;
        }
        let mut evals = [Fr::from_u64(0); 6];
        for (i, e) in evals.iter_mut().enumerate() {
            let o = 432 + 32 * i;
            *e = Fr::from_bytes_be(&b[o..o + 32]).ok_or(format!("scalar {i} is not below r"))?;
        }
        Ok(Proof { points, evals })
    }

    pub fn to_evm_bytes(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(PROOF_EVM_BYTES);
        for p in &self.points {
            v.extend_from_slice(&g1_to_eip2537(p));
        }
        for e in &self.evals {
            v.extend_from_slice(&e.to_bytes_be());
        }
        v
    }
}
