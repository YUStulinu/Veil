//! Fiat-Shamir with Keccak-256. The state starts as
//! keccak256("veil-plonk-v1" || verifying key); each challenge is
//! h = keccak256(state || absorbed bytes), state = h, challenge = h mod r.
//! G1 points are absorbed in their 128-byte EIP-2537 encoding, scalars as
//! 32 big-endian bytes: the same bytes the Solidity verifier hashes.

use crate::curve::{g1_to_eip2537, G1Affine};
use crate::ff::Fr;
use crate::hash::keccak256;

pub struct Transcript {
    state: [u8; 32],
    pending: Vec<u8>,
}

impl Transcript {
    pub fn new(vk_bytes: &[u8]) -> Transcript {
        let mut d = b"veil-plonk-v1".to_vec();
        d.extend_from_slice(vk_bytes);
        Transcript { state: keccak256(&d), pending: vec![] }
    }

    pub fn state(&self) -> [u8; 32] {
        self.state
    }

    pub fn absorb_fr(&mut self, x: &Fr) {
        self.pending.extend_from_slice(&x.to_bytes_be());
    }

    pub fn absorb_g1(&mut self, p: &G1Affine) {
        self.pending.extend_from_slice(&g1_to_eip2537(p));
    }

    pub fn challenge(&mut self) -> Fr {
        let mut d = self.state.to_vec();
        d.append(&mut self.pending);
        self.state = keccak256(&d);
        Fr::from_bytes_be_mod(&self.state)
    }
}
