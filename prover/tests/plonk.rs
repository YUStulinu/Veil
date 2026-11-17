//! PLONK: completeness, and soundness tests: tampered proofs, wrong public
//! inputs, witnesses that break a gate or a copy constraint, invalid points,
//! a proof replayed against another circuit, and random bit flips.

use veil::curve::g1_generator;
use veil::ff::{Field, Fr};
use veil::kzg::Srs;
use veil::plonk::circuit::{synthetic, Circuit, Hint};
use veil::plonk::{preprocess, prove, verify, Proof, ProvingKey};
use veil::rng::Rng;

fn srs() -> Srs {
    Srs::insecure_from_tau(Fr::from_u64(0xdead_beef_1234), 300)
}

/// knows x, y with x y = p and x + y = s (p, s public)
fn factor_circuit() -> Circuit {
    let mut c = Circuit::new("factor");
    let pubs = c.input("ps", true, 2);
    let xy = c.input("xy", false, 2);
    let prod = c.mul(xy[0], xy[1]);
    let sum = c.add(xy[0], xy[1]);
    c.assert_equal(prod, pubs[0]);
    c.assert_equal(sum, pubs[1]);
    c
}

fn setup(c: &Circuit) -> ProvingKey {
    preprocess(c, &srs()).unwrap()
}

fn fr(v: &[u64]) -> Vec<Fr> {
    v.iter().map(|x| Fr::from_u64(*x)).collect()
}

#[test]
fn honest_proofs_verify_and_are_624_bytes() {
    let c = factor_circuit();
    let pk = setup(&c);
    let mut rng = Rng::from_seed([1; 32]);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut rng);
    assert_eq!(verify(&pk.vk, &public, &proof), Ok(()));
    let bytes = proof.to_bytes();
    assert_eq!(bytes.len(), 624);
    assert_eq!(Proof::from_bytes(&bytes).unwrap(), proof);
    // zero knowledge needs fresh blinding: the same statement gives a different proof
    let (proof2, _) = prove(&pk, &public, &w, &mut rng);
    assert_ne!(proof2.to_bytes(), bytes);
    assert_eq!(verify(&pk.vk, &public, &proof2), Ok(()));
    // the verifying key survives serialization
    let vk2 = veil::plonk::VerifyingKey::from_bytes(&pk.vk.to_bytes()).unwrap();
    assert_eq!(vk2, pk.vk);
}

#[test]
fn larger_circuits_verify() {
    for gates in [30usize, 120, 250] {
        let (c, inputs) = synthetic(gates);
        let pk = setup(&c);
        let w = c.witness(&inputs).unwrap();
        let (proof, _) = prove(&pk, &c.public_values(&w), &w, &mut Rng::from_seed([gates as u8; 32]));
        assert_eq!(verify(&pk.vk, &c.public_values(&w), &proof), Ok(()), "{gates} gates");
    }
}

#[test]
fn every_tampered_element_is_rejected() {
    let c = factor_circuit();
    let pk = setup(&c);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_seed([2; 32]));
    let g = g1_generator();
    for i in 0..9 {
        let mut p = proof.clone();
        p.points[i] = p.points[i].to_projective().add_affine(&g).to_affine();
        assert!(verify(&pk.vk, &public, &p).is_err(), "point {i}");
        let mut p = proof.clone();
        p.points[i] = proof.points[(i + 1) % 9]; // swapped with its neighbour
        assert!(verify(&pk.vk, &public, &p).is_err(), "point {i} swapped");
    }
    for i in 0..6 {
        let mut p = proof.clone();
        p.evals[i] += Fr::ONE;
        assert!(verify(&pk.vk, &public, &p).is_err(), "evaluation {i}");
    }
}

#[test]
fn wrong_public_inputs_are_rejected() {
    let c = factor_circuit();
    let pk = setup(&c);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_seed([3; 32]));
    assert!(verify(&pk.vk, &fr(&[78, 18]), &proof).is_err());
    assert!(verify(&pk.vk, &fr(&[77, 19]), &proof).is_err());
    assert!(verify(&pk.vk, &fr(&[77]), &proof).is_err());
    // claiming other public inputs than the witness has: the PI rows break
    let (proof, _) = prove(&pk, &fr(&[78, 18]), &w, &mut Rng::from_seed([4; 32]));
    assert!(verify(&pk.vk, &fr(&[78, 18]), &proof).is_err());
}

#[test]
fn a_witness_that_breaks_a_gate_is_rejected() {
    let c = factor_circuit();
    let pk = setup(&c);
    // 6 * 13 = 78, not 77: the multiplication gate fails, nothing else does
    let w = c.witness_unchecked(&fr(&[77, 19, 6, 13])).unwrap();
    assert!(c.first_failing_gate(&w).is_some());
    for seed in 0..5 {
        let (proof, _) = prove(&pk, &fr(&[77, 19]), &w, &mut Rng::from_seed([seed; 32]));
        assert!(verify(&pk.vk, &fr(&[77, 19]), &proof).is_err());
    }
}

#[test]
fn a_witness_that_breaks_a_copy_constraint_is_rejected() {
    // prod = x y; out = prod + 1 (public). The cheat feeds the second gate
    // another value of prod: each gate holds on its own, the wiring does not.
    let mut c = Circuit::new("copy");
    let out = c.input("out", true, 1)[0];
    let xy = c.input("xy", false, 2);
    let prod = c.mul(xy[0], xy[1]);
    let plus1 = c.eval(Fr::ZERO, Fr::ONE, Fr::ZERO, Fr::ONE, prod, 0);
    c.assert_equal(plus1, out);
    let mut pk = setup(&c);
    let mut w = c.witness(&fr(&[43, 6, 7])).unwrap();
    assert!(c.first_failing_gate(&w).is_none());
    // a fresh variable, prod' = 99, used in the "+1" gate; its output becomes 100
    let fake = w.len() as u32;
    w.push(Fr::from_u64(99));
    let row = 1 + 1; // after the public row and the multiplication
    assert_eq!(pk.wires[0][row], prod);
    pk.wires[0][row] = fake;
    w[plus1 as usize] = Fr::from_u64(100);
    let (proof, _) = prove(&pk, &fr(&[100]), &w, &mut Rng::from_seed([5; 32]));
    assert!(verify(&pk.vk, &fr(&[100]), &proof).is_err());
}

#[test]
fn invalid_points_are_rejected() {
    let c = factor_circuit();
    let pk = setup(&c);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_seed([6; 32]));
    // a point on the curve outside the order-r subgroup
    let mut x = veil::ff::Fq::from_u64(1);
    let bad = loop {
        if let Some(y) = (x.square() * x + veil::ff::Fq::from_u64(4)).sqrt() {
            let p = veil::curve::G1Affine { x, y, infinity: false };
            if !p.is_torsion_free() {
                break p;
            }
        }
        x += veil::ff::Fq::ONE;
    };
    let mut p = proof.clone();
    p.points[7] = bad;
    assert_eq!(verify(&pk.vk, &public, &p), Err(veil::plonk::verifier::VerifyError::BadPoint));
    let mut bytes = proof.to_bytes();
    bytes[7 * 48..8 * 48].copy_from_slice(&veil::curve::g1_to_compressed(&bad));
    assert!(Proof::from_bytes(&bytes).is_err());
    // a scalar >= r
    let mut bytes = proof.to_bytes();
    bytes[432..464].fill(0xff);
    assert!(Proof::from_bytes(&bytes).is_err());
}

#[test]
fn a_proof_for_one_circuit_fails_for_another() {
    let c = factor_circuit();
    let pk = setup(&c);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_seed([7; 32]));
    // the same shape with the sum gate turned into a product
    let mut c2 = c.clone();
    c2.gates[1].q_m = Fr::ONE;
    c2.gates[1].q_l = Fr::ZERO;
    c2.gates[1].q_r = Fr::ZERO;
    let pk2 = setup(&c2);
    assert!(verify(&pk2.vk, &public, &proof).is_err());
}

#[test]
fn random_bit_flips_never_verify() {
    let c = factor_circuit();
    let pk = setup(&c);
    let w = c.witness(&fr(&[77, 18, 7, 11])).unwrap();
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_seed([8; 32]));
    let bytes = proof.to_bytes();
    let mut rng = Rng::from_seed([9; 32]);
    let (mut parse_fail, mut verify_fail) = (0, 0);
    for _ in 0..200 {
        let mut b = bytes.clone();
        let bit = (rng.u64() % (8 * b.len() as u64)) as usize;
        b[bit / 8] ^= 1 << (bit % 8);
        match Proof::from_bytes(&b) {
            Err(_) => parse_fail += 1,
            Ok(p) => {
                assert!(verify(&pk.vk, &public, &p).is_err(), "bit {bit} flipped and still verified");
                verify_fail += 1;
            }
        }
    }
    eprintln!("bit flips: {parse_fail} rejected when parsing, {verify_fail} by the verifier");
}

#[test]
fn identity_commitment_is_handled() {
    // a circuit whose selectors are zero in some column gives identity commitments
    let mut c = Circuit::new("id");
    let x = c.input("x", true, 1)[0];
    let inv = c.new_var(Hint::Inv(x));
    c.assert_gate(Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ZERO, -Fr::ONE, x, inv, 0);
    let pk = setup(&c);
    assert!(pk.vk.commitments[3].infinity, "q_O is zero");
    let w = c.witness(&fr(&[5])).unwrap();
    let (proof, _) = prove(&pk, &fr(&[5]), &w, &mut Rng::from_seed([10; 32]));
    assert_eq!(verify(&pk.vk, &fr(&[5]), &proof), Ok(()));
}
