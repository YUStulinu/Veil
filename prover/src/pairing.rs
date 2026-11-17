//! The optimal ate pairing of BLS12-381, written for clarity over speed:
//! affine Miller loop on the twist, lines mapped into Fq12 through the
//! untwisting (x', y') -> (x' / w^2, y' / w^3), and a final exponentiation
//! whose hard part is a plain square-and-multiply. The verifier needs two
//! pairings, so this costs milliseconds, not seconds.

use crate::curve::{G1Affine, G2Affine};
use crate::ff::{Field, Fq};
use crate::tower::{Fq12, Fq2, Fq6};

/// |x| for the curve parameter x = -0xd201000000010000
const X: u64 = 0xd201000000010000;

/// (q^4 - q^2 + 1) / r
const HARD: [u64; 20] = [0xe516c3f438e3ba79, 0xfa9912aae208ccf1, 0x905ce937335d5b68, 0xc71a2629b0dea236, 0x83774940996754c8, 0x21d160aeb6a1e799, 0x2ed0b283ed237db4, 0x915c97f36c6f1821, 0x67f17fcbde783765, 0x2378b9039096d1b7, 0x7988f8761bdc51dc, 0x2076995003fc77a1, 0x827eca0ba621315b, 0xe5a72bce8d63cb9f, 0xf68f7764c28b6f8a, 0x2f230063cf081517, 0x94506632528d6a9a, 0xd3cde88eeb996ca3, 0xc0bd38c3195c899e, 0x000f686b3d807d01];

/// the line with slope `l` through `t`, at p, times w^3:
/// (l x_t - y_t) + (-l x_p) w^2 + y_p w^3
fn line(l: Fq2, t: &G2Affine, p: &G1Affine) -> Fq12 {
    Fq12 {
        c0: Fq6 { c0: l * t.x - t.y, c1: l.mul_by_fq(-p.x), c2: Fq2::ZERO },
        c1: Fq6 { c0: Fq2::ZERO, c1: Fq2::new(p.y, Fq::ZERO), c2: Fq2::ZERO },
    }
}

/// the product of the Miller loops of all pairs (identity points are skipped)
pub fn miller_loop(pairs: &[(G1Affine, G2Affine)]) -> Fq12 {
    let pairs: Vec<_> = pairs.iter().filter(|(p, q)| !p.infinity && !q.infinity).collect();
    let mut ts: Vec<G2Affine> = pairs.iter().map(|(_, q)| *q).collect();
    let mut f = Fq12::ONE;
    let top = 63 - X.leading_zeros();
    for bit in (0..top).rev() {
        f = f.square();
        for (k, (p, _)) in pairs.iter().enumerate() {
            let t = ts[k];
            let l = (t.x.square() * Fq2::from_u64(3)) * t.y.double().inverse().unwrap();
            f *= line(l, &t, p);
            let x3 = l.square() - t.x.double();
            ts[k] = G2Affine { x: x3, y: l * (t.x - x3) - t.y, infinity: false };
        }
        if (X >> bit) & 1 == 1 {
            for (k, (p, q)) in pairs.iter().enumerate() {
                let t = ts[k];
                // T = jQ with 1 < j < r here, so T != +-Q
                let l = (q.y - t.y) * (q.x - t.x).inverse().unwrap();
                f *= line(l, &t, p);
                let x3 = l.square() - t.x - q.x;
                ts[k] = G2Affine { x: x3, y: l * (t.x - x3) - t.y, infinity: false };
            }
        }
    }
    // x < 0
    f.conjugate()
}

pub fn final_exponentiation(f: &Fq12) -> Fq12 {
    let f1 = f.conjugate() * f.inverse().expect("Miller loop output is never zero"); // f^(q^6 - 1)
    let f2 = f1.frobenius2() * f1; // ^(q^2 + 1)
    f2.pow(&HARD)
}

pub fn pairing(p: &G1Affine, q: &G2Affine) -> Fq12 {
    final_exponentiation(&miller_loop(&[(*p, *q)]))
}

/// prod e(p_i, q_i) == 1
pub fn pairing_check(pairs: &[(G1Affine, G2Affine)]) -> bool {
    final_exponentiation(&miller_loop(pairs)) == Fq12::ONE
}

impl Fq2 {
    pub fn from_u64(v: u64) -> Fq2 {
        Fq2::new(Fq::from_u64(v), Fq::ZERO)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::{g1_generator, g2_generator};
    use crate::ff::Fr;

    #[test]
    fn bilinear_and_non_degenerate() {
        let g = g1_generator();
        let h = g2_generator();
        let e = pairing(&g, &h);
        assert_ne!(e, Fq12::ONE);
        assert_eq!(e.pow(&Fr::MODULUS), Fq12::ONE);
        let a = Fr::from_u64(0x1234_5678_9abc);
        let b = Fr::from_u64(0xfedc_ba98);
        let ga = g.to_projective().mul(&a).to_affine();
        let hb = h.to_projective().mul(&b).to_affine();
        assert_eq!(pairing(&ga, &hb), e.pow(&(a * b).to_limbs()));
        // e(aG, bH) e(-abG, H) = 1
        let gab = g.to_projective().mul(&(a * b)).neg().to_affine();
        assert!(pairing_check(&[(ga, hb), (gab, h)]));
        assert!(!pairing_check(&[(ga, hb), (g, h)]));
    }
}
