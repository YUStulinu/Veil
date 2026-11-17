//! G1 (over Fq) and G2 (over Fq2) of BLS12-381: y^2 = x^3 + b, Jacobian
//! coordinates, and the ZCash serialization (compressed 48 / 96 bytes).

use crate::ff::{Field, Fq, Fr};
use crate::tower::Fq2;

pub trait CurveField: Field {
    const B: Self;
}

impl CurveField for Fq {
    const B: Fq = Fq([0xaa270000000cfff3, 0x53cc0032fc34000a, 0x478fe97a6b0a807f, 0xb1d37ebee6ba24d7, 0x8ec9733bbf78ab2f, 0x09d645513d83de7e]);
}

impl CurveField for Fq2 {
    const B: Fq2 = Fq2 { c0: <Fq as CurveField>::B, c1: <Fq as CurveField>::B };
}

#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct Affine<F> {
    pub x: F,
    pub y: F,
    pub infinity: bool,
}

/// Jacobian: (X, Y, Z) is (X/Z^2, Y/Z^3); Z = 0 is the point at infinity
#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct Projective<F> {
    pub x: F,
    pub y: F,
    pub z: F,
}

pub type G1 = Projective<Fq>;
pub type G2 = Projective<Fq2>;
pub type G1Affine = Affine<Fq>;
pub type G2Affine = Affine<Fq2>;

impl<F: CurveField> PartialEq for Affine<F> {
    fn eq(&self, o: &Self) -> bool {
        (self.infinity && o.infinity) || (!self.infinity && !o.infinity && self.x == o.x && self.y == o.y)
    }
}
impl<F: CurveField> Eq for Affine<F> {}

impl<F: CurveField> Affine<F> {
    pub const IDENTITY: Self = Affine { x: F::ZERO, y: F::ZERO, infinity: true };

    pub fn is_on_curve(&self) -> bool {
        self.infinity || self.y.square() == self.x.square() * self.x + F::B
    }

    pub fn to_projective(&self) -> Projective<F> {
        if self.infinity {
            Projective::IDENTITY
        } else {
            Projective { x: self.x, y: self.y, z: F::ONE }
        }
    }

    pub fn neg(&self) -> Self {
        Affine { x: self.x, y: -self.y, infinity: self.infinity }
    }

    /// in the prime-order subgroup: r P = 0
    pub fn is_torsion_free(&self) -> bool {
        self.to_projective().mul_limbs(&Fr::MODULUS).is_identity()
    }
}

impl<F: CurveField> PartialEq for Projective<F> {
    fn eq(&self, o: &Self) -> bool {
        match (self.is_identity(), o.is_identity()) {
            (true, true) => true,
            (false, false) => {
                let z1 = self.z.square();
                let z2 = o.z.square();
                self.x * z2 == o.x * z1 && self.y * z2 * o.z == o.y * z1 * self.z
            }
            _ => false,
        }
    }
}
impl<F: CurveField> Eq for Projective<F> {}

impl<F: CurveField> Projective<F> {
    pub const IDENTITY: Self = Projective { x: F::ONE, y: F::ONE, z: F::ZERO };

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    pub fn to_affine(&self) -> Affine<F> {
        match self.z.inverse() {
            None => Affine::IDENTITY,
            Some(zi) => {
                let zi2 = zi.square();
                Affine { x: self.x * zi2, y: self.y * zi2 * zi, infinity: false }
            }
        }
    }

    /// dbl-2009-l
    pub fn double(&self) -> Self {
        if self.is_identity() {
            return *self;
        }
        let a = self.x.square();
        let b = self.y.square();
        let c = b.square();
        let d = ((self.x + b).square() - a - c).double();
        let e = a.double() + a;
        let f = e.square();
        let x3 = f - d.double();
        let y3 = e * (d - x3) - c.double().double().double();
        let z3 = (self.y * self.z).double();
        Projective { x: x3, y: y3, z: z3 }
    }

    /// add-2007-bl
    pub fn add(&self, o: &Self) -> Self {
        if self.is_identity() {
            return *o;
        }
        if o.is_identity() {
            return *self;
        }
        let z1z1 = self.z.square();
        let z2z2 = o.z.square();
        let u1 = self.x * z2z2;
        let u2 = o.x * z1z1;
        let s1 = self.y * o.z * z2z2;
        let s2 = o.y * self.z * z1z1;
        let h = u2 - u1;
        let r = (s2 - s1).double();
        if h.is_zero() {
            return if r.is_zero() { self.double() } else { Self::IDENTITY };
        }
        let i = h.double().square();
        let j = h * i;
        let v = u1 * i;
        let x3 = r.square() - j - v.double();
        let y3 = r * (v - x3) - (s1 * j).double();
        let z3 = ((self.z + o.z).square() - z1z1 - z2z2) * h;
        Projective { x: x3, y: y3, z: z3 }
    }

    /// madd-2007-bl: plus an affine point
    pub fn add_affine(&self, o: &Affine<F>) -> Self {
        if o.infinity {
            return *self;
        }
        if self.is_identity() {
            return o.to_projective();
        }
        let z1z1 = self.z.square();
        let u2 = o.x * z1z1;
        let s2 = o.y * self.z * z1z1;
        let h = u2 - self.x;
        let r = (s2 - self.y).double();
        if h.is_zero() {
            return if r.is_zero() { self.double() } else { Self::IDENTITY };
        }
        let hh = h.square();
        let i = hh.double().double();
        let j = h * i;
        let v = self.x * i;
        let x3 = r.square() - j - v.double();
        let y3 = r * (v - x3) - (self.y * j).double();
        let z3 = (self.z + h).square() - z1z1 - hh;
        Projective { x: x3, y: y3, z: z3 }
    }

    pub fn neg(&self) -> Self {
        Projective { x: self.x, y: -self.y, z: self.z }
    }

    /// times an integer given as little-endian limbs
    pub fn mul_limbs(&self, k: &[u64]) -> Self {
        let mut acc = Self::IDENTITY;
        for &limb in k.iter().rev() {
            for i in (0..64).rev() {
                acc = acc.double();
                if (limb >> i) & 1 == 1 {
                    acc = acc.add(self);
                }
            }
        }
        acc
    }

    pub fn mul(&self, k: &Fr) -> Self {
        self.mul_limbs(&k.to_limbs())
    }

    /// many points to affine with one inversion
    pub fn batch_to_affine(v: &[Self]) -> Vec<Affine<F>> {
        let mut acc = F::ONE;
        let mut prefix = Vec::with_capacity(v.len());
        for p in v {
            prefix.push(acc);
            if !p.is_identity() {
                acc *= p.z;
            }
        }
        let mut inv = acc.inverse().unwrap();
        let mut out = vec![Affine::IDENTITY; v.len()];
        for i in (0..v.len()).rev() {
            if v[i].is_identity() {
                continue;
            }
            let zi = inv * prefix[i];
            inv *= v[i].z;
            let zi2 = zi.square();
            out[i] = Affine { x: v[i].x * zi2, y: v[i].y * zi2 * zi, infinity: false };
        }
        out
    }
}

pub fn g1_generator() -> G1Affine {
    Affine {
        x: Fq::from_str_mod("0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb").unwrap(),
        y: Fq::from_str_mod("0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1").unwrap(),
        infinity: false,
    }
}

pub fn g2_generator() -> G2Affine {
    let f = |s| Fq::from_str_mod(s).unwrap();
    Affine {
        x: Fq2::new(
            f("0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8"),
            f("0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e"),
        ),
        y: Fq2::new(
            f("0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801"),
            f("0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be"),
        ),
        infinity: false,
    }
}

// ---- serialization (the ZCash format used by Ethereum and the test vectors) ----

#[derive(Debug, PartialEq)]
pub enum DecodeError {
    Length,
    Flags,
    NotCanonical,
    NotOnCurve,
    NotInSubgroup,
}

pub fn g1_to_compressed(p: &G1Affine) -> [u8; 48] {
    let mut out = [0u8; 48];
    if p.infinity {
        out[0] = 0xc0;
        return out;
    }
    out.copy_from_slice(&p.x.to_bytes_be());
    out[0] |= 0x80;
    if p.y.is_lexicographically_largest() {
        out[0] |= 0x20;
    }
    out
}

pub fn g1_from_compressed(b: &[u8]) -> Result<G1Affine, DecodeError> {
    if b.len() != 48 {
        return Err(DecodeError::Length);
    }
    let (compressed, infinity, sign) = (b[0] & 0x80 != 0, b[0] & 0x40 != 0, b[0] & 0x20 != 0);
    if !compressed {
        return Err(DecodeError::Flags);
    }
    let mut xb = [0u8; 48];
    xb.copy_from_slice(b);
    xb[0] &= 0x1f;
    if infinity {
        return if !sign && xb.iter().all(|&v| v == 0) { Ok(G1Affine::IDENTITY) } else { Err(DecodeError::Flags) };
    }
    let x = Fq::from_bytes_be(&xb).ok_or(DecodeError::NotCanonical)?;
    let mut y = (x.square() * x + Fq::B).sqrt().ok_or(DecodeError::NotOnCurve)?;
    if y.is_lexicographically_largest() != sign {
        y = -y;
    }
    let p = Affine { x, y, infinity: false };
    if !p.is_torsion_free() {
        return Err(DecodeError::NotInSubgroup);
    }
    Ok(p)
}

pub fn g2_to_compressed(p: &G2Affine) -> [u8; 96] {
    let mut out = [0u8; 96];
    if p.infinity {
        out[0] = 0xc0;
        return out;
    }
    out[..48].copy_from_slice(&p.x.c1.to_bytes_be());
    out[48..].copy_from_slice(&p.x.c0.to_bytes_be());
    out[0] |= 0x80;
    if p.y.is_lexicographically_largest() {
        out[0] |= 0x20;
    }
    out
}

pub fn g2_from_compressed(b: &[u8]) -> Result<G2Affine, DecodeError> {
    if b.len() != 96 {
        return Err(DecodeError::Length);
    }
    let (compressed, infinity, sign) = (b[0] & 0x80 != 0, b[0] & 0x40 != 0, b[0] & 0x20 != 0);
    if !compressed {
        return Err(DecodeError::Flags);
    }
    let mut xb = [0u8; 96];
    xb.copy_from_slice(b);
    xb[0] &= 0x1f;
    if infinity {
        return if !sign && xb.iter().all(|&v| v == 0) { Ok(G2Affine::IDENTITY) } else { Err(DecodeError::Flags) };
    }
    let c1 = Fq::from_bytes_be(&xb[..48]).ok_or(DecodeError::NotCanonical)?;
    let c0 = Fq::from_bytes_be(&xb[48..]).ok_or(DecodeError::NotCanonical)?;
    let x = Fq2::new(c0, c1);
    let mut y = (x.square() * x + Fq2::B).sqrt().ok_or(DecodeError::NotOnCurve)?;
    if y.is_lexicographically_largest() != sign {
        y = -y;
    }
    let p = Affine { x, y, infinity: false };
    if !p.is_torsion_free() {
        return Err(DecodeError::NotInSubgroup);
    }
    Ok(p)
}

/// G1 as 128 bytes, each coordinate a 64-byte big-endian word (EIP-2537);
/// the identity is all zeros. The transcript and the Solidity verifier use it.
pub fn g1_to_eip2537(p: &G1Affine) -> [u8; 128] {
    let mut out = [0u8; 128];
    if !p.infinity {
        out[16..64].copy_from_slice(&p.x.to_bytes_be());
        out[80..128].copy_from_slice(&p.y.to_bytes_be());
    }
    out
}

pub fn g2_to_eip2537(p: &G2Affine) -> [u8; 256] {
    let mut out = [0u8; 256];
    if !p.infinity {
        out[16..64].copy_from_slice(&p.x.c0.to_bytes_be());
        out[80..128].copy_from_slice(&p.x.c1.to_bytes_be());
        out[144..192].copy_from_slice(&p.y.c0.to_bytes_be());
        out[208..256].copy_from_slice(&p.y.c1.to_bytes_be());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn b_constant() {
        assert_eq!(Fq::B, Fq::from_u64(4));
        assert_eq!(Fq2::B, Fq2::new(Fq::from_u64(4), Fq::from_u64(4)));
    }

    #[test]
    fn generators_are_on_curve_with_order_r() {
        assert!(g1_generator().is_on_curve() && g1_generator().is_torsion_free());
        assert!(g2_generator().is_on_curve() && g2_generator().is_torsion_free());
    }

    #[test]
    fn group_law() {
        let g = g1_generator().to_projective();
        let a = Fr::from_u64(123456789);
        let b = Fr::from_u64(987654321);
        assert_eq!(g.mul(&a).add(&g.mul(&b)), g.mul(&(a + b)));
        assert_eq!(g.mul(&a).mul(&b), g.mul(&(a * b)));
        assert_eq!(g.add(&g), g.double());
        assert!(g.add(&g.neg()).is_identity());
        assert_eq!(g.double().add_affine(&g1_generator()), g.mul(&Fr::from_u64(3)));
        let h = g2_generator().to_projective();
        assert_eq!(h.mul(&a).add(&h.mul(&b)), h.mul(&(a + b)));
    }

    #[test]
    fn compression_round_trip() {
        let g = g1_generator().to_projective();
        for k in 1..6u64 {
            let p = g.mul(&Fr::from_u64(k * 7919)).to_affine();
            assert_eq!(g1_from_compressed(&g1_to_compressed(&p)).unwrap(), p);
            assert_eq!(g1_from_compressed(&g1_to_compressed(&p.neg())).unwrap(), p.neg());
        }
        let h = g2_generator().to_projective().mul(&Fr::from_u64(31337)).to_affine();
        assert_eq!(g2_from_compressed(&g2_to_compressed(&h)).unwrap(), h);
        assert_eq!(g1_from_compressed(&g1_to_compressed(&G1Affine::IDENTITY)).unwrap(), G1Affine::IDENTITY);
        // the well-known compressed generator
        assert_eq!(&g1_to_compressed(&g1_generator())[..4], &[0x97, 0xf1, 0xd3, 0xa7]);
    }
}
