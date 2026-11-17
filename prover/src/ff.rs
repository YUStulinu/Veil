//! Prime fields in Montgomery form: Fq (the 381-bit base field of BLS12-381)
//! and Fr (the 255-bit scalar field). Little-endian 64-bit limbs.

use std::fmt;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// What the polynomial, curve and pairing code needs from a field.
pub trait Field:
    Copy
    + Clone
    + PartialEq
    + Eq
    + fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
{
    const ZERO: Self;
    const ONE: Self;
    fn is_zero(&self) -> bool {
        *self == Self::ZERO
    }
    fn square(&self) -> Self {
        *self * *self
    }
    fn double(&self) -> Self {
        *self + *self
    }
    fn inverse(&self) -> Option<Self>;
    /// self^e, e as little-endian 64-bit limbs
    fn pow(&self, e: &[u64]) -> Self {
        let mut r = Self::ONE;
        for &limb in e.iter().rev() {
            for i in (0..64).rev() {
                r = r.square();
                if (limb >> i) & 1 == 1 {
                    r *= *self;
                }
            }
        }
        r
    }
}

#[inline(always)]
pub(crate) fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let t = a as u128 + (b as u128) * (c as u128) + carry as u128;
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
pub(crate) fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = a as u128 + b as u128 + carry as u128;
    (t as u64, (t >> 64) as u64)
}

#[inline(always)]
pub(crate) fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let t = (a as u128).wrapping_sub(b as u128 + (borrow >> 63) as u128);
    (t as u64, (t >> 64) as u64)
}

/// a >= b for little-endian limbs
#[inline(always)]
pub(crate) fn geq<const N: usize>(a: &[u64; N], b: &[u64; N]) -> bool {
    for i in (0..N).rev() {
        if a[i] != b[i] {
            return a[i] > b[i];
        }
    }
    true
}

#[inline(always)]
pub(crate) fn sub_limbs<const N: usize>(a: &[u64; N], b: &[u64; N]) -> [u64; N] {
    let mut r = [0u64; N];
    let mut borrow = 0;
    for i in 0..N {
        (r[i], borrow) = sbb(a[i], b[i], borrow);
    }
    r
}

/// Montgomery multiplication (CIOS). Needs p < 2^(64N - 1).
#[inline(always)]
pub(crate) fn mont_mul<const N: usize>(a: &[u64; N], b: &[u64; N], p: &[u64; N], inv: u64) -> [u64; N] {
    let mut t = [0u64; N];
    let mut t_n = 0u64;
    for i in 0..N {
        let mut carry = 0;
        for j in 0..N {
            (t[j], carry) = mac(t[j], a[j], b[i], carry);
        }
        let (s, t_n1) = adc(t_n, carry, 0);
        t_n = s;
        let m = t[0].wrapping_mul(inv);
        let (_, mut carry) = mac(t[0], m, p[0], 0);
        for j in 1..N {
            (t[j - 1], carry) = mac(t[j], m, p[j], carry);
        }
        (t[N - 1], carry) = adc(t_n, carry, 0);
        t_n = t_n1 + carry;
    }
    if t_n != 0 || geq(&t, p) {
        t = sub_limbs(&t, p);
    }
    t
}

macro_rules! prime_field {
    ($name:ident, $n:expr, $bytes:expr, p = $p:expr, r = $r:expr, r2 = $r2:expr, inv = $inv:expr) => {
        #[derive(Copy, Clone, PartialEq, Eq, Hash, Default)]
        #[repr(transparent)]
        pub struct $name(pub(crate) [u64; $n]);

        impl $name {
            pub const MODULUS: [u64; $n] = $p;
            pub const R: [u64; $n] = $r;
            pub const R2: [u64; $n] = $r2;
            pub const INV: u64 = $inv;
            pub const BYTES: usize = $bytes;

            /// from a canonical integer (little-endian limbs); None if >= p
            pub fn from_limbs(l: [u64; $n]) -> Option<Self> {
                if geq(&l, &Self::MODULUS) {
                    return None;
                }
                Some($name(mont_mul(&l, &Self::R2, &Self::MODULUS, Self::INV)))
            }

            pub fn from_u64(v: u64) -> Self {
                let mut l = [0u64; $n];
                l[0] = v;
                Self::from_limbs(l).unwrap()
            }

            pub fn from_i64(v: i64) -> Self {
                if v < 0 {
                    -Self::from_u64(v.unsigned_abs())
                } else {
                    Self::from_u64(v as u64)
                }
            }

            /// the canonical integer
            pub fn to_limbs(&self) -> [u64; $n] {
                let mut one = [0u64; $n];
                one[0] = 1;
                mont_mul(&self.0, &one, &Self::MODULUS, Self::INV)
            }

            /// big-endian bytes of the canonical integer
            pub fn to_bytes_be(&self) -> [u8; $bytes] {
                let l = self.to_limbs();
                let mut out = [0u8; $bytes];
                for i in 0..$bytes {
                    out[$bytes - 1 - i] = (l[i / 8] >> (8 * (i % 8))) as u8;
                }
                out
            }

            /// big-endian bytes; None if the integer is >= p
            pub fn from_bytes_be(b: &[u8]) -> Option<Self> {
                assert_eq!(b.len(), $bytes);
                let mut l = [0u64; $n];
                for i in 0..$bytes {
                    l[i / 8] |= (b[$bytes - 1 - i] as u64) << (8 * (i % 8));
                }
                Self::from_limbs(l)
            }

            /// any number of big-endian bytes, reduced mod p
            pub fn from_bytes_be_mod(b: &[u8]) -> Self {
                let base = Self::from_u64(256);
                let mut acc = Self::ZERO;
                for &x in b {
                    acc = acc * base + Self::from_u64(x as u64);
                }
                acc
            }

            /// decimal or 0x-hex string, reduced mod p; a leading '-' negates
            pub fn from_str_mod(s: &str) -> Option<Self> {
                let (neg, s) = match s.strip_prefix('-') {
                    Some(rest) => (true, rest),
                    None => (false, s),
                };
                let mut acc = Self::ZERO;
                if let Some(h) = s.strip_prefix("0x") {
                    let sixteen = Self::from_u64(16);
                    for c in h.chars() {
                        acc = acc * sixteen + Self::from_u64(c.to_digit(16)? as u64);
                    }
                } else {
                    if s.is_empty() {
                        return None;
                    }
                    let ten = Self::from_u64(10);
                    for c in s.chars() {
                        acc = acc * ten + Self::from_u64(c.to_digit(10)? as u64);
                    }
                }
                Some(if neg { -acc } else { acc })
            }

            /// the "sign" used by point compression: the integer is > (p-1)/2
            pub fn is_lexicographically_largest(&self) -> bool {
                let v = self.to_limbs();
                let neg = (-*self).to_limbs();
                !geq(&neg, &v) // v > p - v
            }
        }

        impl Field for $name {
            const ZERO: Self = $name([0; $n]);
            const ONE: Self = $name($r);
            fn inverse(&self) -> Option<Self> {
                if self.is_zero() {
                    return None;
                }
                let mut e = Self::MODULUS;
                e[0] -= 2;
                Some(self.pow(&e))
            }
            #[inline(always)]
            fn square(&self) -> Self {
                *self * *self
            }
        }

        impl Add for $name {
            type Output = Self;
            #[inline(always)]
            fn add(self, o: Self) -> Self {
                let mut r = [0u64; $n];
                let mut carry = 0;
                for i in 0..$n {
                    (r[i], carry) = adc(self.0[i], o.0[i], carry);
                }
                if geq(&r, &Self::MODULUS) {
                    r = sub_limbs(&r, &Self::MODULUS);
                }
                $name(r)
            }
        }

        impl Sub for $name {
            type Output = Self;
            #[inline(always)]
            fn sub(self, o: Self) -> Self {
                let mut r = [0u64; $n];
                let mut borrow = 0;
                for i in 0..$n {
                    (r[i], borrow) = sbb(self.0[i], o.0[i], borrow);
                }
                if borrow != 0 {
                    let mut carry = 0;
                    for i in 0..$n {
                        (r[i], carry) = adc(r[i], Self::MODULUS[i], carry);
                    }
                }
                $name(r)
            }
        }

        impl Mul for $name {
            type Output = Self;
            #[inline(always)]
            fn mul(self, o: Self) -> Self {
                $name(mont_mul(&self.0, &o.0, &Self::MODULUS, Self::INV))
            }
        }

        impl Neg for $name {
            type Output = Self;
            #[inline(always)]
            fn neg(self) -> Self {
                Self::ZERO - self
            }
        }

        impl AddAssign for $name {
            #[inline(always)]
            fn add_assign(&mut self, o: Self) {
                *self = *self + o;
            }
        }
        impl SubAssign for $name {
            #[inline(always)]
            fn sub_assign(&mut self, o: Self) {
                *self = *self - o;
            }
        }
        impl MulAssign for $name {
            #[inline(always)]
            fn mul_assign(&mut self, o: Self) {
                *self = *self * o;
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "0x")?;
                for b in self.to_bytes_be() {
                    write!(f, "{:02x}", b)?;
                }
                Ok(())
            }
        }
    };
}

prime_field!(Fq, 6, 48,
    p = [0xb9feffffffffaaab, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624, 0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a],
    r = [0x760900000002fffd, 0xebf4000bc40c0002, 0x5f48985753c758ba, 0x77ce585370525745, 0x5c071a97a256ec6d, 0x15f65ec3fa80e493],
    r2 = [0xf4df1f341c341746, 0x0a76e6a609d104f1, 0x8de5476c4c95b6d5, 0x67eb88a9939d83c0, 0x9a793e85b519952d, 0x11988fe592cae3aa],
    inv = 0x89f3fffcfffcfffd);

prime_field!(Fr, 4, 32,
    p = [0xffffffff00000001, 0x53bda402fffe5bfe, 0x3339d80809a1d805, 0x73eda753299d7d48],
    r = [0x00000001fffffffe, 0x5884b7fa00034802, 0x998c4fefecbc4ff5, 0x1824b159acc5056f],
    r2 = [0xc999e990f3f29c6d, 0x2b6cedcb87925c23, 0x05d314967254398f, 0x0748d9d99f59ff11],
    inv = 0xfffffffeffffffff);

impl Fq {
    /// a square root, if one exists (p = 3 mod 4)
    pub fn sqrt(&self) -> Option<Self> {
        const E: [u64; 6] = [0xee7fbfffffffeaab, 0x07aaffffac54ffff, 0xd9cc34a83dac3d89, 0xd91dd2e13ce144af, 0x92c6e9ed90d2eb35, 0x0680447a8e5ff9a6];
        let s = self.pow(&E);
        (s.square() == *self).then_some(s)
    }
}

impl Fr {
    /// 2-adicity of r - 1
    pub const S: u32 = 32;

    /// a primitive 2^k-th root of unity (k <= 32), from 7^((r-1)/2^32)
    pub fn root_of_unity(k: u32) -> Fr {
        assert!(k <= Self::S, "the field has no 2^{k}-th root of unity");
        let mut w = Fr::from_str_mod("0x16a2a19edfe81f20d09b681922c813b4b63683508c2280b93829971f439f0d2b").unwrap();
        for _ in k..Self::S {
            w = w.square();
        }
        w
    }

    /// inverts every nonzero element with one inversion (zeros stay zero)
    pub fn batch_inverse(v: &mut [Fr]) {
        let mut acc = Fr::ONE;
        let mut prefix = Vec::with_capacity(v.len());
        for x in v.iter() {
            prefix.push(acc);
            if !x.is_zero() {
                acc *= *x;
            }
        }
        let mut inv = acc.inverse().unwrap();
        for (x, p) in v.iter_mut().zip(prefix).rev() {
            if !x.is_zero() {
                let t = inv * *x;
                *x = inv * p;
                inv = t;
            }
        }
    }

    /// a uniformly random element from a byte source
    pub fn random(rng: &mut crate::rng::Rng) -> Fr {
        Fr::from_bytes_be_mod(&rng.bytes::<64>())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_agrees_with_integers() {
        // (p-1)^2 = 1, small products, inverses
        let m1 = -Fq::ONE;
        assert_eq!(m1 * m1, Fq::ONE);
        assert_eq!(Fq::from_u64(6) * Fq::from_u64(7), Fq::from_u64(42));
        assert_eq!(Fr::from_u64(3) - Fr::from_u64(5), -Fr::from_u64(2));
        let x = Fq::from_str_mod("123456789012345678901234567890").unwrap();
        assert_eq!(x * x.inverse().unwrap(), Fq::ONE);
        let y = Fr::from_str_mod("0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000000").unwrap();
        assert_eq!(y, -Fr::ONE);
        assert_eq!(Fr::from_str_mod("-5").unwrap(), -Fr::from_u64(5));
    }

    #[test]
    fn byte_round_trip_and_range() {
        let x = Fq::from_str_mod("0xdeadbeef1234").unwrap();
        assert_eq!(Fq::from_bytes_be(&x.to_bytes_be()), Some(x));
        let mut p = [0u8; 48];
        p.copy_from_slice(&(-Fq::ONE).to_bytes_be());
        p[47] += 1; // exactly p: not canonical
        assert_eq!(Fq::from_bytes_be(&p), None);
    }

    #[test]
    fn sqrt_and_roots_of_unity() {
        let x = Fq::from_u64(1234567);
        let s = x.square().sqrt().unwrap();
        assert!(s == x || s == -x);
        assert!(Fq::from_u64(1234567).sqrt().map_or(true, |s| s.square() == Fq::from_u64(1234567)));
        let w = Fr::root_of_unity(32);
        let mut e = [0u64; 4];
        e[0] = 1 << 31;
        assert_ne!(w.pow(&e), Fr::ONE);
        e[0] = 1 << 32;
        assert_eq!(w.pow(&e), Fr::ONE);
    }

    #[test]
    fn batch_inverse_matches() {
        let mut v: Vec<Fr> = (0..10).map(Fr::from_u64).collect();
        Fr::batch_inverse(&mut v);
        assert!(v[0].is_zero());
        for i in 1..10 {
            assert_eq!(v[i] * Fr::from_u64(i as u64), Fr::ONE);
        }
    }
}
