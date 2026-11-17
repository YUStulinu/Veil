//! The extension tower of BLS12-381:
//! Fq2 = Fq[u]/(u^2 + 1), Fq6 = Fq2[v]/(v^3 - xi) with xi = u + 1,
//! Fq12 = Fq6[w]/(w^2 - v).

use crate::ff::{Field, Fq};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::sync::OnceLock;

macro_rules! ops {
    ($t:ident, $($f:ident),+) => {
        impl Add for $t { type Output = Self; #[inline] fn add(self, o: Self) -> Self { $t { $($f: self.$f + o.$f),+ } } }
        impl Sub for $t { type Output = Self; #[inline] fn sub(self, o: Self) -> Self { $t { $($f: self.$f - o.$f),+ } } }
        impl Neg for $t { type Output = Self; #[inline] fn neg(self) -> Self { $t { $($f: -self.$f),+ } } }
        impl AddAssign for $t { #[inline] fn add_assign(&mut self, o: Self) { *self = *self + o; } }
        impl SubAssign for $t { #[inline] fn sub_assign(&mut self, o: Self) { *self = *self - o; } }
        impl MulAssign for $t { #[inline] fn mul_assign(&mut self, o: Self) { *self = *self * o; } }
    };
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}
ops!(Fq2, c0, c1);

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Fq2 {
        Fq2 { c0, c1 }
    }
    pub fn conjugate(&self) -> Fq2 {
        Fq2 { c0: self.c0, c1: -self.c1 }
    }
    /// times xi = u + 1
    #[inline]
    pub fn mul_by_nonresidue(&self) -> Fq2 {
        Fq2 { c0: self.c0 - self.c1, c1: self.c0 + self.c1 }
    }
    pub fn mul_by_fq(&self, k: Fq) -> Fq2 {
        Fq2 { c0: self.c0 * k, c1: self.c1 * k }
    }
    /// the compression "sign": c1 decides unless it is zero
    pub fn is_lexicographically_largest(&self) -> bool {
        if self.c1.is_zero() {
            self.c0.is_lexicographically_largest()
        } else {
            self.c1.is_lexicographically_largest()
        }
    }
    /// a square root (Adj & Rodriguez-Henriquez, algorithm 9, q = 3 mod 4)
    pub fn sqrt(&self) -> Option<Fq2> {
        if self.is_zero() {
            return Some(*self);
        }
        const P_MINUS_3_DIV_4: [u64; 6] = [0xee7fbfffffffeaaa, 0x07aaffffac54ffff, 0xd9cc34a83dac3d89, 0xd91dd2e13ce144af, 0x92c6e9ed90d2eb35, 0x0680447a8e5ff9a6];
        const P_MINUS_1_DIV_2: [u64; 6] = [0xdcff7fffffffd555, 0x0f55ffff58a9ffff, 0xb39869507b587b12, 0xb23ba5c279c2895f, 0x258dd3db21a5d66b, 0x0d0088f51cbff34d];
        let a1 = self.pow(&P_MINUS_3_DIV_4);
        let alpha = a1.square() * *self;
        let a0 = alpha.conjugate() * alpha;
        if a0 == -Fq2::ONE {
            return None;
        }
        let x0 = a1 * *self;
        let x = if alpha == -Fq2::ONE {
            Fq2 { c0: -x0.c1, c1: x0.c0 } // u * x0
        } else {
            (alpha + Fq2::ONE).pow(&P_MINUS_1_DIV_2) * x0
        };
        (x.square() == *self).then_some(x)
    }
}

impl Mul for Fq2 {
    type Output = Self;
    #[inline]
    fn mul(self, o: Self) -> Self {
        let aa = self.c0 * o.c0;
        let bb = self.c1 * o.c1;
        let s = (self.c0 + self.c1) * (o.c0 + o.c1);
        Fq2 { c0: aa - bb, c1: s - aa - bb }
    }
}

impl Field for Fq2 {
    const ZERO: Self = Fq2 { c0: Fq::ZERO, c1: Fq::ZERO };
    const ONE: Self = Fq2 { c0: Fq::ONE, c1: Fq::ZERO };
    fn square(&self) -> Self {
        let a = (self.c0 + self.c1) * (self.c0 - self.c1);
        let b = self.c0 * self.c1;
        Fq2 { c0: a, c1: b.double() }
    }
    fn inverse(&self) -> Option<Self> {
        let t = (self.c0.square() + self.c1.square()).inverse()?;
        Some(Fq2 { c0: self.c0 * t, c1: -(self.c1 * t) })
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fq6 {
    pub c0: Fq2,
    pub c1: Fq2,
    pub c2: Fq2,
}
ops!(Fq6, c0, c1, c2);

impl Fq6 {
    /// times v
    #[inline]
    pub fn mul_by_nonresidue(&self) -> Fq6 {
        Fq6 { c0: self.c2.mul_by_nonresidue(), c1: self.c0, c2: self.c1 }
    }
}

impl Mul for Fq6 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let (a0, a1, a2) = (self.c0, self.c1, self.c2);
        let (b0, b1, b2) = (o.c0, o.c1, o.c2);
        let aa = a0 * b0;
        let bb = a1 * b1;
        let cc = a2 * b2;
        let t1 = ((a1 + a2) * (b1 + b2) - bb - cc).mul_by_nonresidue() + aa;
        let t2 = (a0 + a1) * (b0 + b1) - aa - bb + cc.mul_by_nonresidue();
        let t3 = (a0 + a2) * (b0 + b2) - aa - cc + bb;
        Fq6 { c0: t1, c1: t2, c2: t3 }
    }
}

impl Field for Fq6 {
    const ZERO: Self = Fq6 { c0: Fq2::ZERO, c1: Fq2::ZERO, c2: Fq2::ZERO };
    const ONE: Self = Fq6 { c0: Fq2::ONE, c1: Fq2::ZERO, c2: Fq2::ZERO };
    fn inverse(&self) -> Option<Self> {
        let (a0, a1, a2) = (self.c0, self.c1, self.c2);
        let c0 = a0.square() - (a1 * a2).mul_by_nonresidue();
        let c1 = a2.square().mul_by_nonresidue() - a0 * a1;
        let c2 = a1.square() - a0 * a2;
        let t = (a0 * c0 + (a2 * c1).mul_by_nonresidue() + (a1 * c2).mul_by_nonresidue()).inverse()?;
        Some(Fq6 { c0: c0 * t, c1: c1 * t, c2: c2 * t })
    }
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct Fq12 {
    pub c0: Fq6,
    pub c1: Fq6,
}
ops!(Fq12, c0, c1);

impl Mul for Fq12 {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let aa = self.c0 * o.c0;
        let bb = self.c1 * o.c1;
        let c1 = (self.c0 + self.c1) * (o.c0 + o.c1) - aa - bb;
        Fq12 { c0: aa + bb.mul_by_nonresidue(), c1 }
    }
}

impl Field for Fq12 {
    const ZERO: Self = Fq12 { c0: Fq6::ZERO, c1: Fq6::ZERO };
    const ONE: Self = Fq12 { c0: Fq6::ONE, c1: Fq6::ZERO };
    fn inverse(&self) -> Option<Self> {
        let t = (self.c0.square() - self.c1.square().mul_by_nonresidue()).inverse()?;
        Some(Fq12 { c0: self.c0 * t, c1: -(self.c1 * t) })
    }
}

/// gamma[i] = xi^(i (q^2 - 1) / 6): w^(i q^2) = gamma[i] w^i
fn frobenius2_coeffs() -> &'static [Fq2; 6] {
    static C: OnceLock<[Fq2; 6]> = OnceLock::new();
    C.get_or_init(|| {
        const E: [u64; 12] = [0xdbc7000004bd97b4, 0xea2791da3e5eb271, 0x2e5cb340905834d4, 0xe67542fcdfbd9e8b, 0x86dd1646bd6d9ab6, 0x84e6baef6baeba14, 0x7e32e188819427d5, 0x62c65cd4d924f712, 0x667b9a6188c5174d, 0x507a18262d12b673, 0xe1f8697c705d30fc, 0x0070b3f0c975e54b];
        let g = Fq2::new(Fq::ONE, Fq::ONE).pow(&E);
        let mut c = [Fq2::ONE; 6];
        for i in 1..6 {
            c[i] = c[i - 1] * g;
        }
        c
    })
}

impl Fq12 {
    /// f^(q^6): the inverse for elements of the cyclotomic subgroup
    pub fn conjugate(&self) -> Fq12 {
        Fq12 { c0: self.c0, c1: -self.c1 }
    }

    /// f^(q^2). The coefficients of w^0..w^5 are c0.c0, c1.c0, c0.c1, c1.c1, c0.c2, c1.c2;
    /// conjugating twice is the identity on Fq2, so only the gammas apply.
    pub fn frobenius2(&self) -> Fq12 {
        let g = frobenius2_coeffs();
        Fq12 {
            c0: Fq6 { c0: self.c0.c0, c1: self.c0.c1 * g[2], c2: self.c0.c2 * g[4] },
            c1: Fq6 { c0: self.c1.c0 * g[1], c1: self.c1.c1 * g[3], c2: self.c1.c2 * g[5] },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fq2(a: u64, b: u64) -> Fq2 {
        Fq2::new(Fq::from_u64(a), Fq::from_u64(b))
    }

    fn sample12(k: u64) -> Fq12 {
        let s6 = |o: u64| Fq6 { c0: fq2(k + o, 3 * k + o), c1: fq2(5 + o, 7 * k), c2: fq2(11 * k + o, 13) };
        Fq12 { c0: s6(1), c1: s6(17) }
    }

    #[test]
    fn inverses_and_distributivity() {
        let a = sample12(3);
        let b = sample12(8);
        let c = sample12(21);
        assert_eq!(a * a.inverse().unwrap(), Fq12::ONE);
        assert_eq!(a * (b + c), a * b + a * c);
        assert_eq!((a * b) * c, a * (b * c));
        assert_eq!(a.square(), a * a);
        let x = fq2(9, 4);
        assert_eq!(x * x.inverse().unwrap(), Fq2::ONE);
        let y = Fq6 { c0: x, c1: fq2(1, 2), c2: fq2(3, 7) };
        assert_eq!(y * y.inverse().unwrap(), Fq6::ONE);
    }

    #[test]
    fn frobenius2_is_the_power() {
        let a = sample12(5);
        // q^2 as limbs: compute a^q then ^q
        let q = Fq::MODULUS;
        assert_eq!(a.frobenius2(), a.pow(&q).pow(&q));
    }

    #[test]
    fn fq2_sqrt() {
        for k in 1..20 {
            let x = fq2(k, 3 * k + 1);
            let s = x.square().sqrt().unwrap();
            assert!(s == x || s == -x);
        }
    }
}
