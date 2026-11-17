//! ChaCha20 as a random generator, seeded from the operating system (or a
//! given seed for reproducible tests). The prover's blinding factors come
//! from here, so they must not be predictable.

pub struct Rng {
    key: [u32; 8],
    counter: u64,
    buf: [u8; 64],
    pos: usize,
}

fn quarter(s: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize) {
    s[a] = s[a].wrapping_add(s[b]);
    s[d] = (s[d] ^ s[a]).rotate_left(16);
    s[c] = s[c].wrapping_add(s[d]);
    s[b] = (s[b] ^ s[c]).rotate_left(12);
    s[a] = s[a].wrapping_add(s[b]);
    s[d] = (s[d] ^ s[a]).rotate_left(8);
    s[c] = s[c].wrapping_add(s[d]);
    s[b] = (s[b] ^ s[c]).rotate_left(7);
}

/// one ChaCha20 block (RFC 8439 with a 64-bit counter and zero nonce)
pub fn chacha20_block(key: &[u32; 8], counter: u64, nonce: [u32; 2]) -> [u8; 64] {
    let mut init = [0u32; 16];
    init[..4].copy_from_slice(&[0x61707865, 0x3320646e, 0x79622d32, 0x6b206574]);
    init[4..12].copy_from_slice(key);
    init[12] = counter as u32;
    init[13] = (counter >> 32) as u32;
    init[14] = nonce[0];
    init[15] = nonce[1];
    let mut s = init;
    for _ in 0..10 {
        quarter(&mut s, 0, 4, 8, 12);
        quarter(&mut s, 1, 5, 9, 13);
        quarter(&mut s, 2, 6, 10, 14);
        quarter(&mut s, 3, 7, 11, 15);
        quarter(&mut s, 0, 5, 10, 15);
        quarter(&mut s, 1, 6, 11, 12);
        quarter(&mut s, 2, 7, 8, 13);
        quarter(&mut s, 3, 4, 9, 14);
    }
    let mut out = [0u8; 64];
    for i in 0..16 {
        out[4 * i..4 * i + 4].copy_from_slice(&s[i].wrapping_add(init[i]).to_le_bytes());
    }
    out
}

impl Rng {
    pub fn from_seed(seed: [u8; 32]) -> Rng {
        let mut key = [0u32; 8];
        for i in 0..8 {
            key[i] = u32::from_le_bytes(seed[4 * i..4 * i + 4].try_into().unwrap());
        }
        Rng { key, counter: 0, buf: [0; 64], pos: 64 }
    }

    /// seeded from the operating system's generator
    pub fn from_os() -> Rng {
        Rng::from_seed(os_entropy())
    }

    pub fn fill(&mut self, out: &mut [u8]) {
        for b in out {
            if self.pos == 64 {
                self.buf = chacha20_block(&self.key, self.counter, [0, 0]);
                self.counter += 1;
                self.pos = 0;
            }
            *b = self.buf[self.pos];
            self.pos += 1;
        }
    }

    pub fn bytes<const N: usize>(&mut self) -> [u8; N] {
        let mut b = [0u8; N];
        self.fill(&mut b);
        b
    }

    pub fn u64(&mut self) -> u64 {
        u64::from_le_bytes(self.bytes::<8>())
    }
}

#[cfg(windows)]
fn os_entropy() -> [u8; 32] {
    #[link(name = "bcrypt")]
    extern "system" {
        fn BCryptGenRandom(alg: *mut u8, buf: *mut u8, len: u32, flags: u32) -> i32;
    }
    let mut b = [0u8; 32];
    // BCRYPT_USE_SYSTEM_PREFERRED_RNG = 2
    let st = unsafe { BCryptGenRandom(std::ptr::null_mut(), b.as_mut_ptr(), 32, 2) };
    assert_eq!(st, 0, "BCryptGenRandom failed");
    b
}

#[cfg(all(unix, not(target_arch = "wasm32")))]
fn os_entropy() -> [u8; 32] {
    use std::io::Read;
    let mut b = [0u8; 32];
    std::fs::File::open("/dev/urandom").and_then(|mut f| f.read_exact(&mut b)).expect("/dev/urandom");
    b
}

#[cfg(target_arch = "wasm32")]
fn os_entropy() -> [u8; 32] {
    // the browser passes its entropy in (crypto.getRandomValues) before proving
    crate::wasm::take_seed()
}

#[cfg(test)]
mod tests {
    #[test]
    fn rfc8439_block() {
        // RFC 8439 2.3.2: key 00..1f, counter 1, nonce 000000090000004a00000000
        let mut key = [0u32; 8];
        for i in 0..8 {
            key[i] = u32::from_le_bytes([4 * i as u8, 4 * i as u8 + 1, 4 * i as u8 + 2, 4 * i as u8 + 3]);
        }
        // our layout is counter(64) + nonce(64); the RFC's is counter(32) + nonce(96)
        let counter = 1u64 | (0x09000000u64 << 32);
        let out = super::chacha20_block(&key, counter, [0x4a000000, 0]);
        assert_eq!(&out[..8], &[0x10, 0xf1, 0xe7, 0xe4, 0xd1, 0x3b, 0x59, 0x15]);
    }
}
