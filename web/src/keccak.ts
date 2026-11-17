// Keccak-256 (Ethereum's: padding 0x01), on 64-bit BigInt lanes.

const MASK = (1n << 64n) - 1n;
const RC = [
  0x0000000000000001n, 0x0000000000008082n, 0x800000000000808an, 0x8000000080008000n, 0x000000000000808bn, 0x0000000080000001n,
  0x8000000080008081n, 0x8000000000008009n, 0x000000000000008an, 0x0000000000000088n, 0x0000000080008009n, 0x000000008000000an,
  0x000000008000808bn, 0x800000000000008bn, 0x8000000000008089n, 0x8000000000008003n, 0x8000000000008002n, 0x8000000000000080n,
  0x000000000000800an, 0x800000008000000an, 0x8000000080008081n, 0x8000000000008080n, 0x0000000080000001n, 0x8000000080008008n,
];
const ROT = [0, 1, 62, 28, 27, 36, 44, 6, 55, 20, 3, 10, 43, 25, 39, 41, 45, 15, 21, 8, 18, 2, 61, 56, 14].map(BigInt);
const rotl = (x: bigint, n: bigint) => (n === 0n ? x : ((x << n) | (x >> (64n - n))) & MASK);

function keccakF(a: bigint[]) {
  for (const rc of RC) {
    const c = [0, 1, 2, 3, 4].map((x) => a[x] ^ a[x + 5] ^ a[x + 10] ^ a[x + 15] ^ a[x + 20]);
    for (let x = 0; x < 5; x++) {
      const d = c[(x + 4) % 5] ^ rotl(c[(x + 1) % 5], 1n);
      for (let y = 0; y < 25; y += 5) a[x + y] ^= d;
    }
    const b = new Array<bigint>(25);
    for (let x = 0; x < 5; x++) for (let y = 0; y < 5; y++) b[y + 5 * ((2 * x + 3 * y) % 5)] = rotl(a[x + 5 * y], ROT[x + 5 * y]);
    for (let x = 0; x < 5; x++) for (let y = 0; y < 5; y++) a[x + 5 * y] = b[x + 5 * y] ^ (~b[((x + 1) % 5) + 5 * y] & MASK & b[((x + 2) % 5) + 5 * y]);
    a[0] ^= rc;
  }
}

export function keccak256(data: Uint8Array): Uint8Array {
  const rate = 136;
  const len = Math.ceil((data.length + 1) / rate) * rate;
  const msg = new Uint8Array(len);
  msg.set(data);
  msg[data.length] ^= 0x01;
  msg[len - 1] ^= 0x80;
  const st = new Array<bigint>(25).fill(0n);
  for (let o = 0; o < len; o += rate) {
    for (let i = 0; i < rate / 8; i++) {
      let lane = 0n;
      for (let k = 7; k >= 0; k--) lane = (lane << 8n) | BigInt(msg[o + 8 * i + k]);
      st[i] ^= lane;
    }
    keccakF(st);
  }
  const out = new Uint8Array(32);
  for (let i = 0; i < 4; i++) for (let k = 0; k < 8; k++) out[8 * i + k] = Number((st[i] >> BigInt(8 * k)) & 0xffn);
  return out;
}
