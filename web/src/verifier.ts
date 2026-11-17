// The PLONK verifier in TypeScript: the same transcript and check as the
// Rust and Solidity verifiers, on its own BLS12-381 code.
import {
  R, G1_GEN, G1_ZERO, mod, inv, bytesToBig, bigToBytes,
  g1FromCompressed, g1FromEip2537, g1ToEip2537, g2FromEip2537, g1add, g1mul, g1neg, toAff, toJac, pairingCheck,
} from "./bls.ts";
import type { G1, G1Aff, G2Aff } from "./bls.ts";
import { keccak256 } from "./keccak.ts";

export interface VerifyingKey {
  bytes: Uint8Array;
  n: number;
  numPublic: number;
  k1: bigint;
  k2: bigint;
  omega: bigint;
  commitments: G1Aff[]; // qM qL qR qO qC S1 S2 S3
  g2: [G2Aff, G2Aff]; // [1]_2, [tau]_2
}

export interface Proof {
  points: G1Aff[]; // a b c z tlo tmid thi W W'
  evals: bigint[]; // a b c s1 s2 z(zeta omega)
}

export const PROOF_BYTES = 624;

export function parseVk(b: Uint8Array): VerifyingKey {
  if (b.length !== 7 + 8 + 96 + 8 * 128 + 512 || new TextDecoder().decode(b.slice(0, 7)) !== "VEILVK1") throw new Error("not a verifying key");
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength);
  const fr = (o: number) => bytesToBig(b.slice(o, o + 32));
  const commitments = [];
  for (let i = 0; i < 8; i++) commitments.push(g1FromEip2537(b.slice(111 + 128 * i, 239 + 128 * i)));
  const g2: [G2Aff, G2Aff] = [g2FromEip2537(b.slice(1135, 1391)), g2FromEip2537(b.slice(1391, 1647))];
  return { bytes: b, n: dv.getUint32(7), numPublic: dv.getUint32(11), k1: fr(15), k2: fr(47), omega: fr(79), commitments, g2 };
}

/** checks every point (curve, subgroup) and every scalar (< r) */
export function parseProof(b: Uint8Array): Proof {
  if (b.length !== PROOF_BYTES) throw new Error(`a proof has ${PROOF_BYTES} bytes`);
  const points = [];
  for (let i = 0; i < 9; i++) points.push(g1FromCompressed(b.slice(48 * i, 48 * i + 48)));
  const evals = [];
  for (let i = 0; i < 6; i++) {
    const e = bytesToBig(b.slice(432 + 32 * i, 464 + 32 * i));
    if (e >= R) throw new Error("scalar not below r");
    evals.push(e);
  }
  return { points, evals };
}

class Transcript {
  state: Uint8Array;
  pending: Uint8Array[] = [];
  constructor(vk: Uint8Array) {
    this.state = keccak256(concat([new TextEncoder().encode("veil-plonk-v1"), vk]));
  }
  fr(x: bigint) {
    this.pending.push(bigToBytes(x, 32));
  }
  g1(p: G1Aff) {
    this.pending.push(g1ToEip2537(p));
  }
  challenge(): bigint {
    this.state = keccak256(concat([this.state, ...this.pending]));
    this.pending = [];
    return mod(bytesToBig(this.state), R);
  }
}

function concat(parts: Uint8Array[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let o = 0;
  for (const p of parts) {
    out.set(p, o);
    o += p.length;
  }
  return out;
}

const fmul = (...xs: bigint[]) => xs.reduce((a, b) => (a * b) % R, 1n);
const fpow = (x: bigint, e: bigint) => {
  let r = 1n;
  while (e > 0n) {
    if (e & 1n) r = (r * x) % R;
    x = (x * x) % R;
    e >>= 1n;
  }
  return r;
};

export function verify(vk: VerifyingKey, pub: bigint[], proof: Proof): boolean {
  if (pub.length !== vk.numPublic || pub.some((x) => x < 0n || x >= R)) return false;
  const tr = new Transcript(vk.bytes);
  pub.forEach((x) => tr.fr(x));
  const p = proof.points;
  p.slice(0, 3).forEach((c) => tr.g1(c));
  const beta = tr.challenge();
  const gamma = tr.challenge();
  tr.g1(p[3]);
  const alpha = tr.challenge();
  p.slice(4, 7).forEach((c) => tr.g1(c));
  const zeta = tr.challenge();
  proof.evals.forEach((e) => tr.fr(e));
  const v = tr.challenge();
  tr.g1(p[7]);
  tr.g1(p[8]);
  const u = tr.challenge();

  const [ea, eb, ec, es1, es2, ezw] = proof.evals;
  const n = BigInt(vk.n);
  const zn = fpow(zeta, n);
  const zh = mod(zn - 1n, R);
  const nInv = inv(n, R);
  // L_i(zeta) = omega^i Z_H(zeta) / (n (zeta - omega^i))
  let w = 1n;
  let l0 = 0n;
  let pi = 0n;
  for (let i = 0; i < Math.max(vk.numPublic, 1); i++) {
    const d = mod(zeta - w, R);
    if (d === 0n) return false;
    const li = fmul(w, zh, nInv, inv(d, R));
    if (i === 0) l0 = li;
    if (i < vk.numPublic) pi = mod(pi - pub[i] * li, R);
    w = fmul(w, vk.omega);
  }
  const a2 = fmul(alpha, alpha);
  const permS = fmul(alpha, ea + beta * es1 + gamma, eb + beta * es2 + gamma, ezw);
  const r0 = mod(pi - fmul(permS, ec + gamma) - fmul(a2, l0), R);
  const bz = fmul(beta, zeta);
  const permZ = mod(fmul(alpha, ea + bz + gamma, eb + bz * vk.k1 + gamma, ec + bz * vk.k2 + gamma) + fmul(a2, l0), R);
  const [v2, v3, v4, v5] = [fmul(v, v), fmul(v, v, v), fmul(v, v, v, v), fmul(v, v, v, v, v)];
  const e = mod(-r0 + v * ea + v2 * eb + v3 * ec + v4 * es1 + v5 * es2 + u * ezw, R);
  const [qm, ql, qr, qo, qc, s1, s2, s3] = vk.commitments;
  const terms: [G1Aff, bigint][] = [
    [qm, fmul(ea, eb)], [ql, ea], [qr, eb], [qo, ec], [qc, 1n],
    [p[3], mod(permZ + u, R)], [s3, mod(-permS * beta, R)],
    [p[4], mod(-zh, R)], [p[5], mod(-zh * zn, R)], [p[6], mod(-zh * zn * zn, R)],
    [p[0], v], [p[1], v2], [p[2], v3], [s1, v4], [s2, v5],
    [G1_GEN, mod(-e, R)], [p[7], zeta], [p[8], fmul(u, zeta, vk.omega)],
  ];
  let rhs: G1 = G1_ZERO;
  for (const [pt, k] of terms) rhs = g1add(rhs, g1mul(toJac(pt), mod(k, R)));
  const lhs = g1add(toJac(p[7]), g1mul(toJac(p[8]), u));
  return pairingCheck([
    [toAff(lhs), vk.g2[1]],
    [g1neg(toAff(rhs)), vk.g2[0]],
  ]);
}

/** parse and verify; false for anything malformed */
export function verifyBytes(vk: VerifyingKey, pub: bigint[], proof: Uint8Array): boolean {
  try {
    return verify(vk, pub, parseProof(proof));
  } catch {
    return false;
  }
}
