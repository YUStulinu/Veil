// The TypeScript side: BLS12-381 basics, Keccak, the verifier on a proof
// from the Rust CLI, and the Wasm prover's proof checked by the TS verifier.
import { test } from "node:test";
import assert from "node:assert/strict";
import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { join } from "node:path";
import { G1_GEN, R, g1FromCompressed, g1mul, toAff, toJac, pairing, pairingCheck, g1neg, f12pow, f12eq } from "../src/bls.ts";
import { keccak256 } from "../src/keccak.ts";
import { parseVk, verifyBytes } from "../src/verifier.ts";
import { loadProver } from "../src/prover.ts";

const ROOT = join(import.meta.dirname, "..", "..");
const BUILD = join(ROOT, "build");
const PUBLIC = join(import.meta.dirname, "..", "public");
const hex = (b: Uint8Array) => Buffer.from(b).toString("hex");
/** timings for bench/bench.py */
const timings: Record<string, number> = {};
const record = (k: string, ms: number) => {
  timings[k] = Math.round(ms);
  if (existsSync(BUILD)) writeFileSync(join(BUILD, "web_timings.json"), JSON.stringify(timings, null, 1));
};

test("keccak256 vectors", () => {
  assert.equal(hex(keccak256(new Uint8Array())), "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470");
  assert.equal(hex(keccak256(new TextEncoder().encode("abc"))), "4e03657aea45a94fc7d47ba826c8d667c0d1e6e33a64a036ec44f58fa12d6c45");
  const long = new Uint8Array(300).fill(7); // more than two blocks
  assert.equal(keccak256(long).length, 32);
});

test("G1 decompression of the generator and the pairing's bilinearity", () => {
  const g = g1FromCompressed(Buffer.from("97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb", "hex"));
  assert.equal(g.x, G1_GEN.x);
  assert.equal(g.y, G1_GEN.y);
  const vk = parseVk(readFileSync(join(PUBLIC, "vk.bin")));
  const h = vk.g2[0];
  const a = 0x123456789n, b = 0xfedcban;
  const ga = toAff(g1mul(toJac(G1_GEN), a));
  const gab = toAff(g1mul(toJac(G1_GEN), (a * b) % R));
  // e(abG, H) = e(aG, H)^b, and e(aG, H) e(-aG, H) = 1
  assert.ok(f12eq(pairing(gab, h), f12pow(pairing(ga, h), b)));
  assert.ok(!f12eq(pairing(ga, h), pairing(G1_GEN, h)));
  assert.ok(pairingCheck([[ga, h], [g1neg(ga), h]]));
  assert.ok(!pairingCheck([[ga, h], [g1neg(G1_GEN), h]]));
});

const cli = existsSync(join(BUILD, "proof.bin"));
test("the TS verifier accepts the Rust CLI's Sudoku proof and rejects changes", { skip: !cli && "no build/: run scripts/sudoku.sh" }, () => {
  const vk = parseVk(readFileSync(join(BUILD, "vk.bin")));
  const proof = new Uint8Array(readFileSync(join(BUILD, "proof.bin")));
  const pub = readFileSync(join(BUILD, "public.txt"), "utf8").trim().split(/\s+/).map(BigInt);
  const t0 = performance.now();
  assert.ok(verifyBytes(vk, pub, proof));
  const ms = performance.now() - t0;
  console.log(`TS verification: ${ms.toFixed(0)} ms`);
  record("ts_verify_ms", ms);
  for (const i of [0, 3, 8]) {
    const p = proof.slice();
    p[48 * i + 47] ^= 1;
    assert.equal(verifyBytes(vk, pub, p), false, `point ${i}`);
  }
  const q = proof.slice();
  q[432 + 31] ^= 1;
  assert.equal(verifyBytes(vk, pub, q), false, "an evaluation");
  const other = pub.slice();
  other[0] = other[0] === 0n ? 1n : 0n;
  assert.equal(verifyBytes(vk, other, proof), false, "another puzzle");
});

const wasmReady = existsSync(join(PUBLIC, "veil.wasm"));
test("the Wasm prover's proof verifies in TS; a wrong solution is refused", { skip: !wasmReady && "no web/public/veil.wasm" }, async () => {
  const prove = await loadProver(readFileSync(join(PUBLIC, "veil.wasm")));
  const circuit = readFileSync(join(PUBLIC, "sudoku.circuit"), "utf8");
  const srs = new Uint8Array(readFileSync(join(PUBLIC, "srs.bin")));
  const inputs = readFileSync(join(ROOT, "lang", "examples", "sudoku.inputs"), "utf8");
  const r = prove(circuit, srs, inputs);
  assert.ok(r.ok, r.error);
  console.log(`Wasm proving: ${r.ms.toFixed(0)} ms`);
  record("wasm_prove_ms", r.ms);
  const vk = parseVk(readFileSync(join(PUBLIC, "vk.bin")));
  assert.equal(hex(r.vk!), hex(vk.bytes), "the Wasm prover derives the published verifying key");
  assert.ok(verifyBytes(vk, r.pub!, r.proof!));
  const bad = inputs.replace(/solution: (\d)/, (_, d) => `solution: ${d === "5" ? 4 : 5}`);
  const r2 = prove(circuit, srs, bad);
  assert.equal(r2.ok, false);
});
