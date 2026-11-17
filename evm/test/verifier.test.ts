// The generated Solidity verifier in an EVM (EthereumJS, Prague: EIP-2537):
// a real Sudoku proof from the Rust prover verifies, tampered ones do not,
// and the gas is measured. Needs build/ from `scripts/sudoku.sh`.
import { test } from "node:test";
import assert from "node:assert/strict";
import { readFileSync, writeFileSync, existsSync } from "node:fs";
import { createRequire } from "node:module";
import { join } from "node:path";
import { createEVM } from "@ethereumjs/evm";
import { Common, Hardfork, Mainnet } from "@ethereumjs/common";
import { keccak_256 } from "@noble/hashes/sha3.js";

const require = createRequire(import.meta.url);
const solc = require("solc");
const BUILD = process.env.VEIL_BUILD ?? join(import.meta.dirname, "..", "..", "build");
const R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001n;

function compile(source: string): { bytecode: Uint8Array; size: number } {
  const input = {
    language: "Solidity",
    sources: { "Verifier.sol": { content: source } },
    settings: {
      optimizer: { enabled: true, runs: 200 },
      viaIR: true,
      evmVersion: "prague",
      outputSelection: { "*": { "*": ["evm.bytecode.object", "evm.deployedBytecode.object"] } },
    },
  };
  const out = JSON.parse(solc.compile(JSON.stringify(input)));
  const errors = (out.errors ?? []).filter((e: any) => e.severity === "error");
  if (errors.length) throw new Error(errors.map((e: any) => e.formattedMessage).join("\n"));
  const c = out.contracts["Verifier.sol"].VeilVerifier;
  return { bytecode: Buffer.from(c.evm.bytecode.object, "hex"), size: c.evm.deployedBytecode.object.length / 2 };
}

const word = (x: bigint) => Buffer.from(x.toString(16).padStart(64, "0"), "hex");

function calldata(proof: Uint8Array, pub: bigint[]): Uint8Array {
  const sel = keccak_256(new TextEncoder().encode("verify(bytes,uint256[])")).slice(0, 4);
  const padded = Buffer.concat([Buffer.from(proof), Buffer.alloc((32 - (proof.length % 32)) % 32)]);
  const arrOffset = 64n + 32n + BigInt(padded.length);
  return Buffer.concat([
    Buffer.from(sel),
    word(64n),
    word(arrOffset),
    word(BigInt(proof.length)),
    padded,
    word(BigInt(pub.length)),
    ...pub.map(word),
  ]);
}

async function setup() {
  const sol = readFileSync(join(BUILD, "Verifier.sol"), "utf8");
  const { bytecode, size } = compile(sol);
  const evm = await createEVM({ common: new Common({ chain: Mainnet, hardfork: Hardfork.Prague }) });
  const dep = await evm.runCall({ data: bytecode, gasLimit: 30_000_000n });
  assert.equal(dep.execResult.exceptionError, undefined, "deployment failed");
  const addr = dep.createdAddress!;
  const call = async (proof: Uint8Array, pub: bigint[]) => {
    const r = await evm.runCall({ to: addr, data: calldata(proof, pub), gasLimit: 30_000_000n });
    if (r.execResult.exceptionError) return { ok: false, gas: r.execResult.executionGasUsed, reverted: true };
    const ret = r.execResult.returnValue;
    return { ok: ret.length === 32 && ret[31] === 1, gas: r.execResult.executionGasUsed, reverted: false };
  };
  return { call, size, deployGas: dep.execResult.executionGasUsed };
}

const ready = existsSync(join(BUILD, "Verifier.sol"));

test("the Solidity verifier accepts the Rust prover's Sudoku proof and rejects tampering", { skip: !ready && "no build/: run scripts/sudoku.sh" }, async () => {
  const { call, size, deployGas } = await setup();
  const proof = new Uint8Array(readFileSync(join(BUILD, "proof.evm")));
  const pub = readFileSync(join(BUILD, "public.txt"), "utf8").trim().split(/\s+/).map(BigInt);
  assert.equal(proof.length, 1344);

  const good = await call(proof, pub);
  assert.ok(good.ok, "a valid proof must verify");
  // the calldata cost on top of the execution: 16 gas per nonzero byte, 4 per zero byte
  const data = calldata(proof, pub);
  const calldataGas = data.reduce((g, b) => g + (b === 0 ? 4 : 16), 0);
  const report = { execution_gas: Number(good.gas), calldata_gas: calldataGas, transaction_gas: Number(good.gas) + calldataGas + 21000, calldata_bytes: data.length, deployed_bytes: size, deploy_gas: Number(deployGas) };
  console.log(JSON.stringify(report));
  writeFileSync(join(BUILD, "evm_gas.json"), JSON.stringify(report, null, 1));

  // every proof element changed
  for (let i = 0; i < 9; i++) {
    const p = proof.slice();
    p[128 * i + 127] ^= 1; // the last byte of y: the point leaves the curve
    assert.equal((await call(p, pub)).ok, false, `point ${i} off the curve`);
    const q = proof.slice();
    q.set(proof.slice(128 * ((i + 1) % 9), 128 * ((i + 1) % 9) + 128), 128 * i); // a valid point, the wrong one
    assert.equal((await call(q, pub)).ok, false, `point ${i} swapped`);
  }
  for (let i = 0; i < 6; i++) {
    const p = proof.slice();
    const o = 1152 + 32 * i;
    const v = (BigInt("0x" + Buffer.from(p.slice(o, o + 32)).toString("hex")) + 1n) % R;
    p.set(word(v), o);
    assert.equal((await call(p, pub)).ok, false, `evaluation ${i}`);
  }
  // public inputs: another puzzle, a value >= r, too few
  const other = pub.slice();
  other[2] = other[2] === 0n ? 4n : 0n;
  assert.equal((await call(proof, other)).ok, false, "another puzzle");
  const big = pub.slice();
  big[0] = big[0] + R;
  assert.equal((await call(proof, big)).ok, false, "a public input >= r");
  assert.equal((await call(proof, pub.slice(1))).ok, false, "too few public inputs");
  assert.equal((await call(proof.slice(0, 1343), pub)).ok, false, "a short proof");
});
