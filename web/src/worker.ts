// The heavy work off the page's thread: proving (Rust in Wasm) and
// verifying (TypeScript).
import { loadProver } from "./prover.ts";
import { parseVk, verifyBytes } from "./verifier.ts";
import type { VerifyingKey } from "./verifier.ts";

let prove: Awaited<ReturnType<typeof loadProver>> | null = null;
let circuit = "";
let srs = new Uint8Array();
let vk: VerifyingKey | null = null;

async function init() {
  if (prove) return;
  const get = (f: string) => fetch(new URL(`../public/${f}`, import.meta.url));
  const [w, c, s, k] = await Promise.all([get("veil.wasm"), get("sudoku.circuit"), get("srs.bin"), get("vk.bin")]);
  prove = await loadProver(await w.arrayBuffer());
  circuit = await c.text();
  srs = new Uint8Array(await s.arrayBuffer());
  vk = parseVk(new Uint8Array(await k.arrayBuffer()));
}

self.onmessage = async (ev: MessageEvent) => {
  const m = ev.data;
  try {
    await init();
    if (m.type === "prove") {
      const r = prove!(circuit, srs, m.inputs);
      self.postMessage({ type: "proved", ...r, pub: r.pub?.map(String) });
    } else if (m.type === "verify") {
      const t0 = performance.now();
      const ok = verifyBytes(vk!, m.pub.map(BigInt), m.proof);
      self.postMessage({ type: "verified", ok, ms: performance.now() - t0 });
    }
  } catch (e) {
    self.postMessage({ type: "error", error: String(e) });
  }
};
