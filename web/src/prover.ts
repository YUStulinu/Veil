// The Rust prover compiled to WebAssembly (prover/src/wasm.rs), for the
// browser and for Node.
import { bytesToBig } from "./bls.ts";

export interface ProveResult {
  ok: boolean;
  error?: string;
  proof?: Uint8Array;
  vk?: Uint8Array;
  pub?: bigint[];
  ms: number;
}

const VK_BYTES = 7 + 8 + 96 + 8 * 128 + 512;

export async function loadProver(wasm: BufferSource) {
  const { instance } = await WebAssembly.instantiate(wasm, { env: { veil_now_ms: () => performance.now() } });
  const ex = instance.exports as any;
  const mem = () => new Uint8Array(ex.memory.buffer);
  const put = (b: Uint8Array) => {
    const p = ex.veil_alloc(b.length);
    mem().set(b, p);
    return p;
  };
  return (circuit: string, srs: Uint8Array, inputs: string): ProveResult => {
    const enc = new TextEncoder();
    const c = enc.encode(circuit);
    const i = enc.encode(inputs);
    const seed = crypto.getRandomValues(new Uint8Array(32));
    const ps = put(seed);
    ex.veil_seed(ps);
    ex.veil_free(ps, 32);
    const [pc, pv, pi] = [put(c), put(srs), put(i)];
    const t0 = performance.now();
    const out = ex.veil_prove(pc, c.length, pv, srs.length, pi, i.length);
    const ms = performance.now() - t0;
    ex.veil_free(pc, c.length);
    ex.veil_free(pv, srs.length);
    ex.veil_free(pi, i.length);
    const dv = new DataView(ex.memory.buffer, out, 8);
    const [status, len] = [dv.getUint32(0, true), dv.getUint32(4, true)];
    const body = mem().slice(out + 8, out + 8 + len);
    ex.veil_free(out, 8 + len);
    if (status !== 0) return { ok: false, error: new TextDecoder().decode(body), ms };
    const pub = [];
    for (let o = 624 + VK_BYTES; o < body.length; o += 32) pub.push(bytesToBig(body.slice(o, o + 32)));
    return { ok: true, proof: body.slice(0, 624), vk: body.slice(624, 624 + VK_BYTES), pub, ms };
  };
}
