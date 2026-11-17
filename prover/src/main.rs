//! veil: the command line.
//!
//!   veil prove   --circuit C --inputs I --srs S [--out proof.bin] [--public public.txt] [--vk vk.bin] [--evm proof.evm]
//!   veil verify  --vk vk.bin --proof proof.bin --public public.txt
//!   veil vk      --circuit C --srs S --out vk.bin
//!   veil export-sol --vk vk.bin --out Verifier.sol
//!   veil export-srs --srs S --count N --out srs.bin
//!   veil bench   --log-n K [--gpu] [--reps R]        one JSON line
//!   veil kernels --log-n K [--gpu]                    MSM and NTT timings, one JSON line
//!
//! An SRS S is `ceremony:<trusted_setup.txt>` (the Ethereum KZG ceremony,
//! 4096 powers), `insecure:<powers>` (a known tau: benchmarks only), or a
//! file written by export-srs.

use std::collections::HashMap;
use std::fs;
use std::time::Instant;
use veil::ff::Fr;
use veil::kzg::Srs;
use veil::plonk::circuit::{parse_canonical, parse_inputs, synthetic, to_decimal_plain, Circuit};
use veil::plonk::{preprocess, prove, verify, Proof, VerifyingKey};
use veil::rng::Rng;

fn args() -> (String, HashMap<String, String>) {
    let mut it = std::env::args().skip(1);
    let cmd = it.next().unwrap_or_default();
    let rest: Vec<String> = it.collect();
    let mut m = HashMap::new();
    let mut i = 0;
    while i < rest.len() {
        let k = rest[i].trim_start_matches("--").to_string();
        if i + 1 < rest.len() && !rest[i + 1].starts_with("--") {
            m.insert(k, rest[i + 1].clone());
            i += 2;
        } else {
            m.insert(k, String::new());
            i += 1;
        }
    }
    (cmd, m)
}

fn fail(msg: &str) -> ! {
    eprintln!("veil: {msg}");
    std::process::exit(1)
}

fn need<'a>(m: &'a HashMap<String, String>, k: &str) -> &'a str {
    m.get(k).map(String::as_str).unwrap_or_else(|| fail(&format!("missing --{k}")))
}

fn read(path: &str) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|e| fail(&format!("{path}: {e}")))
}

fn read_text(path: &str) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| fail(&format!("{path}: {e}")))
}

/// The insecure SRS uses this fixed tau, so its files can be cached by the
/// benchmarks. Anyone can forge proofs against it.
const INSECURE_TAU: u64 = 0x5eed_0f7a0;

fn load_srs(spec: &str) -> Srs {
    if let Some(path) = spec.strip_prefix("ceremony:") {
        Srs::from_ceremony_file(&read_text(path)).unwrap_or_else(|e| fail(&e))
    } else if let Some(n) = spec.strip_prefix("insecure:") {
        let n: usize = n.parse().unwrap_or_else(|_| fail("insecure:<powers>"));
        let cache = std::env::temp_dir().join(format!("veil_insecure_srs_{n}.bin"));
        if let Ok(s) = fs::read(&cache).map_err(|e| e.to_string()).and_then(|b| Srs::from_bytes(&b, "insecure-local")) {
            return s;
        }
        let s = Srs::insecure_from_tau(Fr::from_u64(INSECURE_TAU), n);
        let _ = fs::write(&cache, s.to_bytes(n));
        s
    } else {
        Srs::from_bytes(&read(spec), "file").unwrap_or_else(|e| fail(&e))
    }
}

fn read_circuit(path: &str) -> Circuit {
    Circuit::parse(&read_text(path)).unwrap_or_else(|e| fail(&format!("{path}: {e}")))
}

fn parse_public(text: &str) -> Vec<Fr> {
    text.split_whitespace().map(|s| parse_canonical(s).unwrap_or_else(|| fail(&format!("public value {s} is not an integer below r")))).collect()
}

fn read_vk(path: &str) -> VerifyingKey {
    VerifyingKey::from_bytes(&read(path)).unwrap_or_else(|e| fail(&e))
}

fn main() {
    let (cmd, a) = args();
    if a.contains_key("gpu") {
        #[cfg(feature = "cuda")]
        veil::gpu::enable();
        #[cfg(not(feature = "cuda"))]
        fail("built without CUDA (cargo build --features cuda)");
    }
    match cmd.as_str() {
        "prove" => {
            let c = read_circuit(need(&a, "circuit"));
            let inputs = parse_inputs(&c, &read_text(need(&a, "inputs"))).unwrap_or_else(|e| fail(&e));
            let w = c.witness(&inputs).unwrap_or_else(|e| fail(&format!("the inputs do not satisfy the circuit: {e:?}")));
            let srs = load_srs(need(&a, "srs"));
            let t0 = Instant::now();
            let pk = preprocess(&c, &srs).unwrap_or_else(|e| fail(&e));
            let t_setup = t0.elapsed().as_secs_f64();
            let public = c.public_values(&w);
            let t0 = Instant::now();
            let (proof, tm) = prove(&pk, &public, &w, &mut Rng::from_os());
            let t_prove = t0.elapsed().as_secs_f64();
            let bytes = proof.to_bytes();
            fs::write(a.get("out").map_or("proof.bin", String::as_str), &bytes).unwrap();
            let pub_text: Vec<String> = public.iter().map(to_decimal_plain).collect();
            fs::write(a.get("public").map_or("public.txt", String::as_str), pub_text.join(" ") + "\n").unwrap();
            if let Some(p) = a.get("vk") {
                fs::write(p, pk.vk.to_bytes()).unwrap();
            }
            if let Some(p) = a.get("evm") {
                fs::write(p, proof.to_evm_bytes()).unwrap();
            }
            eprintln!(
                "{}: {} rows, n = {}, SRS {}, preprocessing {:.3} s, proving {:.3} s (rounds {:?}), proof {} bytes",
                c.name,
                c.rows(),
                pk.n,
                srs.origin,
                t_setup,
                t_prove,
                tm.rounds.map(|x| (x * 1000.0).round() / 1000.0),
                bytes.len()
            );
        }
        "verify" => {
            let vk = read_vk(need(&a, "vk"));
            let public = parse_public(&read_text(need(&a, "public")));
            let t0 = Instant::now();
            let r = Proof::from_bytes(&read(need(&a, "proof"))).and_then(|p| verify(&vk, &public, &p).map_err(|e| format!("{e:?}")));
            match r {
                Ok(()) => println!("valid ({:.1} ms)", t0.elapsed().as_secs_f64() * 1e3),
                Err(e) => {
                    println!("invalid: {e}");
                    std::process::exit(2);
                }
            }
        }
        "vk" => {
            let c = read_circuit(need(&a, "circuit"));
            let pk = preprocess(&c, &load_srs(need(&a, "srs"))).unwrap_or_else(|e| fail(&e));
            fs::write(need(&a, "out"), pk.vk.to_bytes()).unwrap();
        }
        "export-srs" => {
            let srs = load_srs(need(&a, "srs"));
            let n: usize = need(&a, "count").parse().unwrap_or_else(|_| fail("--count"));
            fs::write(need(&a, "out"), srs.to_bytes(n.min(srs.g1.len()))).unwrap();
        }
        "export-sol" => {
            let vk = read_vk(need(&a, "vk"));
            fs::write(need(&a, "out"), veil::solidity::verifier_contract(&vk)).unwrap();
        }
        "bench" => bench(&a),
        "kernels" => kernels(&a),
        _ => fail("commands: prove, verify, vk, export-srs, export-sol, bench, kernels"),
    }
}

/// a synthetic circuit filling 2^k rows: preprocessing, proving, verifying
fn bench(a: &HashMap<String, String>) {
    let k: u32 = need(a, "log-n").parse().unwrap_or_else(|_| fail("--log-n"));
    let reps: usize = a.get("reps").and_then(|r| r.parse().ok()).unwrap_or(1);
    let n = 1usize << k;
    let (c, inputs) = synthetic(n - 2 - 1); // 2 public rows, the chain, the equality gate
    let srs = load_srs(&format!("insecure:{}", n + 6));
    let w = c.witness(&inputs).unwrap();
    let public = c.public_values(&w);
    let t0 = Instant::now();
    let pk = preprocess(&c, &srs).unwrap();
    let t_setup = t0.elapsed().as_secs_f64();
    assert_eq!(pk.n, n);
    let (mut best, mut rounds, mut proof) = (f64::MAX, [0.0; 5], None);
    let mut rng = Rng::from_os();
    for _ in 0..reps {
        let t0 = Instant::now();
        let (p, tm) = prove(&pk, &public, &w, &mut rng);
        let t = t0.elapsed().as_secs_f64();
        if t < best {
            best = t;
            rounds = tm.rounds;
        }
        proof = Some(p);
    }
    let proof = proof.unwrap();
    let mut tv = f64::MAX;
    for _ in 0..3 {
        let t0 = Instant::now();
        verify(&pk.vk, &public, &proof).unwrap();
        tv = tv.min(t0.elapsed().as_secs_f64());
    }
    println!(
        "{{\"log_n\": {k}, \"rows\": {}, \"gpu\": {}, \"preprocess_s\": {t_setup:.4}, \"prove_s\": {best:.4}, \"rounds_s\": [{}], \"verify_ms\": {:.3}, \"proof_bytes\": {}, \"threads\": {}}}",
        c.rows(),
        a.contains_key("gpu"),
        rounds.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(", "),
        tv * 1e3,
        proof.to_bytes().len(),
        std::thread::available_parallelism().map_or(1, |t| t.get())
    );
}

/// an MSM of size 2^k over SRS points and an NTT of size 2^k
fn kernels(a: &HashMap<String, String>) {
    let k: u32 = need(a, "log-n").parse().unwrap_or_else(|_| fail("--log-n"));
    let n = 1usize << k;
    let srs = load_srs(&format!("insecure:{}", n.max(64)));
    let mut rng = Rng::from_seed([42; 32]);
    let scalars: Vec<Fr> = (0..n).map(|_| Fr::random(&mut rng)).collect();
    let time = |f: &mut dyn FnMut()| {
        let mut best = f64::MAX;
        for _ in 0..3 {
            let t0 = Instant::now();
            f();
            best = best.min(t0.elapsed().as_secs_f64());
        }
        best
    };
    // with --gpu the CUDA kernels directly, at every size
    let gpu = a.contains_key("gpu");
    let msm = |p: &[veil::curve::G1Affine], s: &[Fr]| {
        #[cfg(feature = "cuda")]
        if gpu {
            return veil::gpu::msm_g1(p, s);
        }
        veil::msm::msm_cpu(p, s)
    };
    let ntt = |v: &mut [Fr], w: Fr| {
        #[cfg(feature = "cuda")]
        if gpu {
            return veil::gpu::ntt(v, w);
        }
        veil::poly::ntt_cpu(v, w)
    };
    let mut msm_out = None;
    let t_msm = time(&mut || msm_out = Some(msm(&srs.g1[..n], &scalars)));
    let w = Fr::root_of_unity(k);
    let mut v = scalars.clone();
    let t_ntt = time(&mut || ntt(&mut v, w));
    // with --gpu, the CPU results as a check
    let check = if gpu {
        let cpu = veil::msm::msm_cpu(&srs.g1[..n], &scalars);
        let mut v2 = scalars.clone();
        veil::poly::ntt_cpu(&mut v2, w);
        let mut v3 = scalars.clone();
        ntt(&mut v3, w);
        format!("{}", cpu == msm_out.unwrap() && v2 == v3)
    } else {
        "null".into()
    };
    println!("{{\"log_n\": {k}, \"gpu\": {}, \"msm_s\": {t_msm:.5}, \"ntt_s\": {t_ntt:.5}, \"gpu_matches_cpu\": {check}}}", a.contains_key("gpu"));
}
