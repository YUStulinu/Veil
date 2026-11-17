//! Official test vectors: ethereum/bls12-381-tests (point deserialization)
//! and the EIP-4844 KZG vectors from c-kzg-4844, against the Ethereum KZG
//! ceremony's trusted setup. `scripts/fetch_vectors.py` downloads them into
//! data/; without them the tests are skipped, unless VEIL_REQUIRE_VECTORS=1.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use veil::curve::{g1_from_compressed, g2_from_compressed};
use veil::eip4844 as k;
use veil::hash::unhex;
use veil::kzg::Srs;

fn data_dir() -> Option<PathBuf> {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../data");
    if d.join("trusted_setup.txt").exists() {
        Some(d)
    } else {
        assert!(std::env::var("VEIL_REQUIRE_VECTORS").is_err(), "test vectors missing: run scripts/fetch_vectors.py");
        eprintln!("skipped: no test vectors in data/ (scripts/fetch_vectors.py)");
        None
    }
}

fn srs() -> &'static Srs {
    static S: OnceLock<Srs> = OnceLock::new();
    S.get_or_init(|| Srs::from_ceremony_file(&fs::read_to_string(data_dir().unwrap().join("trusted_setup.txt")).unwrap()).unwrap())
}

/// The vectors' YAML is flat: `key: 'value'`, `key: []`, `key:` + `- 'item'` lines,
/// and `output:` the same way. Values keep their quotes stripped; lists become vectors.
fn parse_yaml(text: &str) -> BTreeMap<String, Vec<String>> {
    let mut m: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut cur = String::new();
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t == "input:" {
            continue;
        }
        if let Some(item) = t.strip_prefix("- ") {
            m.get_mut(&cur).unwrap().push(item.trim_matches('\'').to_string());
        } else if let Some((key, v)) = t.split_once(':') {
            cur = key.to_string();
            let v = v.trim();
            let vals = if v.is_empty() || v == "[]" { vec![] } else { vec![v.trim_matches('\'').to_string()] };
            m.insert(cur.clone(), vals);
        }
    }
    m
}

fn cases(name: &str) -> Vec<(String, BTreeMap<String, Vec<String>>)> {
    let dir = data_dir().unwrap().join("ckzg/tests").join(name).join("kzg-mainnet");
    let mut out = vec![];
    for e in fs::read_dir(dir).unwrap() {
        let p = e.unwrap().path();
        out.push((p.file_name().unwrap().to_string_lossy().into(), parse_yaml(&fs::read_to_string(p.join("data.yaml")).unwrap())));
    }
    assert!(!out.is_empty());
    out
}

fn b(m: &BTreeMap<String, Vec<String>>, key: &str) -> Vec<u8> {
    unhex(&m[key][0]).unwrap_or_default()
}

/// compares a result with the expected output ("null" = invalid input)
fn check<T: std::fmt::Debug + PartialEq>(name: &str, got: Result<T, k::InvalidInput>, want: Option<T>) -> bool {
    let ok = match (&got, &want) {
        (Err(_), None) => true,
        (Ok(g), Some(w)) => g == w,
        _ => false,
    };
    if !ok {
        eprintln!("{name}: got {got:?}, expected {want:?}");
    }
    ok
}

fn out_bool(m: &BTreeMap<String, Vec<String>>) -> Option<bool> {
    m["output"].first().filter(|v| *v != "null").map(|v| v == "true")
}

fn out_bytes(m: &BTreeMap<String, Vec<String>>) -> Option<Vec<u8>> {
    m["output"].first().filter(|v| *v != "null").map(|v| unhex(v).unwrap())
}

#[test]
fn bls_deserialization_g1_g2() {
    let Some(d) = data_dir() else { return };
    let mut n = 0;
    for (dir, g2) in [("deserialization_G1", false), ("deserialization_G2", true)] {
        for e in fs::read_dir(d.join("bls").join(dir)).unwrap() {
            let p = e.unwrap().path();
            let s = fs::read_to_string(&p).unwrap();
            let hexstr = s.split('"').nth(5).unwrap();
            let want = s.contains("\"output\": true");
            let bytes = unhex(hexstr).unwrap_or_default();
            let got = if g2 { g2_from_compressed(&bytes).is_ok() } else { g1_from_compressed(&bytes).is_ok() };
            assert_eq!(got, want, "{}", p.display());
            n += 1;
        }
    }
    eprintln!("deserialization: {n} vectors");
}

#[test]
fn kzg_blob_to_commitment() {
    if data_dir().is_none() {
        return;
    }
    let all = cases("blob_to_kzg_commitment");
    let bad: Vec<_> = all.iter().filter(|(n, m)| !check(n, k::blob_to_kzg_commitment(srs(), &b(m, "blob")).map(|c| c.to_vec()), out_bytes(m))).collect();
    assert!(bad.is_empty(), "{} of {} failed", bad.len(), all.len());
}

#[test]
fn kzg_compute_proof() {
    if data_dir().is_none() {
        return;
    }
    let all = cases("compute_kzg_proof");
    let mut bad = 0;
    for (n, m) in &all {
        let got = k::compute_kzg_proof(srs(), &b(m, "blob"), &b(m, "z")).map(|(p, y)| vec![p.to_vec(), y.to_vec()]);
        let want = if m["output"].len() == 2 { Some(m["output"].iter().map(|v| unhex(v).unwrap()).collect()) } else { None };
        bad += !check(n, got, want) as usize;
    }
    assert_eq!(bad, 0, "of {}", all.len());
}

#[test]
fn kzg_verify_proof() {
    if data_dir().is_none() {
        return;
    }
    let all = cases("verify_kzg_proof");
    let bad = all.iter().filter(|(n, m)| !check(n, k::verify_kzg_proof(srs(), &b(m, "commitment"), &b(m, "z"), &b(m, "y"), &b(m, "proof")), out_bool(m))).count();
    assert_eq!(bad, 0, "of {}", all.len());
}

#[test]
fn kzg_blob_proofs() {
    if data_dir().is_none() {
        return;
    }
    let mut total = 0;
    let mut bad = 0;
    for (n, m) in cases("compute_blob_kzg_proof") {
        bad += !check(&n, k::compute_blob_kzg_proof(srs(), &b(&m, "blob"), &b(&m, "commitment")).map(|p| p.to_vec()), out_bytes(&m)) as usize;
        total += 1;
    }
    for (n, m) in cases("verify_blob_kzg_proof") {
        bad += !check(&n, k::verify_blob_kzg_proof(srs(), &b(&m, "blob"), &b(&m, "commitment"), &b(&m, "proof")), out_bool(&m)) as usize;
        total += 1;
    }
    for (n, m) in cases("verify_blob_kzg_proof_batch") {
        let (bl, cs, ps) = (&m["blobs"], &m["commitments"], &m["proofs"]);
        let got = if bl.len() != cs.len() || cs.len() != ps.len() {
            Err(k::InvalidInput)
        } else {
            // each one alone: the same answer as the random linear combination
            let mut all_ok = Ok(true);
            for i in 0..bl.len() {
                let r = k::verify_blob_kzg_proof(srs(), &unhex(&bl[i]).unwrap_or_default(), &unhex(&cs[i]).unwrap_or_default(), &unhex(&ps[i]).unwrap_or_default());
                match r {
                    Err(e) => {
                        all_ok = Err(e);
                        break;
                    }
                    Ok(false) => all_ok = all_ok.map(|_| false),
                    Ok(true) => {}
                }
            }
            all_ok
        };
        bad += !check(&n, got, out_bool(&m)) as usize;
        total += 1;
    }
    for (n, m) in cases("compute_challenge") {
        let want = out_bytes(&m);
        let got = Ok(k::compute_challenge(&b(&m, "blob"), &b(&m, "commitment")).to_bytes_be().to_vec());
        bad += !check(&n, got, want) as usize;
        total += 1;
    }
    assert_eq!(bad, 0, "of {total}");
    eprintln!("blob proofs: {total} vectors");
}
