//! The browser's entry points (wasm32-unknown-unknown, no bindings crate):
//! the page allocates buffers, writes the circuit, the SRS and the inputs,
//! and calls `veil_prove`. The result buffer is
//!   status (u32: 0 = ok, 1 = bad input, 2 = unsatisfied) | length (u32) | bytes
//! where bytes = proof (624) | verifying key | public inputs (32 bytes each),
//! or an error message.

use crate::kzg::Srs;
use crate::plonk::circuit::{parse_inputs, Circuit};
use crate::plonk::{preprocess, prove};
use crate::rng::Rng;
use std::cell::RefCell;

thread_local! {
    static SEED: RefCell<Option<[u8; 32]>> = const { RefCell::new(None) };
}

#[link(wasm_import_module = "env")]
extern "C" {
    /// milliseconds, from the page (performance.now)
    fn veil_now_ms() -> f64;
}

pub fn now_ms() -> f64 {
    unsafe { veil_now_ms() }
}

/// the 32 bytes of entropy the page wrote with crypto.getRandomValues
pub(crate) fn take_seed() -> [u8; 32] {
    SEED.with(|s| s.borrow_mut().take()).expect("call veil_seed with fresh randomness before proving")
}

#[no_mangle]
pub extern "C" fn veil_alloc(len: usize) -> *mut u8 {
    let mut v = vec![0u8; len];
    let p = v.as_mut_ptr();
    std::mem::forget(v);
    p
}

/// # Safety
/// `p` must come from `veil_alloc(len)`.
#[no_mangle]
pub unsafe extern "C" fn veil_free(p: *mut u8, len: usize) {
    drop(Vec::from_raw_parts(p, len, len));
}

/// # Safety
/// `p` points to 32 readable bytes.
#[no_mangle]
pub unsafe extern "C" fn veil_seed(p: *const u8) {
    let mut s = [0u8; 32];
    s.copy_from_slice(std::slice::from_raw_parts(p, 32));
    SEED.with(|x| *x.borrow_mut() = Some(s));
}

fn result(status: u32, body: &[u8]) -> *mut u8 {
    let mut out = Vec::with_capacity(8 + body.len());
    out.extend_from_slice(&status.to_le_bytes());
    out.extend_from_slice(&(body.len() as u32).to_le_bytes());
    out.extend_from_slice(body);
    let p = veil_alloc(out.len());
    unsafe { std::ptr::copy_nonoverlapping(out.as_ptr(), p, out.len()) };
    p
}

/// # Safety
/// The pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn veil_prove(circuit: *const u8, circuit_len: usize, srs: *const u8, srs_len: usize, inputs: *const u8, inputs_len: usize) -> *mut u8 {
    let text = |p: *const u8, n: usize| String::from_utf8_lossy(std::slice::from_raw_parts(p, n)).into_owned();
    let c = match Circuit::parse(&text(circuit, circuit_len)) {
        Ok(c) => c,
        Err(e) => return result(1, e.as_bytes()),
    };
    let srs = match Srs::from_bytes(std::slice::from_raw_parts(srs, srs_len), "file") {
        Ok(s) => s,
        Err(e) => return result(1, e.as_bytes()),
    };
    let ins = match parse_inputs(&c, &text(inputs, inputs_len)) {
        Ok(v) => v,
        Err(e) => return result(1, e.as_bytes()),
    };
    let w = match c.witness(&ins) {
        Ok(w) => w,
        Err(e) => return result(2, format!("{e:?}").as_bytes()),
    };
    let pk = match preprocess(&c, &srs) {
        Ok(pk) => pk,
        Err(e) => return result(1, e.as_bytes()),
    };
    let public = c.public_values(&w);
    let (proof, _) = prove(&pk, &public, &w, &mut Rng::from_os());
    let mut body = proof.to_bytes();
    body.extend_from_slice(&pk.vk.to_bytes());
    for x in &public {
        body.extend_from_slice(&x.to_bytes_be());
    }
    result(0, &body)
}
