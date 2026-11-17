//! A PLONK constraint system: gates q_M a b + q_L a + q_R b + q_O c + q_C = 0
//! over variables, public inputs, and the witness program that computes
//! every variable from the inputs. The Haskell compiler (`veilc`) writes it
//! as text; tests and benchmarks build it in Rust.
//!
//! Text format:
//!   veil-circuit 1
//!   name <name>
//!   vars <count>                  variable 0 is always 0
//!   input <name> <public|private> <length>      (repeated)
//!   public <count> <var>...       the public variables, in order
//!   hint <out> input <k>          k indexes the inputs flattened in order
//!   hint <out> eval <qM> <qL> <qR> <qC> <a> <b>  out = qM a b + qL a + qR b + qC
//!   hint <out> inv <a>            out = 1/a (0 if a = 0)
//!   gate <qM> <qL> <qR> <qO> <qC> <a> <b> <c>
//! Field constants are decimal, may be negative, and are reduced mod r.

use crate::ff::{Field, Fr};

pub type Var = u32;

#[derive(Clone, Debug, PartialEq)]
pub struct Gate {
    pub q_m: Fr,
    pub q_l: Fr,
    pub q_r: Fr,
    pub q_o: Fr,
    pub q_c: Fr,
    pub a: Var,
    pub b: Var,
    pub c: Var,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Hint {
    Input(usize),
    Eval { q_m: Fr, q_l: Fr, q_r: Fr, q_c: Fr, a: Var, b: Var },
    Inv(Var),
}

#[derive(Clone, Debug)]
pub struct InputDecl {
    pub name: String,
    pub public: bool,
    pub len: usize,
}

#[derive(Clone, Debug, Default)]
pub struct Circuit {
    pub name: String,
    pub vars: usize,
    pub inputs: Vec<InputDecl>,
    pub public: Vec<Var>,
    pub hints: Vec<(Var, Hint)>,
    pub gates: Vec<Gate>,
}

#[derive(Debug)]
pub enum WitnessError {
    Inputs(String),
    /// the inputs do not satisfy the circuit: the first failing gate
    Unsatisfied(usize),
}

impl Circuit {
    pub fn new(name: &str) -> Circuit {
        Circuit { name: name.into(), vars: 1, ..Default::default() }
    }

    /// rows of the PLONK table: one per public input, then the gates
    pub fn rows(&self) -> usize {
        self.public.len() + self.gates.len()
    }

    pub fn total_inputs(&self) -> usize {
        self.inputs.iter().map(|i| i.len).sum()
    }

    // ---- building (tests, benchmarks) ----

    pub fn new_var(&mut self, h: Hint) -> Var {
        let v = self.vars as Var;
        self.vars += 1;
        self.hints.push((v, h));
        v
    }

    pub fn input(&mut self, name: &str, public: bool, len: usize) -> Vec<Var> {
        let base = self.total_inputs();
        self.inputs.push(InputDecl { name: name.into(), public, len });
        let vs: Vec<Var> = (0..len).map(|k| self.new_var(Hint::Input(base + k))).collect();
        if public {
            self.public.extend(&vs);
        }
        vs
    }

    /// out = qM a b + qL a + qR b + qC, with its gate
    pub fn eval(&mut self, q_m: Fr, q_l: Fr, q_r: Fr, q_c: Fr, a: Var, b: Var) -> Var {
        let out = self.new_var(Hint::Eval { q_m, q_l, q_r, q_c, a, b });
        self.gates.push(Gate { q_m, q_l, q_r, q_o: -Fr::ONE, q_c, a, b, c: out });
        out
    }

    pub fn mul(&mut self, a: Var, b: Var) -> Var {
        self.eval(Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ZERO, a, b)
    }

    pub fn add(&mut self, a: Var, b: Var) -> Var {
        self.eval(Fr::ZERO, Fr::ONE, Fr::ONE, Fr::ZERO, a, b)
    }

    pub fn assert_gate(&mut self, q_m: Fr, q_l: Fr, q_r: Fr, q_o: Fr, q_c: Fr, a: Var, b: Var, c: Var) {
        self.gates.push(Gate { q_m, q_l, q_r, q_o, q_c, a, b, c });
    }

    pub fn assert_equal(&mut self, a: Var, b: Var) {
        self.assert_gate(Fr::ZERO, Fr::ONE, -Fr::ONE, Fr::ZERO, Fr::ZERO, a, b, 0);
    }

    // ---- witness ----

    /// all variable values from the flattened inputs; checks every gate
    pub fn witness(&self, inputs: &[Fr]) -> Result<Vec<Fr>, WitnessError> {
        let w = self.witness_unchecked(inputs)?;
        match self.first_failing_gate(&w) {
            Some(i) => Err(WitnessError::Unsatisfied(i)),
            None => Ok(w),
        }
    }

    pub fn witness_unchecked(&self, inputs: &[Fr]) -> Result<Vec<Fr>, WitnessError> {
        if inputs.len() != self.total_inputs() {
            return Err(WitnessError::Inputs(format!("{} inputs given, the circuit takes {}", inputs.len(), self.total_inputs())));
        }
        let mut w = vec![Fr::ZERO; self.vars];
        for (out, h) in &self.hints {
            w[*out as usize] = match h {
                Hint::Input(k) => inputs[*k],
                Hint::Eval { q_m, q_l, q_r, q_c, a, b } => {
                    let (x, y) = (w[*a as usize], w[*b as usize]);
                    *q_m * x * y + *q_l * x + *q_r * y + *q_c
                }
                Hint::Inv(a) => w[*a as usize].inverse().unwrap_or(Fr::ZERO),
            };
        }
        Ok(w)
    }

    pub fn first_failing_gate(&self, w: &[Fr]) -> Option<usize> {
        self.gates.iter().position(|g| {
            let (a, b, c) = (w[g.a as usize], w[g.b as usize], w[g.c as usize]);
            !(g.q_m * a * b + g.q_l * a + g.q_r * b + g.q_o * c + g.q_c).is_zero()
        })
    }

    /// the public values, in order
    pub fn public_values(&self, w: &[Fr]) -> Vec<Fr> {
        self.public.iter().map(|v| w[*v as usize]).collect()
    }

    // ---- text format ----

    pub fn parse(text: &str) -> Result<Circuit, String> {
        let mut c = Circuit::default();
        let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty() && !l.trim_start().starts_with('#'));
        match lines.next() {
            Some((_, l)) if l.trim() == "veil-circuit 1" => {}
            _ => return Err("not a veil-circuit 1 file".into()),
        }
        for (no, line) in lines {
            let err = |m: &str| format!("line {}: {m}: {line}", no + 1);
            let t: Vec<&str> = line.split_whitespace().collect();
            let fr = |s: &str| Fr::from_str_mod(s).ok_or_else(|| err("bad field constant"));
            let var = |s: &str| -> Result<Var, String> {
                let v: Var = s.parse().map_err(|_| err("bad variable"))?;
                Ok(v)
            };
            match (t[0], t.len()) {
                ("name", 2) => c.name = t[1].into(),
                ("vars", 2) => c.vars = t[1].parse().map_err(|_| err("bad count"))?,
                ("input", 4) => c.inputs.push(InputDecl {
                    name: t[1].into(),
                    public: match t[2] {
                        "public" => true,
                        "private" => false,
                        _ => return Err(err("public or private")),
                    },
                    len: t[3].parse().map_err(|_| err("bad length"))?,
                }),
                ("public", _) => {
                    let n: usize = t.get(1).and_then(|s| s.parse().ok()).ok_or_else(|| err("bad count"))?;
                    if t.len() != n + 2 {
                        return Err(err("public count mismatch"));
                    }
                    c.public = t[2..].iter().map(|s| var(s)).collect::<Result<_, _>>()?;
                }
                ("hint", 4) if t[2] == "input" => c.hints.push((var(t[1])?, Hint::Input(t[3].parse().map_err(|_| err("bad index"))?))),
                ("hint", 4) if t[2] == "inv" => c.hints.push((var(t[1])?, Hint::Inv(var(t[3])?))),
                ("hint", 9) if t[2] == "eval" => c.hints.push((
                    var(t[1])?,
                    Hint::Eval { q_m: fr(t[3])?, q_l: fr(t[4])?, q_r: fr(t[5])?, q_c: fr(t[6])?, a: var(t[7])?, b: var(t[8])? },
                )),
                ("gate", 9) => c.gates.push(Gate {
                    q_m: fr(t[1])?,
                    q_l: fr(t[2])?,
                    q_r: fr(t[3])?,
                    q_o: fr(t[4])?,
                    q_c: fr(t[5])?,
                    a: var(t[6])?,
                    b: var(t[7])?,
                    c: var(t[8])?,
                }),
                _ => return Err(err("unknown line")),
            }
        }
        // every variable must be in range, and defined before use
        let n = c.vars as Var;
        let mut defined = vec![false; c.vars];
        if c.vars == 0 {
            return Err("no variables".into());
        }
        defined[0] = true;
        for (out, h) in &c.hints {
            let uses: &[Var] = match h {
                Hint::Input(k) => {
                    if *k >= c.total_inputs() {
                        return Err(format!("hint for {out}: input {k} out of range"));
                    }
                    &[]
                }
                Hint::Eval { a, b, .. } => &[*a, *b],
                Hint::Inv(a) => std::slice::from_ref(a),
            };
            for u in uses {
                if *u >= n || !defined[*u as usize] {
                    return Err(format!("hint for {out} uses {u} before it is defined"));
                }
            }
            if *out >= n || *out == 0 || defined[*out as usize] {
                return Err(format!("variable {out} defined twice or out of range"));
            }
            defined[*out as usize] = true;
        }
        if let Some(v) = defined.iter().position(|d| !d) {
            return Err(format!("variable {v} is never defined"));
        }
        for g in &c.gates {
            if g.a >= n || g.b >= n || g.c >= n {
                return Err("gate variable out of range".into());
            }
        }
        if c.public.iter().any(|v| *v >= n) {
            return Err("public variable out of range".into());
        }
        Ok(c)
    }

    pub fn to_text(&self) -> String {
        let s = |x: &Fr| format!("{}", to_decimal(x));
        let mut out = format!("veil-circuit 1\nname {}\nvars {}\n", self.name, self.vars);
        for i in &self.inputs {
            out += &format!("input {} {} {}\n", i.name, if i.public { "public" } else { "private" }, i.len);
        }
        out += &format!("public {}", self.public.len());
        for v in &self.public {
            out += &format!(" {v}");
        }
        out += "\n";
        for (o, h) in &self.hints {
            out += &match h {
                Hint::Input(k) => format!("hint {o} input {k}\n"),
                Hint::Inv(a) => format!("hint {o} inv {a}\n"),
                Hint::Eval { q_m, q_l, q_r, q_c, a, b } => format!("hint {o} eval {} {} {} {} {a} {b}\n", s(q_m), s(q_l), s(q_r), s(q_c)),
            };
        }
        for g in &self.gates {
            out += &format!("gate {} {} {} {} {} {} {} {}\n", s(&g.q_m), s(&g.q_l), s(&g.q_r), s(&g.q_o), s(&g.q_c), g.a, g.b, g.c);
        }
        out
    }
}

/// A benchmark circuit with `gates` gates: a chain of squarings and
/// additions from a private seed, ending in a public output.
/// Returns the circuit and its inputs [x (public), out (public), y (private)].
pub fn synthetic(gates: usize) -> (Circuit, Vec<Fr>) {
    assert!(gates >= 2);
    let mut c = Circuit::new(&format!("synthetic-{gates}"));
    let pubs = c.input("pub", true, 2);
    let y = c.input("seed", false, 1)[0];
    let (x, out) = (pubs[0], pubs[1]);
    let mut acc = c.add(x, y);
    for i in 0..gates - 2 {
        acc = if i % 2 == 0 { c.mul(acc, acc) } else { c.eval(Fr::ZERO, Fr::ONE, Fr::from_u64(3), Fr::from_u64(i as u64), acc, y) };
    }
    c.assert_equal(acc, out);
    let (xv, yv) = (Fr::from_u64(3), Fr::from_u64(0x5eed));
    let w = c.witness_unchecked(&[xv, Fr::ZERO, yv]).unwrap();
    let outv = w[acc as usize];
    (c, vec![xv, outv, yv])
}

/// decimal, small negatives as "-k" for readability
pub fn to_decimal(x: &Fr) -> String {
    let neg = -*x;
    if neg.to_limbs()[1..].iter().all(|l| *l == 0) && !x.to_limbs()[1..].iter().all(|l| *l == 0) {
        format!("-{}", to_decimal_plain(&neg))
    } else {
        to_decimal_plain(x)
    }
}

/// a public input: decimal digits of an integer below r (no reduction, so
/// one statement has one encoding)
pub fn parse_canonical(s: &str) -> Option<Fr> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v = Fr::from_str_mod(s)?;
    (to_decimal_plain(&v) == s.trim_start_matches('0').to_string() || (v == Fr::ZERO && s.bytes().all(|b| b == b'0'))).then_some(v)
}

pub fn to_decimal_plain(v: &Fr) -> String {
    let mut l = v.to_limbs();
    let mut digits = vec![];
    loop {
        let mut rem: u128 = 0;
        for i in (0..4).rev() {
            let cur = (rem << 64) | l[i] as u128;
            l[i] = (cur / 10) as u64;
            rem = cur % 10;
        }
        digits.push(b'0' + rem as u8);
        if l.iter().all(|x| *x == 0) {
            break;
        }
    }
    digits.reverse();
    String::from_utf8(digits).unwrap()
}

/// inputs file: one line per input array, "name: v1 v2 ...", in declaration order
pub fn parse_inputs(c: &Circuit, text: &str) -> Result<Vec<Fr>, String> {
    let mut by_name = std::collections::HashMap::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let (name, vals) = line.split_once(':').ok_or("expected 'name: values'")?;
        let vals: Vec<Fr> = vals.split_whitespace().map(|v| Fr::from_str_mod(v).ok_or(format!("bad value {v}"))).collect::<Result<_, _>>()?;
        by_name.insert(name.trim().to_string(), vals);
    }
    let mut out = vec![];
    for d in &c.inputs {
        let v = by_name.get(&d.name).ok_or(format!("missing input {}", d.name))?;
        if v.len() != d.len {
            return Err(format!("input {} has {} values, expected {}", d.name, v.len(), d.len));
        }
        out.extend(v);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn text_round_trip_and_witness() {
        let mut c = Circuit::new("t");
        let x = c.input("x", true, 1)[0];
        let y = c.input("y", false, 1)[0];
        let xy = c.mul(x, y);
        let inv = c.new_var(Hint::Inv(y));
        c.assert_gate(Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ZERO, -Fr::ONE, y, inv, 0);
        c.assert_gate(Fr::ZERO, Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::from_i64(-12), xy, 0, 0);
        let c2 = Circuit::parse(&c.to_text()).unwrap();
        assert_eq!(c2.to_text(), c.to_text());
        assert!(c2.witness(&[Fr::from_u64(3), Fr::from_u64(4)]).is_ok());
        assert!(matches!(c2.witness(&[Fr::from_u64(3), Fr::from_u64(5)]), Err(WitnessError::Unsatisfied(2))));
        assert!(matches!(c2.witness(&[Fr::from_u64(0), Fr::from_u64(0)]), Err(WitnessError::Unsatisfied(1))));
        assert_eq!(to_decimal(&Fr::from_i64(-7)), "-7");
        assert_eq!(to_decimal(&Fr::from_u64(1234567)), "1234567");
    }

    #[test]
    fn rejects_bad_files() {
        assert!(Circuit::parse("veil-circuit 1\nvars 2\nhint 1 inv 1\n").is_err());
        assert!(Circuit::parse("veil-circuit 1\nvars 3\nhint 1 inv 0\n").is_err());
        assert!(Circuit::parse("veil-circuit 1\nvars 1\ngate 1 0 0 0 0 0 0 5\n").is_err());
    }
}
