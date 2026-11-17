#!/bin/sh
# The Sudoku end to end: compile the circuit (Haskell), prove and verify it
# with the Ethereum ceremony's SRS (Rust), write the Solidity verifier, and
# the files of the browser demo (the prover in Wasm, the SRS, the key).
set -e
cd "$(dirname "$0")/.."
mkdir -p build lang/build web/public
(cd lang && ghc -O1 -isrc -iapp -outputdir build -o build/veilc app/Main.hs >/dev/null)
lang/build/veilc lang/examples/sudoku.veil -o build/sudoku.circuit
(cd prover && cargo build --release)
V=prover/target/release/veil
SRS=ceremony:data/trusted_setup.txt
$V prove --circuit build/sudoku.circuit --inputs lang/examples/sudoku.inputs --srs $SRS \
  --out build/proof.bin --public build/public.txt --vk build/vk.bin --evm build/proof.evm
$V verify --vk build/vk.bin --proof build/proof.bin --public build/public.txt
$V export-sol --vk build/vk.bin --out build/Verifier.sol
cp build/sudoku.circuit build/vk.bin web/public/
$V export-srs --srs $SRS --count 2054 --out web/public/srs.bin   # n + 6 powers for n = 2048
if rustup target list --installed | grep -q wasm32-unknown-unknown; then
  (cd prover && cargo rustc --lib --release --target wasm32-unknown-unknown --crate-type cdylib -- -C debuginfo=0 -C strip=debuginfo)
  cp prover/target/wasm32-unknown-unknown/release/veil.wasm web/public/
fi
