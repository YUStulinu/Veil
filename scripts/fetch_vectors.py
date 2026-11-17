"""Downloads the official test data into data/ (data only, no code is run):

- ethereum/bls12-381-tests v0.1.2 (JSON): point deserialization vectors
- c-kzg-4844's EIP-4844 KZG vectors (a sparse checkout of tests/ at a pinned commit)
- the Ethereum KZG ceremony's trusted setup (trusted_setup.txt, same commit)

    python scripts/fetch_vectors.py
"""

import hashlib
import os
import subprocess
import tarfile
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
CKZG = "a751e2d92abcd63b01d1400789789f5ea289a8b0"
FILES = {
    "bls.tgz": ("https://github.com/ethereum/bls12-381-tests/releases/download/v0.1.2/bls_tests_json.tar.gz",
                "0dc1dfabfa7b3e0f8e8029f9cb2bb8733597526fa085b18f160d2c4eb0b448bc"),
    "trusted_setup.txt": (f"https://raw.githubusercontent.com/ethereum/c-kzg-4844/{CKZG}/src/trusted_setup.txt",
                          "d39b9f2d047cc9dca2de58f264b6a09448ccd34db967881a6713eacacf0f26b7"),
}
KZG_TESTS = ["blob_to_kzg_commitment", "compute_kzg_proof", "verify_kzg_proof", "compute_blob_kzg_proof",
             "verify_blob_kzg_proof", "verify_blob_kzg_proof_batch", "compute_challenge"]


def main():
    os.makedirs(DATA, exist_ok=True)
    for name, (url, sha) in FILES.items():
        path = os.path.join(DATA, name)
        if not os.path.exists(path):
            print("downloading", url)
            urllib.request.urlretrieve(url, path)
        got = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if got != sha:
            raise SystemExit(f"{name}: sha256 {got}, expected {sha}")
    bls = os.path.join(DATA, "bls")
    if not os.path.isdir(bls):
        with tarfile.open(os.path.join(DATA, "bls.tgz")) as t:
            t.extractall(bls, filter="data")
    ckzg = os.path.join(DATA, "ckzg")
    if not os.path.isdir(ckzg):
        git = lambda *a, **k: subprocess.run(["git", *a], check=True, **k)
        git("clone", "-q", "--filter=blob:none", "--sparse", "--no-checkout", "https://github.com/ethereum/c-kzg-4844.git", ckzg)
        git("-C", ckzg, "sparse-checkout", "set", *[f"tests/{t}" for t in KZG_TESTS])
        git("-C", ckzg, "checkout", "-q", CKZG)
    print("test vectors in", DATA)


if __name__ == "__main__":
    main()
