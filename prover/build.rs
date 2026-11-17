// With the `cuda` feature: compile cuda/veil.cu with nvcc into a static
// library and link it with the static CUDA runtime.
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=cuda/veil.cu");
    if std::env::var("CARGO_FEATURE_CUDA").is_err() {
        return;
    }
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let windows = std::env::var("CARGO_CFG_TARGET_OS").unwrap() == "windows";
    let lib = out.join(if windows { "veilcuda.lib" } else { "libveilcuda.a" });
    let arch = std::env::var("VEIL_CUDA_ARCH").unwrap_or_else(|_| "sm_75".into());
    let mut cmd = Command::new("nvcc");
    cmd.args(["-O3", "-std=c++17", "--lib", &format!("-arch={arch}"), "-o"]).arg(&lib).arg("cuda/veil.cu");
    if windows {
        cmd.args(["-Xcompiler", "/MD,/Zc:preprocessor"]);
        // nvcc wants cl.exe; outside a developer prompt, ask vswhere where it is
        if let Some(cl) = std::env::var("VEIL_CCBIN").ok().or_else(find_cl) {
            cmd.args(["-ccbin", &cl]);
        }
    } else {
        cmd.args(["-Xcompiler", "-fPIC"]);
    }
    let st = cmd.status().expect("nvcc not found (the cuda feature needs the CUDA toolkit)");
    assert!(st.success(), "nvcc failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=veilcuda");
    let cuda = std::env::var("CUDA_PATH").unwrap_or_else(|_| "/usr/local/cuda".into());
    let libdir = PathBuf::from(&cuda).join(if windows { "lib/x64" } else { "lib64" });
    println!("cargo:rustc-link-search=native={}", libdir.display());
    println!("cargo:rustc-link-lib=static=cudart_static");
    if !windows {
        for l in ["dl", "rt", "pthread", "stdc++"] {
            println!("cargo:rustc-link-lib={l}");
        }
    }
}

fn find_cl() -> Option<String> {
    let vswhere = r"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe";
    let out = Command::new(vswhere).args(["-latest", "-products", "*", "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64", "-property", "installationPath"]).output().ok()?;
    let root = PathBuf::from(String::from_utf8_lossy(&out.stdout).trim()).join(r"VC\Tools\MSVC");
    let mut versions: Vec<_> = std::fs::read_dir(&root).ok()?.filter_map(|e| e.ok()).map(|e| e.path()).collect();
    versions.sort();
    let bin = versions.last()?.join(r"bin\Hostx64\x64");
    bin.join("cl.exe").exists().then(|| bin.display().to_string())
}
