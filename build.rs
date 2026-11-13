use std::env;
use std::path::PathBuf;
use std::process::Command;

// Como o build decide entre GPU e CPU
// -----------------------------------
//   MULTIMERGE_CPU_ONLY=1       nao tenta o CUDA: a extensao roda em qualquer
//                               maquina (e e o que vai num wheel publico).
//   MULTIMERGE_CUDA_ARCH=sm_61  arquitetura alvo do nvcc. Sem ela vale o padrao
//                               do toolkit instalado: a mais antiga que ele
//                               ainda aceita, que roda em toda placa mais nova.
//                               (Fixar sm_61 quebrava no CUDA 13, que so
//                               compila de sm_75 em diante.)
//   MULTIMERGE_REQUIRE_CUDA=1   falha o build em vez de cair para CPU.
//
// Sem nvcc, ou se o nvcc falhar (falta o compilador C++ do sistema, toolkit que
// nao aceita a arquitetura pedida), o build segue em modo CPU e avisa. Do lado
// Python, `info_motor()["gpu"]` diz qual dos dois foi compilado.

fn ligado(nome: &str) -> bool {
    matches!(env::var(nome), Ok(v) if !v.is_empty() && v != "0")
}

fn modo_cpu(motivo: &str) {
    if ligado("MULTIMERGE_REQUIRE_CUDA") {
        panic!("MULTIMERGE_REQUIRE_CUDA=1, but CUDA cannot be used: {motivo}");
    }
    println!("cargo:rustc-cfg=cpu_only");
    println!("cargo:warning={motivo}. Building in CPU-only mode (fallback).");
}

fn main() {
    println!("cargo:rustc-check-cfg=cfg(cuda_enabled)");
    println!("cargo:rustc-check-cfg=cfg(cpu_only)");
    for var in ["CUDA_PATH", "CUDA_HOME", "MULTIMERGE_CPU_ONLY", "MULTIMERGE_CUDA_ARCH", "MULTIMERGE_REQUIRE_CUDA"] {
        println!("cargo:rerun-if-env-changed={var}");
    }
    println!("cargo:rerun-if-changed=src/merge_path_kernel.cu");
    println!("cargo:rerun-if-changed=build.rs");

    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();

    // 1. ARM (Mac M1/M2/M3 ou Windows Snapdragon) e macOS nao tem CUDA.
    if target_arch == "aarch64" || target_os == "macos" {
        modo_cpu("Platform without CUDA");
        return;
    }
    if ligado("MULTIMERGE_CPU_ONLY") {
        modo_cpu("MULTIMERGE_CPU_ONLY is set");
        return;
    }

    // 2. Procura o nvcc: no PATH e, depois, nos lugares de costume.
    let Some(nvcc) = find_nvcc() else {
        modo_cpu("CUDA compiler (nvcc) not found");
        return;
    };

    // 3. Compila o kernel. Uma falha aqui nao derruba o build.
    let mut build = cc::Build::new();
    build
        .cuda(true)
        .compiler(&nvcc)
        .file("src/merge_path_kernel.cu")
        .no_default_flags(true)
        .warnings(false)
        .cargo_warnings(false)
        .flag("-allow-unsupported-compiler")
        .define("_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH", None);
    if let Ok(arch) = env::var("MULTIMERGE_CUDA_ARCH") {
        if !arch.is_empty() {
            build.flag(&format!("-arch={arch}"));
        }
    }
    if target_os != "windows" {
        // `no_default_flags` tira o -fPIC, e o objeto vai dentro de uma
        // biblioteca dinamica (a extensao Python). O nvcc nao conhece -fPIC:
        // tem de ser repassado ao compilador do sistema.
        build.flag("-Xcompiler").flag("-fPIC");
    }

    if build.try_compile("merge_kernel").is_err() {
        // A mensagem do nvcc so aparece com `cargo build -vv`; aqui vai a
        // causa de quase sempre.
        let dica = if target_os == "windows" {
            "on Windows it needs cl.exe: use the Visual Studio \"x64 Native Tools\" prompt, or call vcvars64.bat first"
        } else {
            "check g++ and the toolkit version"
        };
        modo_cpu(&format!("nvcc ({}) failed to compile the CUDA kernel ({dica}; details with cargo build -vv)", nvcc.display()));
        return;
    }

    // 4. Descoberta Dinâmica do Caminho da Biblioteca (Cross-Platform)
    if let Some(path) = find_cuda_lib_path(&nvcc) {
        println!("cargo:rustc-link-search=native={}", path.display());
    } else {
        println!("cargo:warning=CUDA compiled, but the library directory was not detected. Linking may fail.");
    }

    println!("cargo:rustc-link-lib=static=cudart_static");
    if target_os == "linux" {
        // O cudart estatico pede estas no Linux.
        for lib in ["dl", "rt", "pthread"] {
            println!("cargo:rustc-link-lib=dylib={lib}");
        }
    }
    println!("cargo:rustc-cfg=cuda_enabled");
}

fn responde(nvcc: &PathBuf) -> bool {
    Command::new(nvcc).arg("--version").output().map(|o| o.status.success()).unwrap_or(false)
}

/// O nvcc do PATH ou, na falta, o de CUDA_PATH / CUDA_HOME / /usr/local/cuda.
fn find_nvcc() -> Option<PathBuf> {
    let direto = PathBuf::from("nvcc");
    if responde(&direto) {
        return Some(direto);
    }
    let exe = if cfg!(windows) { "nvcc.exe" } else { "nvcc" };
    let mut raizes: Vec<PathBuf> = ["CUDA_PATH", "CUDA_HOME"].iter().filter_map(|v| env::var(v).ok()).map(PathBuf::from).collect();
    if !cfg!(windows) {
        raizes.push(PathBuf::from("/usr/local/cuda"));
    }
    raizes.into_iter().map(|r| r.join("bin").join(exe)).find(|c| c.exists() && responde(c))
}

/// Tenta descobrir onde as bibliotecas estáticas do CUDA estão instaladas
fn find_cuda_lib_path(nvcc: &PathBuf) -> Option<PathBuf> {
    let sub: &[&str] = if cfg!(windows) { &["lib", "x64"] } else { &["lib64"] };
    let com_sub = |raiz: PathBuf| {
        let mut p = raiz;
        for s in sub {
            p.push(s);
        }
        p
    };

    // Procura primeiro nas variáveis de ambiente padrão da NVIDIA
    for var in ["CUDA_PATH", "CUDA_HOME"] {
        if let Ok(val) = env::var(var) {
            let path = com_sub(PathBuf::from(val));
            if path.exists() {
                return Some(path);
            }
        }
    }

    // Depois, ao lado do proprio nvcc (<raiz>/bin/nvcc -> <raiz>/lib64).
    if let Some(raiz) = nvcc.parent().and_then(|bin| bin.parent()) {
        let path = com_sub(raiz.to_path_buf());
        if path.exists() {
            return Some(path);
        }
    }

    // Fallback: Se as variáveis não existirem, tenta os caminhos padrão absolutos
    if !cfg!(windows) {
        let default_linux = PathBuf::from("/usr/local/cuda/lib64");
        if default_linux.exists() {
            return Some(default_linux);
        }
    }

    None
}
