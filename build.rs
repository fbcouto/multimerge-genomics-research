use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();

    // 1. Verificação ARM (Mac M1/M2/M3 ou Windows Snapdragon)
    if target_arch == "aarch64" {
        println!("cargo:rustc-cfg=cpu_only");
        println!("cargo:warning=Arquitetura ARM detectada. Compilando no modo CPU-Only (Fallback).");
        return;
    }

    // 2. Avisar ao Cargo para re-avaliar o build se as variáveis do CUDA mudarem
    println!("cargo:rerun-if-env-changed=CUDA_PATH");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    println!("cargo:rerun-if-changed=src/merge_path_kernel.cu");

    // 3. Verificação do Compilador CUDA (nvcc)
    let nvcc_check = Command::new("nvcc").arg("--version").output();

    match nvcc_check {
      Ok(_) => {
            // NVCC encontrado! Vamos compilar o C++ para a GPU.
            cc::Build::new()
                .cuda(true)
                .compiler("nvcc")
                .file("src/merge_path_kernel.cu")
                .no_default_flags(true)
                .warnings(false) 
                .flag("-allow-unsupported-compiler") // <--- VOLTA ESTA LINHA
                .define("_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH", None) // <--- VOLTA ESTA LINHA
                .flag("-arch=sm_61") // Ajuste para a arquitetura Pascal (GTX 1050) ou superior
                .compile("merge_kernel");

            // 4. Descoberta Dinâmica do Caminho da Biblioteca (Cross-Platform)
            if let Some(path) = find_cuda_lib_path() {
                println!("cargo:rustc-link-search=native={}", path.display());
            } else {
                println!("cargo:warning=CUDA compilado, mas o diretório de bibliotecas não foi detectado dinamicamente. A linkagem pode falhar.");
            }

            println!("cargo:rustc-link-lib=static=cudart_static");
            println!("cargo:rustc-cfg=cuda_enabled");
        }
        Err(_) => {
            // Fallback para CPU
            println!("cargo:rustc-cfg=cpu_only");
            println!("cargo:warning=Compilador CUDA (nvcc) não encontrado. Compilando no modo CPU-Only (Fallback).");
        }
    }
}

/// Tenta descobrir onde as bibliotecas estáticas do CUDA estão instaladas
fn find_cuda_lib_path() -> Option<PathBuf> {
    // Procura primeiro nas variáveis de ambiente padrão da NVIDIA
    let env_vars = ["CUDA_PATH", "CUDA_HOME"];
    
    for var in env_vars {
        if let Ok(val) = env::var(var) {
            let mut path = PathBuf::from(val);
            
            // Adiciona a estrutura de pastas correta dependendo do Sistema Operacional
            if cfg!(windows) {
                path.push("lib");
                path.push("x64");
            } else {
                path.push("lib64");
            }
            
            if path.exists() {
                return Some(path);
            }
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