//! Ponte para o kernel CUDA (`merge_path_kernel.cu`) -- e o que responder
//! quando ele nao esta la.
//!
//! As tres operacoes devolvem `Err(motivo)` quando a GPU nao pode ser usada:
//! extensao compilada sem CUDA, maquina sem placa ou sem driver, dado maior
//! que a VRAM. Quem chama decide se cai para a CPU (`lib.rs`, `Dispositivo`).
//! Nada aqui derruba o processo nem devolve lixo.

use std::sync::OnceLock;

/// A extensao foi compilada com o kernel CUDA?
pub const COMPILADA: bool = cfg!(cuda_enabled);

#[cfg(cuda_enabled)]
mod ffi {
    extern "C" {
        pub fn mm_cuda_dispositivos() -> i32;
        pub fn mm_cuda_memoria(livre: *mut usize, total: *mut usize) -> i32;
        pub fn launch_merge_path_kernel(a: *const u64, b: *const u64, out: *mut u64, len_a: usize, len_b: usize) -> i32;
        pub fn launch_bloom_filter_build(kmers: *const u64, num_kmers: usize, bloom: *mut u8, bloom_bytes: usize, num_hashes: i32) -> i32;
        pub fn launch_bloom_filter_query(
            query_kmers: *const u64,
            num_queries: usize,
            bloom: *const u8,
            bloom_bytes: usize,
            num_hashes: i32,
            out_mask: *mut u8,
        ) -> i32;
    }
}

/// O codigo devolvido pelo lancador, em palavras.
#[cfg(cuda_enabled)]
fn motivo(codigo: i32) -> String {
    match codigo {
        -1 => "the input does not fit the kernel's indices".to_string(),
        2 => "not enough GPU memory for this data".to_string(),
        35 => "the NVIDIA driver is older than required".to_string(),
        100 | 101 => "no usable CUDA device".to_string(),
        c => format!("CUDA error {c}"),
    }
}

#[cfg(cuda_enabled)]
fn confere(codigo: i32) -> Result<(), String> {
    if codigo == 0 {
        Ok(())
    } else {
        Err(motivo(codigo))
    }
}

const SEM_CUDA: &str = "extension built without CUDA";
const SEM_PLACA: &str = "no usable CUDA device on this machine";

/// Placas CUDA que este processo enxerga. 0 sem CUDA, sem placa ou sem driver.
pub fn dispositivos() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| {
        #[cfg(cuda_enabled)]
        {
            // SAFETY: funcao sem argumentos; so consulta o runtime do CUDA.
            (unsafe { ffi::mm_cuda_dispositivos() }).max(0) as usize
        }
        #[cfg(not(cuda_enabled))]
        {
            0
        }
    })
}

/// Ha GPU utilizavel agora?
pub fn disponivel() -> bool {
    dispositivos() > 0
}

fn exige_placa() -> Result<(), String> {
    if !COMPILADA {
        Err(SEM_CUDA.to_string())
    } else if !disponivel() {
        Err(SEM_PLACA.to_string())
    } else {
        Ok(())
    }
}

/// `(livre, total)` da placa em uso, em bytes.
pub fn memoria() -> Option<(usize, usize)> {
    if !disponivel() {
        return None;
    }
    #[cfg(cuda_enabled)]
    {
        let (mut livre, mut total) = (0usize, 0usize);
        // SAFETY: os dois ponteiros apontam para variaveis locais vivas.
        if unsafe { ffi::mm_cuda_memoria(&mut livre, &mut total) } == 0 {
            return Some((livre, total));
        }
    }
    None
}

/// Funde dois blocos ja ordenados.
pub fn mesclar(a: &[u64], b: &[u64]) -> Result<Vec<u64>, String> {
    exige_placa()?;
    #[cfg(cuda_enabled)]
    {
        let mut saida = vec![0u64; a.len() + b.len()];
        // SAFETY: `a` e `b` sao lidos nos tamanhos informados e `saida` tem
        // exatamente a soma deles, que e o que o lancador escreve.
        confere(unsafe { ffi::launch_merge_path_kernel(a.as_ptr(), b.as_ptr(), saida.as_mut_ptr(), a.len(), b.len()) })?;
        Ok(saida)
    }
    #[cfg(not(cuda_enabled))]
    {
        let _ = (a, b);
        Err(SEM_CUDA.to_string())
    }
}

/// Constroi o filtro de Bloom de `bytes` bytes. Mesmo resultado, bit a bit,
/// de `bloom::construir`.
pub fn bloom_construir(kmers: &[u64], bytes: usize, num_hashes: i32) -> Result<Vec<u8>, String> {
    exige_placa()?;
    #[cfg(cuda_enabled)]
    {
        let mut filtro = vec![0u8; bytes];
        // SAFETY: `kmers` e lido no tamanho informado; o lancador escreve
        // `bytes` bytes em `filtro`, que tem esse tamanho.
        confere(unsafe { ffi::launch_bloom_filter_build(kmers.as_ptr(), kmers.len(), filtro.as_mut_ptr(), bytes, num_hashes) })?;
        Ok(filtro)
    }
    #[cfg(not(cuda_enabled))]
    {
        let _ = (kmers, bytes, num_hashes);
        Err(SEM_CUDA.to_string())
    }
}

/// Consulta o filtro: um byte por k-mer, 1 = provavelmente presente. Mesmo
/// resultado de `bloom::consultar`.
pub fn bloom_consultar(kmers: &[u64], filtro: &[u8], num_hashes: i32) -> Result<Vec<u8>, String> {
    exige_placa()?;
    #[cfg(cuda_enabled)]
    {
        let mut mascara = vec![0u8; kmers.len()];
        // SAFETY: `kmers` e `filtro` sao lidos nos tamanhos informados; o
        // lancador escreve um byte por k-mer em `mascara`.
        confere(unsafe {
            ffi::launch_bloom_filter_query(kmers.as_ptr(), kmers.len(), filtro.as_ptr(), filtro.len(), num_hashes, mascara.as_mut_ptr())
        })?;
        Ok(mascara)
    }
    #[cfg(not(cuda_enabled))]
    {
        let _ = (kmers, filtro, num_hashes);
        Err(SEM_CUDA.to_string())
    }
}
