//! Filtro de Bloom de k-mers na CPU, bit a bit igual ao do kernel CUDA
//! (`merge_path_kernel.cu`): mesmo hash, mesmas sementes, mesmo arranjo dos
//! bits. E o que roda quando a extensao foi compilada sem CUDA, quando a
//! maquina nao tem placa NVIDIA, ou quando o dado nao cabe na VRAM.
//!
//! Antes deste modulo, sem CUDA o filtro saia todo zerado e a consulta tambem:
//! todo k-mer da amostra era dado como mutacao, sem nenhum aviso.

use rayon::prelude::*;
use std::sync::atomic::{AtomicU32, Ordering};

/// O `gpu_hash_kmer` do kernel.
#[inline]
pub fn hash_kmer(k: u64, semente: u32) -> u64 {
    let mut k = k ^ semente as u64;
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    k ^= k >> 33;
    k
}

/// A semente do i-esimo hash: `i * 0x9E3779B9` em 32 bits, como no kernel.
#[inline]
fn semente(i: i32) -> u32 {
    (i as u32).wrapping_mul(0x9E37_79B9)
}

/// Bits enderecaveis de um filtro de `bytes` bytes. O kernel trabalha em
/// palavras de 32 bits; a sobra (menos de 4 bytes) nao e usada.
#[inline]
pub fn bits_uteis(bytes: usize) -> u64 {
    (bytes / 4) as u64 * 32
}

/// Dobra um k-mer de 128 bits em 64, como a GPU recebe.
#[inline]
pub fn dobra(k: u128) -> u64 {
    ((k >> 64) as u64) ^ (k as u64)
}

/// Constroi o filtro: `bytes` bytes, `num_hashes` bits acesos por k-mer.
pub fn construir(kmers: &[u64], bytes: usize, num_hashes: i32) -> Vec<u8> {
    let bits = bits_uteis(bytes);
    let palavras: Vec<AtomicU32> = (0..bytes / 4).map(|_| AtomicU32::new(0)).collect();
    if bits > 0 {
        kmers.par_iter().for_each(|&k| {
            for i in 0..num_hashes {
                let bit = hash_kmer(k, semente(i)) % bits;
                palavras[(bit / 32) as usize].fetch_or(1u32 << (bit % 32), Ordering::Relaxed);
            }
        });
    }
    // Palavras de 32 bits em little-endian, como a GPU devolve.
    let mut saida = vec![0u8; bytes];
    saida[..palavras.len() * 4]
        .par_chunks_mut(4)
        .with_min_len(1 << 14)
        .zip(palavras.par_iter().with_min_len(1 << 14))
        .for_each(|(destino, p)| destino.copy_from_slice(&p.load(Ordering::Relaxed).to_le_bytes()));
    saida
}

/// Para cada k-mer, 1 se todos os seus bits estao acesos no filtro (o k-mer
/// provavelmente esta na referencia), 0 se algum esta apagado (com certeza nao
/// esta).
pub fn consultar(kmers: &[u64], filtro: &[u8], num_hashes: i32) -> Vec<u8> {
    let bits = bits_uteis(filtro.len());
    if bits == 0 {
        return vec![0u8; kmers.len()];
    }
    kmers
        .par_iter()
        .map(|&k| {
            for i in 0..num_hashes {
                let bit = hash_kmer(k, semente(i)) % bits;
                // O bit `b` da palavra `p` (little-endian) mora no byte
                // `p*4 + b/8`, que e simplesmente `bit / 8`.
                if filtro[(bit / 8) as usize] & (1u8 << (bit % 8)) == 0 {
                    return 0u8;
                }
            }
            1u8
        })
        .collect()
}

#[cfg(test)]
mod testes {
    use super::*;

    #[test]
    fn o_que_entra_e_achado_e_o_resto_quase_nunca() {
        let dentro: Vec<u64> = (0..50_000u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let fora: Vec<u64> = (0..50_000u64).map(|i| (i + 1_000_000).wrapping_mul(0xD6E8_FEB8_6659_FD93)).collect();
        let filtro = construir(&dentro, 1 << 20, 5);
        assert!(consultar(&dentro, &filtro, 5).iter().all(|&m| m == 1), "falso negativo");
        let falsos = consultar(&fora, &filtro, 5).iter().filter(|&&m| m == 1).count();
        assert!(falsos < 50, "falsos positivos demais: {falsos}");
    }

    #[test]
    fn casos_de_borda() {
        assert_eq!(construir(&[1, 2, 3], 0, 5), Vec::<u8>::new());
        assert_eq!(consultar(&[1, 2, 3], &[], 5), vec![0, 0, 0]);
        // Sobra de bytes fora das palavras de 32 bits fica zerada.
        let f = construir(&[7, 8, 9], 10, 3);
        assert_eq!(f.len(), 10);
        assert_eq!(&f[8..], &[0, 0]);
        // Sem hashes nada acende, e a consulta aceita tudo (como no kernel).
        assert!(construir(&[1, 2], 64, 0).iter().all(|&b| b == 0));
        assert_eq!(consultar(&[1, 2], &[0u8; 64], 0), vec![1, 1]);
    }

    #[test]
    fn hash_de_referencia() {
        // Valores fixados para pegar qualquer mudanca na funcao: o filtro
        // gravado por uma versao tem de ser lido pela outra.
        assert_eq!(hash_kmer(0, 0), 0);
        assert_eq!(semente(1), 0x9E37_79B9);
        assert_eq!(semente(3), 0xDAA6_6D2B);
    }
}
