#![allow(dead_code)]
// Ensaios do filtro de Bloom na CPU (src/bloom.rs) contra uma referencia
// escrita aqui do jeito mais ingenuo possivel: um k-mer por vez, palavras de
// 32 bits, sem paralelismo. E a mesma conta do kernel CUDA
// (src/merge_path_kernel.cu); a igualdade com a GPU de verdade e conferida por
// trials/python/gpu_test.py numa maquina que tenha placa.
// Uso: cargo run --release --bin bloom
#[path = "../../../src/bloom.rs"]
mod bloom;

fn hash_ref(mut k: u64, semente: u32) -> u64 {
    k ^= semente as u64;
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}

fn construir_ref(kmers: &[u64], bytes: usize, nh: i32) -> Vec<u8> {
    let mut palavras = vec![0u32; bytes / 4];
    let bits = (palavras.len() as u64) * 32;
    if bits > 0 {
        for &k in kmers {
            for i in 0..nh {
                let bit = hash_ref(k, (i as u32).wrapping_mul(0x9E3779B9)) % bits;
                palavras[(bit / 32) as usize] |= 1u32 << (bit % 32);
            }
        }
    }
    let mut saida: Vec<u8> = palavras.iter().flat_map(|p| p.to_le_bytes()).collect();
    saida.resize(bytes, 0);
    saida
}

fn consultar_ref(kmers: &[u64], filtro: &[u8], nh: i32) -> Vec<u8> {
    let palavras: Vec<u32> = filtro.chunks_exact(4).map(|c| u32::from_le_bytes(c.try_into().unwrap())).collect();
    let bits = (palavras.len() as u64) * 32;
    kmers
        .iter()
        .map(|&k| {
            if bits == 0 {
                return 0u8;
            }
            for i in 0..nh {
                let bit = hash_ref(k, (i as u32).wrapping_mul(0x9E3779B9)) % bits;
                if palavras[(bit / 32) as usize] & (1u32 << (bit % 32)) == 0 {
                    return 0u8;
                }
            }
            1u8
        })
        .collect()
}

/// Gerador simples e reprodutivel (splitmix64).
struct Gerador(u64);
impl Gerador {
    fn prox(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

fn main() {
    let mut g = Gerador(20261113);
    let mut ensaios = 0;

    // Tamanhos de filtro que nao sao potencia de dois nem multiplos de 4, de
    // proposito; quantidades de k-mers de zero a meio milhao; 1 a 10 hashes.
    for &bytes in &[0usize, 3, 4, 10, 4096, 65_537, 1 << 20, 3_000_001] {
        for &n in &[0usize, 1, 1000, 500_000] {
            for &nh in &[1i32, 5, 10] {
                let dentro: Vec<u64> = (0..n).map(|_| g.prox()).collect();
                // Metade da consulta esta no filtro, metade nao.
                let consulta: Vec<u64> = dentro.iter().step_by(2).copied().chain((0..n / 2 + 1).map(|_| g.prox())).collect();

                let filtro = bloom::construir(&dentro, bytes, nh);
                assert_eq!(filtro.len(), bytes, "tamanho do filtro");
                assert!(filtro == construir_ref(&dentro, bytes, nh), "filtro difere da referencia: bytes={bytes} n={n} nh={nh}");

                let mascara = bloom::consultar(&consulta, &filtro, nh);
                assert!(mascara == consultar_ref(&consulta, &filtro, nh), "consulta difere da referencia: bytes={bytes} n={n} nh={nh}");

                // Sem falso negativo: o que entrou e sempre achado.
                if bytes >= 4 {
                    assert!(bloom::consultar(&dentro, &filtro, nh).iter().all(|&m| m == 1), "falso negativo: bytes={bytes} n={n} nh={nh}");
                }
                ensaios += 1;
            }
        }
    }

    // Construir duas vezes da o mesmo filtro (as threads nao mudam o resultado).
    let dentro: Vec<u64> = (0..2_000_000).map(|_| g.prox()).collect();
    let a = bloom::construir(&dentro, 16 << 20, 5);
    let b = bloom::construir(&dentro, 16 << 20, 5);
    assert!(a == b, "o filtro mudou entre duas construcoes");
    assert!(a == construir_ref(&dentro, 16 << 20, 5), "2M de k-mers: difere da referencia");
    // Taxa de falso positivo dentro do esperado: (1 - e^(-k n / m))^k.
    let fora: Vec<u64> = (0..1_000_000).map(|_| g.prox()).collect();
    let falsos = bloom::consultar(&fora, &a, 5).iter().filter(|&&m| m == 1).count() as f64 / fora.len() as f64;
    let m = (16u64 << 20) as f64 * 8.0;
    let esperado = (1.0 - (-5.0 * dentro.len() as f64 / m).exp()).powi(5);
    assert!(falsos < esperado * 2.0 + 1e-5, "falsos positivos {falsos:.2e}, esperado {esperado:.2e}");
    ensaios += 3;

    // Dobra de 128 em 64 bits.
    assert_eq!(bloom::dobra(0), 0);
    assert_eq!(bloom::dobra((7u128 << 64) | 9), 7 ^ 9);
    ensaios += 1;

    println!("{ensaios} Bloom filter tests passed (false positives: {falsos:.2e}, expected {esperado:.2e})");
}
