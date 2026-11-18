use std::fmt;

/// A nossa nova fundação genómica. 
/// Usa 128 bits para permitir K-mers gigantes (até k=63) em Zero-Copy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Kmer128 {
    pub bits: u128,
    pub k: u8, // Guardamos o tamanho do K-mer (ex: 31, 55, 63)
}

impl Kmer128 {
    /// Transforma uma String de DNA (ex: "ATGC...") em zeros e uns ultra-rápidos
    pub fn encode(sequencia: &[u8]) -> Self {
        let mut bits: u128 = 0;
        let k = sequencia.len() as u8;
        
        // Proteção contra K-mers grandes demais para u128
        assert!(k <= 63, "The maximum k-mer size for u128 is 63!");

        for &base in sequencia {
            bits <<= 2; // Empurra os bits 2 casas para a esquerda
            bits |= match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _ => 0b00, // Ignora N's (transforma em A por segurança matemática)
            };
        }
        
        Kmer128 { bits, k }
    }

    /// Faz o inverso: Transforma os bits da máquina de volta em "A, C, G, T"
    pub fn decode(&self) -> String {
        let mut temp_bits = self.bits;
        
        // Criar uma máscara de 2 bits (0b11)
        let mask: u128 = 0b11;

        // Lemos de trás para a frente devido ao shift left na codificação
        let mut letras = Vec::with_capacity(self.k as usize);
        for _ in 0..self.k {
            let base_bits = temp_bits & mask;
            let letra = match base_bits {
                0b00 => 'A',
                0b01 => 'C',
                0b10 => 'G',
                0b11 => 'T',
                _ => unreachable!(),
            };
            letras.push(letra);
            temp_bits >>= 2;
        }
        
        letras.into_iter().rev().collect()
    }
    
    /// Retorna o Prefixo (necessário para a montagem de De Bruijn)
    pub fn prefix(&self) -> u128 {
        self.bits >> 2
    }

    /// Retorna o Sufixo (necessário para a montagem de De Bruijn)
    pub fn suffix(&self) -> u128 {
        let mask = (1 << ((self.k - 1) * 2)) - 1;
        self.bits & mask
    }
}

// =====================================================================
// K-MERS DE UMA SEQUENCIA, UMA BASE POR VEZ
// =====================================================================
//
// Os leitores da montagem (`fasta.rs`, `fastq.rs`) pegavam cada janela de k
// bases, procuravam `N` e `n` nela (duas varreduras) e a codificavam do zero
// com `Kmer128::encode` (k passos): ~3k operacoes por k-mer. Aqui o k-mer
// anda uma base por vez: entra a base nova nos 2 bits de baixo e a mascara
// tira a mais antiga. Os bits sao exatamente os do `encode`, inclusive no
// que ele faz com letra que nao e ACGT nem N (vale A).

/// O codigo de 2 bits de cada byte, como no `encode`; `N_BASE` marca N/n,
/// que tira a janela.
const N_BASE: u8 = 4;
const TABELA: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'N' as usize] = N_BASE;
    t[b'n' as usize] = N_BASE;
    t
};

/// Acrescenta a `saida`, na ordem, o k-mer de cada janela de `k` bases de
/// `seq` que nao tem `N` nem `n` — os mesmos bits que
/// `Kmer128::encode(janela).bits` —, com um passo por base.
pub fn kmers_rolantes(seq: &[u8], k: usize, saida: &mut Vec<u128>) {
    assert!((1..=63).contains(&k), "The k-mer size must be between 1 and 63 (got {k})");
    if seq.len() < k {
        return;
    }
    let mascara: u128 = (1u128 << (2 * k)) - 1;
    let mut bits: u128 = 0;
    let mut seguidas = 0usize; // bases desde o ultimo N
    for &b in seq {
        let c = TABELA[b as usize];
        if c == N_BASE {
            seguidas = 0;
            bits = 0;
            continue;
        }
        bits = ((bits << 2) | c as u128) & mascara;
        seguidas += 1;
        if seguidas >= k {
            saida.push(bits);
        }
    }
}

/// Janelas por tarefa quando uma sequencia longa e dividida entre as threads.
pub const JANELAS_POR_PEDACO: usize = 1 << 20;

/// Como `kmers_rolantes`, mas uma sequencia longa (um cromossomo) e cortada
/// em pedacos de `JANELAS_POR_PEDACO` janelas, feitos em paralelo; cada
/// pedaco le tambem as k-1 bases do seguinte. A ordem da saida e a mesma.
pub fn kmers_da_sequencia(seq: &[u8], k: usize, saida: &mut Vec<u128>) {
    kmers_em_pedacos(seq, k, JANELAS_POR_PEDACO, saida)
}

/// `kmers_da_sequencia` com o tamanho do pedaco escolhido (para os ensaios).
pub fn kmers_em_pedacos(seq: &[u8], k: usize, pedaco: usize, saida: &mut Vec<u128>) {
    use rayon::prelude::*;
    assert!((1..=63).contains(&k), "The k-mer size must be between 1 and 63 (got {k})");
    if seq.len() < k {
        return;
    }
    let janelas = seq.len() - k + 1;
    let pedaco = pedaco.max(1);
    if janelas <= pedaco {
        kmers_rolantes(seq, k, saida);
        return;
    }
    let inicios: Vec<usize> = (0..janelas).step_by(pedaco).collect();
    let partes: Vec<Vec<u128>> = inicios
        .par_iter()
        .map(|&inicio| {
            let fim = (inicio + pedaco).min(janelas);
            let mut v = Vec::with_capacity(fim - inicio);
            kmers_rolantes(&seq[inicio..fim + k - 1], k, &mut v);
            v
        })
        .collect();
    saida.reserve(partes.iter().map(|p| p.len()).sum());
    for p in partes {
        saida.extend_from_slice(&p);
    }
}

// Permite imprimir o Kmer diretamente com {}
impl fmt::Display for Kmer128 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.decode())
    }
}

// =====================================================================
// CAMADA DE COMPATIBILIDADE (LEGADO U64)
// =====================================================================

pub fn encode_kmer(seq: &[u8]) -> u64 {
    let mut kmer = 0;
    for &b in seq {
        kmer <<= 2;
        kmer |= match b {
            b'A' | b'a' => 0b00,
            b'C' | b'c' => 0b01,
            b'G' | b'g' => 0b10,
            b'T' | b't' => 0b11,
            _ => 0b00,
        };
    }
    kmer
}

pub fn decode_kmer(kmer: u64, k: u8) -> String {
    let mut temp = kmer;
    let mut chars = Vec::with_capacity(k as usize);
    for _ in 0..k {
        let base = match temp & 0b11 {
            0b00 => 'A',
            0b01 => 'C',
            0b10 => 'G',
            0b11 => 'T',
            _ => unreachable!(),
        };
        chars.push(base);
        temp >>= 2;
    }
    chars.into_iter().rev().collect()
}

pub fn compute_minhash_signature(kmer: u64, num_hashes: usize) -> Vec<u64> {
    let mut signature = Vec::with_capacity(num_hashes);
    let mut current_hash = kmer;
    for _ in 0..num_hashes {
        current_hash = current_hash.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        signature.push(current_hash);
    }
    signature
}