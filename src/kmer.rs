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