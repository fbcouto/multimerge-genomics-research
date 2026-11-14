use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use crate::kmer::Kmer128;

// Hash rápido e determinístico para as chaves u128 (os nós) e (u128, u128)
// (as arestas). O SipHash padrão do Rust é lento para isso, e com a semente
// aleatória dele a ordem dos contigs empatados mudava de uma execução para
// outra. Os contigs são os mesmos (conferido em 2026-11-14 contra a versão
// anterior; num genoma circular muda só o ponto de partida da volta).
#[derive(Default)]
struct Mistura(u64);

const MULT: u64 = 0x9E37_79B9_7F4A_7C15;

impl Hasher for Mistura {
    #[inline]
    fn finish(&self) -> u64 {
        // finalizador do splitmix64: os bits baixos (o balde) e os altos (a
        // etiqueta do hashbrown) dependem de todos os bits da chave
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ u64::from(b)).wrapping_mul(MULT);
        }
    }
    #[inline]
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(MULT);
    }
    #[inline]
    fn write_u128(&mut self, x: u128) {
        self.write_u64(x as u64);
        self.write_u64((x >> 64) as u64);
    }
    #[inline]
    fn write_usize(&mut self, x: usize) {
        self.write_u64(x as u64);
    }
}

type Rapido = BuildHasherDefault<Mistura>;
type Grafo = HashMap<u128, Vec<u128>, Rapido>;

pub fn construir_grafo_estatisticas(kmers: &[u128], k: u8) -> (usize, usize) {
    let mut grafo: Grafo = HashMap::with_capacity_and_hasher(kmers.len(), Rapido::default());
    for &bits in kmers {
        let kmer = Kmer128 { bits, k };
        let prefixo = kmer.prefix();
        let sufixo = kmer.suffix();
        grafo.entry(prefixo).or_default().push(sufixo);
    }
    let total_nos = grafo.len();
    let total_arestas: usize = grafo.values().map(|v| v.len()).sum();
    (total_nos, total_arestas)
}

pub fn montar_todos_contigs(kmers: &[u128], k: u8) -> Vec<String> {
    println!("🧩 Extracting contigs (u128, K={})...", k);

    let mut grafo: Grafo = HashMap::with_capacity_and_hasher(kmers.len(), Rapido::default());
    let mut in_degree: HashMap<u128, usize, Rapido> = HashMap::with_capacity_and_hasher(kmers.len(), Rapido::default());

    for &bits in kmers {
        let kmer = Kmer128 { bits, k };
        let prefixo = kmer.prefix();
        let sufixo = kmer.suffix();
        grafo.entry(prefixo).or_default().push(sufixo);
        *in_degree.entry(sufixo).or_insert(0) += 1;
        in_degree.entry(prefixo).or_insert(0);
    }

    let mut contigs = Vec::new();
    let mut arestas_visitadas: HashSet<(u128, u128), Rapido> = HashSet::with_capacity_and_hasher(kmers.len(), Rapido::default());
    let mut bolhas_colapsadas = 0;

    for (&u, destinos) in &grafo {
        let out_u = destinos.len();
        let in_u = *in_degree.get(&u).unwrap_or(&0);

        if !(in_u == 1 && out_u == 1) && out_u > 0 {
            for &v in destinos {
                if arestas_visitadas.insert((u, v)) {
                    extrair_caminho_linear(u, v, &grafo, &in_degree, &mut arestas_visitadas, &mut contigs, &mut bolhas_colapsadas, k);
                }
            }
        }
    }

    for (&u, destinos) in &grafo {
        for &v in destinos {
            if arestas_visitadas.insert((u, v)) {
                extrair_caminho_linear(u, v, &grafo, &in_degree, &mut arestas_visitadas, &mut contigs, &mut bolhas_colapsadas, k);
            }
        }
    }

    println!("🫧 Done: {} simple bubble(s) collapsed.", bolhas_colapsadas);
    contigs.sort_by(|a, b| b.len().cmp(&a.len()));
    contigs
}

#[allow(clippy::too_many_arguments)]
fn extrair_caminho_linear(
    inicio_u: u128,
    inicio_v: u128,
    grafo: &Grafo,
    in_degree: &HashMap<u128, usize, Rapido>,
    visitados: &mut HashSet<(u128, u128), Rapido>,
    contigs: &mut Vec<String>,
    bolhas_colapsadas: &mut u32,
    k: u8
) {
    let mut caminho_atual = vec![inicio_u, inicio_v];
    let mut atual = inicio_v;

    while let Some(proximos) = grafo.get(&atual) {
        let out_atual = proximos.len();
        let in_atual = *in_degree.get(&atual).unwrap_or(&0);

        if in_atual == 1 && out_atual == 1 {
            let proximo_no = proximos[0];
            if visitados.insert((atual, proximo_no)) {
                caminho_atual.push(proximo_no);
                atual = proximo_no;
                continue;
            } else { break; }
        }
        else if out_atual == 2 {
            let v1 = proximos[0];
            let v2 = proximos[1];

            if let (Some(d1), Some(d2)) = (grafo.get(&v1), grafo.get(&v2)) {
                if d1.len() == 1 && d2.len() == 1 && d1[0] == d2[0] {
                    let no_convergencia = d1[0];
                    visitados.insert((atual, v1));
                    visitados.insert((atual, v2));
                    visitados.insert((v1, no_convergencia));
                    visitados.insert((v2, no_convergencia));

                    caminho_atual.push(v1);
                    caminho_atual.push(no_convergencia);
                    atual = no_convergencia;
                    *bolhas_colapsadas += 1;
                    continue;
                }
            }
            break;
        } else {
            break;
        }
    }

    // O primeiro nó (k-1 bases) inteiro; de cada nó seguinte, só a última
    // base, que são os 2 bits mais baixos (o `encode` empurra para a esquerda).
    // Antes cada nó era decodificado inteiro, com uma String nova, só para
    // pegar a última letra.
    let mut dna = String::with_capacity(caminho_atual.len() + (k as usize));
    let primeiro_no = Kmer128 { bits: caminho_atual[0], k: k - 1 };
    dna.push_str(&primeiro_no.decode());
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    for &no in &caminho_atual[1..] {
        dna.push(BASES[(no & 3) as usize]);
    }
    contigs.push(dna);
}
