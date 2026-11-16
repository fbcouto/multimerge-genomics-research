//! Grafo de montagem para desenhar, como o Bandage (telas 3 e 7).
//!
//! - `unitig_graph`: o grafo de De Bruijn compactado. Cada nó é um unitig
//!   (caminho máximo sem bifurcação), cada aresta uma sobreposição de k−1
//!   bases entre o fim de um unitig e o começo de outro.
//! - `overlap_graph`: os contigs como nós e as sobreposições exatas entre o
//!   fim de um e o começo de outro, com a escolha que o
//!   `scaffolding::fast_hash_scaffold` faz marcada.
//!
//! Os dois devolvem JSON com um leiaute determinístico: o traçado de cada nó
//! (uma "salsicha" cujo comprimento cresce com o log do tamanho), o ponto do
//! meio dele e o tamanho da folha. Componentes que são um ciclo só viram
//! anéis; os outros, Fruchterman-Reingold em vários níveis; os componentes
//! vão lado a lado, dos maiores para os menores.
//!
//! Relação com `debruijn::montar_todos_contigs`: os unitigs daqui são os
//! contigs de lá. Lá há uma regra de "bolha mínima" (dois ramos de um nó só
//! que reconvergem), mas num grafo de De Bruijn com k ≥ 3 e k-mers distintos
//! ela nunca se aplica: um caminho de duas arestas soletra k+1 bases que
//! ficam todas determinadas pelos dois nós das pontas, então não existem
//! dois caminhos diferentes assim. O que muda é só o começo de um contig
//! circular (lá, o primeiro nó que o hash visitar; aqui, o menor k-mer do
//! ciclo). Este módulo não depende de nenhum outro do crate.

use rayon::prelude::*;
use serde::Serialize;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fmt::Write as _;

/// Lado de um segmento do traçado, nas coordenadas do leiaute.
const UNIDADE: f64 = 10.0;
const NENHUM: u32 = u32::MAX;
/// Quantos nós mostrar quando quem chama não diz.
pub const MAX_NOS_PADRAO: usize = 2000;
/// Acima disto, o grafo de sobreposição é recusado (sequências repetitivas
/// demais para a sobreposição pedida).
const LIMITE_SOBREPOSICOES: usize = 200_000;
const BASES: [u8; 4] = *b"ACGT";

// =====================================================================
//  O grafo de De Bruijn compactado
// =====================================================================

/// Um unitig: um caminho de k-mers (arestas do grafo de (k−1)-mers), dado
/// pelos índices do primeiro e do último no vetor ordenado de k-mers.
#[derive(Clone, Debug)]
struct Unitig {
    primeiro: u32,
    ultimo: u32,
    kmers: u32,
    gc: u64,
    circular: bool,
}

struct Compactado<'a> {
    kmers: &'a [u128],
    k: u8,
    mascara: u128,
    desloc: u32,
    /// índice pelos bits mais altos: os k-mers que começam por `t` estão em
    /// `balde[t]..balde[t + 1]` (a busca fica curta e perto na memória)
    balde: Vec<u32>,
    bits_balde: u32,
    /// O k-mer seguinte dentro do mesmo unitig, ou NENHUM.
    suc: Vec<u32>,
    unitigs: Vec<Unitig>,
}

impl<'a> Compactado<'a> {
    fn novo(kmers: &'a [u128], k: u8) -> Compactado<'a> {
        let desloc = 2 * (k as u32 - 1);
        let bits_balde = (2 * k as u32).min(16);
        let total = 2 * k as u32;
        let mut balde = vec![0u32; (1usize << bits_balde) + 1];
        for &x in kmers {
            balde[(x >> (total - bits_balde)) as usize + 1] += 1;
        }
        for t in 0..(1usize << bits_balde) {
            balde[t + 1] += balde[t];
        }
        Compactado { kmers, k, mascara: (1u128 << desloc) - 1, desloc, balde, bits_balde, suc: Vec::new(), unitigs: Vec::new() }
    }

    /// O primeiro índice com k-mer ≥ x (x < 4^k).
    fn limite_inferior(&self, x: u128) -> usize {
        let t = (x >> (2 * self.k as u32 - self.bits_balde)) as usize;
        let (a, b) = (self.balde[t] as usize, self.balde[t + 1] as usize);
        a + self.kmers[a..b].partition_point(|&v| v < x)
    }

    fn contem(&self, x: u128) -> bool {
        let i = self.limite_inferior(x);
        i < self.kmers.len() && self.kmers[i] == x
    }

    /// Os k-mers que saem do (k−1)-mer `no` (no máximo 4, vizinhos no vetor).
    fn saidas(&self, no: u128) -> std::ops::Range<usize> {
        let lo = self.limite_inferior(no << 2);
        let mut hi = lo;
        while hi < self.kmers.len() && hi - lo < 4 && (self.kmers[hi] >> 2) == no {
            hi += 1;
        }
        lo..hi
    }

    /// Quantos k-mers chegam ao (k−1)-mer `no`.
    fn grau_entrada(&self, no: u128) -> usize {
        (0..4u128).filter(|&b| self.contem((b << self.desloc) | no)).count()
    }

    fn no_inicio(&self, u: &Unitig) -> u128 {
        self.kmers[u.primeiro as usize] >> 2
    }

    fn no_fim(&self, u: &Unitig) -> u128 {
        self.kmers[u.ultimo as usize] & self.mascara
    }

    /// Comprimento em bases: a sequência soletrada (k−1 + nº de k-mers); num
    /// unitig circular, a volta (nº de k-mers).
    fn comprimento(&self, u: &Unitig) -> u64 {
        if u.circular {
            u.kmers as u64
        } else {
            u.kmers as u64 + self.k as u64 - 1
        }
    }

    /// A sequência soletrada, como `montar_todos_contigs` a escreve: o
    /// primeiro k-mer inteiro e a última base de cada um dos seguintes (num
    /// circular, as k−1 primeiras bases se repetem no fim).
    fn soletra(&self, u: &Unitig) -> String {
        let k = self.k as usize;
        let mut s = String::with_capacity(u.kmers as usize + k);
        let x = self.kmers[u.primeiro as usize];
        for t in (0..k).rev() {
            s.push(BASES[((x >> (2 * t)) & 3) as usize] as char);
        }
        let mut e = u.primeiro as usize;
        for _ in 1..u.kmers {
            e = self.suc[e] as usize;
            s.push(BASES[(self.kmers[e] & 3) as usize] as char);
        }
        s
    }

    /// As arestas do grafo de unitigs: U → V quando o (k−1)-mer final de U é
    /// o inicial de V. Ordenadas.
    fn arestas(&self) -> Vec<(u32, u32)> {
        let mut por_primeiro: Vec<(u32, u32)> =
            self.unitigs.iter().enumerate().map(|(id, u)| (u.primeiro, id as u32)).collect();
        por_primeiro.sort_unstable();
        let mut saida = Vec::new();
        for (id, u) in self.unitigs.iter().enumerate() {
            for f in self.saidas(self.no_fim(u)) {
                // todo k-mer que sai de um nó de junção (ou que fecha um
                // ciclo) começa um unitig
                if let Ok(p) = por_primeiro.binary_search_by_key(&(f as u32), |&(pr, _)| pr) {
                    saida.push((id as u32, por_primeiro[p].1));
                } else {
                    debug_assert!(false, "k-mer after a unitig end that starts no unitig");
                }
            }
        }
        saida.sort_unstable();
        saida
    }
}

fn eh_gc(b: u128) -> bool {
    b == 1 || b == 2
}

/// Confere k e os k-mers; devolve-os ordenados e sem repetição.
fn prepara_kmers(kmers: &[u128], k: u8) -> Result<Cow<'_, [u128]>, String> {
    if !(3..=63).contains(&k) {
        return Err(format!("k must be between 3 and 63 for the assembly graph (got {k})"));
    }
    if kmers.is_empty() {
        return Err("No k-mers: the input is empty.".to_string());
    }
    let limite = 1u128 << (2 * k as u32);
    if let Some(&x) = kmers.iter().find(|&&x| x >= limite) {
        let bases = (128 - x.leading_zeros()).div_ceil(2);
        return Err(format!(
            "A k-mer has {bases} bases but k = {k}: were these k-mers read with another k? \
             Use the same k that read them (screen 1 writes it in the file name, e.g. _K31)."
        ));
    }
    if kmers.windows(2).all(|w| w[0] < w[1]) {
        Ok(Cow::Borrowed(kmers))
    } else {
        let mut v = kmers.to_vec();
        v.sort_unstable();
        v.dedup();
        Ok(Cow::Owned(v))
    }
}

/// Os unitigs: os caminhos máximos sem bifurcação, achados como em
/// `montar_todos_contigs` (primeiro os que partem de um nó que não é
/// 1-entra-1-sai, depois os ciclos isolados), na ordem final dos ids:
/// comprimento decrescente, e o menor primeiro k-mer nos empates.
fn compacta(kmers: &[u128], k: u8) -> Compactado<'_> {
    let mut c = Compactado::novo(kmers, k);
    let mascara = c.mascara;
    let n = kmers.len();
    // o sucessor de cada k-mer no seu unitig: o único que sai do nó final,
    // quando esse nó tem uma entrada e uma saída
    let suc: Vec<u32> = {
        let cr = &c;
        (0..n)
            .into_par_iter()
            .map(|i| {
                let y = kmers[i] & mascara;
                let r = cr.saidas(y);
                if r.len() == 1 && cr.grau_entrada(y) == 1 {
                    r.start as u32
                } else {
                    NENHUM
                }
            })
            .collect()
    };
    c.suc = suc;
    let mut tem_pred = vec![false; n];
    for &s in &c.suc {
        if s != NENHUM {
            tem_pred[s as usize] = true;
        }
    }
    let mut visto = vec![false; n];
    let mut unitigs = Vec::new();
    let gc_prefixo = |x: u128| -> u64 {
        // as k−1 primeiras bases do k-mer
        let p = x >> 2;
        (0..(k as u32 - 1)).filter(|&t| eh_gc((p >> (2 * t)) & 3)).count() as u64
    };
    // 1. os que começam num nó que não é 1-entra-1-sai
    for i in 0..n {
        if tem_pred[i] {
            continue;
        }
        let (mut e, mut cont, mut gc) = (i, 0u32, gc_prefixo(kmers[i]));
        loop {
            visto[e] = true;
            cont += 1;
            gc += eh_gc(kmers[e] & 3) as u64;
            let s = c.suc[e];
            if s == NENHUM || visto[s as usize] {
                break;
            }
            e = s as usize;
        }
        unitigs.push(Unitig { primeiro: i as u32, ultimo: e as u32, kmers: cont, gc, circular: false });
    }
    // 2. os ciclos isolados, a partir do menor k-mer de cada um
    for i in 0..n {
        if visto[i] {
            continue;
        }
        let (mut e, mut cont, mut gc, mut circular) = (i, 0u32, 0u64, false);
        loop {
            visto[e] = true;
            cont += 1;
            gc += eh_gc(kmers[e] & 3) as u64;
            let s = c.suc[e];
            if s as usize == i {
                circular = true;
                break;
            }
            if s == NENHUM || visto[s as usize] {
                break; // não acontece: defesa
            }
            e = s as usize;
        }
        if !circular {
            gc += gc_prefixo(kmers[i]);
        }
        unitigs.push(Unitig { primeiro: i as u32, ultimo: e as u32, kmers: cont, gc, circular });
    }
    unitigs.sort_by(|a, b| {
        c.comprimento(b).cmp(&c.comprimento(a)).then(kmers[a.primeiro as usize].cmp(&kmers[b.primeiro as usize]))
    });
    c.unitigs = unitigs;
    c
}

/// As sequências dos unitigs, na ordem dos ids de `unitig_graph`. Num
/// unitig circular, as k−1 primeiras bases se repetem no fim, como em
/// `assemble_debruijn`.
#[allow(dead_code)] // os ensaios usam; serve também para exportar os unitigs
pub fn unitig_sequences(kmers: &[u128], k: u8) -> Result<Vec<String>, String> {
    let kmers = prepara_kmers(kmers, k)?;
    let c = compacta(&kmers, k);
    Ok(c.unitigs.iter().map(|u| c.soletra(u)).collect())
}

// ---------------------------------------------------------------------
//  Saída em JSON
// ---------------------------------------------------------------------

#[derive(Serialize)]
struct NoUnitig {
    id: usize,
    length: u64,
    kmers: u64,
    gc: f64,
    in_degree: usize,
    out_degree: usize,
    circular: bool,
    kind: &'static str,
    component: usize,
    x: f64,
    y: f64,
    size: f64,
    path: Vec<[f64; 2]>,
}

#[derive(Serialize)]
struct ArestaUnitig {
    from: usize,
    to: usize,
    overlap: usize,
}

#[derive(Serialize)]
struct ComponenteJson {
    id: usize,
    nodes: usize,
    edges: usize,
    length: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    kmers: Option<u64>,
    circular: bool,
    shown: &'static str,
    shown_nodes: usize,
}

#[derive(Serialize)]
struct LeiauteJson {
    width: f64,
    height: f64,
    unit: f64,
}

#[derive(Serialize)]
struct ResumoUnitig {
    kmers: u64,
    nodes: usize,
    edges: usize,
    components: usize,
    circular_components: usize,
    length: u64,
    n50: u64,
    longest: u64,
    tips: usize,
    bubbles: usize,
    bubble_branches: usize,
    repeats: usize,
    isolated: usize,
    circular_nodes: usize,
    shown_nodes: usize,
    shown_edges: usize,
    shown_components: usize,
    dropped_components: usize,
    dropped_nodes: usize,
    dropped_kmers: u64,
    partial_component: Option<usize>,
}

#[derive(Serialize)]
struct GrafoUnitig {
    graph: &'static str,
    k: u8,
    max_nodes: usize,
    summary: ResumoUnitig,
    text: String,
    nodes: Vec<NoUnitig>,
    edges: Vec<ArestaUnitig>,
    components: Vec<ComponenteJson>,
    layout: LeiauteJson,
}

/// O grafo de De Bruijn compactado de um conjunto de k-mers (u128, 2 bits
/// por base, uma fita só), com leiaute, em JSON.
///
/// Nós: os unitigs (`id`, `length` em bases, `kmers`, `gc`, `in_degree`,
/// `out_degree`, `circular` — um ciclo isolado —, `kind`, `component`, e o
/// desenho: `x`, `y`, `size`, `path`). `kind`: "tip" (ponta morta com menos
/// de 2k bases que sai de uma bifurcação), "bubble" (um de dois ou mais
/// caminhos paralelos entre os mesmos dois nós, cada um com menos de 3k
/// bases), "repeat" (entra por mais de
/// um caminho e sai por mais de um, sem contar as pontas e contando os ramos
/// de uma bolha como um só), "isolated" (sozinho no seu componente) ou "".
/// Arestas: `from`, `to`, `overlap` (k−1). Com mais de `max_nodes` unitigs,
/// ficam os maiores componentes (e, se o último não couber inteiro, a parte
/// dele mais perto do seu maior unitig); o resumo diz o que ficou de fora.
pub fn unitig_graph(kmers: &[u128], k: u8, max_nodes: usize) -> Result<String, String> {
    if max_nodes == 0 {
        return Err("max_nodes must be at least 1".to_string());
    }
    let kmers = prepara_kmers(kmers, k)?;
    let c = compacta(&kmers, k);
    let n = c.unitigs.len();
    let arestas = c.arestas();
    let mut grau_ent = vec![0usize; n];
    let mut grau_sai = vec![0usize; n];
    for &(a, b) in &arestas {
        grau_sai[a as usize] += 1;
        grau_ent[b as usize] += 1;
    }
    let comprimentos: Vec<u64> = c.unitigs.iter().map(|u| c.comprimento(u)).collect();

    // componentes, dos maiores para os menores; comprimento = bases sem
    // repetir as sobreposições: k-mers + (k−1) por unitig sem antecessor
    let k1 = k as u64 - 1;
    let comps = componentes(n, &arestas, |itens| {
        itens
            .iter()
            .map(|&i| c.unitigs[i as usize].kmers as u64 + if grau_ent[i as usize] == 0 { k1 } else { 0 })
            .sum()
    });

    // tipos
    let mut tipo: Vec<&'static str> = vec![""; n];
    // bolhas: ramos paralelos (mesmo nó de partida e de chegada) com menos
    // de 3k bases cada, o tamanho do que um SNP, um erro de leitura ou um
    // indel pequeno abre; dois trechos longos entre cópias de uma repetição
    // também são paralelos, mas não são uma bolha
    let mut chaves: Vec<(u128, u128, u32)> = Vec::new();
    for (id, u) in c.unitigs.iter().enumerate() {
        let (a, b) = (c.no_inicio(u), c.no_fim(u));
        if !u.circular && a != b && comps.tamanho_de(id) > 1 && comprimentos[id] < 3 * k as u64 {
            chaves.push((a, b, id as u32));
        }
    }
    chaves.sort_unstable();
    let mut bolhas = 0;
    let mut i = 0;
    while i < chaves.len() {
        let mut j = i + 1;
        while j < chaves.len() && chaves[j].0 == chaves[i].0 && chaves[j].1 == chaves[i].1 {
            j += 1;
        }
        if j - i >= 2 {
            bolhas += 1;
            for ch in &chaves[i..j] {
                tipo[ch.2 as usize] = "bubble";
            }
        }
        i = j;
    }
    for (id, u) in c.unitigs.iter().enumerate() {
        if comps.tamanho_de(id) == 1 {
            tipo[id] = "isolated";
            continue;
        }
        if tipo[id] == "bubble" {
            continue;
        }
        let (ge, gs) = (grau_ent[id], grau_sai[id]);
        let curto = comprimentos[id] < 2 * k as u64;
        let ponta5 = ge == 0 && gs >= 1 && c.grau_entrada(c.no_fim(u)) >= 2;
        let ponta3 = gs == 0 && ge >= 1 && c.saidas(c.no_inicio(u)).len() >= 2;
        if curto && (ponta5 || ponta3) {
            tipo[id] = "tip";
        }
    }
    // repetição: entra por mais de um caminho e sai por mais de um. Não
    // contam as pontas, e os ramos de uma bolha contam como um só (o trecho
    // entre duas bolhas de SNP tem 2 entradas e 2 saídas e não é repetição).
    let mut antes: Vec<Vec<u32>> = vec![Vec::new(); n];
    let mut depois: Vec<Vec<u32>> = vec![Vec::new(); n];
    for &(a, b) in &arestas {
        depois[a as usize].push(b);
        antes[b as usize].push(a);
    }
    let mut repeticoes = Vec::new();
    for id in 0..n {
        if !tipo[id].is_empty() || grau_ent[id] < 2 || grau_sai[id] < 2 {
            continue;
        }
        let distintos = |lista: &[u32], ponta: &dyn Fn(&Unitig) -> u128| {
            let mut v: Vec<u128> =
                lista.iter().filter(|&&x| tipo[x as usize] != "tip").map(|&x| ponta(&c.unitigs[x as usize])).collect();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        if distintos(&antes[id], &|w| c.no_inicio(w)) >= 2 && distintos(&depois[id], &|w| c.no_fim(w)) >= 2 {
            repeticoes.push(id);
        }
    }
    for id in repeticoes {
        tipo[id] = "repeat";
    }

    // o que mostrar
    let adj = Csr::novo(n, &arestas);
    let escolha = escolhe(&comps.lista, &adj, &comprimentos, max_nodes);

    // leiaute
    let segs: Vec<usize> = comprimentos.iter().map(|&l| segmentos(l, k as u64)).collect();
    let anel: Vec<bool> = c.unitigs.iter().map(|u| u.circular).collect();
    let desenho = desenha(&segs, &anel, &arestas, &escolha.mostrados);

    // JSON
    let mostrado = |id: usize| !desenho.tracado[id].is_empty();
    let mut nos = Vec::new();
    for (id, u) in c.unitigs.iter().enumerate() {
        if !mostrado(id) {
            continue;
        }
        let total = if u.circular { u.kmers as u64 } else { comprimentos[id] };
        nos.push(NoUnitig {
            id,
            length: comprimentos[id],
            kmers: u.kmers as u64,
            gc: arred(u.gc as f64 / total.max(1) as f64, 10_000.0),
            in_degree: grau_ent[id],
            out_degree: grau_sai[id],
            circular: u.circular,
            kind: tipo[id],
            component: comps.de[id] as usize,
            x: arred(desenho.centro[id][0], 10.0),
            y: arred(desenho.centro[id][1], 10.0),
            size: (2 * segs[id]) as f64 * UNIDADE,
            path: desenho.tracado[id].iter().map(|p| [arred(p[0], 10.0), arred(p[1], 10.0)]).collect(),
        });
    }
    let arestas_json: Vec<ArestaUnitig> = arestas
        .iter()
        .filter(|&&(a, b)| mostrado(a as usize) && mostrado(b as usize))
        .map(|&(a, b)| ArestaUnitig { from: a as usize, to: b as usize, overlap: k as usize - 1 })
        .collect();
    let mut componentes_json = Vec::new();
    let mut circ_comps = 0;
    for (ci, itens) in comps.lista.iter().enumerate() {
        let circular = ciclo_unico(itens, &grau_ent, &grau_sai, comps.arestas[ci]);
        circ_comps += circular as usize;
        let mostr = escolha.mostrados[ci].len();
        if mostr == 0 {
            continue;
        }
        componentes_json.push(ComponenteJson {
            id: ci,
            nodes: itens.len(),
            edges: comps.arestas[ci],
            length: comps.comprimento[ci],
            kmers: Some(itens.iter().map(|&i| c.unitigs[i as usize].kmers as u64).sum()),
            circular,
            shown: if mostr == itens.len() { "all" } else { "part" },
            shown_nodes: mostr,
        });
    }
    let mut ordenados = comprimentos.clone();
    ordenados.sort_unstable_by(|a, b| b.cmp(a));
    let n50 = n50(&ordenados);
    let shown_nodes = nos.len();
    let dropped_kmers: u64 =
        c.unitigs.iter().enumerate().filter(|(id, _)| !mostrado(*id)).map(|(_, u)| u.kmers as u64).sum();
    let resumo = ResumoUnitig {
        kmers: kmers.len() as u64,
        nodes: n,
        edges: arestas.len(),
        components: comps.lista.len(),
        circular_components: circ_comps,
        length: comps.comprimento.iter().sum(),
        n50,
        longest: ordenados[0],
        tips: tipo.iter().filter(|&&t| t == "tip").count(),
        bubbles: bolhas,
        bubble_branches: tipo.iter().filter(|&&t| t == "bubble").count(),
        repeats: tipo.iter().filter(|&&t| t == "repeat").count(),
        isolated: tipo.iter().filter(|&&t| t == "isolated").count(),
        circular_nodes: c.unitigs.iter().filter(|u| u.circular).count(),
        shown_nodes,
        shown_edges: arestas_json.len(),
        shown_components: componentes_json.len(),
        dropped_components: escolha.mostrados.iter().filter(|m| m.is_empty()).count(),
        dropped_nodes: n - shown_nodes,
        dropped_kmers,
        partial_component: escolha.parcial,
    };
    let text = texto_unitigs(&resumo, k, max_nodes);
    let g = GrafoUnitig {
        graph: "unitig",
        k,
        max_nodes,
        summary: resumo,
        text,
        nodes: nos,
        edges: arestas_json,
        components: componentes_json,
        layout: LeiauteJson { width: arred(desenho.largura, 10.0), height: arred(desenho.altura, 10.0), unit: UNIDADE },
    };
    serde_json::to_string(&g).map_err(|e| format!("Could not write the graph as JSON: {e}"))
}

fn texto_unitigs(r: &ResumoUnitig, k: u8, max_nodes: usize) -> String {
    let mut t = format!(
        "{} (k = {k}) make {} in {}, joined by {} of {} bases. Unitig N50: {} bases; longest: {} bases. \
         Tips: {}; bubbles: {} ({} branches); repeats: {}; circular unitigs: {}.",
        conta(r.kmers, "k-mer", "k-mers"),
        conta(r.nodes as u64, "unitig", "unitigs"),
        conta(r.components as u64, "component", "components"),
        conta(r.edges as u64, "overlap", "overlaps"),
        k - 1,
        milhar(r.n50),
        milhar(r.longest),
        milhar(r.tips as u64),
        milhar(r.bubbles as u64),
        milhar(r.bubble_branches as u64),
        milhar(r.repeats as u64),
        milhar(r.circular_nodes as u64),
    );
    if r.dropped_nodes > 0 {
        let _ = write!(
            t,
            " Showing {} of the {} unitigs (max_nodes = {}): the {} largest {}",
            milhar(r.shown_nodes as u64),
            milhar(r.nodes as u64),
            milhar(max_nodes as u64),
            milhar(r.shown_components as u64),
            if r.shown_components == 1 { "component" } else { "components" },
        );
        if r.dropped_components > 0 {
            let _ = write!(
                t,
                "; {} left out ({}, {}).",
                conta(r.dropped_components as u64, "smaller component", "smaller components"),
                conta((r.dropped_nodes) as u64, "unitig", "unitigs"),
                conta(r.dropped_kmers, "k-mer", "k-mers"),
            );
        } else {
            t.push('.');
        }
        if let Some(p) = r.partial_component {
            let _ = write!(t, " Component {p} is shown only in part: the unitigs nearest its longest one.");
        }
    }
    t
}

// =====================================================================
//  O grafo de sobreposição dos contigs (tela 7)
// =====================================================================

#[derive(Serialize)]
struct NoContig {
    id: usize,
    length: u64,
    gc: f64,
    in_degree: usize,
    out_degree: usize,
    circular: bool,
    component: usize,
    scaffold: usize,
    x: f64,
    y: f64,
    size: f64,
    path: Vec<[f64; 2]>,
}

#[derive(Serialize)]
struct ArestaContig {
    from: usize,
    to: usize,
    overlap: usize,
    scaffold_rule: bool,
    used: bool,
}

#[derive(Serialize)]
struct AndaimeJson {
    id: usize,
    contigs: Vec<usize>,
    length: u64,
    circular: bool,
}

#[derive(Serialize)]
struct ResumoContig {
    contigs: usize,
    edges: usize,
    scaffold_rule_edges: usize,
    used_edges: usize,
    longer_overlaps: usize,
    components: usize,
    circular_components: usize,
    scaffolds: usize,
    circular_scaffolds: usize,
    unconnected: usize,
    shown_nodes: usize,
    shown_edges: usize,
    shown_components: usize,
    dropped_components: usize,
    dropped_nodes: usize,
    partial_component: Option<usize>,
}

#[derive(Serialize)]
struct GrafoContig {
    graph: &'static str,
    min_overlap: usize,
    max_nodes: usize,
    summary: ResumoContig,
    text: String,
    nodes: Vec<NoContig>,
    edges: Vec<ArestaContig>,
    components: Vec<ComponenteJson>,
    scaffolds: Vec<AndaimeJson>,
    layout: LeiauteJson,
}

/// O grafo de sobreposição dos contigs (tela 7), com no máximo
/// `MAX_NOS_PADRAO` nós desenhados.
#[allow(dead_code)] // a ligação ao Python usa a versão com limite
pub fn overlap_graph(sequences: &[String], min_overlap: usize) -> Result<String, String> {
    overlap_graph_with_limit(sequences, min_overlap, MAX_NOS_PADRAO)
}

/// O grafo de sobreposição dos contigs, em JSON.
///
/// Nós: os contigs, com o id = a posição na lista de entrada (`length`,
/// `gc`, graus, `circular` = o fim sobrepõe o próprio começo, `component`,
/// `scaffold` e o desenho). Arestas: `from` → `to` quando as últimas bases
/// de `from` são as primeiras de `to`, na mesma fita; `overlap` é a maior
/// dessas sobreposições, de pelo menos `min_overlap` bases. `scaffold_rule`:
/// as últimas `min_overlap` bases de `from` são exatamente as primeiras
/// `min_overlap` de `to` — o teste de `build_scaffolds(…, overlap =
/// min_overlap)`; `used`: o encaixe que ele de fato fez (cada contig recebe
/// um encaixe de cada lado, o primeiro que aparecer). `scaffolds` repete a
/// costura do `build_scaffolds`: os contigs em ordem, o comprimento e se
/// fechou em círculo, na mesma ordem que ele devolve.
pub fn overlap_graph_with_limit(sequences: &[String], min_overlap: usize, max_nodes: usize) -> Result<String, String> {
    if sequences.is_empty() {
        return Err("No contigs: the input is empty.".to_string());
    }
    if min_overlap == 0 {
        return Err("min_overlap must be at least 1 base".to_string());
    }
    if max_nodes == 0 {
        return Err("max_nodes must be at least 1".to_string());
    }
    let seqs: Vec<&[u8]> = sequences.iter().map(|s| s.as_bytes()).collect();
    let n = seqs.len();
    let m = min_overlap;
    let sobrep = sobreposicoes(&seqs, m)?;
    let arestas: Vec<(u32, u32)> = sobrep.iter().map(|&(a, b, _)| (a, b)).collect();
    let regra: Vec<bool> = sobrep
        .iter()
        .map(|&(a, b, _)| {
            let (sa, sb) = (seqs[a as usize], seqs[b as usize]);
            sa.len() >= m && sb.len() >= m && sa[sa.len() - m..] == sb[..m]
        })
        .collect();
    let costura = costura_como_scaffolding(&seqs, m);
    let usada: Vec<bool> =
        sobrep.iter().map(|&(a, b, _)| costura.proximo[a as usize] == b && a != b).collect();

    let mut grau_ent = vec![0usize; n];
    let mut grau_sai = vec![0usize; n];
    for &(a, b) in &arestas {
        grau_sai[a as usize] += 1;
        grau_ent[b as usize] += 1;
    }
    let mut andaime_de = vec![0usize; n];
    for (si, a) in costura.andaimes.iter().enumerate() {
        for &c in &a.contigs {
            andaime_de[c] = si;
        }
    }
    // comprimento de um componente: o dos andaimes que estão nele
    let comps = componentes(n, &arestas, |itens| {
        itens
            .iter()
            .filter(|&&i| costura.andaimes[andaime_de[i as usize]].contigs[0] == i as usize)
            .map(|&i| costura.andaimes[andaime_de[i as usize]].comprimento)
            .sum()
    });
    let adj = Csr::novo(n, &arestas);
    let tamanhos: Vec<u64> = seqs.iter().map(|s| s.len() as u64).collect();
    let escolha = escolhe(&comps.lista, &adj, &tamanhos, max_nodes);
    let base = (m as u64 + 1).max(2);
    let segs: Vec<usize> = seqs.iter().map(|s| segmentos(s.len() as u64, base)).collect();
    let auto: Vec<bool> = {
        let mut v = vec![false; n];
        for &(a, b) in &arestas {
            if a == b {
                v[a as usize] = true;
            }
        }
        v
    };
    let anel: Vec<bool> = (0..n).map(|i| auto[i] && comps.tamanho_de(i) == 1).collect();
    let desenho = desenha(&segs, &anel, &arestas, &escolha.mostrados);
    let mostrado = |id: usize| !desenho.tracado[id].is_empty();

    let mut nos = Vec::new();
    for (id, s) in seqs.iter().enumerate() {
        if !mostrado(id) {
            continue;
        }
        let gc = s.iter().filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c')).count();
        nos.push(NoContig {
            id,
            length: s.len() as u64,
            gc: arred(gc as f64 / s.len().max(1) as f64, 10_000.0),
            in_degree: grau_ent[id],
            out_degree: grau_sai[id],
            circular: auto[id],
            component: comps.de[id] as usize,
            scaffold: andaime_de[id],
            x: arred(desenho.centro[id][0], 10.0),
            y: arred(desenho.centro[id][1], 10.0),
            size: (2 * segs[id]) as f64 * UNIDADE,
            path: desenho.tracado[id].iter().map(|p| [arred(p[0], 10.0), arred(p[1], 10.0)]).collect(),
        });
    }
    let mut arestas_json = Vec::new();
    for (e, &(a, b, l)) in sobrep.iter().enumerate() {
        if mostrado(a as usize) && mostrado(b as usize) {
            arestas_json.push(ArestaContig {
                from: a as usize,
                to: b as usize,
                overlap: l as usize,
                scaffold_rule: regra[e],
                used: usada[e],
            });
        }
    }
    let mut componentes_json = Vec::new();
    let mut circ_comps = 0;
    for (ci, itens) in comps.lista.iter().enumerate() {
        let circular = ciclo_unico(itens, &grau_ent, &grau_sai, comps.arestas[ci]);
        circ_comps += circular as usize;
        let mostr = escolha.mostrados[ci].len();
        if mostr == 0 {
            continue;
        }
        componentes_json.push(ComponenteJson {
            id: ci,
            nodes: itens.len(),
            edges: comps.arestas[ci],
            length: comps.comprimento[ci],
            kmers: None,
            circular,
            shown: if mostr == itens.len() { "all" } else { "part" },
            shown_nodes: mostr,
        });
    }
    let andaimes_json: Vec<AndaimeJson> = costura
        .andaimes
        .iter()
        .enumerate()
        .map(|(i, a)| AndaimeJson { id: i, contigs: a.contigs.clone(), length: a.comprimento, circular: a.circular })
        .collect();
    let shown_nodes = nos.len();
    let resumo = ResumoContig {
        contigs: n,
        edges: arestas.len(),
        scaffold_rule_edges: regra.iter().filter(|&&r| r).count(),
        used_edges: usada.iter().filter(|&&u| u).count(),
        longer_overlaps: sobrep.iter().filter(|&&(_, _, l)| l as usize > m).count(),
        components: comps.lista.len(),
        circular_components: circ_comps,
        scaffolds: andaimes_json.len(),
        circular_scaffolds: andaimes_json.iter().filter(|a| a.circular).count(),
        unconnected: (0..n).filter(|&i| grau_ent[i] == 0 && grau_sai[i] == 0).count(),
        shown_nodes,
        shown_edges: arestas_json.len(),
        shown_components: componentes_json.len(),
        dropped_components: escolha.mostrados.iter().filter(|m| m.is_empty()).count(),
        dropped_nodes: n - shown_nodes,
        partial_component: escolha.parcial,
    };
    let text = texto_contigs(&resumo, &sobrep, &regra, &andaimes_json, m, max_nodes);
    let g = GrafoContig {
        graph: "overlap",
        min_overlap: m,
        max_nodes,
        summary: resumo,
        text,
        nodes: nos,
        edges: arestas_json,
        components: componentes_json,
        scaffolds: andaimes_json,
        layout: LeiauteJson { width: arred(desenho.largura, 10.0), height: arred(desenho.altura, 10.0), unit: UNIDADE },
    };
    serde_json::to_string(&g).map_err(|e| format!("Could not write the graph as JSON: {e}"))
}

fn texto_contigs(
    r: &ResumoContig,
    sobrep: &[(u32, u32, u32)],
    regra: &[bool],
    andaimes: &[AndaimeJson],
    m: usize,
    max_nodes: usize,
) -> String {
    let mais_longo = andaimes.iter().map(|a| a.length).max().unwrap_or(0);
    let mut t = format!(
        "{}: {} of at least {m} bases ({} of exactly {m}, the joins build_scaffolds looks for) in {}. \
         build_scaffolds(overlap = {m}) joins them into {} ({} circular); the longest has {} bases.",
        conta(r.contigs as u64, "contig", "contigs"),
        conta(r.edges as u64, "overlap", "overlaps"),
        milhar(r.scaffold_rule_edges as u64),
        conta(r.components as u64, "component", "components"),
        conta(r.scaffolds as u64, "scaffold", "scaffolds"),
        milhar(r.circular_scaffolds as u64),
        milhar(mais_longo),
    );
    if r.scaffold_rule_edges > r.used_edges {
        let _ = write!(
            t,
            " {} of the exact overlaps {} not used: each contig takes one join at each end, the first one found.",
            milhar((r.scaffold_rule_edges - r.used_edges) as u64),
            if r.scaffold_rule_edges - r.used_edges == 1 { "was" } else { "were" },
        );
    }
    // sobreposições que a regra exata não pega: sugerir o tamanho mais comum
    let mut fora: HashMap<u32, usize> = HashMap::new();
    for (e, &(_, _, l)) in sobrep.iter().enumerate() {
        if !regra[e] {
            *fora.entry(l).or_default() += 1;
        }
    }
    if !fora.is_empty() {
        let total: usize = fora.values().sum();
        let (&comum, _) = fora.iter().max_by_key(|&(l, c)| (*c, std::cmp::Reverse(*l))).unwrap();
        let _ = write!(
            t,
            " {} longer than {m} bases: build_scaffolds joins contigs only when they overlap by exactly `overlap` \
             bases, so it does not use {}; the most common of them is {} bases.",
            if total == 1 { "1 overlap is".to_string() } else { format!("{} overlaps are", milhar(total as u64)) },
            if total == 1 { "it" } else { "them" },
            milhar(comum as u64),
        );
    }
    if r.unconnected > 0 {
        let _ = write!(t, " {} no other.", conta(r.unconnected as u64, "contig overlaps", "contigs overlap"));
    }
    if r.dropped_nodes > 0 {
        let _ = write!(
            t,
            " Showing {} of the {} contigs (max_nodes = {}): the {} largest {}; {} left out.",
            milhar(r.shown_nodes as u64),
            milhar(r.contigs as u64),
            milhar(max_nodes as u64),
            milhar(r.shown_components as u64),
            if r.shown_components == 1 { "component" } else { "components" },
            conta(r.dropped_components as u64, "smaller component", "smaller components"),
        );
        if let Some(p) = r.partial_component {
            let _ = write!(t, " Component {p} is shown only in part: the contigs nearest its longest one.");
        }
    }
    t
}

/// Todas as sobreposições exatas fim → começo de pelo menos `m` bases (a
/// maior de cada par), por hash rolante das janelas de `m` bases. Um contig
/// pode sobrepor a si mesmo (círculo), desde que não inteiro.
fn sobreposicoes(seqs: &[&[u8]], m: usize) -> Result<Vec<(u32, u32, u32)>, String> {
    const B: u64 = 0x0000_0100_0000_01B3;
    let n = seqs.len();
    let pot = (1..m).fold(1u64, |p, _| p.wrapping_mul(B));
    let hash = |s: &[u8]| s.iter().fold(0u64, |h, &c| h.wrapping_mul(B).wrapping_add(c as u64));
    let mut prefixos: HashMap<u64, Vec<u32>> = HashMap::new();
    for (j, s) in seqs.iter().enumerate() {
        if s.len() >= m {
            prefixos.entry(hash(&s[..m])).or_default().push(j as u32);
        }
    }
    let mut visto = vec![NENHUM; n];
    let mut saida: Vec<(u32, u32, u32)> = Vec::new();
    for (i, s) in seqs.iter().enumerate() {
        if s.len() < m {
            continue;
        }
        let mut h = hash(&s[..m]);
        for p in 0..=(s.len() - m) {
            if p > 0 {
                h = h.wrapping_sub((s[p - 1] as u64).wrapping_mul(pot)).wrapping_mul(B).wrapping_add(s[p + m - 1] as u64);
            }
            let Some(cands) = prefixos.get(&h) else { continue };
            let l = s.len() - p;
            for &j in cands {
                let ju = j as usize;
                if visto[ju] == i as u32 || (ju == i && p == 0) || l > seqs[ju].len() {
                    continue;
                }
                if s[p..] == seqs[ju][..l] {
                    visto[ju] = i as u32;
                    saida.push((i as u32, j, l as u32));
                    if saida.len() > LIMITE_SOBREPOSICOES {
                        return Err(format!(
                            "More than {} overlaps of at least {m} bases: the contigs are too repetitive for this view. \
                             Raise min_overlap.",
                            milhar(LIMITE_SOBREPOSICOES as u64)
                        ));
                    }
                }
            }
        }
    }
    saida.sort_unstable();
    Ok(saida)
}

struct Andaime {
    contigs: Vec<usize>,
    comprimento: u64,
    circular: bool,
}

struct Costura {
    /// o encaixe que o `fast_hash_scaffold` faz a partir de cada contig
    proximo: Vec<u32>,
    /// na ordem em que ele devolve os andaimes
    andaimes: Vec<Andaime>,
}

/// A costura do `scaffolding::fast_hash_scaffold`, passo a passo, guardando
/// quais contigs entram em cada andaime (ele só devolve as sequências).
fn costura_como_scaffolding(seqs: &[&[u8]], m: usize) -> Costura {
    let n = seqs.len();
    // 1. o começo de cada contig (o último com um dado começo fica)
    let mut prefixos: HashMap<&[u8], usize> = HashMap::with_capacity(n);
    for (i, c) in seqs.iter().enumerate() {
        if c.len() >= m {
            prefixos.insert(&c[..m], i);
        }
    }
    // 2. o fim de um encaixa no começo de outro
    let mut proximo = vec![NENHUM; n];
    let mut entrada = vec![0u32; n];
    for (i, c) in seqs.iter().enumerate() {
        if c.len() >= m {
            if let Some(&j) = prefixos.get(&c[c.len() - m..]) {
                if i != j && entrada[j] == 0 && proximo[i] == NENHUM {
                    proximo[i] = j as u32;
                    entrada[j] += 1;
                }
            }
        }
    }
    // 3. os caminhos a partir de quem não tem antecessor; depois os anéis
    let mut visto = vec![false; n];
    let mut andaimes: Vec<(Vec<u8>, Andaime)> = Vec::new();
    for i in 0..n {
        if entrada[i] != 0 || visto[i] {
            continue;
        }
        let mut cur = i;
        let mut seq = seqs[cur].to_vec();
        let mut lista = vec![cur];
        visto[cur] = true;
        while proximo[cur] != NENHUM {
            let nx = proximo[cur] as usize;
            if visto[nx] {
                break;
            }
            visto[nx] = true;
            seq.extend_from_slice(&seqs[nx][m..]);
            lista.push(nx);
            cur = nx;
        }
        let mut circular = false;
        if seq.len() >= 2 * m && seq[..m] == seq[seq.len() - m..] {
            circular = true;
            seq.truncate(seq.len() - m);
        }
        let comprimento = seq.len() as u64;
        andaimes.push((seq, Andaime { contigs: lista, comprimento, circular }));
    }
    for i in 0..n {
        if visto[i] {
            continue;
        }
        let mut cur = i;
        let mut seq = seqs[cur].to_vec();
        let mut lista = vec![cur];
        visto[cur] = true;
        let mut circular = false;
        while proximo[cur] != NENHUM {
            let nx = proximo[cur] as usize;
            if visto[nx] {
                circular = true;
                seq.truncate(seq.len() - m);
                break;
            }
            visto[nx] = true;
            seq.extend_from_slice(&seqs[nx][m..]);
            lista.push(nx);
            cur = nx;
        }
        let comprimento = seq.len() as u64;
        andaimes.push((seq, Andaime { contigs: lista, comprimento, circular }));
    }
    // a mesma ordenação estável, do maior para o menor
    andaimes.sort_by_key(|a| std::cmp::Reverse(a.0.len()));
    Costura { proximo, andaimes: andaimes.into_iter().map(|(_, a)| a).collect() }
}

// =====================================================================
//  Componentes e o que mostrar
// =====================================================================

struct Componentes {
    /// item → componente
    de: Vec<u32>,
    /// itens de cada componente (ids crescentes), dos maiores para os menores
    lista: Vec<Vec<u32>>,
    comprimento: Vec<u64>,
    arestas: Vec<usize>,
}

impl Componentes {
    fn tamanho_de(&self, item: usize) -> usize {
        self.lista[self.de[item] as usize].len()
    }
}

fn raiz(pai: &mut [u32], mut x: u32) -> u32 {
    while pai[x as usize] != x {
        pai[x as usize] = pai[pai[x as usize] as usize];
        x = pai[x as usize];
    }
    x
}

/// Componentes (conexos sem direção), ordenados por comprimento, depois nº
/// de itens, depois o menor id.
fn componentes(n: usize, arestas: &[(u32, u32)], comprimento: impl Fn(&[u32]) -> u64) -> Componentes {
    let mut pai: Vec<u32> = (0..n as u32).collect();
    for &(a, b) in arestas {
        let (ra, rb) = (raiz(&mut pai, a), raiz(&mut pai, b));
        if ra != rb {
            let (menor, maior) = if ra < rb { (ra, rb) } else { (rb, ra) };
            pai[maior as usize] = menor;
        }
    }
    let mut grupos: HashMap<u32, Vec<u32>> = HashMap::new();
    let mut ordem_raizes = Vec::new();
    for i in 0..n as u32 {
        let r = raiz(&mut pai, i);
        let g = grupos.entry(r).or_default();
        if g.is_empty() {
            ordem_raizes.push(r);
        }
        g.push(i);
    }
    let mut lista: Vec<(u64, Vec<u32>)> = ordem_raizes
        .into_iter()
        .map(|r| {
            let itens = grupos.remove(&r).unwrap();
            (comprimento(&itens), itens)
        })
        .collect();
    lista.sort_by(|a, b| b.0.cmp(&a.0).then(b.1.len().cmp(&a.1.len())).then(a.1[0].cmp(&b.1[0])));
    let mut de = vec![0u32; n];
    for (ci, (_, itens)) in lista.iter().enumerate() {
        for &i in itens {
            de[i as usize] = ci as u32;
        }
    }
    let mut num_arestas = vec![0usize; lista.len()];
    for &(a, _) in arestas {
        num_arestas[de[a as usize] as usize] += 1;
    }
    let comprimento = lista.iter().map(|x| x.0).collect();
    Componentes { de, lista: lista.into_iter().map(|x| x.1).collect(), comprimento, arestas: num_arestas }
}

/// O componente é um ciclo dirigido só: todo item com uma entrada e uma
/// saída, e tantas arestas quanto itens (conexo ele já é).
fn ciclo_unico(itens: &[u32], grau_ent: &[usize], grau_sai: &[usize], arestas: usize) -> bool {
    arestas == itens.len() && itens.iter().all(|&i| grau_ent[i as usize] == 1 && grau_sai[i as usize] == 1)
}

struct Escolha {
    /// itens mostrados de cada componente (vazio = não mostrado)
    mostrados: Vec<Vec<u32>>,
    parcial: Option<usize>,
}

/// Os maiores componentes primeiro, até `max` itens; o que não couber
/// inteiro entra em parte (os itens mais perto do seu item mais comprido,
/// por busca em largura) e encerra a escolha.
fn escolhe(lista: &[Vec<u32>], adj: &Csr, comprimento: &[u64], max: usize) -> Escolha {
    let mut resta = max;
    let mut mostrados = Vec::with_capacity(lista.len());
    let mut parcial = None;
    for (ci, itens) in lista.iter().enumerate() {
        if resta == 0 {
            mostrados.push(Vec::new());
        } else if itens.len() <= resta {
            resta -= itens.len();
            mostrados.push(itens.clone());
        } else {
            // o mais comprido (o menor id nos empates)
            let raiz = *itens.iter().max_by(|&&a, &&b| comprimento[a as usize].cmp(&comprimento[b as usize]).then(b.cmp(&a))).unwrap();
            let mut marca: HashMap<u32, ()> = HashMap::with_capacity(resta);
            let mut fila = std::collections::VecDeque::new();
            let mut escolhidos = vec![raiz];
            marca.insert(raiz, ());
            fila.push_back(raiz);
            'busca: while let Some(v) = fila.pop_front() {
                for &u in adj.vizinhos(v as usize) {
                    if escolhidos.len() >= resta {
                        break 'busca;
                    }
                    if marca.insert(u, ()).is_none() {
                        escolhidos.push(u);
                        fila.push_back(u);
                    }
                }
            }
            escolhidos.sort_unstable();
            resta = 0;
            parcial = Some(ci);
            mostrados.push(escolhidos);
        }
    }
    Escolha { mostrados, parcial }
}

// =====================================================================
//  Leiaute
// =====================================================================

/// Segmentos do traçado de um nó: 1 + log2(comprimento / base), de 1 a 16.
fn segmentos(comprimento: u64, base: u64) -> usize {
    let r = comprimento as f64 / base.max(1) as f64;
    if r < 2.0 {
        1
    } else {
        (1 + r.log2().floor() as usize).min(16)
    }
}

fn arred(x: f64, escala: f64) -> f64 {
    (x * escala).round() / escala + 0.0
}

struct Desenho {
    /// traçado de cada item (vazio se não mostrado)
    tracado: Vec<Vec<[f64; 2]>>,
    /// o ponto do meio do traçado
    centro: Vec<[f64; 2]>,
    largura: f64,
    altura: f64,
}

/// Desenha os componentes mostrados e os põe lado a lado.
fn desenha(segs: &[usize], anel: &[bool], arestas: &[(u32, u32)], mostrados: &[Vec<u32>]) -> Desenho {
    let n = segs.len();
    let mut tracado: Vec<Vec<[f64; 2]>> = vec![Vec::new(); n];
    let mut local = vec![NENHUM; n];
    // arestas de cada componente, em índices locais
    let mut por_comp: Vec<Vec<(u32, u32)>> = vec![Vec::new(); mostrados.len()];
    let mut comp_de = vec![NENHUM; n];
    for (ci, itens) in mostrados.iter().enumerate() {
        for (li, &i) in itens.iter().enumerate() {
            local[i as usize] = li as u32;
            comp_de[i as usize] = ci as u32;
        }
    }
    for &(a, b) in arestas {
        let (ca, cb) = (comp_de[a as usize], comp_de[b as usize]);
        if ca != NENHUM && ca == cb {
            por_comp[ca as usize].push((local[a as usize], local[b as usize]));
        }
    }
    let mut caixas: Vec<(usize, f64, f64)> = Vec::new();
    for (ci, itens) in mostrados.iter().enumerate() {
        if itens.is_empty() {
            continue;
        }
        let s: Vec<usize> = itens.iter().map(|&i| segs[i as usize]).collect();
        let mut pts = if itens.len() == 1 {
            if anel[itens[0] as usize] {
                vec![anel_de_arcos(&s, &[0])]
            } else {
                vec![(0..=2 * s[0]).map(|t| [t as f64 * UNIDADE, 0.0]).collect()]
            }
        } else if let Some(ordem) = ordem_do_ciclo(itens.len(), &por_comp[ci]) {
            let arcos = anel_de_arcos(&s, &ordem);
            reparte_anel(arcos, &s, &ordem)
        } else {
            let mut p = leiaute_forcas(&s, &por_comp[ci], 0x5EED_0000 ^ (itens.len() as u64));
            gira_para_horizontal(&mut p);
            p
        };
        // origem no canto da caixa
        let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for q in pts.iter().flatten() {
            x0 = x0.min(q[0]);
            y0 = y0.min(q[1]);
            x1 = x1.max(q[0]);
            y1 = y1.max(q[1]);
        }
        for q in pts.iter_mut().flatten() {
            q[0] -= x0;
            q[1] -= y0;
        }
        for (li, &i) in itens.iter().enumerate() {
            tracado[i as usize] = std::mem::take(&mut pts[li]);
        }
        caixas.push((ci, x1 - x0, y1 - y0));
    }
    // empacotamento em prateleiras: os mais altos primeiro
    let margem = 3.0 * UNIDADE;
    let area: f64 = caixas.iter().map(|c| (c.1 + margem) * (c.2 + margem)).sum();
    let maior_l = caixas.iter().map(|c| c.1).fold(0.0, f64::max);
    let alvo = maior_l.max(area.sqrt() * 1.2);
    caixas.sort_by(|a, b| b.2.total_cmp(&a.2).then(b.1.total_cmp(&a.1)).then(a.0.cmp(&b.0)));
    let (mut x, mut y, mut alt_linha, mut largura) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let mut deslocamento = vec![[0.0f64; 2]; mostrados.len()];
    for &(ci, w, h) in &caixas {
        if x > 0.0 && x + w > alvo {
            x = 0.0;
            y += alt_linha + margem;
            alt_linha = 0.0;
        }
        deslocamento[ci] = [x, y];
        largura = largura.max(x + w);
        alt_linha = alt_linha.max(h);
        x += w + margem;
    }
    let altura = y + alt_linha;
    let mut centro = vec![[0.0f64; 2]; n];
    for i in 0..n {
        if tracado[i].is_empty() {
            continue;
        }
        let d = deslocamento[comp_de[i] as usize];
        for q in tracado[i].iter_mut() {
            q[0] += d[0];
            q[1] += d[1];
        }
        centro[i] = tracado[i][tracado[i].len() / 2];
    }
    Desenho { tracado, centro, largura, altura }
}

/// Se as arestas formam um ciclo dirigido só, a ordem dele a partir do item 0.
fn ordem_do_ciclo(n: usize, arestas: &[(u32, u32)]) -> Option<Vec<usize>> {
    if arestas.len() != n {
        return None;
    }
    let mut prox = vec![NENHUM; n];
    let mut ent = vec![0usize; n];
    for &(a, b) in arestas {
        if prox[a as usize] != NENHUM {
            return None;
        }
        prox[a as usize] = b;
        ent[b as usize] += 1;
    }
    if ent.iter().any(|&e| e != 1) {
        return None;
    }
    let mut ordem = Vec::with_capacity(n);
    let mut v = 0usize;
    loop {
        ordem.push(v);
        v = prox[v] as usize;
        if v == 0 {
            break;
        }
        if ordem.len() > n {
            return None;
        }
    }
    (ordem.len() == n).then_some(ordem)
}

/// Um anel: um arco por item (comprimento = 2 × segmentos × unidade), com um
/// vão de uma unidade entre um e o seguinte. Devolve os pontos de todos os
/// arcos, na ordem dada, concatenados com marcas de onde cada um termina
/// (usado por `reparte_anel`); para um item só, o próprio traçado.
fn anel_de_arcos(segs: &[usize], ordem: &[usize]) -> Vec<[f64; 2]> {
    let circunf: f64 = ordem.iter().map(|&i| (2 * segs[i]) as f64 * UNIDADE + UNIDADE).sum::<f64>().max(6.0 * UNIDADE);
    let r = circunf / std::f64::consts::TAU;
    let mut pts = Vec::new();
    let mut ang = -std::f64::consts::FRAC_PI_2 + 0.5 * UNIDADE / r;
    for &i in ordem {
        let arco = (2 * segs[i]) as f64 * UNIDADE;
        let npts = pontos_do_arco(segs[i]);
        for t in 0..npts {
            let a = ang + (arco / r) * t as f64 / (npts - 1) as f64;
            pts.push([r + r * a.cos(), r + r * a.sin()]);
        }
        ang += (arco + UNIDADE) / r;
    }
    pts
}

/// Pontos de um arco de 2·segs unidades: um a cada meia unidade (ímpar: o
/// do meio é o centro do nó).
fn pontos_do_arco(segs: usize) -> usize {
    4 * segs + 1
}

/// Separa os pontos de `anel_de_arcos` em um traçado por item (índice local).
fn reparte_anel(pts: Vec<[f64; 2]>, segs: &[usize], ordem: &[usize]) -> Vec<Vec<[f64; 2]>> {
    let mut saida = vec![Vec::new(); segs.len()];
    let mut pos = 0;
    for &i in ordem {
        let npts = pontos_do_arco(segs[i]);
        saida[i] = pts[pos..pos + npts].to_vec();
        pos += npts;
    }
    saida
}

/// Gira o desenho para o eixo principal ficar na horizontal.
fn gira_para_horizontal(p: &mut [Vec<[f64; 2]>]) {
    let (mut sx, mut sy, mut cnt) = (0.0, 0.0, 0.0);
    for q in p.iter().flatten() {
        sx += q[0];
        sy += q[1];
        cnt += 1.0;
    }
    if cnt < 2.0 {
        return;
    }
    let (mx, my) = (sx / cnt, sy / cnt);
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    for q in p.iter().flatten() {
        let (dx, dy) = (q[0] - mx, q[1] - my);
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    let teta = 0.5 * (2.0 * sxy).atan2(sxx - syy);
    let (c, s) = (teta.cos(), teta.sin());
    for q in p.iter_mut().flatten() {
        let (dx, dy) = (q[0] - mx, q[1] - my);
        // rotação por −teta
        q[0] = dx * c + dy * s;
        q[1] = -dx * s + dy * c;
    }
}

/// Leiaute por forças de um componente: cada item vira uma corrente de
/// 2·segmentos + 1 contas (ímpar: a do meio é o centro do nó), cada aresta
/// liga a última conta de um à primeira do outro; Fruchterman-Reingold em
/// vários níveis e, no fim, nenhuma conta a menos de 0,4 unidade de outra.
/// Devolve as contas de cada item.
fn leiaute_forcas(segs: &[usize], arestas: &[(u32, u32)], semente: u64) -> Vec<Vec<[f64; 2]>> {
    let mut base = Vec::with_capacity(segs.len() + 1);
    let mut total = 0u32;
    for &s in segs {
        base.push(total);
        total += 2 * s as u32 + 1;
    }
    let mut ligacoes: Vec<(u32, u32)> = Vec::with_capacity(total as usize + arestas.len());
    for (i, &s) in segs.iter().enumerate() {
        for t in 0..2 * s as u32 {
            ligacoes.push((base[i] + t, base[i] + t + 1));
        }
    }
    for &(a, b) in arestas {
        let fim_a = base[a as usize] + 2 * segs[a as usize] as u32;
        ligacoes.push((fim_a, base[b as usize]));
    }
    let mut pos = fr_multinivel(total as usize, &ligacoes, semente);
    separa(&mut pos, 0.4 * UNIDADE);
    segs.iter().enumerate().map(|(i, &s)| pos[base[i] as usize..base[i] as usize + 2 * s + 1].to_vec()).collect()
}

// ---------------------------------------------------------------------
//  Fruchterman-Reingold em vários níveis (Walshaw 2003; passo adaptativo
//  de Hu 2005)
// ---------------------------------------------------------------------

struct Sorteio(u64);

impl Sorteio {
    fn proximo(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn real(&mut self) -> f64 {
        (self.proximo() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Vizinhança sem direção, sem laços nem repetição.
struct Csr {
    inicio: Vec<u32>,
    viz: Vec<u32>,
}

impl Csr {
    fn novo(n: usize, arestas: &[(u32, u32)]) -> Csr {
        let mut pares: Vec<(u32, u32)> = Vec::with_capacity(arestas.len() * 2);
        for &(a, b) in arestas {
            if a != b {
                pares.push((a, b));
                pares.push((b, a));
            }
        }
        pares.sort_unstable();
        pares.dedup();
        let mut inicio = vec![0u32; n + 1];
        for &(a, _) in &pares {
            inicio[a as usize + 1] += 1;
        }
        for i in 0..n {
            inicio[i + 1] += inicio[i];
        }
        Csr { inicio, viz: pares.into_iter().map(|(_, b)| b).collect() }
    }
    fn n(&self) -> usize {
        self.inicio.len() - 1
    }
    fn vizinhos(&self, v: usize) -> &[u32] {
        &self.viz[self.inicio[v] as usize..self.inicio[v + 1] as usize]
    }
}

/// Emparelha cada vértice com o vizinho livre de menor massa (os mais
/// leves primeiro). Devolve o mapa fino → grosso e o nº de vértices grossos.
fn emparelha(g: &Csr, massa: &[u32]) -> (Vec<u32>, usize) {
    let n = g.n();
    let mut ordem: Vec<u32> = (0..n as u32).collect();
    ordem.sort_by_key(|&v| (massa[v as usize], g.vizinhos(v as usize).len() as u32, v));
    let mut mapa = vec![NENHUM; n];
    let mut nc = 0u32;
    for &v in &ordem {
        let v = v as usize;
        if mapa[v] != NENHUM {
            continue;
        }
        let mut melhor: Option<(u32, u32)> = None;
        for &u in g.vizinhos(v) {
            if mapa[u as usize] == NENHUM {
                let c = (massa[u as usize], u);
                if melhor.is_none_or(|m| c < m) {
                    melhor = Some(c);
                }
            }
        }
        mapa[v] = nc;
        if let Some((_, u)) = melhor {
            mapa[u as usize] = nc;
        }
        nc += 1;
    }
    (mapa, nc as usize)
}

fn fr_multinivel(n: usize, arestas: &[(u32, u32)], semente: u64) -> Vec<[f64; 2]> {
    let mut grafos = vec![Csr::novo(n, arestas)];
    let mut massas = vec![vec![1u32; n]];
    let mut mapas: Vec<Vec<u32>> = Vec::new();
    loop {
        let g = grafos.last().unwrap();
        let nl = g.n();
        if nl <= 40 {
            break;
        }
        let (mapa, nc) = emparelha(g, massas.last().unwrap());
        if nc as f64 > 0.9 * nl as f64 {
            break;
        }
        let mut ar: Vec<(u32, u32)> = Vec::with_capacity(g.viz.len() / 2);
        for v in 0..nl {
            for &u in g.vizinhos(v) {
                let (a, b) = (mapa[v], mapa[u as usize]);
                if a < b {
                    ar.push((a, b));
                }
            }
        }
        let mut mc = vec![0u32; nc];
        for v in 0..nl {
            mc[mapa[v] as usize] += massas.last().unwrap()[v];
        }
        grafos.push(Csr::novo(nc, &ar));
        massas.push(mc);
        mapas.push(mapa);
    }
    let k_de = |l: usize| UNIDADE * (n as f64 / grafos[l].n() as f64).sqrt();
    let mut rng = Sorteio(semente);
    let topo = grafos.len() - 1;
    let nc = grafos[topo].n();
    let lado = k_de(topo) * (nc as f64).sqrt();
    let mut pos: Vec<[f64; 2]> = (0..nc).map(|_| [rng.real() * lado, rng.real() * lado]).collect();
    fr(&mut pos, &grafos[topo], k_de(topo), 300, lado * 0.25);
    for l in (0..topo).rev() {
        let kl = k_de(l);
        let mapa = &mapas[l];
        pos = (0..grafos[l].n())
            .map(|v| {
                let p = pos[mapa[v] as usize];
                let a = rng.real() * std::f64::consts::TAU;
                [p[0] + 0.3 * kl * a.cos(), p[1] + 0.3 * kl * a.sin()]
            })
            .collect();
        fr(&mut pos, &grafos[l], kl, iteracoes_do_nivel(l, grafos[l].n()), 0.75 * kl);
    }
    pos
}

/// Quantas iterações em cada nível: muitas nos grafos pequenos (custa pouco
/// e o desenho fica mais liso), poucas nos grandes (o formato já vem dos
/// níveis grossos; medido em 2026-11-16: 5.000 unitigs de leituras, 30/25
/// iterações 0,83 s de CPU, 15/12 0,36 s, desenho quase igual).
fn iteracoes_do_nivel(nivel: usize, n: usize) -> usize {
    match (nivel, n) {
        (_, 0..=400) => 60,
        (0, 401..=4000) => 40,
        (0, _) => 15,
        (_, 401..=2000) => 30,
        _ => 12,
    }
}

/// Grade de células quadradas (lado = o alcance da força) para achar os
/// vizinhos próximos: os itens ordenados por célula e uma tabela de
/// espalhamento célula → faixa. Só as células ocupadas existem, então um
/// desenho comprido não custa memória nem células maiores.
struct Grade {
    lado: f64,
    x0: f64,
    y0: f64,
    chaves: Vec<(u64, u32)>,
    /// (célula, início, fim) com endereçamento aberto; célula u64::MAX = vazio
    tabela: Vec<(u64, u32, u32)>,
    mascara: usize,
}

fn celula_de(p: &[f64; 2], lado: f64, x0: f64, y0: f64) -> (u64, u64) {
    let f = |v: f64| ((v / lado).floor().clamp(0.0, (1u64 << 30) as f64) as u64) + 1;
    (f(p[0] - x0), f(p[1] - y0))
}

fn chave_celula(cx: u64, cy: u64) -> u64 {
    (cx << 32) | cy
}

fn espalha(x: u64) -> usize {
    let mut z = x.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 29)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    (z ^ (z >> 32)) as usize
}

impl Grade {
    fn nova(pos: &[[f64; 2]], lado: f64) -> Grade {
        let mut g = Grade { lado, x0: 0.0, y0: 0.0, chaves: Vec::new(), tabela: Vec::new(), mascara: 0 };
        g.refaz(pos);
        g
    }
    /// Remonta a grade para as posições novas, aproveitando a memória.
    fn refaz(&mut self, pos: &[[f64; 2]]) {
        self.x0 = pos.iter().map(|p| p[0]).fold(f64::MAX, f64::min);
        self.y0 = pos.iter().map(|p| p[1]).fold(f64::MAX, f64::min);
        let cap = (2 * pos.len()).next_power_of_two().max(16);
        self.mascara = cap - 1;
        self.tabela.clear();
        self.tabela.resize(cap, (u64::MAX, 0, 0));
        let mut chaves = std::mem::take(&mut self.chaves);
        chaves.clear();
        let (lado, x0, y0) = (self.lado, self.x0, self.y0);
        chaves.par_extend(pos.par_iter().enumerate().map(|(i, p)| {
            let (cx, cy) = celula_de(p, lado, x0, y0);
            (chave_celula(cx, cy), i as u32)
        }));
        chaves.par_sort_unstable();
        let mut i = 0;
        while i < chaves.len() {
            let mut j = i + 1;
            while j < chaves.len() && chaves[j].0 == chaves[i].0 {
                j += 1;
            }
            let mut h = espalha(chaves[i].0) & self.mascara;
            while self.tabela[h].0 != u64::MAX {
                h = (h + 1) & self.mascara;
            }
            self.tabela[h] = (chaves[i].0, i as u32, j as u32);
            i = j;
        }
        self.chaves = chaves;
    }
    /// A célula de um ponto (contada a partir de 1, para a vizinha −1 existir).
    fn celula(&self, p: &[f64; 2]) -> (u64, u64) {
        celula_de(p, self.lado, self.x0, self.y0)
    }
    /// Chama `f` para cada item das 9 células em volta de `p` (a dele incluída).
    fn vizinhos(&self, p: &[f64; 2], mut f: impl FnMut(usize)) {
        let (cx, cy) = self.celula(p);
        for y in (cy - 1)..=(cy + 1) {
            for x in (cx - 1)..=(cx + 1) {
                let chave = chave_celula(x, y);
                let mut h = espalha(chave) & self.mascara;
                loop {
                    let e = self.tabela[h];
                    if e.0 == chave {
                        for &(_, j) in &self.chaves[e.1 as usize..e.2 as usize] {
                            f(j as usize);
                        }
                        break;
                    }
                    if e.0 == u64::MAX {
                        break;
                    }
                    h = (h + 1) & self.mascara;
                }
            }
        }
    }
}

/// Direção determinística para dois pontos coincidentes.
fn direcao(i: usize, j: usize) -> [f64; 2] {
    let a = Sorteio(((i as u64) << 32) ^ j as u64 ^ 0xA5A5).real() * std::f64::consts::TAU;
    [a.cos(), a.sin()]
}

/// Iterações de Fruchterman-Reingold: repulsão K²/d (todos os pares até 400
/// pontos; acima, só até 2 K, pela grade, como no FR original), atração d²/K nas arestas;
/// passo limitado e adaptado à energia (Hu 2005). A força de cada ponto é
/// somada em paralelo, cada uma na sua ordem fixa: o resultado não depende
/// do número de threads.
fn fr(pos: &mut [[f64; 2]], g: &Csr, k: f64, iteracoes: usize, passo0: f64) {
    let n = pos.len();
    if n < 2 {
        return;
    }
    let exato = n <= 400;
    let k2 = k * k;
    let corte = 2.0 * k;
    let corte2 = corte * corte;
    let minimo2 = 1e-6 * k2;
    let mut f = vec![[0.0f64; 2]; n];
    let mut passo = passo0;
    let mut energia_ant = f64::INFINITY;
    let mut progresso = 0;
    let mut grade = if exato { None } else { Some(Grade::nova(pos, corte)) };
    for it in 0..iteracoes {
        if it > 0 {
            if let Some(gr) = grade.as_mut() {
                gr.refaz(pos);
            }
        }
        {
            let p: &[[f64; 2]] = pos;
            f.par_iter_mut().enumerate().with_min_len(512).for_each(|(i, fi)| {
                let (xi, yi) = (p[i][0], p[i][1]);
                let mut acc = [0.0f64; 2];
                let mut repele = |j: usize| {
                    if j == i {
                        return;
                    }
                    let (mut dx, mut dy) = (xi - p[j][0], yi - p[j][1]);
                    let mut d2 = dx * dx + dy * dy;
                    if d2 < minimo2 {
                        // pontos coincidentes: uma direção fixa, oposta para o outro
                        let u = direcao(i.min(j), i.max(j));
                        let sinal = if i < j { 1e-3 * k } else { -1e-3 * k };
                        dx = u[0] * sinal;
                        dy = u[1] * sinal;
                        d2 = dx * dx + dy * dy;
                    }
                    if !exato && d2 > corte2 {
                        return;
                    }
                    let s = k2 / d2;
                    acc[0] += dx * s;
                    acc[1] += dy * s;
                };
                match &grade {
                    None => (0..n).for_each(&mut repele),
                    Some(gr) => gr.vizinhos(&p[i], &mut repele),
                }
                for &j in g.vizinhos(i) {
                    let (dx, dy) = (p[j as usize][0] - xi, p[j as usize][1] - yi);
                    let s = (dx * dx + dy * dy).sqrt() / k;
                    acc[0] += dx * s;
                    acc[1] += dy * s;
                }
                *fi = acc;
            });
        }
        let mut energia = 0.0;
        for i in 0..n {
            let m2 = f[i][0] * f[i][0] + f[i][1] * f[i][1];
            energia += m2;
            let m = m2.sqrt();
            if m > 1e-12 {
                let d = m.min(passo) / m;
                pos[i][0] += f[i][0] * d;
                pos[i][1] += f[i][1] * d;
            }
        }
        if energia < energia_ant {
            progresso += 1;
            if progresso >= 5 {
                progresso = 0;
                passo /= 0.9;
            }
        } else {
            progresso = 0;
            passo *= 0.9;
        }
        energia_ant = energia;
        if passo < 0.005 * k {
            break;
        }
    }
}

/// Afasta os pares de pontos a menos de `dmin` (algumas rodadas).
fn separa(pos: &mut [[f64; 2]], dmin: f64) {
    let n = pos.len();
    if n < 2 {
        return;
    }
    let alvo2 = (dmin * 0.999) * (dmin * 0.999);
    let mut grade = Grade::nova(pos, dmin);
    for rodada in 0..30 {
        if rodada > 0 {
            grade.refaz(pos);
        }
        let mut mexeu = false;
        let mut perto: Vec<usize> = Vec::new();
        for i in 0..n {
            let p = pos[i];
            perto.clear();
            grade.vizinhos(&p, |j| {
                if j > i {
                    perto.push(j);
                }
            });
            for &j in &perto {
                let (dx, dy) = (pos[j][0] - pos[i][0], pos[j][1] - pos[i][1]);
                let d2 = dx * dx + dy * dy;
                if d2 >= alvo2 {
                    continue;
                }
                let d = d2.sqrt();
                // a direção unitária de i para j e quanto cada um anda
                let (ux, uy, empurra) = if d < 1e-9 {
                    let u = direcao(i, j);
                    (u[0], u[1], dmin / 2.0)
                } else {
                    (dx / d, dy / d, (dmin - d) / 2.0)
                };
                pos[i][0] -= ux * empurra;
                pos[i][1] -= uy * empurra;
                pos[j][0] += ux * empurra;
                pos[j][1] += uy * empurra;
                mexeu = true;
            }
        }
        if !mexeu {
            break;
        }
    }
}

// =====================================================================
//  Auxiliares
// =====================================================================

fn n50(desc: &[u64]) -> u64 {
    let total: u64 = desc.iter().sum();
    let mut acum = 0u64;
    for &l in desc {
        acum += l;
        if 2 * acum >= total {
            return l;
        }
    }
    0
}

fn milhar(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn conta(n: u64, um: &str, varios: &str) -> String {
    format!("{} {}", milhar(n), if n == 1 { um } else { varios })
}
