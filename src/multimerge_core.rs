// COPIA de adaptive-parallel-multimerge-sort/src/multimerge.rs.
// O original e a fonte da verdade e guarda os testes do motor; nao edite
// esta copia -- altere la, rode os testes la e copie de novo.

// =====================================================================
//  ADAPTIVE PARALLEL MULTIMERGE SORT -- motor Rust
//
//  Dois motores sobre a mesma deteccao de runs e a mesma arvore de merge:
//
//    multi_merge_sort          NAO estavel. E o padrao.
//    multi_merge_sort_stable   estavel. Para quando a ordem relativa de chaves
//                              iguais importa (SGBD, ordenacao em cascata).
//
//  Cada um tem as variantes `_by` (comparador) e `_by_key` (extrator de chave).
//
//  sort(arr)
//    |- n pequeno .......................... sort da std
//    |- ESCUDO (janelas espalhadas) ........ caotico -> caminho caotico
//    '- DETECCAO de runs (monoide, paralela)
//         |- runs curtos consecutivos ...... viram UM bloco, ordenado pela std
//         |- runs descendentes ............. revertidos no lugar
//         |- runs ja em ordem global ....... fundidos sem mover dado
//         '- o que sobrar .................. arvore de merge k-vias
//
//  ONDE OS DOIS MOTORES DIFEREM
//
//  1. Run descendente com chaves repetidas. Um trecho 9 9 7 7 3 e descendente
//     (nao-crescente), mas inverte-lo troca a ordem dos iguais. O instavel
//     inverte e pronto. O estavel inverte e depois desinverte cada bloco de
//     iguais (`restore_ties`), o que devolve a ordem original entre eles. Os
//     dois resolvem o caso em O(n); antes o motor quebrava o trecho em um run
//     por plato e pagava a arvore inteira.
//  2. Bloco sem estrutura. O estavel usa o `sort` da std (driftsort), o
//     instavel usa o `sort_unstable` (ipnsort).
//  3. Dado caotico. O estavel ordena em blocos e funde pela arvore; o instavel
//     delega ao quicksort paralelo, que trabalha no lugar e nao aloca buffer.
//
//  O motor nao tem `unsafe`. `T: Copy` e o que permite ao borrow checker
//  provar o merge bidirecional: as duas threads leem as fontes por referencia
//  compartilhada e escrevem em metades disjuntas do destino.
// =====================================================================

use rayon::prelude::*;
use std::cmp::Ordering::{self, Equal, Greater, Less};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

/// Fan-out of the k-way merge.
pub const KWAY_FANOUT: usize = 8;

// =====================================================================
//  API PUBLICA
// =====================================================================

/// Sorts the slice. **Not stable**: equal elements may be reordered.
///
/// This is the default engine. It detects pre-existing order (ascending runs,
/// descending runs, blocks already in global order) and only pays for the
/// disorder that is actually there. Use [`multi_merge_sort_stable`] when the
/// relative order of equal keys must be preserved.
///
/// Allocates a buffer of `arr.len()` elements only when runs have to be
/// merged. Chaotic input is sorted in place.
pub fn multi_merge_sort<T: Ord + Copy + Send + Sync>(arr: &mut [T]) {
    sort_with(arr, &T::cmp, false, &Tuning::for_type::<T>(), None);
}

/// [`multi_merge_sort`] with a comparator.
pub fn multi_merge_sort_by<T, F>(arr: &mut [T], compare: F)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with(arr, &compare, false, &Tuning::for_type::<T>(), None);
}

/// [`multi_merge_sort`] with a key extractor.
pub fn multi_merge_sort_by_key<T, K, F>(arr: &mut [T], key: F)
where
    T: Copy + Send + Sync,
    K: Ord,
    F: Fn(&T) -> K + Sync,
{
    sort_with(arr, &|a: &T, b: &T| key(a).cmp(&key(b)), false, &Tuning::for_type::<T>(), None);
}

/// Sorts the slice, **preserving the relative order of equal elements**.
///
/// Same adaptive machinery as [`multi_merge_sort`]. Always allocates a buffer
/// of `arr.len()` elements unless the input turns out to be a single run.
pub fn multi_merge_sort_stable<T: Ord + Copy + Send + Sync>(arr: &mut [T]) {
    sort_with(arr, &T::cmp, true, &Tuning::for_type::<T>(), None);
}

/// [`multi_merge_sort_stable`] with a comparator.
pub fn multi_merge_sort_stable_by<T, F>(arr: &mut [T], compare: F)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with(arr, &compare, true, &Tuning::for_type::<T>(), None);
}

/// [`multi_merge_sort_stable`] with a key extractor.
pub fn multi_merge_sort_stable_by_key<T, K, F>(arr: &mut [T], key: F)
where
    T: Copy + Send + Sync,
    K: Ord,
    F: Fn(&T) -> K + Sync,
{
    sort_with(arr, &|a: &T, b: &T| key(a).cmp(&key(b)), true, &Tuning::for_type::<T>(), None);
}

/// [`multi_merge_sort_by`] reusing the caller's working memory.
///
/// `scratch` is grown to `arr.len()` on first use and kept between calls. A
/// server that sorts repeatedly pays for the buffer once: no allocation and no
/// fresh pages to fault in on later calls. Its contents are unspecified
/// afterwards. Pass `T::cmp` as the comparator for the natural order.
pub fn multi_merge_sort_by_with_buffer<T, F>(arr: &mut [T], scratch: &mut Vec<T>, compare: F)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with(arr, &compare, false, &Tuning::for_type::<T>(), Some(scratch));
}

/// [`multi_merge_sort_stable_by`] reusing the caller's working memory. See
/// [`multi_merge_sort_by_with_buffer`].
pub fn multi_merge_sort_stable_by_with_buffer<T, F>(
    arr: &mut [T],
    scratch: &mut Vec<T>,
    compare: F,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with(arr, &compare, true, &Tuning::for_type::<T>(), Some(scratch));
}

/// Run decomposition of `arr` as signed lengths: positive for a non-decreasing
/// run, negative for a non-increasing one.
///
/// A negative run may contain repeated keys. Reversing it sorts it, but only
/// [`multi_merge_sort_stable`] restores the original order among the equals.
pub fn detect_global_trend<T: Ord + Sync>(arr: &[T]) -> Vec<i64> {
    let t = Tuning::for_type::<T>();
    detect_runs(arr, &T::cmp, t.detect_block, t.min_run, usize::MAX)
        .unwrap_or_default()
        .iter()
        .map(|r| if r.is_desc() { -(r.len() as i64) } else { r.len() as i64 })
        .collect()
}

// =====================================================================
//  AJUSTES
//
//  Um unico lugar para os limiares. `pub(crate)` para que os testes forcem
//  recursao profunda com dados pequenos, e para que uma sonda de medicao
//  varra valores sem recompilar o motor.
// =====================================================================

// `Structured` e `Chaotic` so sao construidos pelos testes e pela sonda de
// medicao: forcam um caminho para exercita-lo ou cronometra-lo isolado.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Route {
    /// O escudo decide.
    Auto,
    /// Forca deteccao de runs + arvore.
    Structured,
    /// Forca o caminho caotico.
    Chaotic,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Tuning {
    /// Ate aqui o sort da std resolve sozinho.
    pub small_sort: usize,
    /// Bloco da deteccao de runs (com 1 elemento de sobreposicao).
    pub detect_block: usize,
    /// Run mais curto que isto nao vale um no da arvore: runs curtos
    /// consecutivos sao agrupados e ordenados como um bloco.
    pub min_run: usize,
    /// A arvore so entra quando, depois de encadear, os runs tem em media
    /// pelo menos este comprimento.
    ///
    /// Em u64 a arvore sobre runs naturais e o caminho de blocos empatam de
    /// 256 ate 16384 elementos por run, e os blocos ganham de 8% a 15% entre
    /// 32 e 64. Com chave cara (`&str`) a arvore ganha com folga mesmo em run
    /// curto, porque reordenar o bloco repete comparacoes que os runs ja
    /// tinham resolvido: ordenado com 1% de trocas (run medio ~50), 367 ms
    /// contra 520; 64 runs de 31 mil strings, 175 contra 262. O valor fica
    /// baixo para proteger a chave cara, que e onde o erro custa mais.
    pub struct_min_avg_run: usize,
    /// O mesmo limiar quando a comparacao e barata. Ali o custo e mover
    /// memoria, nao comparar, e os blocos empatam com a arvore de 256 para
    /// cima e ganham abaixo: ordenado com 1% de trocas em 10M u64, 128 ms nos
    /// blocos contra 146 na arvore.
    pub struct_min_avg_run_cheap: usize,
    /// Piso do numero de segmentos aceitos, para arrays pequenos.
    pub min_segs: usize,
    /// Folga do aborto antecipado da deteccao.
    pub abort_slack: usize,
    /// Tamanho do bloco ordenado pela std (caminho caotico estavel e teto do
    /// agrupamento de runs curtos).
    pub chunk: usize,
    /// Grao paralelo do merge: acima disto corta por co-rank / multisequencia.
    pub par_grain: usize,
    /// Acima disto um pedaco de merge binario roda bidirecional (2 tarefas).
    pub bidir_min: usize,
    /// No da arvore binaria abaixo disto recursa sem `rayon::join`.
    pub seq_node: usize,
    /// Subarvore ate este tamanho usa merge binario; acima, k-vias.
    pub kway_min_elems: usize,
    /// Tile do merge k-vias, em elementos.
    pub tile: usize,
    /// Escudo: acima deste percentual de inversoes de direcao os runs sao
    /// curtos demais para a deteccao valer a passada.
    pub chaos_percent: usize,
    /// Escudo: acima deste percentual o dado e aleatorio de fato, nao apenas
    /// de runs curtos. So interessa ao motor instavel.
    pub random_percent: usize,
    /// Comparacao "barata": ate este custo medio, em nanossegundos, medido na
    /// amostra do escudo.
    pub cheap_ns: usize,
    /// Fixa o veredito de custo em vez de medir. Para os testes.
    pub force_cheap: Option<bool>,
    pub route: Route,
}

impl Tuning {
    pub(crate) fn for_type<T>() -> Self {
        let e = std::mem::size_of::<T>().max(1);
        Tuning {
            small_sort: 4096,
            detect_block: 32_768,
            min_run: 16,
            struct_min_avg_run: 48,
            struct_min_avg_run_cheap: 256,
            min_segs: 8,
            abort_slack: 64,
            chunk: 65_536,
            par_grain: 65_536,
            bidir_min: 4096,
            seq_node: 2048,
            kway_min_elems: ((2usize << 20) / e).max(2),
            tile: (32_768 / (8 * e)).clamp(128, 2048),
            chaos_percent: 9,
            random_percent: 35,
            cheap_ns: 10,
            force_cheap: None,
            route: Route::Auto,
        }
    }

    /// Quantos segmentos um array de `n` elementos pode ter para a arvore
    /// ainda compensar.
    fn allowed_segs(&self, n: usize, cheap: bool) -> usize {
        let avg = if cheap { self.struct_min_avg_run_cheap } else { self.struct_min_avg_run };
        (n / avg.max(1)).max(self.min_segs)
    }
}

// =====================================================================
//  ENTRADA
// =====================================================================

pub(crate) fn sort_with<T, F>(
    arr: &mut [T],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let n = arr.len();
    if n < 2 {
        return;
    }
    if n <= t.small_sort {
        leaf_sort(arr, cmp, stable);
        return;
    }
    let probe = match t.route {
        Route::Chaotic => Probe { rate: 100, cheap: false },
        Route::Structured => Probe { rate: 0, cheap: false },
        Route::Auto => sample(arr, cmp, t.cheap_ns),
    };
    let rate = probe.rate;
    let cheap = t.force_cheap.unwrap_or(probe.cheap);
    if rate > t.chaos_percent {
        let kind = if rate > t.random_percent { Chaos::Random } else { Chaos::ShortRuns };
        sort_chaotic(arr, cmp, stable, t, scratch, kind, cheap);
        return;
    }
    // ATALHO DO INSTAVEL: dado globalmente ascendente, com desordem espalhada
    // e comparacao barata, e o melhor caso do quicksort paralelo -- ordenado
    // com 1% de trocas em 10M u64 sai em 55 ms, contra 113-128 de qualquer
    // caminho de merge. Exige que haja desordem na amostra (senao a deteccao
    // resolve em uma passada) e que a comparacao seja barata (com `&str` o
    // mesmo dado leva 733 ms nele contra 367 da arvore).
    if !stable
        && t.route == Route::Auto
        && cheap
        && rate >= 1
        && globally_ascending(arr, cmp)
    {
        arr.par_sort_unstable_by(|a, b| cmp(a, b));
        return;
    }
    // A deteccao desiste sozinha quando o dado se mostra fragmentado demais
    // para a arvore: o custo de ter perguntado fica limitado a uma fracao da
    // varredura, em vez da passada inteira.
    let limit = match t.route {
        Route::Auto => 2 * t.allowed_segs(n, cheap) + t.abort_slack,
        _ => usize::MAX,
    };
    match detect_runs(arr, cmp, t.detect_block, t.min_run, limit) {
        Some(runs) => sort_structured(arr, &runs, cmp, stable, t, scratch, cheap),
        None => sort_chaotic(arr, cmp, stable, t, scratch, Chaos::Fragmented, cheap),
    }
}

#[inline]
fn leaf_sort<T, F>(s: &mut [T], cmp: &F, stable: bool)
where
    F: Fn(&T, &T) -> Ordering,
{
    if stable {
        s.sort_by(|a, b| cmp(a, b));
    } else {
        s.sort_unstable_by(|a, b| cmp(a, b));
    }
}

// =====================================================================
//  FASE 0: ESCUDO
//
//  Oito janelas espalhadas pelo array inteiro, nao uma janela no meio: um
//  ponto unico deixava 0,004% do dado decidir o caminho dos outros 99,996%.
//
//  Empates NAO contam como direcao. 4 4 3 3 2 2 e um run descendente perfeito,
//  e contar `<=` como "subindo" o classificava como caos puro.
//
//  O veredito deixou de ser critico: errar para "estruturado" custa uma
//  passada de deteccao, porque runs curtos sao agrupados e ordenados em bloco
//  de qualquer forma; errar para "caotico" custa a mesma passada, porque os
//  blocos ja ordenados sao fundidos sem mover dado.
// =====================================================================

/// O que a amostra do escudo diz sobre o dado.
struct Probe {
    /// Percentual de comparacoes em que a direcao inverteu. Aleatorio de alta
    /// cardinalidade da ~66; de duas chaves, ~50; um dente de serra de
    /// comprimento L da 200/L. Ordenado, invertido e platos dao zero.
    rate: usize,
    /// A comparacao e barata? Medido no relogio, sobre as mesmas 512
    /// comparacoes: um inteiro custa 1 a 3 ns, uma string fora de cache passa
    /// de 50. O erro so tem uma direcao possivel -- uma preempcao no meio da
    /// amostra faz o barato parecer caro -- e essa e a direcao conservadora:
    /// o motor deixa de usar um atalho, nunca usa um que nao devia.
    cheap: bool,
}

fn sample<T, F>(arr: &[T], cmp: &F, cheap_ns: usize) -> Probe
where
    F: Fn(&T, &T) -> Ordering,
{
    const WINDOWS: usize = 8;
    const SPAN: usize = 64;
    let n = arr.len();
    if n < 3 {
        return Probe { rate: 0, cheap: false };
    }
    let clock = std::time::Instant::now();
    let span = SPAN.min(n - 1);
    let last_start = n - 1 - span;
    let (mut seen, mut changes) = (0usize, 0usize);
    for w in 0..WINDOWS {
        let start = last_start / (WINDOWS - 1) * w;
        let mut dir = Equal;
        for i in start..start + span {
            let o = cmp(&arr[i], &arr[i + 1]);
            if o != Equal {
                if dir != Equal && o != dir {
                    changes += 1;
                }
                dir = o;
            }
        }
        seen += span;
        if last_start == 0 {
            break;
        }
    }
    let nanos = clock.elapsed().as_nanos() as usize;
    Probe {
        rate: changes * 100 / seen.max(1),
        cheap: nanos <= cheap_ns * seen,
    }
}

/// O array e globalmente ascendente, a menos de elementos fora do lugar?
///
/// 256 amostras igualmente espacadas, em 16 grupos; a MEDIANA de cada grupo
/// tem que nao decrescer. A mediana ignora elementos isolados fora do lugar
/// (ordenado com 1% de trocas passa), e runs intercalados reprovam: a cada
/// troca de run o nivel cai e fica la embaixo.
///
/// Serve para uma decisao so: se o motor instavel pode entregar um dado
/// fragmentado ao quicksort. Ver `sort_chaotic`.
fn globally_ascending<T, F>(arr: &[T], cmp: &F) -> bool
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    const GROUPS: usize = 16;
    const PER: usize = 16;
    let n = arr.len();
    if n < GROUPS * PER * 4 {
        return false;
    }
    let step = n / (GROUPS * PER);
    let mut prev: Option<T> = None;
    for g in 0..GROUPS {
        let mut s = [arr[0]; PER];
        for (k, slot) in s.iter_mut().enumerate() {
            *slot = arr[(g * PER + k) * step];
        }
        s.sort_unstable_by(|a, b| cmp(a, b));
        let med = s[PER / 2];
        if let Some(p) = prev {
            if cmp(&p, &med) == Greater {
                return false;
            }
        }
        prev = Some(med);
    }
    true
}

// =====================================================================
//  FASE 1: DETECCAO DE RUNS
//
//  Um run e descrito por (tipo, comprimento) empacotado em 64 bits.
//
//    ASC        nao-decrescente
//    FLAT       todos iguais -- neutro: serve tanto a um vizinho ascendente
//               quanto a um descendente
//    DESC       estritamente decrescente
//    DESC_TIES  nao-crescente COM chaves repetidas
//
//  A direcao de um run e decidida pelo PRIMEIRO par estritamente diferente,
//  nao pelo primeiro par. Comecar por `<=` classificava todo plato como
//  ascendente, e um array invertido com chaves repetidas nunca era visto como
//  invertido: virava um run por plato.
// =====================================================================

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Run(u64);

const ASC: u64 = 0;
const FLAT: u64 = 1;
const DESC: u64 = 2;
const DESC_TIES: u64 = 3;

impl Run {
    #[inline]
    fn new(kind: u64, len: usize) -> Run {
        Run((kind << 62) | len as u64)
    }
    #[inline]
    fn kind(self) -> u64 {
        self.0 >> 62
    }
    #[inline]
    fn len(self) -> usize {
        (self.0 & ((1u64 << 62) - 1)) as usize
    }
    #[inline]
    fn is_desc(self) -> bool {
        self.kind() >= DESC
    }
}

/// Decompoe `arr` em runs, acrescentando a `out`.
///
/// Devolve uma estimativa de quantos SEGMENTOS este trecho custaria a arvore:
/// um trecho de runs curtos conta uma vez (vira um bloco ordenado), e um run
/// longo que encadeia com o anterior nao conta (vai ser fundido sem mover
/// dado). E essa contagem que alimenta o aborto antecipado da deteccao.
fn scan_runs<T, F>(arr: &[T], cmp: &F, out: &mut Vec<Run>, min_run: usize) -> usize
where
    F: Fn(&T, &T) -> Ordering,
{
    let n = arr.len();
    let first = out.len();
    let mut head = 0usize;
    let mut count = SegCount { segs: 0, have_prev: false, prev_long: false, prev_max: 0, in_short: false };
    while head < n {
        // Plato inicial: anda ate o primeiro par estritamente diferente.
        let mut tail = head + 1;
        let mut ord = Equal;
        while tail < n {
            ord = cmp(&arr[tail - 1], &arr[tail]);
            if ord != Equal {
                break;
            }
            tail += 1;
        }
        if tail == n {
            out.push(Run::new(FLAT, n - head));
            count.add(arr, cmp, head, n, false, min_run);
            break;
        }
        let plateau = tail - head > 1;
        tail += 1;
        if ord == Less {
            while tail < n && cmp(&arr[tail - 1], &arr[tail]) != Greater {
                tail += 1;
            }
            out.push(Run::new(ASC, tail - head));
            count.add(arr, cmp, head, tail, false, min_run);
        } else {
            let mut ties = plateau;
            while tail < n {
                match cmp(&arr[tail - 1], &arr[tail]) {
                    Less => break,
                    Equal => ties = true,
                    Greater => {}
                }
                tail += 1;
            }
            // ELEMENTO ROUBADO. A varredura gulosa entrega ao run ascendente
            // anterior o primeiro elemento deste run descendente -- que e o
            // MAXIMO dele. Em ... 3 9 | 8 7 6 5 4 o 9 fica a esquerda, e ai o
            // run ascendente termina em 9 e nunca encadeia com o descendente,
            // cujo minimo e 4. Devolver o 9 deixa 3 <= 4 e os dois viram um
            // run so. Sem isto o dente de serra em faixas disjuntas ficava com
            // um run a cada troca de direcao (medido: 1,25x de ganho contra
            // 11x com a devolucao).
            let mut start = head;
            if out.len() > first {
                let last = out.len() - 1;
                let prev = out[last];
                if prev.kind() == ASC
                    && prev.len() >= 2
                    && cmp(&arr[head - 2], &arr[tail - 1]) != Greater
                {
                    let l = prev.len() - 1;
                    out[last] = Run::new(if l == 1 { FLAT } else { ASC }, l);
                    start = head - 1;
                    count.prev_max = head - 2;
                }
            }
            out.push(Run::new(if ties { DESC_TIES } else { DESC }, tail - start));
            count.add(arr, cmp, start, tail, true, min_run);
        }
        head = tail;
    }
    count.segs
}

/// Contabilidade de segmentos de um bloco, para o aborto antecipado.
struct SegCount {
    segs: usize,
    have_prev: bool,
    /// O run anterior era longo?
    prev_long: bool,
    /// Indice do MAIOR elemento do run anterior.
    prev_max: usize,
    /// Estamos dentro de um trecho de runs curtos (ja contado)?
    in_short: bool,
}

impl SegCount {
    /// Registra o run `arr[start..end]`.
    fn add<T, F>(&mut self, arr: &[T], cmp: &F, start: usize, end: usize, desc: bool, min_run: usize)
    where
        F: Fn(&T, &T) -> Ordering,
    {
        let (min_i, max_i) = if desc { (end - 1, start) } else { (start, end - 1) };
        if end - start >= min_run {
            // O primeiro run do bloco nao conta: nao da para saber daqui se
            // ele encadeia com o bloco anterior.
            let chained = !self.have_prev
                || (self.prev_long && cmp(&arr[self.prev_max], &arr[min_i]) != Greater);
            if !chained {
                self.segs += 1;
            }
            self.prev_long = true;
            self.in_short = false;
        } else {
            if !self.in_short {
                self.segs += 1;
                self.in_short = true;
            }
            self.prev_long = false;
        }
        self.have_prev = true;
        self.prev_max = max_i;
    }
}

/// Costura dois runs que COMPARTILHAM um elemento: o ultimo de `l` e o
/// primeiro de `r`. Opera so sobre o metadata, sem ler o array -- e a operacao
/// do monoide. Devolve o run da esquerda e, se sobrar, o da direita.
///
/// A associatividade vale a menos de decomposicoes equivalentes: agrupamentos
/// diferentes podem cortar um elemento isolado em lugares diferentes, e todas
/// as respostas descrevem o mesmo array corretamente.
fn join_runs(l: Run, r: Run) -> (Run, Option<Run>) {
    let (lk, ll, rk, rl) = (l.kind(), l.len(), r.kind(), r.len());
    let kind = match (lk, rk) {
        (FLAT, FLAT) => Some(FLAT),
        (ASC, ASC) | (ASC, FLAT) | (FLAT, ASC) => Some(ASC),
        // FLAT ao lado de um descendente: entra no run. Se o plato tem mais
        // de um elemento, o run passa a ter chaves repetidas.
        (FLAT, _) => Some(if ll > 1 { DESC_TIES } else { rk }),
        (_, FLAT) => Some(if rl > 1 { DESC_TIES } else { lk }),
        // Direcoes opostas: nao ha como unir.
        (ASC, _) | (_, ASC) => None,
        _ => Some(if lk == DESC_TIES || rk == DESC_TIES { DESC_TIES } else { DESC }),
    };
    match kind {
        Some(k) => (Run::new(k, ll + rl - 1), None),
        None => {
            // O elemento compartilhado fica com a esquerda.
            let rest = match rl - 1 {
                0 => None,
                1 => Some(Run::new(FLAT, 1)),
                m => Some(Run::new(rk, m)),
            };
            (l, rest)
        }
    }
}

/// Costura, no lugar, o run em `mark - 1` com o run em `mark`.
fn stitch_at(v: &mut Vec<Run>, mark: usize) {
    if mark == 0 || mark >= v.len() {
        return;
    }
    let (l, r) = join_runs(v[mark - 1], v[mark]);
    v[mark - 1] = l;
    match r {
        Some(r) => v[mark] = r,
        None => {
            v.remove(mark);
        }
    }
}

fn stitch(mut left: Vec<Run>, right: Vec<Run>) -> Vec<Run> {
    if left.is_empty() {
        return right;
    }
    if right.is_empty() {
        return left;
    }
    let last = left.len() - 1;
    let (l, r) = join_runs(left[last], right[0]);
    left[last] = l;
    if let Some(r) = r {
        left.push(r);
    }
    left.extend_from_slice(&right[1..]);
    left
}

/// Deteccao paralela. Blocos com 1 elemento de sobreposicao, reduzidos pelo
/// monoide. O `fold` do Rayon da a cada tarefa UM vetor que acumula os blocos
/// consecutivos dela, em vez de um vetor por bloco.
///
/// A versao anterior tinha dois niveis (macro 32768 / micro 512). Medido em
/// dado ordenado, um nivel e dois niveis empatam, e ambos empatam com um
/// `is_sorted` paralelo cru: a varredura e sequencial e a largura de banda
/// manda. O nivel micro so custava -- um `Vec` a cada 512 elementos.
///
/// `seg_limit` e o aborto antecipado: quando a soma dos segmentos estimados
/// passa do limite, os blocos restantes nao sao varridos e a funcao devolve
/// `None`. Dado quase ordenado com 1% de trocas, por exemplo, estoura o
/// limite depois de varrer uma fracao minima do array.
fn detect_runs<T, F>(
    arr: &[T],
    cmp: &F,
    block: usize,
    min_run: usize,
    seg_limit: usize,
) -> Option<Vec<Run>>
where
    T: Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let n = arr.len();
    if n == 0 {
        return Some(Vec::new());
    }
    if n == 1 {
        return Some(vec![Run::new(FLAT, 1)]);
    }
    let block = block.max(2);
    let step = block - 1;
    let blocks = (n - 2) / step + 1;
    let segs = AtomicUsize::new(0);
    let runs = (0..blocks)
        .into_par_iter()
        .fold(Vec::new, |mut acc: Vec<Run>, i| {
            if segs.load(Relaxed) > seg_limit {
                return acc;
            }
            let start = i * step;
            let end = (start + block).min(n);
            let mark = acc.len();
            let s = scan_runs(&arr[start..end], cmp, &mut acc, min_run);
            segs.fetch_add(s, Relaxed);
            stitch_at(&mut acc, mark);
            acc
        })
        .reduce(Vec::new, stitch);
    if segs.load(Relaxed) > seg_limit {
        None
    } else {
        Some(runs)
    }
}

// =====================================================================
//  FASE 2: NORMALIZACAO
//
//  Transforma a lista de runs em segmentos ASCENDENTES, no lugar:
//
//    Keep      ja ascendente
//    Reverse   descendente -> inverte (e desfaz os empates, se estavel)
//    Sort      varios runs curtos seguidos -> ordena o trecho de uma vez
//
//  Depois disso a arvore so enxerga runs ascendentes, e a direcao deixa de
//  ser assunto dela.
// =====================================================================

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Op {
    Keep,
    Reverse { ties: bool },
    Sort,
}

struct Plan {
    lens: Vec<usize>,
    ops: Vec<Op>,
    /// Elementos que cairam em segmentos `Sort`.
    sorted_elems: usize,
}

struct Pending {
    len: usize,
    count: usize,
    first: Run,
}

impl Plan {
    fn natural(&mut self, r: Run) {
        self.lens.push(r.len());
        self.ops.push(if r.is_desc() {
            Op::Reverse { ties: r.kind() == DESC_TIES }
        } else {
            Op::Keep
        });
    }

    fn flush(&mut self, p: &mut Pending) {
        match p.count {
            0 => {}
            // Um run curto isolado ja e monotono: nao precisa de sort.
            1 => self.natural(p.first),
            _ => {
                self.lens.push(p.len);
                self.ops.push(Op::Sort);
                self.sorted_elems += p.len;
            }
        }
        p.len = 0;
        p.count = 0;
    }
}

/// Agrupa runs curtos consecutivos em blocos a ordenar.
///
/// Fundir run a run so compensa quando o run e longo. Medido em dente de
/// serra de comprimento L: com L=16 ordenar em bloco ganha da arvore, com L=32
/// empata, com L>=64 a arvore ganha. Um trecho de runs curtos vira um bloco de
/// ate `chunk` elementos -- e um array aleatorio que escapou do escudo vira
/// simplesmente uma sequencia de blocos, sem penhasco.
fn plan_segments(runs: &[Run], min_run: usize, chunk: usize) -> Plan {
    let mut plan = Plan {
        lens: Vec::with_capacity(runs.len()),
        ops: Vec::with_capacity(runs.len()),
        sorted_elems: 0,
    };
    let mut pend = Pending { len: 0, count: 0, first: Run::new(FLAT, 0) };
    for &r in runs {
        let len = r.len();
        if len < min_run {
            if pend.count > 0 && pend.len + len > chunk {
                plan.flush(&mut pend);
            }
            if pend.count == 0 {
                pend.first = r;
            }
            pend.len += len;
            pend.count += 1;
        } else {
            plan.flush(&mut pend);
            plan.natural(r);
        }
    }
    plan.flush(&mut pend);
    plan
}

/// Aplica as operacoes em paralelo, cada uma na sua fatia.
fn apply_ops<T, F>(arr: &mut [T], plan: &Plan, cmp: &F, stable: bool)
where
    T: Send,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let mut jobs: Vec<(&mut [T], Op)> = Vec::new();
    let mut rest: &mut [T] = arr;
    for (&len, &op) in plan.lens.iter().zip(&plan.ops) {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(len);
        if op != Op::Keep {
            jobs.push((head, op));
        }
        rest = tail;
    }
    jobs.into_par_iter().for_each(|(s, op)| match op {
        Op::Reverse { ties } => {
            parallel_reverse(s);
            if ties && stable {
                restore_ties(s, cmp);
            }
        }
        Op::Sort => leaf_sort(s, cmp, stable),
        Op::Keep => {}
    });
}

fn parallel_reverse<T: Send>(arr: &mut [T]) {
    let n = arr.len();
    if n <= 100_000 {
        arr.reverse();
        return;
    }
    let (left, right) = arr.split_at_mut(n / 2);
    left.par_iter_mut()
        .zip(right.par_iter_mut().rev())
        .for_each(|(a, b)| std::mem::swap(a, b));
}

/// `s` esta ordenado, mas cada bloco de chaves iguais ficou na ordem INVERSA
/// da original (acabou de sair de um `reverse`). Desinverte cada bloco.
///
/// Como `s` ja esta ordenado, as bordas do bloco que contem o meio saem por
/// busca binaria, e as duas metades seguem em paralelo.
fn restore_ties<T, F>(s: &mut [T], cmp: &F)
where
    T: Send,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let n = s.len();
    if n < 2 {
        return;
    }
    if n <= 16_384 {
        let mut i = 0usize;
        while i < n {
            let mut j = i + 1;
            while j < n && cmp(&s[i], &s[j]) == Equal {
                j += 1;
            }
            if j - i > 1 {
                s[i..j].reverse();
            }
            i = j;
        }
        return;
    }
    let mid = n / 2;
    let lo = s[..mid].partition_point(|x| cmp(x, &s[mid]) == Less);
    let hi = mid + s[mid..].partition_point(|x| cmp(x, &s[mid]) != Greater);
    let (left, rest) = s.split_at_mut(lo);
    let (group, right) = rest.split_at_mut(hi - lo);
    rayon::join(
        || restore_ties(left, cmp),
        || {
            rayon::join(
                || {
                    if group.len() > 1 {
                        parallel_reverse(group)
                    }
                },
                || restore_ties(right, cmp),
            )
        },
    );
}

/// Funde, so nos offsets, segmentos ascendentes adjacentes que ja estao em
/// ordem global (ultimo de um <= primeiro do seguinte). Uma comparacao por
/// fronteira, nenhum elemento movido.
///
/// Antes, o atalho "ja em ordem" vivia DENTRO da arvore e era um memcpy por
/// nivel: mil dentes em faixas disjuntas pagavam log(1000) passadas de copia.
/// Medido em 10M u64: 112 ms -> 9,9 ms.
fn coalesce<T, F>(arr: &[T], offsets: Vec<usize>, cmp: &F) -> Vec<usize>
where
    T: Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let s = offsets.len() - 1;
    if s <= 1 {
        return offsets;
    }
    let is_break = |i: usize| cmp(&arr[offsets[i] - 1], &arr[offsets[i]]) == Greater;
    let breaks: Vec<bool> = if s > 4096 {
        (1..s).into_par_iter().map(is_break).collect()
    } else {
        (1..s).map(is_break).collect()
    };
    let mut out = Vec::with_capacity(s + 1);
    out.push(offsets[0]);
    for i in 1..s {
        if breaks[i - 1] {
            out.push(offsets[i]);
        }
    }
    out.push(offsets[s]);
    out
}

fn prefix_sums(lens: &[usize]) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(lens.len() + 1);
    let mut off = 0usize;
    offsets.push(0);
    for &l in lens {
        off += l;
        offsets.push(off);
    }
    offsets
}

fn sort_structured<T, F>(
    arr: &mut [T],
    runs: &[Run],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
    cheap: bool,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let n = arr.len();
    let plan = plan_segments(runs, t.min_run, t.chunk);
    let offsets = prefix_sums(&plan.lens);

    if t.route == Route::Auto {
        // A decisao exata, com o metadata inteiro na mao: quantos runs sobram
        // DEPOIS de encadear? Mil dentes em faixas disjuntas sao mil runs e
        // um unico run encadeado; mil dentes intercalados sao mil runs mesmo.
        let allowed = t.allowed_segs(n, cheap);
        let fragmented = plan.lens.len() > allowed
            && runs_after_coalesce(arr, &plan, &offsets, cmp) > allowed;
        if fragmented {
            sort_chaotic(arr, cmp, stable, t, scratch, Chaos::Fragmented, cheap);
            return;
        }
        // O escudo foi enganado e a maior parte do dado nao tem run que se
        // aproveite. O estavel segue: os blocos `Sort` do plano ja sao o
        // caminho caotico dele. O instavel escolhe pelo comprimento medio do
        // run bruto -- aleatorio da ~2,4, dente de serra curto da o proprio
        // comprimento do dente.
        if !stable && plan.sorted_elems * 2 > n {
            let kind = if runs.len() * 6 > n { Chaos::Random } else { Chaos::ShortRuns };
            sort_chaotic(arr, cmp, stable, t, scratch, kind, cheap);
            return;
        }
    }

    apply_ops(arr, &plan, cmp, stable);
    let offsets = coalesce(arr, offsets, cmp);
    merge_segments(arr, &offsets, cmp, t, scratch);
}

/// Quantos runs sobrariam depois de normalizar e encadear -- calculado ANTES
/// de mover qualquer elemento, so com os extremos de cada segmento. Um
/// segmento `Sort` nao tem extremos conhecidos e conta como quebra.
fn runs_after_coalesce<T, F>(arr: &[T], plan: &Plan, offsets: &[usize], cmp: &F) -> usize
where
    T: Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let s = plan.lens.len();
    if s <= 1 {
        return s;
    }
    let is_break = |i: usize| {
        let (prev, next) = (plan.ops[i - 1], plan.ops[i]);
        if prev == Op::Sort || next == Op::Sort {
            return true;
        }
        let max_prev = match prev {
            Op::Reverse { .. } => &arr[offsets[i - 1]],
            _ => &arr[offsets[i] - 1],
        };
        let min_next = match next {
            Op::Reverse { .. } => &arr[offsets[i + 1] - 1],
            _ => &arr[offsets[i]],
        };
        cmp(max_prev, min_next) == Greater
    };
    let breaks = if s > 4096 {
        (1..s).into_par_iter().filter(|&i| is_break(i)).count()
    } else {
        (1..s).filter(|&i| is_break(i)).count()
    };
    breaks + 1
}

// =====================================================================
//  CAMINHO CAOTICO
// =====================================================================

/// Por que o dado veio parar no caminho caotico.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Chaos {
    /// Aleatorio de fato: a direcao inverte em mais de um terco da amostra.
    Random,
    /// Runs curtos, mas com ordem local (dente de serra de 16 a 60).
    ShortRuns,
    /// Runs de verdade, so que em numero grande demais para a arvore.
    Fragmented,
}

fn sort_chaotic<T, F>(
    arr: &mut [T],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
    kind: Chaos,
    cheap: bool,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    // INSTAVEL: o quicksort paralelo trabalha no lugar e sem buffer, mas so
    // em dois casos.
    //
    //  - Dado aleatorio. Contra ordenar blocos e fundir, em 10M u64: 124 ms
    //    contra 138-154; 52 contra 64-82 (64 chaves); 18 contra 50-59
    //    (2 chaves).
    //  - Dado globalmente ascendente com elementos fora do lugar E comparacao
    //    barata (ver o atalho em `sort_with`).
    //
    // Fora disso ele NAO e seguro. Runs intercalados sao um padrao ruim para
    // a escolha de pivo, e com chave cara o estrago e de ordem de grandeza:
    // 8 runs de 250 mil strings levam 1942 ms nele contra 109 ms da arvore;
    // 64 runs, 1121 ms contra 175. E em dente de serra curto de u64 os blocos
    // ja ganham (L=32: 116 ms contra 148).
    if !stable {
        let quick = match kind {
            Chaos::Random => true,
            Chaos::ShortRuns => false,
            Chaos::Fragmented => cheap && globally_ascending(arr, cmp),
        };
        if quick {
            arr.par_sort_unstable_by(|a, b| cmp(a, b));
            return;
        }
    }
    // BLOCOS: ordenados pelo `sort` da std (driftsort, que aproveita os runs
    // de dentro do bloco) e fundidos pela arvore. E o caminho do motor
    // estavel, e o do instavel quando o quicksort nao serve.
    //
    // Contra o `par_sort` do Rayon, em 10M u64: 175 ms contra 191 (aleatorio),
    // 102 contra 185 (1024 chaves), 89 contra 190 (64 chaves), 78 contra 154
    // (2 chaves). O bloco de 65536 venceu ou empatou em toda a varredura
    // (16K, 64K, 256K, n/8, n/4) -- o mesmo valor que o motor C++ ja usava.
    let chunk = t.chunk.max(1);
    let n = arr.len();
    if n <= chunk {
        leaf_sort(arr, cmp, true);
        return;
    }
    arr.par_chunks_mut(chunk).for_each(|c| leaf_sort(c, cmp, true));
    let mut offsets: Vec<usize> = (0..n).step_by(chunk).collect();
    offsets.push(n);
    let offsets = coalesce(arr, offsets, cmp);
    merge_segments(arr, &offsets, cmp, t, scratch);
}

// =====================================================================
//  FASE 3: ARVORE DE MERGE
//
//  Entrada: `offsets` delimita runs ASCENDENTES em `arr`.
//
//  O BUFFER
//
//  Sem buffer do chamador, ele nasce como COPIA PARALELA do array. Isso troca
//  duas coisas da versao anterior -- o preenchimento sequencial
//  (`vec![arr[0]; n]`, uma thread escrevendo n elementos que seriam
//  sobrescritos) e a copia das folhas de paridade impar -- por uma passada
//  paralela so. Com os dois lados iguais na largada, toda folha ja esta onde
//  precisa estar.
//
//  A medicao que tinha concluido que inicializar o buffer "nao custa nada"
//  (457,4 contra 458,5 ms) estava confundida: quando `arr[0]` e o inteiro
//  zero, `vec![arr[0]; n]` vira `calloc` e o sistema entrega paginas zeradas
//  sob demanda. Com `arr[0] != 0` o preenchimento custa ~14 ms por 80 MB
//  (64 runs em 10M u64: 64,0 ms com zero na frente, 77,5 sem). Dado real
//  quase nunca comeca com zero.
//
//  Com buffer do chamador (`*_with_buffer`), ja inicializado e do tamanho
//  certo, nao ha nada a criar nem pagina nova a tocar: so as folhas de
//  paridade impar copiam o proprio trecho.
// =====================================================================

fn merge_segments<T, F>(
    arr: &mut [T],
    offsets: &[usize],
    cmp: &F,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    if offsets.len() <= 2 {
        return;
    }
    let n = arr.len();
    match scratch {
        Some(buf) if buf.len() >= n => {
            tree_kway(arr, &mut buf[..n], offsets, cmp, t, false, true);
        }
        Some(buf) => {
            arr.par_iter().copied().collect_into_vec(buf);
            tree_kway(arr, &mut buf[..], offsets, cmp, t, false, false);
        }
        None => {
            let mut buffer: Vec<T> = Vec::with_capacity(n);
            arr.par_iter().copied().collect_into_vec(&mut buffer);
            tree_kway(arr, &mut buffer, offsets, cmp, t, false, false);
        }
    }
}

/// Recursao K-aria sobre os runs. Hibrida: k-vias nos niveis que vao a DRAM,
/// binaria abaixo de `kway_min_elems`, onde o trecho ja mora em cache.
///
/// Contrato: o dado de entrada esta em `v`. Na saida o trecho ordenado esta
/// em `buf` se `into_buf`, senao em `v`. Se `copy_leaves` e falso, `buf` ja
/// nasceu como copia de `v` e a folha nao tem nada a fazer; se e verdadeiro,
/// `buf` tem lixo valido e a folha de paridade impar copia o proprio trecho.
fn tree_kway<T, F>(v: &mut [T], buf: &mut [T], offsets: &[usize], cmp: &F, t: &Tuning, into_buf: bool, copy_leaves: bool)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let runs = offsets.len() - 1;
    if runs <= 1 {
        if into_buf && copy_leaves {
            par_copy(buf, v, t.par_grain);
        }
        return;
    }
    let base = offsets[0];
    let total = offsets[runs] - base;
    if runs == 2 || total <= t.kway_min_elems {
        tree_binary(v, buf, offsets, cmp, t, into_buf, copy_leaves);
        return;
    }

    // Cortes balanceados por ELEMENTO, nao por contagem de runs.
    let g = KWAY_FANOUT.min(runs);
    let mut cm = [0usize; KWAY_FANOUT + 1];
    cm[g] = runs;
    for i in 1..g {
        let target = base + total / g * i + total % g * i / g;
        let idx = offsets.partition_point(|&o| o < target);
        cm[i] = idx.clamp(cm[i - 1] + 1, runs - (g - i));
    }

    {
        let mut parts: Vec<(&mut [T], &mut [T], usize, usize)> = Vec::with_capacity(g);
        let mut vv: &mut [T] = &mut *v;
        let mut bb: &mut [T] = &mut *buf;
        let mut prev = 0usize;
        for i in 0..g {
            let e = offsets[cm[i + 1]] - base;
            let (vl, vr) = std::mem::take(&mut vv).split_at_mut(e - prev);
            let (bl, br) = std::mem::take(&mut bb).split_at_mut(e - prev);
            parts.push((vl, bl, cm[i], cm[i + 1]));
            vv = vr;
            bb = br;
            prev = e;
        }
        parts.into_par_iter().for_each(|(sv, sb, lo, hi)| {
            tree_kway(sv, sb, &offsets[lo..=hi], cmp, t, !into_buf, copy_leaves);
        });
    }

    let (src, dst): (&[T], &mut [T]) = if into_buf { (&*v, buf) } else { (&*buf, v) };

    // Grupos adjacentes ja em ordem viram UM stream.
    let mut streams: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
    let mut ns = 0usize;
    let mut gs = 0usize;
    for i in 0..g {
        let e = offsets[cm[i + 1]] - base;
        let join_next = i + 1 < g && cmp(&src[e - 1], &src[e]) != Greater;
        if !join_next {
            streams[ns] = &src[gs..e];
            ns += 1;
            gs = e;
        }
    }
    kway_merge(&streams[..ns], dst, cmp, t);
}

/// Arvore binaria, com split balanceado por elementos.
fn tree_binary<T, F>(v: &mut [T], buf: &mut [T], offsets: &[usize], cmp: &F, t: &Tuning, into_buf: bool, copy_leaves: bool)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let runs = offsets.len() - 1;
    if runs <= 1 {
        if into_buf && copy_leaves {
            par_copy(buf, v, t.par_grain);
        }
        return;
    }
    let base = offsets[0];
    let total = offsets[runs] - base;
    let target = base + total / 2;
    let split = offsets.partition_point(|&o| o < target).clamp(1, runs - 1);
    let mid = offsets[split] - base;
    {
        let (vl, vr) = v.split_at_mut(mid);
        let (bl, br) = buf.split_at_mut(mid);
        let (lo, ro) = (&offsets[..=split], &offsets[split..]);
        if total > t.seq_node {
            rayon::join(
                || tree_binary(vl, bl, lo, cmp, t, !into_buf, copy_leaves),
                || tree_binary(vr, br, ro, cmp, t, !into_buf, copy_leaves),
            );
        } else {
            tree_binary(vl, bl, lo, cmp, t, !into_buf, copy_leaves);
            tree_binary(vr, br, ro, cmp, t, !into_buf, copy_leaves);
        }
    }
    let (src, dst): (&[T], &mut [T]) = if into_buf { (&*v, buf) } else { (&*buf, v) };
    let (a, b) = src.split_at(mid);
    merge_two(a, b, dst, cmp, t);
}

// =====================================================================
//  MERGE DE DOIS: CO-RANK + BIDIRECIONAL
// =====================================================================

/// Maior `i` tal que `a[..i]` e `b[..k - i]` sao os `k` menores, com empate
/// para `a`. E o corte que mantem o merge estavel.
fn co_rank<T, F>(k: usize, a: &[T], b: &[T], cmp: &F) -> (usize, usize)
where
    F: Fn(&T, &T) -> Ordering,
{
    let (m, n) = (a.len(), b.len());
    let mut lo = k.saturating_sub(n);
    let mut hi = k.min(m);
    while lo < hi {
        let i = lo + (hi - lo + 1) / 2;
        let j = k - i;
        if j == n || cmp(&a[i - 1], &b[j]) != Greater {
            lo = i;
        } else {
            hi = i - 1;
        }
    }
    (lo, k - lo)
}

/// Merge para FRENTE: emite os `dest.len()` menores. Empate -> A.
#[inline]
fn merge_front<T, F>(a: &[T], b: &[T], dest: &mut [T], cmp: &F)
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    let n = dest.len();
    let (mut i, mut j, mut k) = (0usize, 0usize, 0usize);
    while k < n && i < a.len() && j < b.len() {
        if cmp(&b[j], &a[i]) == Less {
            dest[k] = b[j];
            j += 1;
        } else {
            dest[k] = a[i];
            i += 1;
        }
        k += 1;
    }
    if k < n {
        let rest = n - k;
        if i < a.len() {
            dest[k..].copy_from_slice(&a[i..i + rest]);
        } else {
            dest[k..].copy_from_slice(&b[j..j + rest]);
        }
    }
}

/// Merge para TRAS: emite os `dest.len()` maiores. Empate -> B.
///
/// A regra de desempate e a COMPLEMENTAR da de `merge_front`, e e isso que faz
/// o teorema valer: os k menores (empate -> A) e os total-k maiores
/// (empate -> B) particionam o multiconjunto para QUALQUER k. Nenhum ponto de
/// corte precisa ser calculado -- cada lado so conta as proprias saidas.
#[inline]
fn merge_back<T, F>(a: &[T], b: &[T], dest: &mut [T], cmp: &F)
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    let (mut qa, mut qb, mut k) = (a.len(), b.len(), dest.len());
    while k > 0 && qa > 0 && qb > 0 {
        k -= 1;
        if cmp(&b[qb - 1], &a[qa - 1]) == Less {
            qa -= 1;
            dest[k] = a[qa];
        } else {
            qb -= 1;
            dest[k] = b[qb];
        }
    }
    if k > 0 {
        if qa > 0 {
            dest[..k].copy_from_slice(&a[qa - k..qa]);
        } else {
            dest[..k].copy_from_slice(&b[qb - k..qb]);
        }
    }
}

/// Merge paralelo de dois runs, de fontes imutaveis. Corta por co-rank ate o
/// grao; em cada pedaco tenta primeiro os atalhos de faixas disjuntas e so
/// entao funde -- com uma thread vindo da frente e outra do fim.
fn merge_two<T, F>(a: &[T], b: &[T], dest: &mut [T], cmp: &F, t: &Tuning)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let total = dest.len();
    if total > t.par_grain && total >= 2 {
        let k = total / 2;
        let (i, j) = co_rank(k, a, b, cmp);
        let (al, ar) = a.split_at(i);
        let (bl, br) = b.split_at(j);
        let (dl, dr) = dest.split_at_mut(k);
        rayon::join(
            || merge_two(al, bl, dl, cmp, t),
            || merge_two(ar, br, dr, cmp, t),
        );
        return;
    }
    if a.is_empty() {
        dest.copy_from_slice(b);
        return;
    }
    if b.is_empty() {
        dest.copy_from_slice(a);
        return;
    }
    // Faixas disjuntas: e uma concatenacao.
    if cmp(&a[a.len() - 1], &b[0]) != Greater {
        let (d1, d2) = dest.split_at_mut(a.len());
        d1.copy_from_slice(a);
        d2.copy_from_slice(b);
        return;
    }
    if cmp(&b[b.len() - 1], &a[0]) == Less {
        let (d1, d2) = dest.split_at_mut(b.len());
        d1.copy_from_slice(b);
        d2.copy_from_slice(a);
        return;
    }
    if total > t.bidir_min {
        let (df, db) = dest.split_at_mut(total / 2);
        rayon::join(
            || merge_front(a, b, df, cmp),
            || merge_back(a, b, db, cmp),
        );
    } else {
        merge_front(a, b, dest, cmp);
    }
}

fn par_copy<T: Copy + Send + Sync>(dest: &mut [T], src: &[T], grain: usize) {
    if dest.len() <= grain.saturating_mul(4) {
        dest.copy_from_slice(src);
        return;
    }
    let g = grain.max(1);
    dest.par_chunks_mut(g)
        .zip(src.par_chunks(g))
        .for_each(|(d, s)| d.copy_from_slice(s));
}

// =====================================================================
//  MERGE K-VIAS COM BLOCAGEM DE CACHE
//
//  O merge binario paga log2(runs) passadas sobre a memoria, e a maquina
//  satura a banda: o custo e o NUMERO DE PASSADAS. Um merge de K vias paga
//  logK(runs).
//
//  NAO usa loser tree: medido, custa ~12x mais por elemento que o merge
//  simples, porque `head[w]` com `w` em tempo de execucao nao e promovido a
//  registrador. Cada tile e cortado por `multiseq_partition` e fundido por uma
//  arvore BINARIA que cabe em L1 -- o kernel rapido continua o mesmo.
// =====================================================================

/// Acha `out[i]` com soma igual a `rank`, tal que os prefixos sao exatamente
/// os `rank` menores elementos. Empates sao distribuidos em ordem de stream,
/// o que mantem a estabilidade e concorda com o desempate de `merge_front`.
///
/// As buscas binarias rodam so dentro de `[lo, hi)`. Todo pivo sai de dentro
/// do intervalo vivo, entao fica estritamente entre os pivos anteriores, e
/// nenhum elemento fora do intervalo pode mudar de lado. O resultado e
/// identico ao da busca no run inteiro; contado em 8 runs, as comparacoes por
/// elemento do merge caem de 6,24 para 5,22 (o ideal e 4,00).
fn multiseq_partition<T, F>(seq: &[&[T]], rank: usize, out: &mut [usize], cmp: &F)
where
    F: Fn(&T, &T) -> Ordering,
{
    let ns = seq.len();
    let mut lo = [0usize; KWAY_FANOUT];
    let mut hi = [0usize; KWAY_FANOUT];
    let mut plt = [0usize; KWAY_FANOUT];
    let mut ple = [0usize; KWAY_FANOUT];
    for i in 0..ns {
        hi[i] = seq[i].len();
    }

    loop {
        // Pivo: o meio do intervalo mais largo.
        let mut m = usize::MAX;
        let mut w = 0usize;
        for i in 0..ns {
            if hi[i] - lo[i] > w {
                w = hi[i] - lo[i];
                m = i;
            }
        }
        if m == usize::MAX {
            out[..ns].copy_from_slice(&lo[..ns]);
            return;
        }
        let p = &seq[m][lo[m] + w / 2];

        let (mut lt, mut le) = (0usize, 0usize);
        for i in 0..ns {
            let s = seq[i];
            plt[i] = lo[i] + s[lo[i]..hi[i]].partition_point(|x| cmp(x, p) == Less);
            ple[i] = plt[i] + s[plt[i]..hi[i]].partition_point(|x| cmp(x, p) != Greater);
            lt += plt[i];
            le += ple[i];
        }

        if lt <= rank && rank <= le {
            let mut need = rank - lt;
            for i in 0..ns {
                let take = (ple[i] - plt[i]).min(need);
                out[i] = plt[i] + take;
                need -= take;
            }
            return;
        }
        if lt > rank {
            hi[..ns].copy_from_slice(&plt[..ns]);
        } else {
            lo[..ns].copy_from_slice(&ple[..ns]);
        }
    }
}

/// Um nivel da arvore binaria do tile: funde pares de segmentos contiguos de
/// `src` em `dst`. Devolve os novos segmentos (offset, len) e quantos sao.
fn merge_level<T, F>(
    src: &[T],
    segs: &[(usize, usize)],
    dst: &mut [T],
    cmp: &F,
) -> ([(usize, usize); KWAY_FANOUT], usize)
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    let mut out = [(0usize, 0usize); KWAY_FANOUT];
    let mut cnt = 0usize;
    let mut off = 0usize;
    let mut i = 0usize;
    while i < segs.len() {
        let (o1, l1) = segs[i];
        let l = if i + 1 == segs.len() {
            dst[off..off + l1].copy_from_slice(&src[o1..o1 + l1]);
            l1
        } else {
            let (_, l2) = segs[i + 1];
            let (left, right) = src[o1..o1 + l1 + l2].split_at(l1);
            merge_front(left, right, &mut dst[off..off + l1 + l2], cmp);
            l1 + l2
        };
        out[cnt] = (off, l);
        off += l;
        cnt += 1;
        i += 2;
    }
    (out, cnt)
}

/// Merge de 3..=K runs com blocagem de cache: a DRAM ve uma leitura e uma
/// escrita, em vez de log2(K) passadas.
fn kway_tiled<T, F>(runs: &[&[T]], dest: &mut [T], cmp: &F, tile: usize)
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    let ns = runs.len();
    let total = dest.len();
    if ns == 1 {
        dest.copy_from_slice(&runs[0][..total]);
        return;
    }
    if ns == 2 {
        merge_front(runs[0], runs[1], dest, cmp);
        return;
    }
    let filler = match runs.iter().find(|r| !r.is_empty()) {
        Some(r) => r[0],
        None => return,
    };
    let tile = tile.max(1);
    let mut scratch = vec![filler; 2 * tile];
    let (s1, s2) = scratch.split_at_mut(tile);

    let mut cur = [0usize; KWAY_FANOUT];
    let mut done = 0usize;

    while done < total {
        let want = tile.min(total - done);

        // Janela restrita: um tile consome no maximo `want` de qualquer run.
        let mut win: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
        for i in 0..ns {
            let rem = &runs[i][cur[i]..];
            win[i] = &rem[..rem.len().min(want)];
        }
        let mut cut = [0usize; KWAY_FANOUT];
        multiseq_partition(&win[..ns], want, &mut cut[..ns], cmp);

        // Nivel 1: das janelas para s1.
        let mut segs = [(0usize, 0usize); KWAY_FANOUT];
        let mut cnt = 0usize;
        let mut off = 0usize;
        let mut i = 0usize;
        while i < ns {
            let l = if i + 1 == ns {
                let l = cut[i];
                s1[off..off + l].copy_from_slice(&win[i][..l]);
                l
            } else {
                let l = cut[i] + cut[i + 1];
                merge_front(
                    &win[i][..cut[i]],
                    &win[i + 1][..cut[i + 1]],
                    &mut s1[off..off + l],
                    cmp,
                );
                l
            };
            segs[cnt] = (off, l);
            off += l;
            cnt += 1;
            i += 2;
        }

        // Niveis intermediarios em s1/s2 (residentes em L1); o ULTIMO nivel
        // escreve direto no destino -- a unica escrita que chega na DRAM.
        let out = &mut dest[done..done + want];
        if cnt > 2 {
            let (segs2, _) = merge_level(&s1[..want], &segs[..cnt], s2, cmp);
            let (o1, l1) = segs2[0];
            let (_, l2) = segs2[1];
            let (left, right) = s2[o1..o1 + l1 + l2].split_at(l1);
            merge_front(left, right, out, cmp);
        } else {
            let (o1, l1) = segs[0];
            let (_, l2) = segs[1];
            let (left, right) = s1[o1..o1 + l1 + l2].split_at(l1);
            merge_front(left, right, out, cmp);
        }

        for i in 0..ns {
            cur[i] += cut[i];
        }
        done += want;
    }
}

/// Merge k-vias paralelo: corta por rank e recursa.
///
/// Streams vazios sao descartados a cada nivel. Depois de um corte e comum
/// sobrar um ou dois streams com dado num pedaco; sem o descarte eles
/// atravessavam os tres niveis do tile como se fossem oito. Com UM stream e
/// copia, com DOIS e o merge binario por co-rank -- que escala com as threads,
/// onde antes o caso de dois runs ficava preso em duas.
fn kway_merge<T, F>(runs: &[&[T]], dest: &mut [T], cmp: &F, t: &Tuning)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    let mut live: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
    let mut ns = 0usize;
    for r in runs {
        if !r.is_empty() {
            live[ns] = r;
            ns += 1;
        }
    }
    let runs = &live[..ns];
    let total = dest.len();
    match ns {
        0 => return,
        1 => {
            par_copy(dest, runs[0], t.par_grain);
            return;
        }
        2 => {
            merge_two(runs[0], runs[1], dest, cmp, t);
            return;
        }
        _ => {}
    }
    if total <= t.par_grain || total < 2 {
        kway_tiled(runs, dest, cmp, t.tile);
        return;
    }

    let mut cut = [0usize; KWAY_FANOUT];
    multiseq_partition(runs, total / 2, &mut cut[..ns], cmp);

    let mut left: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
    let mut right: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
    let mut ltot = 0usize;
    for i in 0..ns {
        let (l, r) = runs[i].split_at(cut[i]);
        left[i] = l;
        right[i] = r;
        ltot += cut[i];
    }
    let (dl, dr) = dest.split_at_mut(ltot);
    rayon::join(
        || kway_merge(&left[..ns], dl, cmp, t),
        || kway_merge(&right[..ns], dr, cmp, t),
    );
}
