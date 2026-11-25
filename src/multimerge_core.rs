// COPIA de adaptive-parallel-multimerge-sort/src/multimerge.rs.
// O original e a fonte da verdade e guarda os testes do motor; nao edite
// esta copia -- altere la, rode os testes la e copie de novo.

// =====================================================================
//  ADAPTIVE PARALLEL MULTIMERGE SORT -- Rust engine
//
//  Two engines over the same run detection and the same merge tree:
//
//    multi_merge_sort          NOT stable. The default.
//    multi_merge_sort_stable   stable. For when the relative order of equal
//                              keys matters (DBMS, cascaded sorting).
//
//  Each one has the variants `_by` (comparator) and `_by_key` (key extractor).
//
//  sort(arr)
//    |- small n ............................ std sort
//    |- SHIELD (scattered windows) ......... chaotic -> chaotic path
//    '- run DETECTION (monoid, parallel)
//         |- consecutive short runs ........ become ONE block, sorted by the std
//         |- descending runs ............... reversed in place
//         |- runs already in global order .. merged without moving data
//         '- whatever is left .............. k-way merge tree
//
//  WHERE THE TWO ENGINES DIFFER
//
//  1. Descending run with repeated keys. A stretch 9 9 7 7 3 is descending
//     (non-increasing), but reversing it swaps the order of the equals. The
//     unstable engine reverses and is done. The stable one reverses and then
//     un-reverses each block of equals (`restore_ties`), which gives back the
//     original order among them. Both solve the case in O(n); before, the
//     engine broke the stretch into one run per plateau and paid for the
//     whole tree.
//  2. Block without structure. The stable engine uses the std `sort`
//     (driftsort), the unstable one uses `sort_unstable` (ipnsort).
//  3. Chaotic data. The stable engine sorts in blocks and merges through the
//     tree; the unstable one delegates to the parallel quicksort, which works
//     in place and allocates no buffer.
//
//  RADIX KEY
//
//  The `_radix`, `_by_radix_key` and `_by_prefix_key` entries tell the engine
//  what the BITS of the key are. With that, random data no longer goes
//  through comparison: one radix pass over the high bits and each bucket
//  solved inside the cache (see "CHAOTIC PATH WITH A KEY"). This holds for
//  both engines; the rest of the path -- shield, runs, tree -- does not change.
//
//  The engine has no `unsafe`. `T: Copy` is what lets the borrow checker
//  prove the bidirectional merge: the two threads read the sources through a
//  shared reference and write to disjoint halves of the destination.
// =====================================================================

use rayon::prelude::*;
use std::cmp::Ordering::{self, Equal, Greater, Less};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

/// Fan-out of the k-way merge.
pub const KWAY_FANOUT: usize = 8;

// =====================================================================
//  PUBLIC API
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

/// Integer types whose order is the order of an unsigned integer of up to 64
/// bits. Implemented for the primitive integers.
///
/// A type that implements it can use [`multi_merge_sort_radix`] and
/// [`multi_merge_sort_stable_radix`], where chaotic input is sorted by radix
/// instead of by comparison.
pub trait RadixKey: Ord + Copy + Send + Sync {
    /// Order-preserving map to `u64`: `a < b` exactly when
    /// `a.to_radix() < b.to_radix()`.
    fn to_radix(self) -> u64;
    /// Inverse of [`RadixKey::to_radix`].
    fn from_radix(bits: u64) -> Self;
}

macro_rules! radix_key_unsigned {
    ($($t:ty),*) => {$(
        impl RadixKey for $t {
            #[inline(always)]
            fn to_radix(self) -> u64 {
                self as u64
            }
            #[inline(always)]
            fn from_radix(bits: u64) -> Self {
                bits as $t
            }
        }
    )*};
}
radix_key_unsigned!(u8, u16, u32, u64, usize);

macro_rules! radix_key_signed {
    ($($t:ty => $u:ty),*) => {$(
        impl RadixKey for $t {
            #[inline(always)]
            fn to_radix(self) -> u64 {
                ((self as $u) ^ (1 << (<$u>::BITS - 1))) as u64
            }
            #[inline(always)]
            fn from_radix(bits: u64) -> Self {
                ((bits as $u) ^ (1 << (<$u>::BITS - 1))) as $t
            }
        }
    )*};
}
radix_key_signed!(i8 => u8, i16 => u16, i32 => u32, i64 => u64, isize => usize);

/// [`multi_merge_sort`] for integer keys. **Not stable.**
///
/// Same engine and same result. The difference is the chaotic path: random
/// input is sorted by a parallel radix pass plus in-cache bucket sorts instead
/// of the parallel quicksort. That path allocates a buffer of `arr.len()`
/// elements, except when the keys span 16 bits or less (counting sort, no
/// buffer). If the buffer cannot be allocated, or the keys are too skewed for
/// radix to pay off, the quicksort takes over.
pub fn multi_merge_sort_radix<T: RadixKey>(arr: &mut [T]) {
    sort_with_key(arr, &T::cmp, false, &Tuning::for_type::<T>(), None, &Bare);
}

/// [`multi_merge_sort_stable`] for integer keys. See
/// [`multi_merge_sort_radix`].
pub fn multi_merge_sort_stable_radix<T: RadixKey>(arr: &mut [T]) {
    sort_with_key(arr, &T::cmp, true, &Tuning::for_type::<T>(), None, &Bare);
}

/// Sorts by a `u64` key that decides the order completely: elements with equal
/// keys are equal. **Not stable.** Chaotic input is sorted by radix on the key.
pub fn multi_merge_sort_by_radix_key<T, K>(arr: &mut [T], key: K)
where
    T: Copy + Send + Sync,
    K: Fn(&T) -> u64 + Sync,
{
    let cmp = |a: &T, b: &T| key(a).cmp(&key(b));
    sort_with_key(arr, &cmp, false, &Tuning::for_type::<T>(), None, &Exact(&key));
}

/// [`multi_merge_sort_by_radix_key`], preserving the relative order of
/// elements with equal keys.
pub fn multi_merge_sort_stable_by_radix_key<T, K>(arr: &mut [T], key: K)
where
    T: Copy + Send + Sync,
    K: Fn(&T) -> u64 + Sync,
{
    let cmp = |a: &T, b: &T| key(a).cmp(&key(b));
    sort_with_key(arr, &cmp, true, &Tuning::for_type::<T>(), None, &Exact(&key));
}

/// Sorts with `compare`, using `prefix` as a radix key for chaotic input.
/// **Not stable.**
///
/// `prefix` must agree with `compare` wherever it can tell two elements apart:
/// `prefix(a) < prefix(b)` implies `compare(a, b) == Less`. Elements with the
/// same prefix are ordered by `compare`. This is the entry for keys longer
/// than 64 bits, such as strings: the prefix is the first 8 bytes, big-endian.
pub fn multi_merge_sort_by_prefix_key<T, K, F>(arr: &mut [T], prefix: K, compare: F)
where
    T: Copy + Send + Sync,
    K: Fn(&T) -> u64 + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with_key(arr, &compare, false, &Tuning::for_type::<T>(), None, &Prefix(&prefix));
}

/// [`multi_merge_sort_by_prefix_key`], preserving the relative order of equal
/// elements.
pub fn multi_merge_sort_stable_by_prefix_key<T, K, F>(arr: &mut [T], prefix: K, compare: F)
where
    T: Copy + Send + Sync,
    K: Fn(&T) -> u64 + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    sort_with_key(arr, &compare, true, &Tuning::for_type::<T>(), None, &Prefix(&prefix));
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
//  TUNING
//
//  A single place for the thresholds. `pub(crate)` so that the tests can force
//  deep recursion with small data, and so that a measurement probe can sweep
//  values without recompiling the engine.
// =====================================================================

// `Structured` and `Chaotic` are only built by the tests and by the measurement
// probe: they force a path in order to exercise it or time it in isolation.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Route {
    /// The shield decides.
    Auto,
    /// Forces run detection + tree.
    Structured,
    /// Forces the chaotic path.
    Chaotic,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Tuning {
    /// Up to this size the std sort handles it alone.
    pub small_sort: usize,
    /// Block of the run detection (with 1 element of overlap).
    pub detect_block: usize,
    /// A run shorter than this is not worth a tree node: consecutive short
    /// runs are grouped and sorted as one block.
    pub min_run: usize,
    /// The tree only comes in when, after chaining, the runs have on average
    /// at least this length.
    ///
    /// On u64 the tree over natural runs and the block path tie from 256
    /// up to 16384 elements per run, and the blocks win by 8% to 15% between
    /// 32 and 64. With an expensive key (`&str`) the tree wins comfortably even
    /// on a short run, because re-sorting the block repeats comparisons the runs
    /// had already resolved: sorted with 1% swaps (average run ~50), 367 ms
    /// vs 520; 64 runs of 31K strings, 175 vs 262. The value stays low
    /// to protect the expensive key, which is where a mistake costs the most.
    pub struct_min_avg_run: usize,
    /// The same threshold when the comparison is cheap. There the cost is
    /// moving memory, not comparing, and the blocks tie with the tree from 256
    /// up and win below: sorted with 1% swaps on 10M u64, 128 ms with the
    /// blocks vs 146 with the tree.
    pub struct_min_avg_run_cheap: usize,
    /// Floor on the number of accepted segments, for small arrays.
    pub min_segs: usize,
    /// Slack of the detection's early abort.
    pub abort_slack: usize,
    /// Size of the block sorted by the std (stable chaotic path, and cap on
    /// the grouping of short runs).
    pub chunk: usize,
    /// Parallel grain of the merge: above this it cuts by co-rank / multisequence.
    pub par_grain: usize,
    /// Above this a binary merge piece runs bidirectionally (2 tasks).
    pub bidir_min: usize,
    /// A binary-tree node below this recurses without `rayon::join`.
    pub seq_node: usize,
    /// A subtree up to this size uses the binary merge; above it, k-way.
    pub kway_min_elems: usize,
    /// Tile of the k-way merge, in elements.
    pub tile: usize,
    /// Shield: above this percentage of direction reversals the runs are too
    /// short for detection to be worth the pass.
    pub chaos_percent: usize,
    /// Shield: above this percentage the data is truly random, not just made
    /// of short runs. Only matters to the unstable engine.
    pub random_percent: usize,
    /// "Cheap" comparison: up to this average cost, in nanoseconds, measured
    /// on the shield's sample.
    pub cheap_ns: usize,
    /// Pins the cost verdict instead of measuring it. For the tests.
    pub force_cheap: Option<bool>,
    /// Radix: below this size the chaotic path of the UNSTABLE engine goes by
    /// comparison. Measured with the engine (minimum of 60, 8 threads): on i32
    /// radix loses up to 65K (0.87 vs 0.64 ms) and wins from 131K
    /// on (1.04 vs 1.30); on u64 it ties at 262K and wins at 1M.
    pub radix_min: usize,
    /// The same floor for the STABLE engine, which is lower: there the
    /// competitor is the stable sort in blocks, and radix already wins at 32K
    /// (0.59 vs 0.79 ms on i32; 0.54 vs 0.80 on u64) and ties at 16K.
    pub radix_min_stable: usize,
    /// Radix: the UNSTABLE engine's floor on random data with 1 thread. `radix_min` was
    /// measured with 8 threads, against the parallel quicksort; without parallelism radix wins
    /// from 8K elements on (1.3x to 1.9x at 8K; 1.4x to 3.1x from 16K to 262K).
    pub radix_min_serial: usize,
    /// Radix: the same floor with more than one thread, in elements PER THREAD, capped
    /// at `radix_min`. The parallel quicksort scales with the threads and a small radix
    /// does not. Measured: u64 starts winning at 16-32K with 2 threads, 32-65K with 4, 131K
    /// with 8 and 262-524K with 24; a 16-byte pair, at 16-65K, 131K, 262K and 524K.
    pub radix_min_thread: usize,
    /// Radix: bits of the digit of the MSD pass. 8 beat 10 (more buckets cost
    /// more in the scatter than they save in the leaves).
    pub radix_digit: u32,
    /// Radix: a "bare" key whose range fits in these bits is solved by
    /// counting, without a buffer.
    pub radix_count_bits: u32,
    /// Radix: the bucket is sorted by LSD when up to these bits remain and the
    /// comparison is cheap; above that, by the std sort. On i32 (24 remain)
    /// LSD wins; on full-range u64 (56 remain) it loses.
    pub radix_lsd_cheap: u32,
    /// The same threshold when the comparison is expensive: there LSD pays off
    /// even with more passes (6-byte record, 39 remain: 411 ms vs 741 on 30M).
    pub radix_lsd_costly: u32,
    /// Radix: elements per piece of the parallel histogram and scatter.
    pub radix_piece: usize,
    /// Radix: a bucket above this, and above 1/`radix_big_div` of the stretch,
    /// goes through the routine again instead of going to a sequential leaf.
    pub radix_big: usize,
    pub radix_big_div: usize,
    /// Radix: when the leaf is NOT a short LSD (more bits remain than
    /// `radix_lsd_cheap`: wide key, text), a bucket above this size gets
    /// a second MSD cut before the leaf. Measured with the engine, 8
    /// threads, cut at 64K vs no cut: 6-byte text in u64,
    /// 30M, 356 -> 320 ms; random u64, 100M, 1190 -> 1073 (stable 1796 ->
    /// 1211). With a short LSD leaf the cut HURTS (i32, 100M: 562 -> 606),
    /// which is why the rule looks at the leaf.
    pub radix_recut: usize,
    /// Radix: a bucket with more than 1/`radix_skew` of the data reveals a
    /// skewed key; nothing is moved and the comparison path takes over.
    pub radix_skew: usize,
    /// Radix: a bucket smaller than this goes straight to the std sort.
    pub radix_leaf_min: usize,
    /// Radix: with a radix key and a cheap comparison, the merge tree is only used when at most
    /// this many runs remain after chaining. Measured from 320K to 5M elements, 1 to 24 threads:
    /// with 8 runs the tree wins or ties; with 48 radix already wins or ties in every case.
    pub radix_max_runs: usize,
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
            radix_min: if e <= 4 { 131_072 } else { 524_288 },
            // Was 32_768. Measured again, against the stable block: at 8K radix wins
            // by 1.2x to 1.7x with 1 to 8 threads, at 16K by 1.1x to 2.5x.
            radix_min_stable: 8_192,
            radix_min_serial: 8_192,
            radix_min_thread: if e <= 8 { 16_384 } else { 32_768 },
            radix_digit: 8,
            radix_count_bits: 16,
            radix_lsd_cheap: 24,
            radix_lsd_costly: 40,
            radix_piece: 65_536,
            radix_big: 2 << 20,
            radix_big_div: 16,
            radix_recut: 65_536,
            radix_skew: 4,
            radix_leaf_min: 128,
            radix_max_runs: 32,
            route: Route::Auto,
        }
    }

    /// Below how many elements the chaotic path goes by comparison. It depends on who
    /// radix's competitor is on that data.
    fn radix_floor(&self, stable: bool, kind: Chaos) -> usize {
        let threads = rayon::current_num_threads().max(1);
        // Sorted and merged blocks: the stable engine's path, and the unstable engine's when
        // the quicksort is not suitable (short runs, fragmented data). Radix beats them from
        // early on. With many threads the floor rises up to the old value: with 24 threads
        // radix loses at 8K elements (0.7x) and ties at 16K.
        if stable || kind != Chaos::Random {
            return self.radix_min_stable.max(1024 * threads).min(4 * self.radix_min_stable);
        }
        // Parallel quicksort.
        if threads == 1 {
            self.radix_min_serial
        } else {
            self.radix_min.min(self.radix_min_thread.saturating_mul(threads))
        }
    }

    /// How many segments an array of `n` elements may have for the tree to
    /// still pay off.
    fn allowed_segs(&self, n: usize, cheap: bool) -> usize {
        let avg = if cheap { self.struct_min_avg_run_cheap } else { self.struct_min_avg_run };
        (n / avg.max(1)).max(self.min_segs)
    }
}

// =====================================================================
//  ENTRY
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
    sort_with_key(arr, cmp, stable, t, scratch, &NoKey);
}

/// `sort_with` with the key strategy of the chaotic path: `NoKey` (comparison
/// only), `Bare`, `Exact` or `Prefix`.
pub(crate) fn sort_with_key<T, F, K>(
    arr: &mut [T],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
    key: &K,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
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
    // With a bare key the element is an integer and the comparison is cheap by definition:
    // the clock does not decide.
    let cheap = t.force_cheap.unwrap_or(probe.cheap || K::BARE);
    if rate > t.chaos_percent {
        let kind = if rate > t.random_percent { Chaos::Random } else { Chaos::ShortRuns };
        sort_chaotic(arr, cmp, stable, t, scratch, kind, cheap, key);
        return;
    }
    // UNSTABLE SHORTCUT: globally ascending data, with scattered disorder
    // and a cheap comparison, is the best case of the parallel quicksort --
    // sorted with 1% swaps on 10M u64 finishes in 55 ms, vs 113-128 for any
    // merge path. It requires disorder in the sample (otherwise detection
    // solves it in one pass) and a cheap comparison (with `&str` the
    // same data takes 733 ms in it vs 367 in the tree).
    if !stable
        && t.route == Route::Auto
        && cheap
        && rate >= 1
        && globally_ascending(arr, cmp)
    {
        arr.par_sort_unstable_by(|a, b| cmp(a, b));
        return;
    }
    // Detection gives up on its own when the data turns out too fragmented
    // for the tree: the cost of having asked stays limited to a fraction of
    // the scan, instead of the whole pass.
    let limit = match t.route {
        Route::Auto => 2 * t.allowed_segs(n, cheap) + t.abort_slack,
        _ => usize::MAX,
    };
    match detect_runs(arr, cmp, t.detect_block, t.min_run, limit) {
        Some(runs) => sort_structured(arr, &runs, cmp, stable, t, scratch, cheap, key),
        None => sort_chaotic(arr, cmp, stable, t, scratch, Chaos::Fragmented, cheap, key),
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
//  PHASE 0: SHIELD
//
//  Eight windows spread across the whole array, not one window in the middle:
//  a single spot let 0.004% of the data decide the path of the other 99.996%.
//
//  Ties do NOT count as direction. 4 4 3 3 2 2 is a perfect descending run,
//  and counting `<=` as "going up" classified it as pure chaos.
//
//  The verdict is no longer critical: erring toward "structured" costs one
//  detection pass, because short runs are grouped and sorted as a block
//  anyway; erring toward "chaotic" costs the same pass, because the blocks
//  that are already sorted are merged without moving data.
// =====================================================================

/// What the shield's sample says about the data.
struct Probe {
    /// Percentage of comparisons in which the direction reversed. Random data of
    /// high cardinality gives ~66; with two keys, ~50; a sawtooth of
    /// length L gives 200/L. Sorted, reversed and plateaus give zero.
    rate: usize,
    /// Is the comparison cheap? Measured on the clock, over the same 512
    /// comparisons: an integer costs 1 to 3 ns, an out-of-cache string goes over
    /// 50. The error has only one possible direction -- a preemption in the middle
    /// of the sample makes cheap look expensive -- and that is the conservative
    /// direction: the engine forgoes a shortcut, it never uses one it should not.
    cheap: bool,
}

fn sample<T, F>(arr: &[T], cmp: &F, cheap_ns: usize) -> Probe
where
    T: Copy,
    F: Fn(&T, &T) -> Ordering,
{
    const WINDOWS: usize = 8;
    const SPAN: usize = 64;
    let n = arr.len();
    if n < 3 {
        return Probe { rate: 0, cheap: false };
    }
    let span = SPAN.min(n - 1);
    let last_start = n - 1 - span;
    // The windows are read once BEFORE the clock starts. The clock measures the cost of the
    // comparison, and until now it also measured the cache miss of bringing each window from
    // memory, which is a cost of reading the array: with u64 (one-instruction comparison) the
    // 512 comparisons took 4.7 us cold vs 1.3 warm, the limit is 5.1, and the verdict came out
    // "expensive" in 19 of 60 calls on 100M (25 of 60 on 10M on an i7). With the wrong verdict
    // the unstable engine loses the quicksort shortcut: 100M nearly sorted, ~735 ms in 7 of 12
    // calls vs ~295 in the others. What the comparison fetches OUTSIDE the element (the bytes
    // of a string) is still cold and still counts: &str of 12 to 27 letters gives 30 to 60 us
    // in both forms, "expensive" in 60 of 60.
    for w in 0..WINDOWS {
        let start = last_start / (WINDOWS - 1) * w;
        for x in &arr[start..=start + span] {
            std::hint::black_box(*x);
        }
        if last_start == 0 {
            break;
        }
    }
    let clock = std::time::Instant::now();
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

/// Is the array globally ascending, apart from out-of-place elements?
///
/// 256 equally spaced samples, in 16 groups; the MEDIAN of each group
/// must not decrease. The median ignores isolated out-of-place elements
/// (sorted with 1% swaps passes), and interleaved runs fail: at each
/// change of run the level drops and stays down there.
///
/// It serves a single decision: whether the unstable engine may hand
/// fragmented data to the quicksort. See `sort_chaotic`.
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
    // A random position inside each step, instead of the start of the step. With a fixed
    // step the sample falls in phase with periodic data: on 100M u64 with 64 interleaved runs
    // there are exactly 4 samples per run, always at the same fractions of the run, the 16
    // medians come out nearly equal and do not decrease, and the function answers "ascending"
    // for the worst shape the quicksort can receive (762 ms vs 382 for the tree). Until
    // 2026-11-20 the question was only asked for data with thousands of short runs; since data
    // with more than `radix_max_runs` long runs also goes to the chaotic path, the question
    // comes up with dozens of runs, and that is where the phase effect shows up.
    let mut seed = 0xA5A5_1F3D_9E37_79B9u64 ^ n as u64;
    let mut prev: Option<T> = None;
    for g in 0..GROUPS {
        let mut s = [arr[0]; PER];
        for (k, slot) in s.iter_mut().enumerate() {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            *slot = arr[(g * PER + k) * step + (seed % step as u64) as usize];
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
//  PHASE 1: RUN DETECTION
//
//  A run is described by (kind, length) packed into 64 bits.
//
//    ASC        non-decreasing
//    FLAT       all equal -- neutral: serves an ascending neighbor as well
//               as a descending one
//    DESC       strictly decreasing
//    DESC_TIES  non-increasing WITH repeated keys
//
//  The direction of a run is decided by the FIRST strictly different pair,
//  not by the first pair. Starting with `<=` classified every plateau as
//  ascending, and a reversed array with repeated keys was never seen as
//  reversed: it became one run per plateau.
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

/// Decomposes `arr` into runs, appending to `out`.
///
/// Returns an estimate of how many SEGMENTS this stretch would cost the tree:
/// a stretch of short runs counts once (it becomes a sorted block), and a long
/// run that chains with the previous one does not count (it will be merged
/// without moving data). It is this count that feeds the early abort of detection.
fn scan_runs<T, F>(arr: &[T], cmp: &F, out: &mut Vec<Run>, min_run: usize) -> usize
where
    F: Fn(&T, &T) -> Ordering,
{
    let n = arr.len();
    let first = out.len();
    let mut head = 0usize;
    let mut count = SegCount { segs: 0, have_prev: false, prev_long: false, prev_max: 0, in_short: false };
    while head < n {
        // Initial plateau: walk up to the first strictly different pair.
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
            // STOLEN ELEMENT. The greedy scan hands the previous ascending run
            // the first element of this descending run -- which is its
            // MAXIMUM. In ... 3 9 | 8 7 6 5 4 the 9 stays on the left, and then
            // the ascending run ends at 9 and never chains with the descending
            // one, whose minimum is 4. Giving the 9 back leaves 3 <= 4 and the
            // two become a single run. Without this the sawtooth over disjoint
            // ranges ended up with one run per change of direction (measured:
            // 1.25x gain vs 11x with the element given back).
            //
            // The theft can take MORE than one element: the ascending scan
            // accepts ties, so in ... 3 9 9 | 8 7 6 both 9s stay on the
            // left. Giving back only one left the ascending run ending at 9 and
            // the whole array went to a merge (measured on int32, 30M, 1000
            // teeth: 13.6 ms without the duplicate, 66 ms with ONE duplicate at
            // the top of a tooth). The whole final plateau is given back; with
            // two or more equals the descending run now has a repeated key.
            let mut start = head;
            if out.len() > first {
                let last = out.len() - 1;
                let prev = out[last];
                if prev.kind() == ASC && prev.len() >= 2 {
                    let mut k = 1;
                    while k < prev.len() - 1 && cmp(&arr[head - 1 - k], &arr[head - 1]) == Equal {
                        k += 1;
                    }
                    if cmp(&arr[head - 1 - k], &arr[tail - 1]) != Greater {
                        let l = prev.len() - k;
                        out[last] = Run::new(if l == 1 { FLAT } else { ASC }, l);
                        start = head - k;
                        count.prev_max = head - 1 - k;
                        if k > 1 {
                            ties = true;
                        }
                    }
                }
            }
            out.push(Run::new(if ties { DESC_TIES } else { DESC }, tail - start));
            count.add(arr, cmp, start, tail, true, min_run);
        }
        head = tail;
    }
    count.segs
}

/// Segment accounting of a block, for the early abort.
struct SegCount {
    segs: usize,
    have_prev: bool,
    /// Was the previous run long?
    prev_long: bool,
    /// Index of the LARGEST element of the previous run.
    prev_max: usize,
    /// Are we inside a stretch of short runs (already counted)?
    in_short: bool,
}

impl SegCount {
    /// Records the run `arr[start..end]`.
    fn add<T, F>(&mut self, arr: &[T], cmp: &F, start: usize, end: usize, desc: bool, min_run: usize)
    where
        F: Fn(&T, &T) -> Ordering,
    {
        let (min_i, max_i) = if desc { (end - 1, start) } else { (start, end - 1) };
        if end - start >= min_run {
            // The first run of the block does not count: from here there is no
            // way to know whether it chains with the previous block.
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

/// Stitches two runs that SHARE an element: the last one of `l` and the
/// first one of `r`. Works on the metadata only, without reading the array --
/// it is the monoid operation. Returns the left run and, if any is left, the right one.
///
/// Associativity holds up to equivalent decompositions: different groupings
/// may cut an isolated element at different places, and all the answers
/// describe the same array correctly.
fn join_runs(l: Run, r: Run) -> (Run, Option<Run>) {
    let (lk, ll, rk, rl) = (l.kind(), l.len(), r.kind(), r.len());
    let kind = match (lk, rk) {
        (FLAT, FLAT) => Some(FLAT),
        (ASC, ASC) | (ASC, FLAT) | (FLAT, ASC) => Some(ASC),
        // FLAT next to a descending run: it joins the run. If the plateau has
        // more than one element, the run now has repeated keys.
        (FLAT, _) => Some(if ll > 1 { DESC_TIES } else { rk }),
        (_, FLAT) => Some(if rl > 1 { DESC_TIES } else { lk }),
        // Opposite directions: there is no way to join them.
        (ASC, _) | (_, ASC) => None,
        _ => Some(if lk == DESC_TIES || rk == DESC_TIES { DESC_TIES } else { DESC }),
    };
    match kind {
        Some(k) => (Run::new(k, ll + rl - 1), None),
        None => {
            // The shared element stays with the left one.
            let rest = match rl - 1 {
                0 => None,
                1 => Some(Run::new(FLAT, 1)),
                m => Some(Run::new(rk, m)),
            };
            (l, rest)
        }
    }
}

/// Stitches, in place, the run at `mark - 1` with the run at `mark`.
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

/// Parallel detection. Blocks with 1 element of overlap, reduced by the
/// monoid. Rayon's `fold` gives each task ONE vector that accumulates its
/// consecutive blocks, instead of one vector per block.
///
/// The previous version had two levels (macro 32768 / micro 512). Measured on
/// sorted data, one level and two levels tie, and both tie with a raw parallel
/// `is_sorted`: the scan is sequential and bandwidth is what counts. The micro
/// level only cost -- a `Vec` every 512 elements.
///
/// `seg_limit` is the early abort: when the sum of the estimated segments
/// goes past the limit, the remaining blocks are not scanned and the function
/// returns `None`. Nearly sorted data with 1% swaps, for example, blows the
/// limit after scanning a tiny fraction of the array.
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
//  PHASE 2: NORMALIZATION
//
//  Turns the list of runs into ASCENDING segments, in place:
//
//    Keep      already ascending
//    Reverse   descending -> reverse (and un-reverse the ties, if stable)
//    Sort      several short runs in a row -> sort the stretch in one go
//
//  After that the tree only sees ascending runs, and direction is no longer
//  its concern.
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
    /// Elements that fell into `Sort` segments.
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
            // An isolated short run is already monotonic: it needs no sort.
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

/// Groups consecutive short runs into blocks to be sorted.
///
/// Merging run by run only pays off when the run is long. Measured on a
/// sawtooth of length L: with L=16 sorting as a block beats the tree, with L=32
/// it ties, with L>=64 the tree wins. A stretch of short runs becomes a block of
/// up to `chunk` elements -- and a random array that slipped past the shield
/// simply becomes a sequence of blocks, with no cliff.
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

/// Applies the operations in parallel, each one on its own slice.
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

/// `s` is sorted, but each block of equal keys ended up in the REVERSE of the
/// original order (it just came out of a `reverse`). Un-reverses each block.
///
/// Since `s` is already sorted, the edges of the block that contains the middle
/// come out of a binary search, and the two halves proceed in parallel.
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

/// Merges, in the offsets only, adjacent ascending segments that are already in
/// global order (last of one <= first of the next). One comparison per
/// boundary, no element moved.
///
/// Before, the "already in order" shortcut lived INSIDE the tree and was a memcpy
/// per level: a thousand teeth over disjoint ranges paid log(1000) copy passes.
/// Measured on 10M u64: 112 ms -> 9.9 ms.
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

#[allow(clippy::too_many_arguments)]
fn sort_structured<T, F, K>(
    arr: &mut [T],
    runs: &[Run],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    scratch: Option<&mut Vec<T>>,
    cheap: bool,
    key: &K,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    let n = arr.len();
    let plan = plan_segments(runs, t.min_run, t.chunk);
    let offsets = prefix_sums(&plan.lens);

    if t.route == Route::Auto {
        // The exact decision, with all the metadata at hand: how many runs remain
        // AFTER chaining? A thousand teeth over disjoint ranges are a thousand runs and
        // one single chained run; a thousand interleaved teeth really are a thousand runs.
        let allowed = t.allowed_segs(n, cheap);
        let fragmented = plan.lens.len() > allowed
            && runs_after_coalesce(arr, &plan, &offsets, cmp) > allowed;
        if fragmented {
            sort_chaotic(arr, cmp, stable, t, scratch, Chaos::Fragmented, cheap, key);
            return;
        }
        // With a radix key and a cheap comparison, merging run by run only beats radix when
        // the runs are few: each level of the tree moves the whole array, and radix solves it
        // in two passes whatever the order. Only from 4 scatter pieces on (262K elements): at
        // 20K and 80K the tree fits in cache and wins (up to 2x with 4 to 8 threads).
        if K::ENABLED
            && cheap
            && n >= 4 * t.radix_piece
            && plan.lens.len() > t.radix_max_runs
            && runs_after_coalesce(arr, &plan, &offsets, cmp) > t.radix_max_runs
        {
            sort_chaotic(arr, cmp, stable, t, scratch, Chaos::Fragmented, cheap, key);
            return;
        }
        // The shield was fooled and most of the data has no run worth using.
        // The stable engine carries on: the plan's `Sort` blocks already are
        // its chaotic path. The unstable engine chooses by the average length
        // of the raw run -- random gives ~2.4, a short sawtooth gives the
        // tooth length itself.
        if !stable && plan.sorted_elems * 2 > n {
            let kind = if runs.len() * 6 > n { Chaos::Random } else { Chaos::ShortRuns };
            sort_chaotic(arr, cmp, stable, t, scratch, kind, cheap, key);
            return;
        }
    }

    apply_ops(arr, &plan, cmp, stable);
    let offsets = coalesce(arr, offsets, cmp);
    merge_segments(arr, &offsets, cmp, t, scratch);
}

/// How many runs would remain after normalizing and chaining -- computed BEFORE
/// moving any element, using only the extremes of each segment. A
/// `Sort` segment has no known extremes and counts as a break.
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
//  CHAOTIC PATH
// =====================================================================

/// Why the data ended up on the chaotic path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Chaos {
    /// Truly random: the direction reverses in more than a third of the sample.
    Random,
    /// Short runs, but with local order (sawtooth of 16 to 60).
    ShortRuns,
    /// Real runs, only too many of them for the tree.
    Fragmented,
}

#[allow(clippy::too_many_arguments)]
fn sort_chaotic<T, F, K>(
    arr: &mut [T],
    cmp: &F,
    stable: bool,
    t: &Tuning,
    mut scratch: Option<&mut Vec<T>>,
    kind: Chaos,
    cheap: bool,
    key: &K,
) where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    // RADIX: truly random data and a known key. Holds for both
    // engines. If the routine gives up (skewed key, no memory for the
    // buffer) nothing was moved and the comparison path below proceeds unchanged.
    // Short runs and fragmented data still go by comparison: there is local
    // order there that radix does not exploit.
    // The quicksort's best case: fragmented but globally ascending data, with a
    // cheap comparison. It beats radix there, so it is checked first.
    let ascending_noise = !stable && cheap && kind == Chaos::Fragmented && globally_ascending(arr, cmp);
    // With a radix key and a cheap comparison, short runs and fragmented data also go to
    // radix: it costs almost the same in any order, and comparison blocks take twice as long.
    let radix_fits = kind == Chaos::Random || (cheap && !ascending_noise);
    // The floor depends on who radix's competitor is on this data: see `radix_floor`.
    let radix_min = t.radix_floor(stable, kind);
    if K::ENABLED && radix_fits && arr.len() >= radix_min {
        // In the stable engine the comparison leaf is the stable `sort`, more expensive
        // than `sort_unstable`: LSD pays off with more passes even with a
        // cheap comparison (6-byte text in u64, 30M: 384 ms with an LSD
        // leaf vs 554 with a comparison leaf).
        let rx = Radix {
            cmp,
            key,
            stable,
            t,
            lsd_max: if cheap && !stable { t.radix_lsd_cheap } else { t.radix_lsd_costly },
        };
        if radix_sort(arr, &rx, scratch.as_mut().map(|v| &mut **v)) {
            return;
        }
    }
    // UNSTABLE: the parallel quicksort works in place and without a buffer, but
    // only in two cases.
    //
    //  - Random data. Against sorting blocks and merging, on 10M u64: 124 ms
    //    vs 138-154; 52 vs 64-82 (64 keys); 18 vs 50-59
    //    (2 keys).
    //  - Globally ascending data with out-of-place elements AND a cheap
    //    comparison (see the shortcut in `sort_with`).
    //
    // Outside of that it is NOT safe. Interleaved runs are a bad pattern for
    // pivot selection, and with an expensive key the damage is an order of
    // magnitude: 8 runs of 250K strings take 1942 ms in it vs 109 ms in the
    // tree; 64 runs, 1121 ms vs 175. And on a short u64 sawtooth the blocks
    // already win (L=32: 116 ms vs 148).
    if !stable {
        let quick = match kind {
            Chaos::Random => true,
            Chaos::ShortRuns => false,
            Chaos::Fragmented => ascending_noise,
        };
        if quick {
            arr.par_sort_unstable_by(|a, b| cmp(a, b));
            return;
        }
    }
    // BLOCKS: sorted by the std `sort` (driftsort, which exploits the runs
    // inside the block) and merged by the tree. It is the path of the stable
    // engine, and of the unstable one when the quicksort is not suitable.
    //
    // Against Rayon's `par_sort`, on 10M u64: 175 ms vs 191 (random),
    // 102 vs 185 (1024 keys), 89 vs 190 (64 keys), 78 vs 154
    // (2 keys). The block of 65536 won or tied across the whole sweep
    // (16K, 64K, 256K, n/8, n/4) -- the same value the C++ engine already used.
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
//  CHAOTIC PATH WITH A KEY: RADIX
//
//  When the caller says what the bits of the key are, random data does not
//  need comparison:
//
//    1. SAMPLE (for free): 1024 scattered elements. If one bucket of the high
//       bits would get more than 1/4 of them, the key is skewed and the
//       routine gives up before scanning the array.
//    2. RANGE: one parallel scan (AND and OR of the keys) tells which bits
//       vary. A "bare" key with up to 16 varying bits: COUNTING -- parallel
//       histogram and parallel fill, without a buffer.
//    3. HISTOGRAM of the 8 highest bits that vary. If a bucket would exceed
//       1/4 of the data (the sample was wrong), it gives up: nothing was moved yet.
//    4. Only then is the buffer allocated and the array SCATTERED into 256
//       buckets, in parallel: the destination is sliced beforehand into
//       disjoint ranges [bucket][piece], each piece writes into its own.
//    5. Each bucket, in parallel and inside the cache: sequential LSD radix
//       if few bits remain, std sort if many remain.
//    6. A bucket that is too big goes through the routine again.
//
//  Measured on 30M / 100M random, 8 threads, buffer allocated on the clock,
//  against the parallel quicksort (prototype of 2026-11-17):
//  i32 344 -> 152 ms / 1274 -> 528; u64 406 -> 278 / 1463 -> 1084; 64 values
//  on i32 132 -> 23; skewed (90% in 0..1000) ties.
//
//  Why MSD + leaf and not pure LSD: LSD scans the whole array twice per
//  key byte (16 scans on u64, where it loses to the quicksort: 540
//  vs 372 ms). The single MSD pass leaves the rest in cache and does not
//  depend on the width of the key.
//
//  STABILITY. The scatter walks the pieces in order and each bucket receives
//  first what came from piece 0, then from 1...: it preserves the order of
//  the equals. LSD is stable by construction. The comparison leaf uses the
//  stable `sort` when the engine is the stable one. Counting only exists for
//  a bare key, where equals are indistinguishable.
//
//  The routine returns `false` WITHOUT HAVING MOVED ANYTHING when it does not
//  apply; the caller goes on along the comparison path.
// =====================================================================

/// How the chaotic path sees the key of an element.
pub(crate) trait Keyer<T>: Sync {
    /// Is there a key? `NoKey` turns radix off at compile time.
    const ENABLED: bool;
    /// Equal keys imply elements equal under the comparison. Enables the LSD
    /// leaf. `false` for a prefix: there ties are always broken by comparison.
    const EXACT: bool;
    /// The element is the key itself and is rebuilt from the bits. Enables
    /// counting.
    const BARE: bool;
    fn bits(&self, x: &T) -> u64;
    fn rebuild(&self, bits: u64) -> T;
}

/// No key: comparison only. It is what the entries without `radix` in the name use.
pub(crate) struct NoKey;
impl<T> Keyer<T> for NoKey {
    const ENABLED: bool = false;
    const EXACT: bool = false;
    const BARE: bool = false;
    #[inline(always)]
    fn bits(&self, _: &T) -> u64 {
        0
    }
    fn rebuild(&self, _: u64) -> T {
        unreachable!("NoKey does not rebuild an element")
    }
}

/// The element is an integer: it is the key.
pub(crate) struct Bare;
impl<T: RadixKey> Keyer<T> for Bare {
    const ENABLED: bool = true;
    const EXACT: bool = true;
    const BARE: bool = true;
    #[inline(always)]
    fn bits(&self, x: &T) -> u64 {
        x.to_radix()
    }
    #[inline(always)]
    fn rebuild(&self, bits: u64) -> T {
        T::from_radix(bits)
    }
}

/// The extracted key decides the whole order; the element may carry more data.
pub(crate) struct Exact<K>(pub K);
impl<T, K: Fn(&T) -> u64 + Sync> Keyer<T> for Exact<K> {
    const ENABLED: bool = true;
    const EXACT: bool = true;
    const BARE: bool = false;
    #[inline(always)]
    fn bits(&self, x: &T) -> u64 {
        (self.0)(x)
    }
    fn rebuild(&self, _: u64) -> T {
        unreachable!("an extracted key does not rebuild an element")
    }
}

/// The extracted key is only the start of the order; the comparator breaks ties.
pub(crate) struct Prefix<K>(pub K);
impl<T, K: Fn(&T) -> u64 + Sync> Keyer<T> for Prefix<K> {
    const ENABLED: bool = true;
    const EXACT: bool = false;
    const BARE: bool = false;
    #[inline(always)]
    fn bits(&self, x: &T) -> u64 {
        (self.0)(x)
    }
    fn rebuild(&self, _: u64) -> T {
        unreachable!("a prefix does not rebuild an element")
    }
}

/// What the routine carries from one level to the next.
struct Radix<'a, F, K> {
    cmp: &'a F,
    key: &'a K,
    stable: bool,
    t: &'a Tuning,
    /// LSD leaf: up to how many remaining bits.
    lsd_max: u32,
}

/// Comparison sort of a large stretch, when radix is not suitable for it.
fn radix_fallback<T, F>(a: &mut [T], cmp: &F, stable: bool)
where
    T: Send,
    F: Fn(&T, &T) -> Ordering + Sync,
{
    if stable {
        a.par_sort_by(|x, y| cmp(x, y));
    } else {
        a.par_sort_unstable_by(|x, y| cmp(x, y));
    }
}

/// Sorts `arr` by radix. Returns `false`, without having moved any element,
/// when the key is skewed, when the prefix tells nothing apart or when
/// there is no memory for the buffer.
fn radix_sort<T, F, K>(arr: &mut [T], rx: &Radix<F, K>, scratch: Option<&mut Vec<T>>) -> bool
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    let n = arr.len();
    if n < 2 {
        return true;
    }
    let t = rx.t;
    // A skewed key (one bucket with more than 1/4 of the data) makes radix give up, because
    // the scatter does not divide the work. In the unstable engine the alternative is the
    // quicksort, which handles repeated keys in passing. In the STABLE one the alternative is
    // the blocks sorted by comparison and merged, and there are two cases in which radix
    // finishes in ONE pass and beats them:
    //
    //  a single key holds more than half of the data: it is done in a three-way split. From 4
    //          scatter pieces on (at 300K: 1.2x to 2.9x).
    //  one digit covers all the bits that vary: the scatter already leaves everything sorted.
    //          From 8 pieces on (two and three keys: 0.8x to 1.1x at 300K; 1.0x to 1.9x at 1M).
    //
    // Outside of them (heavy tail, concentration in a wide range) radix would need several
    // levels, each moving almost everything, and it loses: 0.44x to 0.92x. It still gives up.
    let one_key = K::EXACT && rx.stable && n >= 4 * t.radix_piece;
    let one_digit = one_key && n >= 8 * t.radix_piece;
    let skew = radix_skewed(arr, rx.key, t);
    match skew {
        Skew::No => {}
        Skew::OneKey if one_key => {}
        Skew::OneDigit if one_digit => {}
        _ => return false,
    }
    // The sample saw a dominant key: counting it comes before any other
    // scan. If the count contradicts the sample, giving up cost one read pass.
    let mut pivot = None;
    if skew == Skew::OneKey {
        match radix_pivot_plan(arr, rx.key, t) {
            Some(p) if p.eq == n => return true, // a single key: nothing to sort
            Some(p) => pivot = Some(p),
            None => return false,
        }
    }
    let mut plan = None;
    if pivot.is_none() {
        let (and, or) = radix_range(arr, rx.key);
        let diff = and ^ or;
        if diff == 0 {
            // All keys equal. With an exact key there is nothing to sort; with a
            // prefix, it tells nothing apart and the comparison takes over.
            return K::EXACT;
        }
        let hi = 64 - diff.leading_zeros();
        if K::BARE && hi <= t.radix_count_bits {
            radix_count(arr, rx.key, or, hi);
            return true;
        }
        // The histogram comes BEFORE the buffer: if the key is skewed, giving up
        // did not even cost the allocation.
        plan = radix_plan(arr, rx.key, t, hi, one_digit);
        // The sample did not see it, but the exact histogram showed concentration: it is worth
        // checking whether it is a single key, also before the buffer.
        if plan.is_none() && one_key {
            pivot = radix_pivot_plan(arr, rx.key, t);
        }
        if plan.is_none() && pivot.is_none() {
            return false;
        }
    }
    let mut owned: Vec<T> = Vec::new();
    let buf: &mut [T] = match scratch {
        Some(v) if v.len() >= n => &mut v[..n],
        Some(v) => {
            v.clear();
            if v.try_reserve_exact(n).is_err() {
                return false;
            }
            arr.par_iter().copied().collect_into_vec(v);
            &mut v[..]
        }
        None => {
            if owned.try_reserve_exact(n).is_err() {
                return false;
            }
            arr.par_iter().copied().collect_into_vec(&mut owned);
            &mut owned[..]
        }
    };
    match (&plan, &pivot) {
        (Some(p), _) => {
            radix_scatter(arr, buf, p, rx.key);
            radix_buckets(arr, buf, p, rx, 0);
        }
        (None, Some(pv)) => radix_pivot(arr, buf, pv, rx, 0),
        (None, None) => unreachable!("without a plan, radix has already given up"),
    }
    true
}

/// What the sample says about the key.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Skew {
    /// No bucket of the high bits exceeds 1/`radix_skew` of the sample.
    No,
    /// Skewed, but one digit covers all the bits that vary.
    OneDigit,
    /// Skewed, and a single key holds at least half of the sample.
    OneKey,
    /// Skewed and spread over more than one digit, with no dominant key.
    Wide,
}

/// Is the key skewed? Looks at 1024 scattered elements: if their range is
/// wide and one bucket of the high bits gets more than 1/`radix_skew` of the
/// sample, the scatter would not divide the work.
fn radix_skewed<T, K>(a: &[T], key: &K, t: &Tuning) -> Skew
where
    K: Keyer<T>,
{
    let n = a.len();
    let m = 1024.min(n);
    let step = n / m;
    // A random position inside each step, instead of the start of the step. With a fixed step,
    // periodic data fools the sample: 512 runs of the same length in 640K elements put half of
    // the samples on the SMALLEST value of each run, and the key looks skewed when it is not.
    let mut seed = 0x9E37_79B9_7F4A_7C15u64 ^ n as u64;
    let mut keys: Vec<u64> = (0..m)
        .map(|i| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            key.bits(&a[i * step + (seed % step as u64) as usize])
        })
        .collect();
    let (mut and, mut or) = (u64::MAX, 0u64);
    for &b in &keys {
        and &= b;
        or |= b;
    }
    let hi = 64 - (and ^ or).leading_zeros();
    if K::BARE && hi <= t.radix_count_bits {
        return Skew::No; // counting candidate; the exact scan decides
    }
    if hi == 0 {
        return Skew::OneKey; // the whole sample has the same key
    }
    let b = t.radix_digit.clamp(1, hi.min(16));
    let shift = hi - b;
    let mask = (1u64 << b) - 1;
    let mut h = vec![0usize; 1usize << b];
    for &k in &keys {
        h[((k >> shift) & mask) as usize] += 1;
    }
    if !h.iter().any(|&c| c * t.radix_skew > m) {
        return Skew::No;
    }
    // Skewed. Does radix still finish in one pass?
    if shift == 0 {
        return Skew::OneDigit;
    }
    keys.sort_unstable();
    let (mut best, mut i) = (0usize, 0usize);
    while i < m {
        let mut j = i + 1;
        while j < m && keys[j] == keys[i] {
            j += 1;
        }
        best = best.max(j - i);
        i = j;
    }
    if best * 2 >= m {
        Skew::OneKey
    } else {
        Skew::Wide
    }
}

/// AND and OR of all the keys: the bits in which they differ are the ones that vary.
fn radix_range<T, K>(a: &[T], key: &K) -> (u64, u64)
where
    T: Sync,
    K: Keyer<T>,
{
    a.par_chunks(1 << 16)
        .map(|c| {
            c.iter().fold((u64::MAX, 0u64), |(and, or), x| {
                let b = key.bits(x);
                (and & b, or | b)
            })
        })
        .reduce(|| (u64::MAX, 0u64), |l, r| (l.0 & r.0, l.1 | r.1))
}

/// Bare key with a small range: counts and fills, without a buffer. `or` is the OR
/// of all the keys and `hi` how many low bits vary.
fn radix_count<T, K>(a: &mut [T], key: &K, or: u64, hi: u32)
where
    T: Copy + Send + Sync,
    K: Keyer<T>,
{
    let n = a.len();
    let nb = 1usize << hi;
    let low = (nb - 1) as u64;
    let base = or & !low; // the upper bits are equal in all the keys
    let pieces = (n >> 20).clamp(1, 16);
    let chunk = (n + pieces - 1) / pieces;
    let hists: Vec<Vec<usize>> = a
        .par_chunks(chunk)
        .map(|c| {
            let mut h = vec![0usize; nb];
            for x in c {
                h[(key.bits(x) & low) as usize] += 1;
            }
            h
        })
        .collect();
    let mut parts: Vec<(&mut [T], T)> = Vec::new();
    let mut rest: &mut [T] = a;
    for d in 0..nb {
        let c: usize = hists.iter().map(|h| h[d]).sum();
        if c > 0 {
            let (s, r) = rest.split_at_mut(c);
            parts.push((s, key.rebuild(base | d as u64)));
            rest = r;
        }
    }
    parts.into_par_iter().for_each(|(s, v)| {
        if s.len() > (1 << 17) {
            s.par_chunks_mut(1 << 16).for_each(|c| c.fill(v));
        } else {
            s.fill(v);
        }
    });
}

/// Histogram of the MSD pass over a stretch.
struct RadixPlan {
    /// Bits below the digit: what is left for the leaf.
    shift: u32,
    mask: u64,
    nb: usize,
    /// Elements per piece.
    chunk: usize,
    /// One histogram per piece.
    hists: Vec<Vec<usize>>,
    /// Bucket boundaries in the destination (`nb + 1` values).
    bounds: Vec<usize>,
}

/// Counts the MSD digit -- the `radix_digit` highest bits among the `hi` that
/// vary. `None` if a bucket would get more than 1/`radix_skew` of the stretch.
///
/// With `one_digit`, the dominant bucket is accepted when the digit covers all the bits that
/// vary: after the scatter nothing is left to sort.
fn radix_plan<T, K>(a: &[T], key: &K, t: &Tuning, hi: u32, one_digit: bool) -> Option<RadixPlan>
where
    T: Sync,
    K: Keyer<T>,
{
    let n = a.len();
    let b = t.radix_digit.clamp(1, hi.min(16));
    let shift = hi - b;
    let nb = 1usize << b;
    let mask = (nb - 1) as u64;
    let pieces = (n / t.radix_piece.max(1)).clamp(1, 128);
    let chunk = (n + pieces - 1) / pieces;
    let hists: Vec<Vec<usize>> = a
        .par_chunks(chunk)
        .map(|c| {
            let mut h = vec![0usize; nb];
            for x in c {
                h[((key.bits(x) >> shift) & mask) as usize] += 1;
            }
            h
        })
        .collect();
    let mut bounds = vec![0usize; nb + 1];
    for d in 0..nb {
        bounds[d + 1] = bounds[d] + hists.iter().map(|h| h[d]).sum::<usize>();
    }
    let limit = n / t.radix_skew.max(1);
    if !(one_digit && shift == 0) && (0..nb).any(|d| bounds[d + 1] - bounds[d] > limit) {
        return None;
    }
    Some(RadixPlan { shift, mask, nb, chunk, hists, bounds })
}

/// Scatters `src` into `dst` according to the plan, in parallel. The destination is
/// sliced beforehand: bucket by bucket and, inside the bucket, piece by piece. Each
/// piece receives its own `nb` ranges and writes through them -- without `unsafe`, and
/// in the order that keeps the equals as they came.
fn radix_scatter<T, K>(src: &[T], dst: &mut [T], p: &RadixPlan, key: &K)
where
    T: Copy + Send + Sync,
    K: Keyer<T>,
{
    let mut slots: Vec<Vec<&mut [T]>> =
        (0..p.hists.len()).map(|_| Vec::with_capacity(p.nb)).collect();
    let mut rest: &mut [T] = dst;
    for d in 0..p.nb {
        for (piece, h) in p.hists.iter().enumerate() {
            let (a, r) = rest.split_at_mut(h[d]);
            slots[piece].push(a);
            rest = r;
        }
    }
    src.par_chunks(p.chunk).zip(slots.into_par_iter()).for_each(|(c, mut cur)| {
        for x in c {
            let d = ((key.bits(x) >> p.shift) & p.mask) as usize;
            let s = std::mem::take(&mut cur[d]);
            if let Some((head, tail)) = s.split_first_mut() {
                *head = *x;
                cur[d] = tail;
            }
        }
    });
}

/// Sorts each bucket. The data is in `buf`; the result ends up in `arr`.
fn radix_buckets<T, F, K>(arr: &mut [T], buf: &mut [T], p: &RadixPlan, rx: &Radix<F, K>, depth: u32)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    let t = rx.t;
    let big = (arr.len() / t.radix_big_div.max(1)).max(t.radix_big);
    let shift = p.shift;
    // Narrow leaf = LSD with few passes: there only a truly big bucket
    // re-enters. Wide leaf (many passes or comparison): a second MSD
    // cut comes out cheaper than sorting a bucket that does not fit in cache.
    let narrow = K::EXACT && shift <= t.radix_lsd_cheap;
    let recut = if narrow { big } else { t.radix_recut.min(big) };
    let mut pairs: Vec<(&mut [T], &mut [T])> = Vec::with_capacity(p.nb);
    let (mut ra, mut rb): (&mut [T], &mut [T]) = (arr, buf);
    for w in p.bounds.windows(2) {
        let len = w[1] - w[0];
        let (a, r) = ra.split_at_mut(len);
        ra = r;
        let (s, r) = rb.split_at_mut(len);
        rb = r;
        if len > 0 {
            pairs.push((a, s));
        }
    }
    pairs.into_par_iter().for_each(|(a, s)| {
        let len = a.len();
        if K::EXACT && shift == 0 {
            // The digit covered all the bits that vary: identical keys.
            a.copy_from_slice(s);
        } else if len > recut && depth < 6 {
            // Big bucket (data with concentration, or a wide leaf): back
            // to `arr` and the routine again, parallel inside, with its
            // own piece of the buffer as scratch.
            par_copy(a, s, t.par_grain);
            radix_bucket(a, s, rx, depth + 1);
        } else if K::EXACT && len >= t.radix_leaf_min && shift <= rx.lsd_max {
            radix_lsd(s, a, rx.key, shift);
        } else {
            a.copy_from_slice(s);
            if len > big {
                radix_fallback(a, rx.cmp, rx.stable);
            } else {
                leaf_sort(a, rx.cmp, rx.stable);
            }
        }
    });
}

/// The routine again on a big bucket. The data is in `a`; `s` is scratch.
fn radix_bucket<T, F, K>(a: &mut [T], s: &mut [T], rx: &Radix<F, K>, depth: u32)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    let (and, or) = radix_range(a, rx.key);
    let diff = and ^ or;
    if diff == 0 {
        if !K::EXACT {
            radix_fallback(a, rx.cmp, rx.stable);
        }
        return;
    }
    let hi = 64 - diff.leading_zeros();
    if K::BARE && hi <= rx.t.radix_count_bits {
        radix_count(a, rx.key, or, hi);
        return;
    }
    // Here radix has already started and the alternative is sorting the whole bucket by
    // comparison: the two one-pass cases hold for a bucket of any size.
    let one_pass = K::EXACT && rx.stable;
    match radix_plan(a, rx.key, rx.t, hi, one_pass) {
        Some(p) => {
            radix_scatter(a, s, &p, rx.key);
            radix_buckets(a, s, &p, rx, depth);
        }
        None => {
            if one_pass && depth < 6 {
                if let Some(pv) = radix_pivot_plan(a, rx.key, rx.t) {
                    radix_pivot(a, s, &pv, rx, depth);
                    return;
                }
            }
            radix_fallback(a, rx.cmp, rx.stable)
        }
    }
}

/// The dominant key of a stretch and the exact count of the three parts.
struct PivotPlan {
    pivot: u64,
    /// Elements per piece.
    chunk: usize,
    /// Per piece: how many are smaller than, equal to and greater than `pivot`.
    hists: Vec<[usize; 3]>,
    lt: usize,
    eq: usize,
}

/// 0: smaller than the dominant key; 1: equal; 2: greater.
#[inline(always)]
fn pivot_side(bits: u64, pivot: u64) -> usize {
    (bits >= pivot) as usize + (bits > pivot) as usize
}

/// Does a single key hold at least half of the stretch? The candidate is the mode of 256
/// samples; the exact count, in parallel, confirms it. Moves nothing.
fn radix_pivot_plan<T, K>(a: &[T], key: &K, t: &Tuning) -> Option<PivotPlan>
where
    T: Sync,
    K: Keyer<T>,
{
    let n = a.len();
    if n < 2 {
        return None;
    }
    let m = 256.min(n);
    let step = n / m;
    let mut seed = 0xD1B5_4A32_D192_ED03u64 ^ n as u64;
    let mut keys: Vec<u64> = (0..m)
        .map(|i| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            key.bits(&a[i * step + (seed % step as u64) as usize])
        })
        .collect();
    keys.sort_unstable();
    let (mut pivot, mut best, mut i) = (keys[0], 0usize, 0usize);
    while i < m {
        let mut j = i + 1;
        while j < m && keys[j] == keys[i] {
            j += 1;
        }
        if j - i > best {
            best = j - i;
            pivot = keys[i];
        }
        i = j;
    }
    // Slack for the sample's error; the count is what decides.
    if best * 5 < m * 2 {
        return None;
    }
    let pieces = (n / t.radix_piece.max(1)).clamp(1, 128);
    let chunk = (n + pieces - 1) / pieces;
    let hists: Vec<[usize; 3]> = a
        .par_chunks(chunk)
        .map(|c| {
            let mut h = [0usize; 3];
            for x in c {
                h[pivot_side(key.bits(x), pivot)] += 1;
            }
            h
        })
        .collect();
    let lt: usize = hists.iter().map(|h| h[0]).sum();
    let eq: usize = hists.iter().map(|h| h[1]).sum();
    // The entry sample requires half; the count accepts from 40% on. With the same cut in
    // both, a key with exactly half of the data would pass the sample and fail here by chance,
    // paying for the scans without any gain (measured: 0.91x to 0.93x at 1M).
    if eq * 5 < n * 2 {
        return None;
    }
    Some(PivotPlan { pivot, chunk, hists, lt, eq })
}

/// Splits the stretch into smaller than, equal to and greater than the dominant key. The
/// data is in `a`; `s` is scratch; the result ends up in `a`.
///
/// The scatter walks the pieces in order, like the MSD pass: the equals come out done and in
/// the order they came. The outer parts, which sum to less than half, go back to the routine.
fn radix_pivot<T, F, K>(a: &mut [T], s: &mut [T], p: &PivotPlan, rx: &Radix<F, K>, depth: u32)
where
    T: Copy + Send + Sync,
    F: Fn(&T, &T) -> Ordering + Sync,
    K: Keyer<T>,
{
    let t = rx.t;
    {
        let mut slots: Vec<Vec<&mut [T]>> =
            (0..p.hists.len()).map(|_| Vec::with_capacity(3)).collect();
        let mut rest: &mut [T] = &mut *s;
        for d in 0..3 {
            for (piece, h) in p.hists.iter().enumerate() {
                let (head, r) = std::mem::take(&mut rest).split_at_mut(h[d]);
                slots[piece].push(head);
                rest = r;
            }
        }
        a.par_chunks(p.chunk).zip(slots.into_par_iter()).for_each(|(c, mut cur)| {
            for x in c {
                let d = pivot_side(rx.key.bits(x), p.pivot);
                let slot = std::mem::take(&mut cur[d]);
                if let Some((head, tail)) = slot.split_first_mut() {
                    *head = *x;
                    cur[d] = tail;
                }
            }
        });
    }
    let (a_lt, a_rest) = a.split_at_mut(p.lt);
    let (a_eq, a_gt) = a_rest.split_at_mut(p.eq);
    let (s_lt, s_rest) = s.split_at_mut(p.lt);
    let (s_eq, s_gt) = s_rest.split_at_mut(p.eq);
    let finish = |x: &mut [T], y: &mut [T]| {
        if x.is_empty() {
            return;
        }
        par_copy(x, y, t.par_grain);
        if x.len() > t.small_sort {
            radix_bucket(x, y, rx, depth + 1);
        } else {
            leaf_sort(x, rx.cmp, rx.stable);
        }
    };
    rayon::join(
        || par_copy(a_eq, s_eq, t.par_grain),
        || rayon::join(|| finish(a_lt, s_lt), || finish(a_gt, s_gt)),
    );
}

/// Sequential LSD radix of a bucket, over the lower `low_bits` of the key
/// (the upper ones are equal across the whole bucket). The data is in `src`; the
/// result ends up in `dst`. A pass over a constant digit is skipped.
fn radix_lsd<T, K>(src: &mut [T], dst: &mut [T], key: &K, low_bits: u32)
where
    T: Copy,
    K: Keyer<T>,
{
    let passes = (low_bits + 7) / 8;
    let mut in_src = true;
    for p in 0..passes {
        let sh = p * 8;
        let (a, b): (&[T], &mut [T]) =
            if in_src { (&*src, &mut *dst) } else { (&*dst, &mut *src) };
        let mut h = [0usize; 256];
        for x in a {
            h[((key.bits(x) >> sh) & 0xFF) as usize] += 1;
        }
        if h.iter().any(|&c| c == a.len()) {
            continue;
        }
        let mut pos = [0usize; 256];
        let mut acc = 0usize;
        for d in 0..256 {
            pos[d] = acc;
            acc += h[d];
        }
        for x in a {
            let d = ((key.bits(x) >> sh) & 0xFF) as usize;
            b[pos[d]] = *x;
            pos[d] += 1;
        }
        in_src = !in_src;
    }
    if in_src {
        dst.copy_from_slice(src);
    }
}

// =====================================================================
//  PHASE 3: MERGE TREE
//
//  Input: `offsets` delimits ASCENDING runs in `arr`.
//
//  THE BUFFER
//
//  Without a caller's buffer, it is born as a PARALLEL COPY of the array. That
//  trades two things from the previous version -- the sequential fill
//  (`vec![arr[0]; n]`, one thread writing n elements that would be
//  overwritten) and the copy of the odd-parity leaves -- for a single
//  parallel pass. With both sides equal at the start, every leaf already is
//  where it needs to be.
//
//  The measurement that had concluded that initializing the buffer "costs
//  nothing" (457.4 vs 458.5 ms) was confounded: when `arr[0]` is the integer
//  zero, `vec![arr[0]; n]` becomes `calloc` and the system hands out zeroed
//  pages on demand. With `arr[0] != 0` the fill costs ~14 ms per 80 MB
//  (64 runs on 10M u64: 64.0 ms with a zero in front, 77.5 without). Real
//  data almost never starts with zero.
//
//  With a caller's buffer (`*_with_buffer`), already initialized and of the
//  right size, there is nothing to create and no new page to touch: only the
//  odd-parity leaves copy their own stretch.
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

/// K-ary recursion over the runs. Hybrid: k-way on the levels that go to DRAM,
/// binary below `kway_min_elems`, where the stretch already lives in cache.
///
/// Contract: the input data is in `v`. On exit the sorted stretch is in
/// `buf` if `into_buf`, otherwise in `v`. If `copy_leaves` is false, `buf` was
/// born as a copy of `v` and the leaf has nothing to do; if it is true, `buf`
/// holds valid garbage and the odd-parity leaf copies its own stretch.
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

    // Cuts balanced by ELEMENT, not by run count.
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

    // Adjacent groups already in order become ONE stream.
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

/// Binary tree, with a split balanced by elements.
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
//  TWO-WAY MERGE: CO-RANK + BIDIRECTIONAL
// =====================================================================

/// Largest `i` such that `a[..i]` and `b[..k - i]` are the `k` smallest, with
/// ties going to `a`. It is the cut that keeps the merge stable.
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

/// FORWARD merge: emits the `dest.len()` smallest. Tie -> A.
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

/// BACKWARD merge: emits the `dest.len()` largest. Tie -> B.
///
/// The tie-breaking rule is the COMPLEMENT of the one in `merge_front`, and that is
/// what makes the theorem hold: the k smallest (tie -> A) and the total-k largest
/// (tie -> B) partition the multiset for ANY k. No cut point needs to be
/// computed -- each side only counts its own outputs.
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

/// Parallel merge of two runs, from immutable sources. Cuts by co-rank down to the
/// grain; in each piece it first tries the disjoint-range shortcuts and only
/// then merges -- with one thread coming from the front and another from the end.
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
    // Disjoint ranges: it is a concatenation.
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
//  K-WAY MERGE WITH CACHE BLOCKING
//
//  The binary merge pays log2(runs) passes over memory, and the machine
//  saturates the bandwidth: the cost is the NUMBER OF PASSES. A K-way merge
//  pays logK(runs).
//
//  It does NOT use a loser tree: measured, it costs ~12x more per element than
//  the simple merge, because `head[w]` with a run-time `w` is not promoted to a
//  register. Each tile is cut by `multiseq_partition` and merged by a BINARY
//  tree that fits in L1 -- the fast kernel stays the same.
// =====================================================================

/// Finds `out[i]` summing to `rank`, such that the prefixes are exactly the
/// `rank` smallest elements. Ties are distributed in stream order, which
/// keeps stability and agrees with the tie-breaking of `merge_front`.
///
/// The binary searches run only inside `[lo, hi)`. Every pivot comes from
/// inside the live interval, so it lies strictly between the previous pivots,
/// and no element outside the interval can change sides. The result is
/// identical to that of the search over the whole run; counted on 8 runs, the
/// comparisons per merge element drop from 6.24 to 5.22 (the ideal is 4.00).
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
        // Pivot: the middle of the widest interval.
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

/// One level of the tile's binary tree: merges pairs of contiguous segments of
/// `src` into `dst`. Returns the new segments (offset, len) and how many there are.
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

/// Merge of 3..=K runs with cache blocking: the DRAM sees one read and one
/// write, instead of log2(K) passes.
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

        // Restricted window: a tile consumes at most `want` from any run.
        let mut win: [&[T]; KWAY_FANOUT] = [&[]; KWAY_FANOUT];
        for i in 0..ns {
            let rem = &runs[i][cur[i]..];
            win[i] = &rem[..rem.len().min(want)];
        }
        let mut cut = [0usize; KWAY_FANOUT];
        multiseq_partition(&win[..ns], want, &mut cut[..ns], cmp);

        // Level 1: from the windows into s1.
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

        // Intermediate levels in s1/s2 (resident in L1); the LAST level
        // writes straight to the destination -- the only write that reaches DRAM.
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

/// Parallel k-way merge: cuts by rank and recurses.
///
/// Empty streams are dropped at each level. After a cut it is common for one
/// or two streams with data to remain in a piece; without dropping them they
/// went through the three levels of the tile as if they were eight. With ONE
/// stream it is a copy, with TWO it is the binary merge by co-rank -- which
/// scales with the threads, where before the two-run case was stuck on two.
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
