#ifndef SUBNETR_CORE_H
#define SUBNETR_CORE_H

#include <vector>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>

namespace subnetr {

static const double NEG_INF = -std::numeric_limits<double>::infinity();

// ---------------------------------------------------------------------------
// Objective functions for the peeling sequence.
//
// A candidate node set with n nodes carries total edge weight `total`.
//
//   OBJ_GENERALIZED (default)  total / n^(2*lambda), evaluated on the log
//                              scale as log(total) - 2*lambda*log(n). The
//                              exponent sets how hard size is penalized:
//                              lambda = 0 maximizes total weight and returns
//                              the whole graph, lambda = 0.5 maximizes average
//                              degree, and lambda approaching 1 maximizes edge
//                              density and collapses onto small cliques.
//   OBJ_DENSITY                total / (n(n-1)/2), pure edge density. Given no
//                              size floor this degenerates to the single
//                              densest pair, so it is only useful with a
//                              substantial min_size.
// ---------------------------------------------------------------------------
enum Objective { OBJ_GENERALIZED = 0, OBJ_DENSITY = 1 };

inline double obj_score(double total, int n, double lambda, int obj) {
  if (n < 2 || total <= 0.0) return NEG_INF;
  const double dn = static_cast<double>(n);
  if (obj == OBJ_DENSITY) return total / (0.5 * dn * (dn - 1.0));
  return std::log(total) - 2.0 * lambda * std::log(dn);
}

// ---------------------------------------------------------------------------
// Sparse symmetric graph in compressed-row form, holding only the edges that
// survived screening.
//
// Screening keeps a high quantile of the edge weights -- typically the top 1
// to 10 percent -- so the graph the algorithm actually works on is sparse by
// construction. A dense representation would spend nearly all of its time
// visiting zeros: at 200 nodes screened at the 99th percentile that is 19900
// matrix entries to reach 199 real edges.
// ---------------------------------------------------------------------------
struct Graph {
  int n = 0;
  std::vector<int>    ptr;  // n + 1 row offsets
  std::vector<int>    idx;  // 2E neighbour ids
  std::vector<double> wgt;  // 2E weights

  void build(int n_, const std::vector<int>& ei, const std::vector<int>& ej,
             const std::vector<double>& ew) {
    n = n_;
    const size_t E = ei.size();
    ptr.assign(n + 1, 0);
    for (size_t e = 0; e < E; ++e) { ++ptr[ei[e] + 1]; ++ptr[ej[e] + 1]; }
    for (int a = 0; a < n; ++a) ptr[a + 1] += ptr[a];
    idx.resize(ptr[n]);
    wgt.resize(ptr[n]);
    std::vector<int> fill(ptr.begin(), ptr.end() - 1);
    for (size_t e = 0; e < E; ++e) {
      const int a = ei[e], b = ej[e];
      idx[fill[a]] = b; wgt[fill[a]] = ew[e]; ++fill[a];
      idx[fill[b]] = a; wgt[fill[b]] = ew[e]; ++fill[b];
    }
  }
};

// ---------------------------------------------------------------------------
// Lazily-updated binary min-heap keyed on weighted degree.
//
// Peeling needs the minimum-degree node n times while degrees keep falling.
// Rather than rescanning every active node -- O(n) per step, O(n^2) overall --
// each degree change pushes a fresh entry and stamps the node with a new
// version. Entries whose version is stale are discarded when popped, so the
// total work is O(E log E) rather than O(n^2).
// ---------------------------------------------------------------------------
struct HeapItem {
  double key;
  int node;
  int version;
};

struct HeapCmp {
  // std::push_heap builds a max-heap, so invert to get the minimum on top.
  // Ties break toward the lower node index, which keeps results deterministic.
  bool operator()(const HeapItem& a, const HeapItem& b) const {
    if (a.key != b.key) return a.key > b.key;
    return a.node > b.node;
  }
};

// ---------------------------------------------------------------------------
// Scratch space reused across peeling rounds and across permutations.
// ---------------------------------------------------------------------------
struct Workspace {
  std::vector<double>   deg;
  std::vector<int>      ver;
  std::vector<char>     removed;
  std::vector<int>      loc;      // global node -> local slot, or -1
  std::vector<int>      order;    // removal order, local slots
  std::vector<int>      active;   // global ids of the nodes still in play
  std::vector<int>      buf;
  std::vector<char>     mark;
  std::vector<HeapItem> heap;

  // Edge-list buffers for rebuilding a permuted graph.
  std::vector<int>    ei, ej;
  std::vector<double> ew;

  void reserve(int N) {
    if (static_cast<int>(deg.size()) < N) {
      deg.resize(N); ver.resize(N); removed.resize(N); loc.resize(N);
      order.resize(N); active.resize(N); buf.resize(N); mark.resize(N);
    }
  }
};

// ---------------------------------------------------------------------------
// One greedy peeling pass over the subgraph induced on `ws.active[0..n)`.
//
// Repeatedly removes the minimum-weighted-degree node, scoring the set that
// remains after each removal, and returns the step whose remaining set scores
// highest. Returns -1 when the subgraph carries no weight at all.
//
// `ws.loc` must already map each active global id to its local slot and every
// other id to -1.
// ---------------------------------------------------------------------------
inline int peel_once(const Graph& g, int n, double lambda, int obj,
                     int min_size, Workspace& ws) {
  if (n < min_size + 1) return -1;

  const int* av = ws.active.data();
  double* deg = ws.deg.data();
  int* ver = ws.ver.data();
  char* removed = ws.removed.data();
  const int* loc = ws.loc.data();

  double total = 0.0;
  ws.heap.clear();
  ws.heap.reserve(static_cast<size_t>(n) * 2);

  for (int a = 0; a < n; ++a) {
    const int u = av[a];
    double s = 0.0;
    for (int p = g.ptr[u]; p < g.ptr[u + 1]; ++p) {
      if (loc[g.idx[p]] >= 0) s += g.wgt[p];
    }
    deg[a] = s;
    ver[a] = 0;
    removed[a] = 0;
    total += s;
    ws.heap.push_back(HeapItem{s, a, 0});
  }
  total *= 0.5;
  if (!(total > 0.0)) return -1;

  std::make_heap(ws.heap.begin(), ws.heap.end(), HeapCmp());

  int    best_step  = -1;
  double best_score = NEG_INF;

  for (int step = 0; step < n; ++step) {
    int jmin = -1;
    while (!ws.heap.empty()) {
      const HeapItem it = ws.heap.front();
      std::pop_heap(ws.heap.begin(), ws.heap.end(), HeapCmp());
      ws.heap.pop_back();
      if (!removed[it.node] && it.version == ver[it.node]) { jmin = it.node; break; }
    }
    if (jmin < 0) return -1;  // cannot happen, but never loop on a bad heap

    ws.order[step] = jmin;
    total -= deg[jmin];
    if (total < 0.0) total = 0.0;  // guard against floating-point drift
    removed[jmin] = 1;

    const int u = av[jmin];
    for (int p = g.ptr[u]; p < g.ptr[u + 1]; ++p) {
      const int b = loc[g.idx[p]];
      if (b < 0 || removed[b]) continue;
      deg[b] -= g.wgt[p];
      if (deg[b] < 0.0) deg[b] = 0.0;
      ws.heap.push_back(HeapItem{deg[b], b, ++ver[b]});
      std::push_heap(ws.heap.begin(), ws.heap.end(), HeapCmp());
    }

    const int n_rem = n - step - 1;
    if (n_rem >= min_size) {
      const double sc = obj_score(total, n_rem, lambda, obj);
      if (sc > best_score) { best_score = sc; best_step = step; }
    }
  }
  return best_step;
}

// ---------------------------------------------------------------------------
// Partition of the node set: a sequence of extracted subnetworks followed by
// one background block holding everything left over.
// ---------------------------------------------------------------------------
struct Partition {
  std::vector<int>    order;    // 0-based node ids, blocks concatenated
  std::vector<int>    sizes;    // block sizes; the last entry is the background
  std::vector<double> n_edge;   // supra-threshold edge count per block
  int n_clusters = 0;           // genuine subnetworks, excluding the background
};

// Count edges with both endpoints inside `block`.
inline double count_internal(const Graph& g, const int* block, int len,
                             Workspace& ws) {
  char* mark = ws.mark.data();
  for (int a = 0; a < len; ++a) mark[block[a]] = 1;
  double c = 0.0;
  for (int a = 0; a < len; ++a) {
    const int u = block[a];
    for (int p = g.ptr[u]; p < g.ptr[u + 1]; ++p) if (mark[g.idx[p]]) c += 1.0;
  }
  for (int a = 0; a < len; ++a) mark[block[a]] = 0;
  return 0.5 * c;
}

// ---------------------------------------------------------------------------
// Iteratively peel subnetworks out of a screened graph.
//
// Each round extracts the densest core; the nodes peeled off in that round
// become the input to the next round. Extraction stops as soon as the
// remainder carries no supra-threshold weight, which is typically after two or
// three rounds, avoiding a long tail of meaningless singleton blocks.
// ---------------------------------------------------------------------------
inline Partition extract_partition(const Graph& g, double lambda, int obj,
                                   int min_size, int max_clusters,
                                   Workspace& ws) {
  const int N = g.n;
  Partition out;
  out.order.reserve(N);

  ws.reserve(N);
  std::fill(ws.mark.begin(), ws.mark.begin() + N, 0);
  for (int i = 0; i < N; ++i) { ws.active[i] = i; ws.loc[i] = i; }
  int n_active = N;

  while (n_active >= min_size + 1 && out.n_clusters < max_clusters) {
    const int n = n_active;
    const int best = peel_once(g, n, lambda, obj, min_size, ws);
    if (best < 0) break;  // remainder carries no weight

    // The cluster is what survives `best + 1` removals, listed densest-first,
    // i.e. in reverse order of removal.
    const size_t off = out.order.size();
    for (int t = n - 1; t > best; --t) out.order.push_back(ws.active[ws.order[t]]);
    const int csize = n - best - 1;
    out.sizes.push_back(csize);
    out.n_edge.push_back(count_internal(g, out.order.data() + off, csize, ws));
    out.n_clusters += 1;

    // The peeled-off nodes feed the next round. Sorting keeps tie-breaking
    // deterministic and memory access ordered.
    for (int t = 0; t <= best; ++t) ws.buf[t] = ws.active[ws.order[t]];
    std::sort(ws.buf.begin(), ws.buf.begin() + best + 1);
    for (int t = 0; t < n; ++t) ws.loc[ws.active[t]] = -1;
    n_active = best + 1;
    for (int t = 0; t < n_active; ++t) {
      ws.active[t] = ws.buf[t];
      ws.loc[ws.buf[t]] = t;
    }
  }

  for (int t = 0; t < n_active; ++t) out.order.push_back(ws.active[t]);
  out.sizes.push_back(n_active);
  out.n_edge.push_back(
    count_internal(g, out.order.data() + out.order.size() - n_active,
                   n_active, ws));

  for (int t = 0; t < n_active; ++t) ws.loc[ws.active[t]] = -1;
  return out;
}

// ---------------------------------------------------------------------------
// Counter-based random number generation.
//
// Each permutation derives its own stream from (base_seed, permutation index)
// via splitmix64, so results do not depend on how permutations are spread
// across threads, and a run split into chunks reproduces a single run exactly.
// ---------------------------------------------------------------------------
struct Rng {
  uint64_t s[4];

  static inline uint64_t splitmix64(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

  explicit Rng(uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) s[i] = splitmix64(x);
  }

  static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

  inline uint64_t next() {  // xoshiro256++
    const uint64_t r = rotl(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;    s[3] = rotl(s[3], 45);
    return r;
  }

  // Lemire's nearly-divisionless bounded generator: uniform on [0, bound).
  inline uint64_t bounded(uint64_t bound) {
    __uint128_t m = static_cast<__uint128_t>(next()) * static_cast<__uint128_t>(bound);
    uint64_t l = static_cast<uint64_t>(m);
    if (l < bound) {
      const uint64_t t = (-bound) % bound;
      while (l < t) {
        m = static_cast<__uint128_t>(next()) * static_cast<__uint128_t>(bound);
        l = static_cast<uint64_t>(m);
      }
    }
    return static_cast<uint64_t>(m >> 64);
  }
};

template <typename T>
inline void shuffle_in_place(T* v, size_t n, Rng& rng) {
  for (size_t i = n - 1; i > 0; --i) {
    const size_t j = static_cast<size_t>(rng.bounded(static_cast<uint64_t>(i) + 1));
    std::swap(v[i], v[j]);
  }
}

// Open-addressing hash set of uint64, used by Floyd's sampling algorithm.
struct U64Set {
  std::vector<uint64_t> slot;
  std::vector<char>     used;
  uint64_t mask = 0;

  void reset(size_t capacity) {
    size_t cap = 8;
    while (cap < capacity * 2) cap <<= 1;
    slot.assign(cap, 0);
    used.assign(cap, 0);
    mask = cap - 1;
  }
  // Returns true if `x` was newly inserted.
  bool insert(uint64_t x) {
    uint64_t h = x * 0x9E3779B97F4A7C15ULL;
    size_t i = static_cast<size_t>((h >> 32) & mask);
    while (used[i]) {
      if (slot[i] == x) return false;
      i = (i + 1) & mask;
    }
    used[i] = 1; slot[i] = x;
    return true;
  }
};

// ---------------------------------------------------------------------------
// Draw the positions of the E surviving edges under a random permutation of
// the M-length edge-weight vector.
//
// Permuting the screened weight vector and keeping its non-zeros is the same,
// in distribution, as choosing E of the M positions uniformly at random and
// dealing the E surviving weights into them in random order. Doing it that way
// costs O(E) instead of O(M): at 400 nodes screened at the 99th percentile
// that is 798 operations rather than 79800.
// ---------------------------------------------------------------------------
inline void sample_positions(uint64_t M, uint64_t E, Rng& rng, U64Set& set,
                             std::vector<uint64_t>& out) {
  out.clear();
  if (E == 0) return;
  if (E * 4 > M) {
    // Dense case: a straight partial shuffle beats Floyd's algorithm.
    std::vector<uint64_t> all(M);
    for (uint64_t t = 0; t < M; ++t) all[t] = t;
    for (uint64_t i = 0; i < E; ++i) {
      const uint64_t j = i + rng.bounded(M - i);
      std::swap(all[i], all[j]);
      out.push_back(all[i]);
    }
    return;
  }
  set.reset(static_cast<size_t>(E));
  out.reserve(static_cast<size_t>(E));
  for (uint64_t j = M - E; j < M; ++j) {
    const uint64_t t = rng.bounded(j + 1);
    if (set.insert(t)) out.push_back(t);
    else { set.insert(j); out.push_back(j); }
  }
}

// Decode a position in the column-major strict lower triangle back to its
// (row, column) node pair. This is the inverse of R's `W[lower.tri(W)]` and of
// the column-major vectorization used throughout the package.
inline void decode_position(uint64_t t, int N, int& i, int& j) {
  const double b = 2.0 * N - 1.0;
  int jj = static_cast<int>((b - std::sqrt(b * b - 8.0 * static_cast<double>(t))) / 2.0);
  if (jj < 0) jj = 0;
  if (jj > N - 2) jj = N - 2;
  auto offset = [&](int c) -> uint64_t {
    return static_cast<uint64_t>(c) * (2ULL * N - 1ULL - c) / 2ULL;
  };
  while (jj > 0 && offset(jj) > t) --jj;
  while (jj < N - 2 && offset(jj + 1) <= t) ++jj;
  j = jj;
  i = static_cast<int>(t - offset(jj)) + jj + 1;
}

}  // namespace subnetr

#endif  // SUBNETR_CORE_H
