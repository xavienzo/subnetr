#include <Rcpp.h>
#include "subnet_core.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace subnetr;

// [[Rcpp::export]]
bool subnetr_has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

// ---------------------------------------------------------------------------
// Extract subnetworks from an observed weighted adjacency matrix.
// `W` must be symmetric; the diagonal is ignored. `r` screens the weights.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List subnet_extract_cpp(NumericMatrix W, double r, double lambda, int obj,
                        int min_size, int max_clusters) {
  const int N = W.nrow();
  if (W.ncol() != N) stop("`W` must be square.");

  std::vector<int> ei, ej;
  std::vector<double> ew;
  for (int j = 0; j < N; ++j) {
    for (int i = j + 1; i < N; ++i) {
      const double x = W(i, j);
      if (x >= r && x != 0.0) { ei.push_back(i); ej.push_back(j); ew.push_back(x); }
    }
  }

  Graph g;
  g.build(N, ei, ej, ew);

  Workspace ws;
  ws.reserve(N);
  Partition p = extract_partition(g, lambda, obj, min_size, max_clusters, ws);

  IntegerVector order(p.order.size());
  for (size_t i = 0; i < p.order.size(); ++i) order[i] = p.order[i] + 1;  // 1-based

  return List::create(
    _["order"]      = order,
    _["sizes"]      = wrap(p.sizes),
    _["n_edge"]     = wrap(p.n_edge),
    _["n_clusters"] = p.n_clusters);
}

// ---------------------------------------------------------------------------
// Permutation null.
//
// For each permutation the screened edge weights are redistributed at random
// over the graph and the same extraction is rerun. We return per-block
// supra-threshold edge counts and block sizes rather than p-values, so that no
// R math library call happens inside a parallel region.
//
// Screening commutes with permutation, so the weight vector is screened once
// up front. Only the surviving weights are then placed, which makes the cost
// of a permutation proportional to the number of edges that survived rather
// than to the number of node pairs.
//
// `perm_offset` lets a run be split into chunks that reproduce a single
// contiguous run exactly, because each permutation's RNG stream is derived
// from (seed, permutation index) alone.
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List subnet_perm_cpp(NumericVector wvec, int N, double r, double lambda,
                     int obj, int min_size, int max_clusters,
                     int n_perm, double seed, int perm_offset, int n_threads) {
  const uint64_t M = static_cast<uint64_t>(wvec.size());
  if (M != static_cast<uint64_t>(N) * (N - 1) / 2)
    stop("`wvec` length does not match `N`.");

  std::vector<double> kept;
  kept.reserve(static_cast<size_t>(M) / 20 + 8);
  for (uint64_t t = 0; t < M; ++t) {
    const double x = wvec[t];
    if (x >= r && x != 0.0) kept.push_back(x);
  }
  const uint64_t E = kept.size();

  const int ncol = max_clusters + 1;  // blocks, plus the background block
  NumericMatrix k_mat(n_perm, ncol);  // supra-threshold edges per block
  NumericMatrix m_mat(n_perm, ncol);  // possible edges per block
  IntegerVector nc_vec(n_perm);       // extracted subnetworks, excl. background

  const uint64_t base = static_cast<uint64_t>(seed);

#ifdef _OPENMP
  if (n_threads > 1) omp_set_num_threads(n_threads);
#else
  (void)n_threads;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    Workspace ws;
    ws.reserve(N);
    Graph g;
    U64Set set;
    std::vector<uint64_t> pos;
    std::vector<double> vals(kept.size());

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 8)
#endif
    for (int b = 0; b < n_perm; ++b) {
      // Derive an independent stream from (base seed, absolute permutation id)
      // so the result is invariant to how permutations are chunked or threaded.
      uint64_t s = base ^ (0x9E3779B97F4A7C15ULL *
                           static_cast<uint64_t>(perm_offset + b + 1));
      Rng rng(Rng::splitmix64(s));

      sample_positions(M, E, rng, set, pos);
      std::copy(kept.begin(), kept.end(), vals.begin());
      if (E > 1) shuffle_in_place(vals.data(), vals.size(), rng);

      ws.ei.resize(E); ws.ej.resize(E); ws.ew.resize(E);
      for (uint64_t e = 0; e < E; ++e) {
        int i, j;
        decode_position(pos[e], N, i, j);
        ws.ei[e] = i; ws.ej[e] = j; ws.ew[e] = vals[e];
      }
      g.build(N, ws.ei, ws.ej, ws.ew);

      Partition p = extract_partition(g, lambda, obj, min_size, max_clusters, ws);

      const int nb = std::min(static_cast<int>(p.sizes.size()), ncol);
      for (int c = 0; c < nb; ++c) {
        const double sz = p.sizes[c];
        k_mat(b, c) = p.n_edge[c];
        m_mat(b, c) = 0.5 * sz * (sz - 1.0);
      }
      nc_vec[b] = p.n_clusters;
    }
  }

  return List::create(_["k"] = k_mat, _["m"] = m_mat, _["nc"] = nc_vec);
}

// ---------------------------------------------------------------------------
// Exposed for testing: decode positions in the vectorized lower triangle back
// to node pairs, so the mapping can be checked against R's lower.tri().
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerMatrix decode_position_cpp(NumericVector t, int N) {
  IntegerMatrix out(t.size(), 2);
  for (int e = 0; e < t.size(); ++e) {
    int i, j;
    decode_position(static_cast<uint64_t>(t[e]), N, i, j);
    out(e, 0) = i + 1;
    out(e, 1) = j + 1;
  }
  return out;
}
