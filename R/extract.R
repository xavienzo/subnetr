#' Extract candidate subnetworks by greedy peeling
#'
#' Partitions the nodes of a weighted graph into a sequence of increasingly
#' less dense blocks. Each round repeatedly deletes the minimum-degree node,
#' scores the set that remains after every deletion, and keeps the
#' best-scoring set as one subnetwork; the deleted nodes are then fed back in
#' for the next round. Everything left once the remainder carries no
#' supra-threshold weight is returned as a single background block.
#'
#' This is the extraction step only. It reports how dense each block is but
#' attaches no significance to it -- use [subnet()] for the full analysis with
#' a calibrated permutation test.
#'
#' @section Choosing an objective:
#' A candidate node set \eqn{S} of `n` nodes carrying total suprathreshold edge
#' weight `w` is scored by
#'
#' \describe{
#'   \item{`"generalized"`}{\eqn{w / n^{2\lambda}}, the default. This is the
#'     adaptive density function \eqn{f(S; \lambda_W) = |W(S)| / |S|^{\lambda_W}}
#'     of Wu et al. (2022), and equivalently the \eqn{\ell_0} graph norm
#'     shrinkage criterion \eqn{\log\|U\|_1 - \lambda_0 \log\|U\|_0} of Chen et
#'     al. (2023), which rewards edge weight within a subnetwork while
#'     penalizing its size.}
#'   \item{`"density"`}{\eqn{w / \binom{n}{2}}, pure edge density, ignoring
#'     `lambda`. With no size floor this degenerates to the single densest pair
#'     of nodes, so use it only with a substantial `min_size`.}
#' }
#'
#' @section Parameterization of lambda:
#' The `lambda` used here equals \eqn{\lambda_0} of Chen et al. (2023) and half
#' of \eqn{\lambda_W} of Wu et al. (2022), whose exponent runs over \[1, 2\].
#' Concretely, `lambda = 0.5` (\eqn{\lambda_W = 1}) gives the degree density
#' \eqn{f_1}, the objective of Charikar (2000), and `lambda = 1`
#' (\eqn{\lambda_W = 2}) gives the area density \eqn{f_2}. Larger values favour
#' smaller and denser subnetworks; `lambda = 0` places every node in one
#' subnetwork. Values between 0.5 and 0.7 are the usual working range, and
#' [tune_subnet()] selects a value from the data when none is supplied.
#'
#' @param W Symmetric numeric matrix of edge weights, typically
#'   \eqn{-\log_{10}} p-values from [edge_stats()]. The diagonal is ignored.
#' @param threshold Screening threshold on the edge weights. Edges below it are
#'   set to zero before peeling.
#' @param lambda Size-penalty exponent in \[0, 1\]; see the section on choosing
#'   an objective. Ignored when `objective = "density"`.
#' @param objective One of `"generalized"` (default) or `"density"`.
#' @param min_size Smallest number of nodes a block may contain. Blocks of two
#'   nodes are single edges and carry almost no evidence, so the default is 3.
#' @param max_clusters Cap on the number of extraction rounds. Peeling returns
#'   the densest block first and each subsequent round works on sparser
#'   leftovers, so the blocks that matter appear early; the cap simply stops
#'   the algorithm from grinding through a long tail of noise blocks. Raising
#'   it costs time and, in practice, does not change which blocks reach
#'   significance.
#'
#' @return An object of class `subnet_partition`, a list with components
#'   \item{blocks}{List of integer vectors of node indices, densest block
#'     first. The final element is the background block and may be empty.}
#'   \item{sizes}{Number of nodes in each block.}
#'   \item{n_edge}{Number of supra-threshold edges inside each block.}
#'   \item{density}{Proportion of possible within-block edges that are
#'     supra-threshold.}
#'   \item{n_clusters}{Number of extracted subnetworks, excluding the
#'     background block.}
#'   \item{order}{The node ordering that produces the block-diagonal
#'     arrangement, i.e. `blocks` concatenated.}
#'   \item{threshold, lambda, objective, min_size}{The settings used.}
#'   \item{pi0}{Proportion of all edges that are supra-threshold.}
#'
#' @references
#' Chen S, Zhang Y, Wu Q, Bi C, Kochunov P, Hong LE (2023). Identifying
#' covariate-related subnetworks for whole-brain connectome analysis.
#' *Biostatistics*, kxad007. \doi{10.1093/biostatistics/kxad007}
#'
#' Wu Q, Huang X, Culbreth AJ, Waltz JA, Hong LE, Chen S (2022). Extracting
#' brain disease-related connectome subgraphs by adaptive dense subgraph
#' discovery. *Biometrics* 78(4), 1566-1578. \doi{10.1111/biom.13537}
#'
#' @seealso [subnet()] for extraction plus inference, [tune_subnet()] for
#'   choosing `lambda` and `threshold`.
#'
#' @examples
#' sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.15,
#'                    seed = 1)
#' part <- subnet_extract(sim$W, threshold = quantile(vech(sim$W), 0.95))
#' part
#'
#' @export
subnet_extract <- function(W, threshold, lambda = 0.6,
                           objective = c("generalized", "density"),
                           min_size = 3L, max_clusters = 25L) {
  W <- check_adjacency(W)
  obj <- match_objective(objective)
  objective <- names(OBJ_CODES)[OBJ_CODES == obj]

  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("`threshold` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0 || lambda > 1) {
    stop("`lambda` must be a single number in [0, 1].", call. = FALSE)
  }
  min_size <- max(2L, as.integer(min_size))
  max_clusters <- max(1L, as.integer(max_clusters))

  res <- subnet_extract_cpp(W, threshold, lambda, obj, min_size, max_clusters)
  new_partition(res, W, threshold, lambda, objective, min_size)
}

new_partition <- function(res, W, threshold, lambda, objective, min_size) {
  sizes <- as.integer(res$sizes)
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  blocks <- Map(function(a, b) if (b >= a) res$order[a:b] else integer(0),
                starts, ends)

  m <- sizes * (sizes - 1) / 2
  density <- ifelse(m > 0, res$n_edge / m, NA_real_)

  wv <- vech(W)
  pi0 <- mean(wv >= threshold)

  structure(
    list(blocks = blocks,
         sizes = sizes,
         n_edge = as.numeric(res$n_edge),
         density = density,
         n_clusters = as.integer(res$n_clusters),
         order = as.integer(res$order),
         threshold = threshold,
         lambda = lambda,
         objective = objective,
         min_size = min_size,
         pi0 = pi0,
         n_nodes = nrow(W)),
    class = "subnet_partition")
}

## Block-wise log p-values for a partition, background block included.
partition_logp <- function(part) {
  m <- part$sizes * (part$sizes - 1) / 2
  block_logp(part$n_edge, m, part$pi0)
}
