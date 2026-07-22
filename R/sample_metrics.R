#' @importFrom stats median quantile rnorm

## Shared draw-based helpers, reused by more than one sample_metric.

#' @noRd
.pairwise_dist_matrix <- function(A, B) {
  aa <- rowSums(A^2)
  bb <- rowSums(B^2)
  ab <- A %*% t(B)
  d2 <- outer(aa, bb, "+") - 2 * ab
  sqrt(pmax(d2, 0))
}

#' @noRd
.mean_offdiag <- function(D) {
  n <- nrow(D)
  if (n < 2) return(0)
  (sum(D) - sum(diag(D))) / (n * (n - 1))
}

#' @noRd
.wasserstein1d <- function(x, y) {
  n <- max(length(x), length(y))
  qx <- stats::quantile(x, probs = seq(0, 1, length.out = n + 1), type = 7, names = FALSE)
  qy <- stats::quantile(y, probs = seq(0, 1, length.out = n + 1), type = 7, names = FALSE)
  mean(abs(qx - qy))
}

#' @title Draw-Based Group-Difference Metrics
#'
#' @description
#' A \code{sample_metric} is a lightweight S3 object describing which
#' draw-based (rather than closed-form Gaussian) distance to compute between
#' two groups' posterior samples. Unlike \code{distance_metric}
#' objects, these don't need \code{mu}/\code{Sigma} at all -- they work
#' directly from posterior draws (e.g. from \code{\link{sample_posterior}}),
#' which is also why they're driven by \code{\link{compute_group_diff_samples}}
#' rather than \code{\link{compute_group_diff}}.
#'
#' @param metric A \code{sample_metric} object.
#' @param draws1,draws2 Numeric matrices of posterior draws for the two
#'   groups (rows = draws, columns = shared IDs, same column order).
#'   Row counts may differ between the two groups.
#' @param ... Unused; included for extensibility.
#'
#' @return \code{evaluate_sample_metric()} returns a single numeric value.
#' @export
evaluate_sample_metric <- function(metric, draws1, draws2, ...) {
  UseMethod("evaluate_sample_metric")
}

#' @title Energy Distance Metric
#' @description The (U-statistic) energy distance between two groups'
#'   posterior draws: \eqn{2\,\overline{\|X-Y\|} - \overline{\|X-X'\|} - \overline{\|Y-Y'\|}},
#'   with within-group terms excluding self-pairs. Unlike
#'   \code{\link{ovl_metric}}, it does not require the joint density to
#'   overlap in every dimension simultaneously, so it does not collapse
#'   toward a degenerate value as the number of shared IDs grows.
#' @return A \code{sample_metric} object.
#' @export
energy_metric <- function() {
  structure(list(), class = c("energy_metric", "sample_metric"))
}

#' @export
evaluate_sample_metric.energy_metric <- function(metric, draws1, draws2, ...) {
  cross   <- .pairwise_dist_matrix(draws1, draws2)
  within1 <- .pairwise_dist_matrix(draws1, draws1)
  within2 <- .pairwise_dist_matrix(draws2, draws2)
  2 * mean(cross) - .mean_offdiag(within1) - .mean_offdiag(within2)
}

#' @title Maximum Mean Discrepancy (MMD) Metric
#' @description Squared MMD between two groups' posterior draws under a
#'   Gaussian (RBF) kernel. When \code{bandwidth} is \code{NULL} (default),
#'   it is set via the median heuristic (the median pairwise distance in the
#'   pooled sample of both groups).
#' @param bandwidth The Gaussian kernel bandwidth, or \code{NULL} for the
#'   median-heuristic default.
#' @return A \code{sample_metric} object.
#' @export
mmd_metric <- function(bandwidth = NULL) {
  structure(list(bandwidth = bandwidth), class = c("mmd_metric", "sample_metric"))
}

#' @export
evaluate_sample_metric.mmd_metric <- function(metric, draws1, draws2, ...) {
  cross   <- .pairwise_dist_matrix(draws1, draws2)
  within1 <- .pairwise_dist_matrix(draws1, draws1)
  within2 <- .pairwise_dist_matrix(draws2, draws2)

  bw <- metric$bandwidth
  if (is.null(bw)) {
    pooled_dist <- .pairwise_dist_matrix(rbind(draws1, draws2), rbind(draws1, draws2))
    bw <- stats::median(pooled_dist[upper.tri(pooled_dist)])
    if (!is.finite(bw) || bw <= 0) bw <- 1
  }

  gaussian_kernel <- function(D) exp(-D^2 / (2 * bw^2))
  .mean_offdiag(gaussian_kernel(within1)) + .mean_offdiag(gaussian_kernel(within2)) -
    2 * mean(gaussian_kernel(cross))
}

#' @title Sliced Wasserstein Distance Metric
#' @description Average of the 1D Wasserstein-1 distance (via sorted order
#'   statistics/interpolated quantiles) between two groups' posterior draws,
#'   projected onto \code{n_projections} random unit directions in ID-space.
#'   Explicitly designed to avoid the curse of dimensionality that afflicts
#'   the general (non-sliced) multivariate optimal-transport distance.
#' @param n_projections Number of random projection directions to average
#'   over. Defaults to \code{50}.
#' @return A \code{sample_metric} object.
#' @export
sliced_wasserstein_metric <- function(n_projections = 50) {
  structure(list(n_projections = n_projections), class = c("sliced_wasserstein_metric", "sample_metric"))
}

#' @export
evaluate_sample_metric.sliced_wasserstein_metric <- function(metric, draws1, draws2, ...) {
  d <- ncol(draws1)
  vals <- vapply(seq_len(metric$n_projections), function(i) {
    v <- stats::rnorm(d)
    v <- v / sqrt(sum(v^2))
    .wasserstein1d(draws1 %*% v, draws2 %*% v)
  }, numeric(1))
  mean(vals)
}

#' @title Compute a Pairwise Draw-Based Group-Difference Metric
#'
#' @description
#' The draw-based counterpart to \code{\link{compute_group_diff}}: for every
#' pair of groups present in \code{sample_distrib}, aligns their posterior
#' draws by shared ID (via the same melt-order invariant documented for
#' \code{\link{compute_multi_diff}}) and evaluates \code{metric} via
#' \code{\link{evaluate_sample_metric}}. All \code{sample_metric}s are
#' symmetric, and the diagonal is fixed at \code{0} by convention (a group's
#' distance to itself), rather than computed from finite-sample draws, which
#' would be a slightly-biased nonzero estimate rather than exactly \code{0}.
#'
#' @param sample_distrib A data frame, typically from
#'   \code{\link{sample_posterior}}, with columns \code{ID}, \code{Group}, \code{Sample}.
#' @param metric A \code{sample_metric} object, e.g. \code{\link{energy_metric}()},
#'   \code{\link{mmd_metric}()}, \code{\link{sliced_wasserstein_metric}()}.
#'
#' @return A square, symmetric matrix of metric values, one row/column per group, diagonal 0.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 300)
#' compute_group_diff_samples(samples, energy_metric())
compute_group_diff_samples <- function(sample_distrib, metric) {
  required_cols <- c("ID", "Group", "Sample")
  if (!all(required_cols %in% names(sample_distrib))) {
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }
  if (!inherits(metric, "sample_metric")) {
    stop("'metric' must be a sample_metric object (e.g. energy_metric(), mmd_metric(), sliced_wasserstein_metric()).")
  }

  groups     <- sort(unique(sample_distrib$Group))
  num_groups <- length(groups)
  diff_matrix <- matrix(0, num_groups, num_groups)
  rownames(diff_matrix) <- colnames(diff_matrix) <- groups
  if (num_groups < 2) return(diff_matrix)

  pairs <- utils::combn(groups, 2, simplify = FALSE)
  for (p in pairs) {
    g1 <- p[1]
    g2 <- p[2]
    ids1 <- unique(sample_distrib$ID[sample_distrib$Group == g1])
    ids2 <- unique(sample_distrib$ID[sample_distrib$Group == g2])
    if (!setequal(ids1, ids2)) {
      stop(paste0(
        "Groups '", g1, "' and '", g2, "' do not share the same set of IDs: ",
        "only in '", g1, "': [", paste(setdiff(ids1, ids2), collapse = ", "), "]; ",
        "only in '", g2, "': [", paste(setdiff(ids2, ids1), collapse = ", "), "]."
      ))
    }
    shared_ids <- sort(ids1)
    draws1 <- extract_draw_matrix(sample_distrib, g1, shared_ids)
    draws2 <- extract_draw_matrix(sample_distrib, g2, shared_ids)

    val <- evaluate_sample_metric(metric, draws1, draws2)
    diff_matrix[g1, g2] <- val
    diff_matrix[g2, g1] <- val
  }
  diff_matrix
}
