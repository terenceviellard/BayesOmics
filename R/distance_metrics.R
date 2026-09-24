#' @importFrom stats pnorm pchisq

## Shared (mu1, mu2, Sigma1, Sigma2) helpers used by more than one metric's
## evaluate_metric() method, so they aren't duplicated inside each method.

#' @noRd
.kl_divergence <- function(mu1, mu2, Sigma1, Sigma2) {
  d <- length(mu1)
  inv2 <- chol_inv_jitter(Sigma2, pen_diag = 1e-6)
  delta <- mu2 - mu1
  logdet1 <- as.numeric(determinant(Sigma1, logarithm = TRUE)$modulus)
  logdet2 <- as.numeric(determinant(Sigma2, logarithm = TRUE)$modulus)
  0.5 * (sum(diag(inv2 %*% Sigma1)) + as.numeric(t(delta) %*% inv2 %*% delta) - d + (logdet2 - logdet1))
}

#' @noRd
.matrix_sqrt_sym <- function(M) {
  eig <- eigen(M, symmetric = TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values)) %*% t(eig$vectors)
}

## Shared helpers behind ovl_metric()'s Monte Carlo/KDE fallback (used when
## Sigma1 is not a scalar multiple of Sigma2, so the exact closed-form ratio
## does not apply -- see evaluate_metric.ovl_metric() below).

#' @noRd
#'
#' @details Silverman's per-dimension bandwidth rule of thumb for a diagonal
#'   (product) Gaussian KDE fit on `draws` (n x d): `h_k = sd(draws[, k]) *
#'   n^(-1/(d+4))`. Floored at a small positive value so a (near-)degenerate
#'   dimension (near-zero variance across draws) does not collapse that
#'   axis's kernel to a point mass and blow up the log-density below.
.kde_bandwidth <- function(draws) {
  n <- nrow(draws)
  d <- ncol(draws)
  h <- apply(draws, 2, stats::sd) * n^(-1 / (d + 4))
  pmax(h, 1e-8)
}

#' @noRd
#'
#' @details Log-density of a diagonal (product) Gaussian KDE fit on
#'   `ref_draws` (n_ref x d, per-dimension bandwidth `bw`), evaluated at
#'   every row of `query` (n_query x d). Computed in log-space -- rather
#'   than multiplying `d` per-dimension kernel values directly -- to avoid
#'   underflow when `d` is large, via a numerically stable log-sum-exp over
#'   reference points (the max-subtraction trick). O(n_query * n_ref * d).
.log_kde_eval <- function(query, ref_draws, bw) {
  n_ref <- nrow(ref_draws)
  d     <- ncol(ref_draws)
  bw_mat_q <- matrix(bw, nrow(query), d, byrow = TRUE)
  bw_mat_r <- matrix(bw, n_ref,        d, byrow = TRUE)
  scaled_query <- query / bw_mat_q
  scaled_ref   <- ref_draws / bw_mat_r
  # sq_dist[i, j] = sum_k ((query[i, k] - ref_draws[j, k]) / bw[k])^2
  sq_dist <- outer(rowSums(scaled_query^2), rowSums(scaled_ref^2), "+") -
    2 * scaled_query %*% t(scaled_ref)
  log_kernel <- -0.5 * d * log(2 * pi) - sum(log(bw)) - 0.5 * sq_dist
  m <- apply(log_kernel, 1, max)
  log_sum <- m + log(rowSums(exp(log_kernel - m)))
  log_sum - log(n_ref)
}

#' @noRd
#'
#' @details Monte Carlo estimate of \eqn{OVL = \int\min(f_1,f_2)\,dx} for two
#'   general (not necessarily proportional-covariance) Gaussians, using the
#'   two-sample identity \eqn{OVL = E_{f_1}[\min(1, f_2/f_1)] = E_{f_2}[\min(1,
#'   f_1/f_2)]} (both exact in the population limit; averaged here for lower
#'   variance). `n_mc` fresh draws are simulated directly from
#'   \eqn{N(\mu_1,\Sigma_1)} and \eqn{N(\mu_2,\Sigma_2)} (no posterior
#'   samples/`sample_posterior()` call needed -- `mu`/`Sigma` fully
#'   characterize each Gaussian already), then `f_1`/`f_2` are estimated by a
#'   diagonal Gaussian KDE on each simulated sample (see `.kde_bandwidth()`/
#'   `.log_kde_eval()`). Noisier (Monte Carlo + KDE bias) and increasingly
#'   unreliable as `d = length(mu1)` grows (the curse of dimensionality for
#'   KDE) -- this is only reached when the exact formula does not apply.
.ovl_mc_fallback <- function(mu1, mu2, Sigma1, Sigma2, n_mc) {
  draws1 <- matrix(mvtnorm::rmvnorm(n_mc, mean = mu1, sigma = Sigma1), nrow = n_mc)
  draws2 <- matrix(mvtnorm::rmvnorm(n_mc, mean = mu2, sigma = Sigma2), nrow = n_mc)

  bw1 <- .kde_bandwidth(draws1)
  bw2 <- .kde_bandwidth(draws2)

  log_f1_at_1 <- .log_kde_eval(draws1, draws1, bw1)
  log_f2_at_1 <- .log_kde_eval(draws1, draws2, bw2)
  log_f1_at_2 <- .log_kde_eval(draws2, draws1, bw1)
  log_f2_at_2 <- .log_kde_eval(draws2, draws2, bw2)

  term1 <- mean(exp(pmin(0, log_f2_at_1 - log_f1_at_1)))
  term2 <- mean(exp(pmin(0, log_f1_at_2 - log_f2_at_2)))
  0.5 * (term1 + term2)
}

#' @title Group-Difference Metrics for BayesOmics Posteriors
#'
#' @description
#' A \code{distance_metric} is a lightweight S3 object describing *which*
#' pairwise divergence/distance to compute between two groups' Gaussian
#' posteriors, without itself holding any data. Pass one to
#' \code{\link{compute_group_diff}} (or, for sample-based metrics,
#' \code{\link{compute_group_diff_samples}}) to get a full group x group
#' matrix of values.
#'
#' Three S3 generics drive the dispatch, each with a default method so a new
#' metric only needs to implement \code{evaluate_metric}:
#' \describe{
#'   \item{\code{evaluate_metric(metric, mu1, mu2, Sigma1, Sigma2, scale1, scale2, same_kernel, ...)}}{
#'     Computes the metric's value for one ordered pair of groups. \code{mu1}/\code{mu2}
#'     are the groups' posterior means (aligned by ID), \code{Sigma1}/\code{Sigma2}
#'     their posterior covariances (from the internal \code{get_sigmak()} helper), \code{scale1}/\code{scale2}
#'     the \code{n_obs + lambda_0} divisors from \code{\link{posterior_mean}},
#'     and \code{same_kernel} whether the two groups share the same \code{kernel_key}.}
#'   \item{\code{requires_shared_kernel(metric)}}{\code{TRUE} if the metric has
#'     no valid formula at all unless the two groups share the same
#'     \code{kernel_key} (no metric here needs this -- \code{\link{ovl_metric}}
#'     instead switches internally to a Monte Carlo/KDE fallback when they
#'     don't, see its own documentation). Defaults to \code{FALSE}.}
#'   \item{\code{is_symmetric_metric(metric)}}{\code{TRUE} if swapping the two
#'     groups never changes the value (\code{\link{kl_metric}} and
#'     \code{\link{mahalanobis_metric}} are not -- the latter is evaluated
#'     under \eqn{\Sigma_1^{-1}} specifically, so it differs from its swapped
#'     counterpart whenever \eqn{\Sigma_1 \neq \Sigma_2}). Defaults to \code{TRUE}.}
#' }
#'
#' @param metric A \code{distance_metric} object.
#' @param mu1,mu2 Named numeric vectors, the two groups' posterior means (same IDs).
#' @param Sigma1,Sigma2 The two groups' posterior covariance matrices (same ID order as \code{mu1}/\code{mu2}).
#' @param scale1,scale2 The two groups' posterior scale (\code{n_obs + lambda_0}); \code{NULL} if unavailable.
#' @param same_kernel Whether the two groups share the same \code{kernel_key}.
#' @param ... Passed through by decorator metrics; unused by concrete metrics.
#'
#' @return \code{evaluate_metric()} returns a single numeric value.
#'   \code{requires_shared_kernel()}/\code{is_symmetric_metric()} return a single logical.
#' @export
evaluate_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  UseMethod("evaluate_metric")
}

#' @rdname evaluate_metric
#' @export
requires_shared_kernel <- function(metric) {
  UseMethod("requires_shared_kernel")
}

#' @export
requires_shared_kernel.distance_metric <- function(metric) FALSE

#' @rdname evaluate_metric
#' @export
is_symmetric_metric <- function(metric) {
  UseMethod("is_symmetric_metric")
}

#' @export
is_symmetric_metric.distance_metric <- function(metric) TRUE

# =============================================================================
# Concrete analytic metrics
# =============================================================================

#' @title Mahalanobis Distance Metric
#' @description The Mahalanobis distance between two groups' posterior means
#'   under \eqn{\Sigma_1}: \eqn{D = \sqrt{(\mu_1-\mu_2)^\top\Sigma_1^{-1}(\mu_1-\mu_2)}}.
#'   This is the same quantity \code{\link{calculate_group_overlaps}} computes
#'   internally before converting it to an overlap coefficient. Asymmetric
#'   whenever \eqn{\Sigma_1 \neq \Sigma_2} -- swapping the two groups uses
#'   \eqn{\Sigma_2^{-1}} instead, generally a different value -- so
#'   \code{\link{compute_group_diff}} computes both directions separately
#'   rather than mirroring one across the diagonal (see
#'   \code{is_symmetric_metric.mahalanobis_metric} below).
#' @return A \code{distance_metric} object.
#' @export
mahalanobis_metric <- function() {
  structure(list(), class = c("mahalanobis_metric", "distance_metric"))
}

#' @export
is_symmetric_metric.mahalanobis_metric <- function(metric) FALSE

#' @export
evaluate_metric.mahalanobis_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  delta <- mu1 - mu2
  inv1 <- chol_inv_jitter(Sigma1, pen_diag = 1e-6)
  sqrt(as.numeric(t(delta) %*% inv1 %*% delta))
}

#' @title Gaussian Overlapping Coefficient (OVL) Metric
#' @description The overlapping coefficient between two groups' posteriors,
#'   as used by \code{\link{calculate_group_overlaps}}. Same call in every
#'   case: when the two groups share the same \code{kernel_key} (so
#'   \eqn{\Sigma_2 = c\Sigma_1} for a scalar \eqn{c =
#'   \mathrm{scale}_1/\mathrm{scale}_2}), the exact closed-form ratio is
#'   used. Otherwise (e.g. \code{pooled = FALSE}, or two groups fit with
#'   genuinely different kernels), there is no closed form -- the
#'   log-likelihood-ratio boundary between the two Gaussians then has a
#'   different weight per axis, not reducible to one (noncentral)
#'   chi-squared variable -- so a Monte Carlo/KDE estimate is used
#'   automatically instead (\code{.ovl_mc_fallback()}): \code{n_mc} draws
#'   simulated directly from each group's \eqn{N(\mu,\Sigma)}, then the
#'   two-sample identity \eqn{OVL = E_{f_1}[\min(1, f_2/f_1)] = E_{f_2}[\min(1,
#'   f_1/f_2)]} with \eqn{f_1}/\eqn{f_2} estimated by a diagonal Gaussian KDE.
#'   Noisier and increasingly unreliable as the number of shared IDs grows
#'   (the curse of dimensionality for KDE) -- only reached when the exact
#'   formula does not apply, so it never affects the pooled case.
#' @param n_mc Number of Monte Carlo draws used by the KDE fallback (see
#'   above); irrelevant (unused) whenever the two groups share a
#'   \code{kernel_key}, since the exact formula applies then. Defaults to
#'   \code{2000}.
#' @return A \code{distance_metric} object.
#' @export
ovl_metric <- function(n_mc = 2000) {
  structure(list(n_mc = n_mc), class = c("ovl_metric", "distance_metric"))
}

#' @export
evaluate_metric.ovl_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  if (!same_kernel) {
    return(.ovl_mc_fallback(mu1, mu2, Sigma1, Sigma2, metric$n_mc))
  }
  if (is.null(scale1) || is.null(scale2)) {
    stop("ovl_metric() requires 'scale1' and 'scale2' (the posterior scale of each group) for the closed-form (same-kernel) case.")
  }
  d <- length(mu1)
  delta <- mu1 - mu2
  c_ratio <- scale1 / scale2

  scale_rel_tol <- 1e-6
  near_equal_scale <- isTRUE(scale1 == scale2) || abs(c_ratio - 1) < scale_rel_tol

  if (near_equal_scale) {
    if (!isTRUE(scale1 == scale2)) {
      warning(sprintf(
        paste0(
          "ovl_metric: nearly identical scales (scale1 = %.6g, scale2 = %.6g, ratio = %.10f); the ",
          "c != 1 formula is numerically unstable this close to c = 1, so they are treated as equal ",
          "using the larger posterior covariance (smaller scale)."
        ),
        scale1, scale2, c_ratio
      ))
    }
    Sigma_eq <- if (scale1 <= scale2) Sigma1 else Sigma2
    inv_eq <- chol_inv_jitter(Sigma_eq, pen_diag = 1e-6)
    D2 <- as.numeric(t(delta) %*% inv_eq %*% delta)
    return(2 * stats::pnorm(-sqrt(D2) / 2))
  }

  inv1 <- chol_inv_jitter(Sigma1, pen_diag = 1e-6)
  D2 <- as.numeric(t(delta) %*% inv1 %*% delta)
  lambda1 <- D2 / (1 - c_ratio)^2
  lambda2 <- c_ratio * D2 / (1 - c_ratio)^2
  t_val <- (D2 - d * (1 - c_ratio) * log(c_ratio)) / (1 - c_ratio)^2
  if (c_ratio < 1) {
    stats::pchisq(c_ratio * t_val, df = d, ncp = lambda1) + 1 - stats::pchisq(t_val, df = d, ncp = lambda2)
  } else {
    stats::pchisq(t_val, df = d, ncp = lambda2) + 1 - stats::pchisq(c_ratio * t_val, df = d, ncp = lambda1)
  }
}

#' @title Kullback-Leibler Divergence Metric
#' @description The (asymmetric) KL divergence between two groups' Gaussian
#'   posteriors, \eqn{\mathrm{KL}(\mathcal N_1\|\mathcal N_2)} (\code{direction = "1to2"})
#'   or \eqn{\mathrm{KL}(\mathcal N_2\|\mathcal N_1)} (\code{"2to1"}).
#' @param direction Which direction of the (asymmetric) divergence to compute.
#' @return A \code{distance_metric} object.
#' @export
kl_metric <- function(direction = c("1to2", "2to1")) {
  direction <- match.arg(direction)
  structure(list(direction = direction), class = c("kl_metric", "distance_metric"))
}

#' @export
is_symmetric_metric.kl_metric <- function(metric) FALSE

#' @export
evaluate_metric.kl_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  if (identical(metric$direction, "2to1")) {
    .kl_divergence(mu2, mu1, Sigma2, Sigma1)
  } else {
    .kl_divergence(mu1, mu2, Sigma1, Sigma2)
  }
}

#' @title Symmetrized Jeffreys Divergence Metric
#' @description \eqn{J = \mathrm{KL}(\mathcal N_1\|\mathcal N_2) + \mathrm{KL}(\mathcal N_2\|\mathcal N_1)},
#'   computed by summing two \code{\link{kl_metric}} evaluations rather than
#'   re-deriving the formula.
#' @return A \code{distance_metric} object.
#' @export
jeffreys_metric <- function() {
  structure(list(), class = c("jeffreys_metric", "distance_metric"))
}

#' @export
evaluate_metric.jeffreys_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  .kl_divergence(mu1, mu2, Sigma1, Sigma2) + .kl_divergence(mu2, mu1, Sigma2, Sigma1)
}

#' @title Bhattacharyya Distance Metric
#' @description The Bhattacharyya distance between two groups' Gaussian posteriors.
#' @return A \code{distance_metric} object.
#' @export
bhattacharyya_metric <- function() {
  structure(list(), class = c("bhattacharyya_metric", "distance_metric"))
}

#' @export
evaluate_metric.bhattacharyya_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  Sigma_avg <- (Sigma1 + Sigma2) / 2
  delta <- mu1 - mu2
  inv_avg <- chol_inv_jitter(Sigma_avg, pen_diag = 1e-6)
  logdet_avg <- as.numeric(determinant(Sigma_avg, logarithm = TRUE)$modulus)
  logdet1 <- as.numeric(determinant(Sigma1, logarithm = TRUE)$modulus)
  logdet2 <- as.numeric(determinant(Sigma2, logarithm = TRUE)$modulus)
  0.125 * as.numeric(t(delta) %*% inv_avg %*% delta) + 0.5 * (logdet_avg - 0.5 * (logdet1 + logdet2))
}

#' @title Hellinger Distance Metric
#' @description \eqn{H = \sqrt{1 - e^{-D_B}}}, derived from the
#'   \code{\link{bhattacharyya_metric}} distance \eqn{D_B}. Bounded in \eqn{[0, 1]}.
#' @return A \code{distance_metric} object.
#' @export
hellinger_metric <- function() {
  structure(list(), class = c("hellinger_metric", "distance_metric"))
}

#' @export
evaluate_metric.hellinger_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  db <- evaluate_metric(bhattacharyya_metric(), mu1, mu2, Sigma1, Sigma2)
  sqrt(1 - exp(-db))
}

#' @title Wasserstein-2 (Frechet) Distance Metric
#' @description The 2-Wasserstein distance between two multivariate Gaussian
#'   posteriors. When the two groups share the same \code{kernel_key} (so
#'   \eqn{\Sigma_2 = c\Sigma_1}), a closed-form shortcut is used that avoids a
#'   matrix square root entirely; the general eigendecomposition-based formula
#'   is used otherwise.
#' @return A \code{distance_metric} object.
#' @export
wasserstein_metric <- function() {
  structure(list(), class = c("wasserstein_metric", "distance_metric"))
}

#' @export
evaluate_metric.wasserstein_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  mean_term <- sum((mu1 - mu2)^2)
  if (same_kernel && !is.null(scale1) && !is.null(scale2)) {
    c_ratio <- scale1 / scale2
    return(sqrt(mean_term + (1 - sqrt(c_ratio))^2 * sum(diag(Sigma1))))
  }
  sqrt1 <- .matrix_sqrt_sym(Sigma1)
  inner <- sqrt1 %*% Sigma2 %*% sqrt1
  sqrt_inner <- .matrix_sqrt_sym(inner)
  sqrt(max(mean_term + sum(diag(Sigma1 + Sigma2 - 2 * sqrt_inner)), 0))
}

#' @title Mean Absolute Difference of Means Metric
#' @description The mean absolute difference between the two groups' posterior
#'   means over shared IDs, evaluated marginally (per ID) and ignoring the
#'   covariance structure entirely: \eqn{\frac{1}{d}\sum_k |\mu_{1k}-\mu_{2k}|}.
#'   \code{Sigma1}/\code{Sigma2} are accepted (for interface consistency with
#'   \code{\link{evaluate_metric}}) but unused.
#' @return A \code{distance_metric} object.
#' @export
mean_diff_metric <- function() {
  structure(list(), class = c("mean_diff_metric", "distance_metric"))
}

#' @export
evaluate_metric.mean_diff_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  mean(abs(mu1 - mu2))
}

# =============================================================================
# Decorators: wrap a base metric to get a variant
# =============================================================================

#' @title Per-Feature-Normalized Metric
#' @description Divides a base metric's value by \eqn{d^{power}} (\eqn{d} =
#'   number of shared IDs), turning a joint quantity into a per-feature rate
#'   that stays roughly constant (rather than growing/collapsing) as \eqn{d}
#'   changes under a constant per-feature effect. Forwards
#'   \code{\link{requires_shared_kernel}}/\code{\link{is_symmetric_metric}} to
#'   the wrapped metric.
#' @param base A \code{distance_metric} object to normalize.
#' @param power The exponent applied to \eqn{d} (e.g. \code{0.5} for \eqn{\sqrt d}, \code{1} for \eqn{d}).
#' @return A \code{distance_metric} object.
#' @export
per_feature_metric <- function(base, power) {
  if (!inherits(base, "distance_metric")) {
    stop("'base' must be a distance_metric object.")
  }
  structure(list(base = base, power = power), class = c("per_feature_metric", "distance_metric"))
}

#' @export
requires_shared_kernel.per_feature_metric <- function(metric) requires_shared_kernel(metric$base)

#' @export
is_symmetric_metric.per_feature_metric <- function(metric) is_symmetric_metric(metric$base)

#' @export
evaluate_metric.per_feature_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  d <- length(mu1)
  base_val <- evaluate_metric(metric$base, mu1, mu2, Sigma1, Sigma2, scale1 = scale1, scale2 = scale2, same_kernel = same_kernel, ...)
  base_val / d^metric$power
}

#' @title Marginal-Averaged Metric
#' @description Averages a base metric's *univariate* value over each shared
#'   ID individually, rather than evaluating it jointly on the full
#'   \eqn{d}-dimensional posterior. This sidesteps metrics like
#'   \code{\link{ovl_metric}} requiring *every* dimension to overlap
#'   simultaneously, at the cost of losing joint (cross-ID correlation)
#'   information. Forwards \code{\link{requires_shared_kernel}}/
#'   \code{\link{is_symmetric_metric}} to the wrapped metric.
#' @param base A \code{distance_metric} object to marginalize.
#' @return A \code{distance_metric} object.
#' @export
marginal_metric <- function(base) {
  if (!inherits(base, "distance_metric")) {
    stop("'base' must be a distance_metric object.")
  }
  structure(list(base = base), class = c("marginal_metric", "distance_metric"))
}

#' @export
requires_shared_kernel.marginal_metric <- function(metric) requires_shared_kernel(metric$base)

#' @export
is_symmetric_metric.marginal_metric <- function(metric) is_symmetric_metric(metric$base)

#' @export
evaluate_metric.marginal_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  ids <- names(mu1)
  vals <- vapply(ids, function(k) {
    evaluate_metric(metric$base, mu1[k], mu2[k], Sigma1[k, k, drop = FALSE], Sigma2[k, k, drop = FALSE],
                     scale1 = scale1, scale2 = scale2, same_kernel = same_kernel, ...)
  }, numeric(1))
  mean(vals)
}

#' @title Total Variation Distance
#' @description Total variation distance between the two groups' posteriors,
#'   \eqn{TV(f,g) = 1 - \int \min(f,g)\,dx}. This is a general identity for
#'   any two densities (not an approximation specific to Gaussians): since
#'   \eqn{\int(f+g)\,dx = 2}, \eqn{\int|f-g|\,dx = 2 - 2\int\min(f,g)\,dx}, so
#'   \eqn{TV = \tfrac12\int|f-g|\,dx = 1 - \int\min(f,g)\,dx} always. Because
#'   \code{\link{ovl_metric}} already computes \eqn{\int\min(f,g)\,dx} in
#'   closed form, TVD is just \code{1 - } that value -- a thin decorator
#'   around it rather than a new closed-form derivation. Forwards
#'   \code{\link{requires_shared_kernel}}/\code{\link{is_symmetric_metric}} to
#'   the wrapped metric.
#' @param base A \code{distance_metric} object whose value is an overlap
#'   coefficient (\code{\link{ovl_metric}} by default). Passing anything else
#'   is meaningless (TVD is only defined via an overlap coefficient), but
#'   left as a parameter so a decorated OVL variant (e.g. \code{marginal_metric(ovl_metric())})
#'   can be turned into its own TVD analogue too.
#' @return A \code{distance_metric} object.
#' @export
tvd_metric <- function(base = ovl_metric()) {
  if (!inherits(base, "distance_metric")) {
    stop("'base' must be a distance_metric object.")
  }
  structure(list(base = base), class = c("tvd_metric", "distance_metric"))
}

#' @export
requires_shared_kernel.tvd_metric <- function(metric) requires_shared_kernel(metric$base)

#' @export
is_symmetric_metric.tvd_metric <- function(metric) is_symmetric_metric(metric$base)

#' @export
evaluate_metric.tvd_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  1 - evaluate_metric(metric$base, mu1, mu2, Sigma1, Sigma2, scale1 = scale1, scale2 = scale2, same_kernel = same_kernel, ...)
}

#' @title Fraction of Individually Significant Features
#' @description For each shared ID, computes the z-score
#'   \eqn{(\mu_{1k}-\mu_{2k})/\sqrt{\Sigma_{1,kk}+\Sigma_{2,kk}}} and returns
#'   the proportion of IDs whose absolute z-score exceeds \code{threshold} --
#'   a nonparametric, per-feature summary that stays a proportion in
#'   \eqn{[0, 1]} regardless of \eqn{d}.
#' @param threshold The absolute z-score threshold above which an ID counts
#'   as individually significant. Defaults to \code{1.96}.
#' @return A \code{distance_metric} object.
#' @export
significant_fraction_metric <- function(threshold = 1.96) {
  structure(list(threshold = threshold), class = c("significant_fraction_metric", "distance_metric"))
}

#' @export
evaluate_metric.significant_fraction_metric <- function(metric, mu1, mu2, Sigma1, Sigma2, scale1 = NULL, scale2 = NULL, same_kernel = FALSE, ...) {
  ids <- names(mu1)
  var1 <- diag(Sigma1)[ids]
  var2 <- diag(Sigma2)[ids]
  z <- (mu1 - mu2) / sqrt(var1 + var2)
  mean(abs(z) > metric$threshold)
}

# =============================================================================
# Driver
# =============================================================================

#' @title Compute a Pairwise Group-Difference Metric
#'
#' @description
#' Generalizes \code{\link{calculate_group_overlaps}} to any
#' \code{distance_metric}: for every pair of groups in \code{results},
#' aligns their posterior means/covariances by shared ID and evaluates
#' \code{metric} via \code{\link{evaluate_metric}}. Every metric requires the
#' two groups being compared to share the same set of IDs (checked via
#' \code{setequal}); metrics for which \code{\link{requires_shared_kernel}}
#' is \code{TRUE} additionally require the two groups to share the same
#' \code{kernel_key}.
#'
#' @param results A list, typically from \code{\link{posterior_mean}},
#'   with elements \code{kernels} and \code{groups} (one entry per group,
#'   each with \code{muk}, \code{id_to_input}, \code{kernel_key}, \code{scale}).
#' @param metric A \code{distance_metric} object, e.g.
#'   \code{\link{ovl_metric}()}, \code{\link{kl_metric}()}, \code{\link{wasserstein_metric}()}.
#' @param max_groups_warn Emit a warning about the \eqn{O(G^2 d^3)} cost of
#'   the underlying matrix inversions if the number of groups exceeds this. Defaults to \code{50}.
#' @param max_dim_warn Emit the same warning if the number of shared IDs \eqn{d}
#'   exceeds this. Defaults to \code{500}.
#'
#' @return A square matrix of metric values, one row/column per group. The
#'   diagonal is each group's self-comparison value (computed via
#'   \code{evaluate_metric}, so it is \code{1} for \code{\link{ovl_metric}}
#'   and \code{0} for distance-like metrics, without special-casing). For
#'   metrics where \code{\link{is_symmetric_metric}} is \code{FALSE} (e.g.
#'   \code{\link{kl_metric}}), the matrix itself is not symmetric.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- posterior_mean(data, kern)
#' compute_group_diff(posterior, mahalanobis_metric())
compute_group_diff <- function(results, metric, max_groups_warn = 50, max_dim_warn = 500) {
  if (!is.list(results) || !all(c("kernels", "groups") %in% names(results))) {
    stop("'results' must be the list returned by posterior_mean() (with 'kernels' and 'groups').")
  }
  if (!inherits(metric, "distance_metric")) {
    stop("'metric' must be a distance_metric object (e.g. ovl_metric(), kl_metric(), wasserstein_metric()).")
  }
  groups      <- results$groups
  group_names <- names(groups)
  num_groups  <- length(group_names)

  required_group_fields <- c("muk", "id_to_input", "kernel_key", "scale")
  missing_fields <- vapply(groups, function(g) !all(required_group_fields %in% names(g)), logical(1))
  if (any(missing_fields)) {
    stop(paste0(
      "Each group in results$groups must contain 'muk', 'id_to_input', 'kernel_key' and 'scale'; ",
      "missing in group(s): ", paste(group_names[missing_fields], collapse = ", "), "."
    ))
  }

  if (num_groups >= 2) {
    d <- length(groups[[1]]$muk)
    if (num_groups > max_groups_warn || d > max_dim_warn) {
      warning(sprintf(
        paste0(
          "compute_group_diff: comparing %d groups with ~%d shared IDs requires ",
          "up to %d matrix inversions of size %dx%d (cost is O(G^2 * d^3)); this may be slow."
        ),
        num_groups, d, choose(num_groups, 2), d, d
      ))
    }
  }

  diff_matrix <- matrix(0, num_groups, num_groups)
  rownames(diff_matrix) <- colnames(diff_matrix) <- group_names

  if (num_groups == 0) return(diff_matrix)

  for (i in seq_len(num_groups)) {
    g <- groups[[i]]
    ids <- names(g$muk)
    Sigma <- get_sigmak(g, results$kernels)[ids, ids, drop = FALSE]
    diff_matrix[i, i] <- evaluate_metric(metric, g$muk[ids], g$muk[ids], Sigma, Sigma,
                                          scale1 = g$scale, scale2 = g$scale, same_kernel = TRUE)
  }

  if (num_groups < 2) return(diff_matrix)

  symmetric           <- is_symmetric_metric(metric)
  needs_shared_kernel  <- requires_shared_kernel(metric)

  for (i in seq_len(num_groups - 1)) {
    for (j in (i + 1):num_groups) {
      group1 <- groups[[i]]
      group2 <- groups[[j]]

      ids1 <- names(group1$muk)
      ids2 <- names(group2$muk)
      if (!setequal(ids1, ids2)) {
        stop(paste0(
          "Groups '", group_names[i], "' and '", group_names[j], "' do not share the same set of IDs: ",
          "only in '", group_names[i], "': [", paste(setdiff(ids1, ids2), collapse = ", "), "]; ",
          "only in '", group_names[j], "': [", paste(setdiff(ids2, ids1), collapse = ", "), "]."
        ))
      }

      same_kernel <- identical(group1$kernel_key, group2$kernel_key)
      if (needs_shared_kernel && !same_kernel) {
        stop(paste0(
          "Groups '", group_names[i], "' and '", group_names[j], "' do not share the same kernel matrix ",
          "(different sets of Input values). This metric requires both groups' posterior ",
          "covariances to derive from the same raw kernel matrix, only scaled differently by 'scale'."
        ))
      }

      mu1    <- group1$muk[ids1]
      mu2    <- group2$muk[ids1]
      Sigma1 <- get_sigmak(group1, results$kernels)[ids1, ids1, drop = FALSE]
      Sigma2 <- get_sigmak(group2, results$kernels)[ids1, ids1, drop = FALSE]

      val_ij <- evaluate_metric(metric, mu1, mu2, Sigma1, Sigma2,
                                 scale1 = group1$scale, scale2 = group2$scale, same_kernel = same_kernel)
      diff_matrix[group_names[i], group_names[j]] <- val_ij

      if (symmetric) {
        diff_matrix[group_names[j], group_names[i]] <- val_ij
      } else {
        val_ji <- evaluate_metric(metric, mu2, mu1, Sigma2, Sigma1,
                                   scale1 = group2$scale, scale2 = group1$scale, same_kernel = same_kernel)
        diff_matrix[group_names[j], group_names[i]] <- val_ji
      }
    }
  }
  diff_matrix
}

## ===========================================================================
## group_diff(): convenience wrapper -- sensible default metric, nicer print.
## ===========================================================================

#' @noRd
#'
#' @details Builds a human-readable call-like label for a \code{distance_metric}
#' object from its class name and fields, e.g. \code{ovl_metric()},
#' \code{per_feature_metric(wasserstein_metric(), power = 0.5)},
#' \code{kl_metric(direction = 1to2)}. Generic (works for any current or
#' future metric via introspection) rather than hardcoded per class: a
#' \code{base} field (decorators) is recursed into first, then every other
#' field is rendered as \code{name = value}.
.metric_label <- function(metric) {
  ctor <- class(metric)[1]
  args <- character(0)
  if (!is.null(metric$base)) {
    args <- c(args, .metric_label(metric$base))
  }
  other <- metric[setdiff(names(metric), "base")]
  if (length(other) > 0) {
    args <- c(args, paste0(names(other), " = ", vapply(other, function(v) paste(format(v), collapse = ", "), character(1))))
  }
  paste0(ctor, "(", paste(args, collapse = ", "), ")")
}

#' @title Compare Two Groups' Posteriors, with a Sensible Default Metric
#'
#' @description
#' A convenience wrapper around \code{\link{compute_group_diff}}. When
#' \code{metric = NULL} (the default), the metric is chosen automatically
#' from the number of shared IDs \code{d}: \code{\link{ovl_metric}()} for
#' \code{d <= id_threshold} (still interpretable, not yet dimension-collapsed),
#' or \code{\link{per_feature_metric}(\link{wasserstein_metric}(), power =
#' 0.5)} otherwise -- a per-feature-normalized distance that keeps growing
#' informatively with \code{d} instead of collapsing toward a degenerate
#' value (see \code{dev/30_package_dev/distance_metrics/distance-metrics.Rmd}
#' for the empirical study behind this choice, and the "Choice of metric"
#' article). Supplying \code{metric} explicitly bypasses this choice
#' entirely and behaves exactly like \code{\link{compute_group_diff}}.
#'
#' The returned matrix carries the metric actually used (as a label) and
#' prints it above the matrix, so an auto-selected default is never silently
#' invisible even if the result is stored and printed again later.
#'
#' @param results A list, typically from \code{\link{posterior_mean}}; see
#'   \code{\link{compute_group_diff}}.
#' @param metric A \code{distance_metric} object, or \code{NULL} (default) to
#'   auto-select based on the number of shared IDs (see \code{id_threshold}).
#' @param id_threshold The number of shared IDs at or below which
#'   \code{\link{ovl_metric}()} is auto-selected (above it,
#'   \code{per_feature_metric(wasserstein_metric(), power = 0.5)} is used
#'   instead). Only used when \code{metric = NULL}. Defaults to \code{6}.
#' @param ... Forwarded to \code{\link{compute_group_diff}} (e.g.
#'   \code{max_groups_warn}, \code{max_dim_warn}).
#'
#' @return A \code{group_diff_result} object: the same matrix
#'   \code{\link{compute_group_diff}} returns (still usable as a plain
#'   matrix, e.g. \code{result["g1", "g2"]}), with an extra
#'   \code{"metric_label"} attribute and its own \code{print} method.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- posterior_mean(data, kern)
#' group_diff(posterior)                # d = 8 > 6 -> per_feature_metric(wasserstein_metric(), power = 0.5)
#' group_diff(posterior, ovl_metric())  # explicit metric, bypasses auto-selection
group_diff <- function(results, metric = NULL, id_threshold = 6, ...) {
  if (!is.list(results) || !all(c("kernels", "groups") %in% names(results)) || length(results$groups) == 0) {
    stop("'results' must be the list returned by posterior_mean() (with 'kernels' and 'groups').")
  }
  if (is.null(metric)) {
    d <- length(results$groups[[1]]$muk)
    metric <- if (d <= id_threshold) ovl_metric() else per_feature_metric(wasserstein_metric(), power = 0.5)
  } else if (!inherits(metric, "distance_metric")) {
    stop("'metric' must be NULL or a distance_metric object (e.g. ovl_metric(), kl_metric(), wasserstein_metric()).")
  }
  mat <- compute_group_diff(results, metric, ...)
  attr(mat, "metric_label") <- .metric_label(metric)
  class(mat) <- c("group_diff_result", class(mat))
  mat
}

#' @title Print a Group-Difference Result
#' @description Prints the metric that was used (see \code{\link{group_diff}}),
#'   then the underlying matrix.
#' @param x A \code{group_diff_result} object, from \code{\link{group_diff}}.
#' @param ... Unused, included for S3 consistency.
#' @return \code{x}, invisibly.
#' @export
print.group_diff_result <- function(x, ...) {
  cat("Metric:", attr(x, "metric_label"), "\n\n")
  plain <- x
  attr(plain, "metric_label") <- NULL
  class(plain) <- setdiff(class(plain), "group_diff_result")
  print(plain, ...)
  invisible(x)
}
