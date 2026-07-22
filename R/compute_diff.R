#' @importFrom stats pnorm pchisq
#' @importFrom dplyr %>%
#' @importFrom rlang .data

#' @title Compute Overlapping Coefficient between Groups
#'
#' @description
#' Computes a symmetric matrix of pairwise overlapping coefficients (OVL) between
#' every pair of groups.
#'
#' For two groups sharing the same underlying kernel matrix (same \code{kernel_key},
#' i.e. the same set of Input values), their posterior covariances are
#' \eqn{\Sigma_1 = \Sigma/(\lambda_0+N_1)} and \eqn{\Sigma_2 = \Sigma/(\lambda_0+N_2)} for
#' a shared raw kernel matrix \eqn{\Sigma}, so \eqn{\Sigma_2 = c\,\Sigma_1} with
#' \eqn{c = \mathrm{scale}_1/\mathrm{scale}_2}. Writing
#' \eqn{D^2 = (\mu_2-\mu_1)^\top\Sigma_1^{-1}(\mu_2-\mu_1)} and \eqn{d} the dimension
#' (number of shared IDs), the exact Gaussian overlap is used:
#'
#' \itemize{
#'   \item \eqn{c = 1}: \deqn{OVL = 2 \Phi\!\left(-D/2\right)}
#'   \item \eqn{c \ne 1}, with
#'     \eqn{\lambda_1 = D^2/(1-c)^2}, \eqn{\lambda_2 = cD^2/(1-c)^2},
#'     \eqn{t = (D^2 - d(1-c)\ln c)/(1-c)^2}:
#'     \deqn{OVL = \begin{cases}
#'       F_{\chi^2_d(\lambda_1)}(ct) + 1 - F_{\chi^2_d(\lambda_2)}(t), & 0<c<1 \\
#'       F_{\chi^2_d(\lambda_2)}(t) + 1 - F_{\chi^2_d(\lambda_1)}(ct), & c>1
#'     \end{cases}}
#' }
#'
#' Two groups with different \code{kernel_key} (different sets of Input values) do not
#' share a common raw kernel matrix \eqn{\Sigma}, so this closed-form ratio does not
#' apply; an error is raised in that case.
#'
#' When \eqn{c} is not exactly 1 but very close to it (relative difference below
#' \code{1e-6}), the \eqn{c \ne 1} formula above becomes numerically unstable
#' (it divides by \eqn{(1-c)^2}). In that case a warning is issued and the
#' \eqn{c = 1} formula is used instead, with \eqn{\Sigma_1} taken from whichever
#' group has the smaller scale (i.e. the larger, more conservative posterior
#' covariance).
#'
#' Cost: each pair of groups requires one matrix inversion of size \eqn{d \times d}
#' (\eqn{d} = number of shared IDs), so the total cost is \eqn{O(G^2 d^3)} for
#' \eqn{G} groups -- e.g. 10 groups with 400 shared IDs already means about 45
#' inversions of 400x400 matrices. A warning is issued if this is likely to be
#' slow (see \code{max_groups_warn}/\code{max_dim_warn}).
#'
#' @param results A list, typically from \code{\link{multi_posterior_mean}},
#'   with elements \code{kernels} and \code{groups} (one entry per group,
#'   each with \code{muk}, \code{id_to_input}, \code{kernel_key}, \code{scale}).
#' @param max_groups_warn Emit a warning about the \eqn{O(G^2 d^3)} cost above
#'   (see Description) if the number of groups exceeds this. Defaults to \code{50}.
#' @param max_dim_warn Emit the same warning if the number of shared IDs \eqn{d}
#'   exceeds this. Defaults to \code{500}.
#'
#' @return A symmetric matrix of OVL coefficients in \eqn{[0, 1]}, with 1 on the diagonal.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- multi_posterior_mean(data, kern)
#' calculate_group_overlaps(posterior)
calculate_group_overlaps <- function(results, max_groups_warn = 50, max_dim_warn = 500) {
  compute_group_diff(results, ovl_metric(), max_groups_warn = max_groups_warn, max_dim_warn = max_dim_warn)
}

## Extracts an (n_draws x n_ids) matrix of posterior draws for one group, with
## columns in the caller-supplied `ids` order. Relies on sample_posterior()'s
## documented melt order (Sample = as.vector(mat), one n_draws x n_ids matrix
## per group) being preserved in `sample_distrib`: filtering once for the
## whole group and splitting locally (rather than filtering once per ID) is a
## single O(n_rows) pass, and the per-ID length check is the closest available
## sanity check on that ordering invariant actually holding.
#' @noRd
extract_draw_matrix <- function(sample_distrib, group, ids) {
  sub   <- sample_distrib %>% dplyr::filter(.data$Group == group, .data$ID %in% ids)
  by_id <- split(sub$Sample, sub$ID)
  lens  <- vapply(by_id[ids], length, integer(1))
  if (length(unique(lens)) > 1) {
    stop(paste0(
      "Group '", group, "' has IDs with different numbers of posterior draws in ",
      "'sample_distrib' (", paste(sort(unique(lens)), collapse = ", "), "); ",
      "compute_multi_diff() requires every ID within a group to contribute the same ",
      "number of draws, in the same draw order, as produced by sample_posterior()."
    ))
  }
  mat <- do.call(cbind, by_id[ids])
  colnames(mat) <- ids
  mat
}

#' @title Compute the Multivariate Distribution of Group Differences
#'
#' @description
#' For every pair of groups present in \code{sample_distrib}, computes the
#' empirical distribution (over posterior draws) of the number of IDs for
#' which group1's posterior draw exceeds group2's -- a multivariate,
#' uncertainty-aware complement to the single scalar overlapping coefficient
#' returned by \code{\link{calculate_group_overlaps}}.
#'
#' This relies on a documented invariant of \code{\link{sample_posterior}}'s
#' output: for a given group, its \code{n} draws are melted from a single
#' \code{n x n_ids} matrix (\code{Sample = as.vector(mat)}), so filtering
#' \code{sample_distrib} by \code{(Group, ID)} recovers draws in a consistent
#' order across every ID in that group, without needing an explicit draw index
#' column. This only holds if \code{sample_distrib}'s row order has not been
#' shuffled since \code{sample_posterior()} produced it, and every ID within a
#' group contributes the same number of draws (checked, and an error raised
#' otherwise; the actual per-draw pairing itself cannot be independently
#' verified from a melted data frame with no explicit draw index).
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}.
#' @param results An optional list, typically from
#'    \code{\link{multi_posterior_mean}}, with elements \code{kernels} and
#'    \code{groups}. If supplied, \code{\link{calculate_group_overlaps}} is
#'    used to attach an exact \code{Overlap_coef} to the result. If \code{NULL}
#'    (default), \code{Overlap_coef} is omitted.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{Diff_proba}}{A tibble with columns \code{Group1}, \code{Group2},
#'       \code{Nb_id} (from 0 to the number of shared IDs), \code{Proba} (the
#'       probability mass at that count) and \code{Cumul_proba} (its cumulative sum).}
#'     \item{\code{Diff_mean}}{A tibble with columns \code{ID}, \code{Group} and
#'       \code{Mean}, the posterior mean of each (ID, Group) pair.}
#'     \item{\code{Overlap_coef}}{A tibble with columns \code{Group1}, \code{Group2}
#'       and \code{Overlap_coef}. Only present when \code{results} is supplied.}
#'   }
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 3, nb_sample = 5, diff_group = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' multi_diff <- compute_multi_diff(samples, results = posterior)
#' multi_diff$Diff_proba
compute_multi_diff <- function(sample_distrib, results = NULL) {
  required_cols <- c("ID", "Group", "Sample")
  if (!all(required_cols %in% names(sample_distrib))) {
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  groups <- sample_distrib$Group %>% unique() %>% sort()
  if (length(groups) < 2) {
    stop("compute_multi_diff() requires at least two groups in 'sample_distrib'.")
  }
  pairs <- utils::combn(groups, 2, simplify = FALSE)

  proba_rows <- lapply(pairs, function(p) {
    g1 <- p[1]
    g2 <- p[2]

    ids1 <- sample_distrib %>% dplyr::filter(.data$Group == g1) %>% dplyr::pull(.data$ID) %>% unique()
    ids2 <- sample_distrib %>% dplyr::filter(.data$Group == g2) %>% dplyr::pull(.data$ID) %>% unique()
    if (!setequal(ids1, ids2)) {
      stop(paste0(
        "Groups '", g1, "' and '", g2, "' do not share the same set of IDs: ",
        "only in '", g1, "': [", paste(setdiff(ids1, ids2), collapse = ", "), "]; ",
        "only in '", g2, "': [", paste(setdiff(ids2, ids1), collapse = ", "), "]."
      ))
    }
    shared_ids <- sort(ids1)

    mat1 <- extract_draw_matrix(sample_distrib, g1, shared_ids)
    mat2 <- extract_draw_matrix(sample_distrib, g2, shared_ids)
    if (nrow(mat1) != nrow(mat2)) {
      stop(paste0(
        "Groups '", g1, "' and '", g2, "' have different numbers of posterior draws (",
        nrow(mat1), " vs ", nrow(mat2), "); compute_multi_diff() requires the same number ",
        "of draws for every group being compared (as produced by a single ",
        "sample_posterior(results, n) call)."
      ))
    }

    n_ids   <- length(shared_ids)
    n_draws <- nrow(mat1)
    counts  <- rowSums(mat1 > mat2)
    tab     <- tabulate(counts + 1, nbins = n_ids + 1)
    proba   <- tab / n_draws
    cumul   <- cumsum(proba)

    tibble::tibble(Group1 = g1, Group2 = g2, Nb_id = 0:n_ids, Proba = proba, Cumul_proba = cumul)
  })
  Diff_proba <- dplyr::bind_rows(proba_rows)

  Diff_mean <- sample_distrib %>%
    dplyr::group_by(.data$ID, .data$Group) %>%
    dplyr::summarise(Mean = mean(.data$Sample), .groups = "drop")

  out <- list(Diff_proba = Diff_proba, Diff_mean = Diff_mean)

  if (!is.null(results)) {
    if (!is.list(results) || !all(c("kernels", "groups") %in% names(results))) {
      stop("'results' must be the list returned by multi_posterior_mean() (with 'kernels' and 'groups'), or NULL.")
    }
    result_groups <- names(results$groups)
    if (!setequal(result_groups, groups)) {
      stop(paste0(
        "'results' groups do not match the groups present in 'sample_distrib': ",
        "only in 'results': [", paste(setdiff(result_groups, groups), collapse = ", "), "]; ",
        "only in 'sample_distrib': [", paste(setdiff(groups, result_groups), collapse = ", "), "]."
      ))
    }
    overlap_mat <- calculate_group_overlaps(results)
    out$Overlap_coef <- dplyr::bind_rows(lapply(pairs, function(p) {
      tibble::tibble(Group1 = p[1], Group2 = p[2], Overlap_coef = overlap_mat[p[1], p[2]])
    }))
  }

  out
}
