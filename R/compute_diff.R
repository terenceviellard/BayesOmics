#' @importFrom stats pnorm pchisq

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
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' calculate_group_overlaps(posterior)
calculate_group_overlaps <- function(results, max_groups_warn = 50, max_dim_warn = 500) {
  if (!is.list(results) || !all(c("kernels", "groups") %in% names(results))) {
    stop("'results' must be the list returned by multi_posterior_mean() (with 'kernels' and 'groups').")
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
    # TODO: no test in tests/testthat/ verifies that this warning actually
    # fires at the documented thresholds (max_groups_warn/max_dim_warn) --
    # only the message text seems covered, not the boundary condition itself
    # (e.g. num_groups == max_groups_warn vs max_groups_warn + 1).
    if (num_groups > max_groups_warn || d > max_dim_warn) {
      warning(sprintf(
        paste0(
          "calculate_group_overlaps: comparing %d groups with ~%d shared IDs requires ",
          "up to %d matrix inversions of size %dx%d (cost is O(G^2 * d^3)); this may be slow."
        ),
        num_groups, d, choose(num_groups, 2), d, d
      ))
    }
  }

  overlap_matrix <- diag(num_groups)
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- group_names

  if (num_groups < 2) return(overlap_matrix)

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

      if (!identical(group1$kernel_key, group2$kernel_key)) {
        stop(paste0(
          "Groups '", group_names[i], "' and '", group_names[j], "' do not share the same kernel matrix ",
          "(different sets of Input values). The closed-form overlap requires both groups' posterior ",
          "covariances to derive from the same raw kernel matrix, only scaled differently by 'scale'."
        ))
      }

      mu1    <- group1$muk[ids1]
      mu2    <- group2$muk[ids1]
      d      <- length(ids1)
      scale1 <- group1$scale
      scale2 <- group2$scale
      delta  <- mu1 - mu2
      c_ratio <- scale1 / scale2

      scale_rel_tol <- 1e-6
      near_equal_scale <- isTRUE(scale1 == scale2) || abs(c_ratio - 1) < scale_rel_tol

      if (near_equal_scale) {
        if (!isTRUE(scale1 == scale2)) {
          warning(sprintf(
            paste0(
              "calculate_group_overlaps: groups '%s' and '%s' have nearly identical ",
              "scales (scale1 = %.6g, scale2 = %.6g, ratio = %.10f); the c != 1 formula ",
              "is numerically unstable this close to c = 1, so they are treated as equal ",
              "using the larger posterior covariance (smaller scale)."
            ),
            group_names[i], group_names[j], scale1, scale2, c_ratio
          ))
        }
        smaller_group <- if (scale1 <= scale2) group1 else group2
        Sigma_eq <- get_sigmak(smaller_group, results$kernels)[ids1, ids1, drop = FALSE]
        inv_eq   <- chol_inv_jitter(Sigma_eq, pen_diag = 1e-6)
        D2       <- as.numeric(t(delta) %*% inv_eq %*% delta)
        OV <- 2 * stats::pnorm(-sqrt(D2) / 2)
      } else {
        Sigma1 <- get_sigmak(group1, results$kernels)[ids1, ids1, drop = FALSE]
        inv1   <- chol_inv_jitter(Sigma1, pen_diag = 1e-6)
        D2     <- as.numeric(t(delta) %*% inv1 %*% delta)
        lambda1 <- D2 / (1 - c_ratio)^2
        lambda2 <- c_ratio * D2 / (1 - c_ratio)^2
        t       <- (D2 - d * (1 - c_ratio) * log(c_ratio)) / (1 - c_ratio)^2
        if (c_ratio < 1) {
          OV <- stats::pchisq(c_ratio * t, df = d, ncp = lambda1) +
            1 - stats::pchisq(t, df = d, ncp = lambda2)
        } else {
          OV <- stats::pchisq(t, df = d, ncp = lambda2) +
            1 - stats::pchisq(c_ratio * t, df = d, ncp = lambda1)
        }
      }

      overlap_matrix[group_names[i], group_names[j]] <- OV
      overlap_matrix[group_names[j], group_names[i]] <- OV
    }
  }
  return(overlap_matrix)
}
