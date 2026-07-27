#' @noRd
# TODO: format_input_key() and get_sigmak() are @noRd (internal) but are
# called directly from tests/testthat/ -- their signature is therefore a de
# facto stable contract without being part of the public API. Either export
# them officially or accept that test breakage on refactor is expected.
format_input_key <- function(x) sprintf("%.17g", x)

#' @noRd
#'
#' @details `group_entry$id_to_input` is either a named numeric vector (one
#'   scalar position per ID -- the legacy/hand-crafted-fixture shape) or an
#'   n_ids x D matrix with `rownames = ID` (the shape `multi_posterior_mean()`
#'   itself now always produces, D >= 1). Both are turned into one
#'   `row_input_key()` string per ID, so re-indexing the shared kernel matrix
#'   works identically either way.
#' @details `group_entry$obs_noise` (added by `multi_posterior_mean()`, default
#'   0) is a fixed, known observation-noise variance to add back to the
#'   kernel-derived covariance before dividing by `scale`. It exists because
#'   HP-fitting workflows built around `optim_hp(..., prior_cov = )` (see
#'   `R/optim_kernel.R`'s `resolve_prior_cov()`) add `prior_cov` as a SEPARATE
#'   additive nugget to the likelihood, so the fitted kernel's own HPs (e.g.
#'   `variance`) only capture whatever correlated structure remains BEYOND
#'   that fixed floor -- not the full marginal noise variance the posterior
#'   formula `Sigma_theta / (N + lambda_0)` assumes. Passing the SAME
#'   `prior_cov` used at fit time as `obs_noise` here reconstructs the correct
#'   total variance. Left at its default of 0, behavior is unchanged from
#'   before this parameter existed (kernel HPs assumed to already capture the
#'   full marginal variance themselves, e.g. via an explicit `NoiseKernel()`
#'   term composed into the kernel).
get_sigmak <- function(group_entry, kernels) {
  ids <- names(group_entry$muk)
  id_to_input <- group_entry$id_to_input
  sub <- if (is.matrix(id_to_input)) id_to_input[ids, , drop = FALSE] else id_to_input[ids]
  input_keys <- row_input_key(sub)
  kern_mat <- kernels[[group_entry$kernel_key]]
  noise <- group_entry$obs_noise
  if (is.null(noise)) noise <- 0
  sigmak <- (kern_mat[input_keys, input_keys, drop = FALSE] + diag(noise, length(ids))) / group_entry$scale
  dimnames(sigmak) <- list(ids, ids)
  sigmak
}

#' Compute Posterior Means for Different Groups
#'
#' This function computes the posterior means for different groups within a dataset.
#' It expects the input data frame to contain a column named 'Group' that identifies
#' the groups, as well as columns 'ID', 'Output', and 'Input' -- and, for a
#' multi-dimensional Input, an additional 'Input_ID' column (see
#' \code{\link{normalize_input_cols}}): one row per (Group, ID, Sample,
#' Input_ID), with 'Output' repeated identically across the rows of one
#' observation. A legacy data frame with just a single 'Input' column (no
#' 'Input_ID') keeps working unchanged.
#'
#' The posterior distribution is given by:
#'
#' \eqn{p(\mathbf{\mu} \mid y_1, \dots, y_N, \Sigma_{\hat{\theta}}) = \mathcal{N}\left(\mathbf{\mu}; \ \dfrac{\lambda_0 \mu_0 + \sum_{n=1}^{N} y_n}{N + \lambda_0}, \dfrac{1}{N + \lambda_0} \Sigma_{\hat{\theta}}\right)}
#'
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @param data A data frame containing the data to be analyzed. Must include columns 'Group', 'ID', 'Output', and 'Input' (plus 'Input_ID' for a multi-dimensional Input).
#' @param kern A kernel object (from the keRnel package) used to compute pairwise covariances.
#' @param mu_0 Prior mean parameter.
#' @param lambda_0 Prior precision parameter.
#' @param obs_noise Fixed, known observation-noise variance to add back into
#'   the posterior covariance (see `get_sigmak()`'s `@details`). Needed
#'   whenever `kern`'s hyperparameters were fit with `optim_hp(..., prior_cov
#'   = )` treating that same value as a separate additive nugget rather than
#'   composing it into `kern` itself (e.g. via a `NoiseKernel()` term) --
#'   otherwise the reported credible interval only reflects uncertainty about
#'   whatever correlated structure the kernel captures BEYOND that fixed
#'   noise floor, which collapses toward zero whenever there is little such
#'   structure to find (the exact mechanism behind BayesOmics's severe
#'   under-coverage on ProteoBayes's own univariate/multivariate paper
#'   scenarios, `dev/benchmark_server/E11_proteobayes_paper_replay_server.R`).
#'   Defaults to 0 (previous behavior: `kern`'s own HPs assumed to already
#'   capture the full marginal variance).
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{kernels}}{A named list of kernel/correlation matrices, one per
#'       distinct set of Input positions found across groups. Each matrix has
#'       \code{dimnames} set to a string key per position (see
#'       \code{row_input_key()}), and is shared by reference across every group
#'       with that same set of positions (no duplication).}
#'     \item{\code{groups}}{A named list (one entry per group) with: \code{muk},
#'       a named vector of posterior means keyed by ID; \code{id_to_input}, an
#'       n_ids x D matrix (\code{rownames = ID}) giving that group's ID -> Input
#'       position mapping; \code{kernel_key}, which entry of \code{kernels} to
#'       use; \code{scale}, the divisor (\code{n_obs + lambda_0}) applied to
#'       that kernel matrix to get the posterior covariance; and
#'       \code{obs_noise} (this call's value, reused by `get_sigmak()`). The
#'       (internal) \code{get_sigmak()} helper reconstructs the actual
#'       (ID-aligned) posterior covariance matrix for a group.}
#'   }
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- multi_posterior_mean(data, kern)
#' posterior$groups[["1"]]$muk
multi_posterior_mean <- function(data, kern, mu_0 = 1, lambda_0 = 1, obs_noise = 0) {
  # FIXME: mu_0 defaults to 1, which is unusual for a Gaussian prior mean
  # (0 would be the conventional uninformative default). Verify this is
  # intentional (e.g. tied to a specific use case) and document the
  # rationale in @details, or change the default to 0.
  # === Initial checks ===
  required_cols <- c("Group", "ID", "Output", "Input")
  if (!all(required_cols %in% names(data))) {
    stop(paste0("The following columns are missing: ", paste(setdiff(required_cols, names(data)), collapse = ", ")))
  }
  data <- normalize_input_cols(data)
  if (!is.numeric(data$Input) || !is.numeric(data$Output) || !is.numeric(mu_0)) {
    stop("The 'Input' and 'Output' columns must be numeric, and mu_0 must be numeric.")
  }
  if (any(!is.finite(data$Input)) || any(!is.finite(data$Output))) {
    stop("'Input' and 'Output' columns must not contain NaN, Inf, or NA values.")
  }
  if (lambda_0 <= 0 || !is.numeric(lambda_0)) {
    stop("lambda_0 must be numeric and strictly positive.")
  }
  if (any(is.na(data$Group)) || any(is.na(data$ID))) {
    stop("The 'Group' and 'ID' columns must not contain NA values.")
  }
  if (!inherits(kern, "kernel")) {
    stop("The 'kern' argument must be a valid kernel object from the keRnel package.")
  }

  # === Convert Group to character if necessary ===
  if (!is.character(data$Group)) {
    data$Group <- as.character(data$Group)
  }

  # === Check for empty groups ===
  group_info <- data %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::summarise(
      n_id = dplyr::n_distinct(.data$ID),
      n_unique_inputs = dplyr::n_distinct(.data$Input),
      .groups = "drop"
    )

  empty_mask <- group_info$n_id == 0 | group_info$n_unique_inputs == 0
  if (any(empty_mask)) {
    empty_groups <- group_info$Group[empty_mask]
    stop(paste0("Empty groups: ", paste(empty_groups, collapse = ", ")))
  }

  # === Calculate df_mu ===
  # Output is repeated identically across the D rows (one per Input_ID) of a
  # single observation, so both the raw sum and the raw row count are
  # inflated D-fold; dividing by n_distinct(Input_ID) recovers the true
  # per-observation sum/count. D = 1 (legacy scalar Input, or already-1-row
  # Input_ID) makes this division a no-op, so the D = 1 case is numerically
  # identical to before this generalization.
  df_mu <- data %>%
    dplyr::group_by(.data$Group, .data$ID) %>%
    dplyr::summarise(
      n_dim = dplyr::n_distinct(.data$Input_ID),
      sum_outputs = sum(.data$Output, na.rm = TRUE) / n_dim,
      lengths = dplyr::n() / n_dim,
      .groups = "drop"
    )
  if (any(abs(df_mu$lengths - round(df_mu$lengths)) > 1e-8)) {
    stop(paste0(
      "Inconsistent number of Input_ID rows across replicate observations ",
      "within at least one (Group, ID): each replicate must have exactly one ",
      "row per Input_ID dimension."
    ))
  }
  df_mu$lengths <- round(df_mu$lengths)
  df_mu$muk_vector <- (lambda_0 * mu_0 + df_mu$sum_outputs) / (df_mu$lengths + lambda_0)

  # === Per-(Group, ID) Input position, deduped across replicate rows (Input
  # is unique per ID within a group, checked below) so muk and the kernel
  # matrix can be re-aligned by ID rather than by position. Deduplicating
  # first means input_matrix_by_id() only sees each (ID, Input_ID) pair once
  # -- a genuine contradiction (the same ID/Input_ID recorded with two
  # different Input values) still surfaces as its "duplicated Input_ID" error.
  id_input <- dplyr::distinct(data, .data$Group, .data$ID, .data$Input_ID, .data$Input)
  id_input_by_group <- split(id_input, id_input$Group)

  # === Pre-compute subsets for group_data ===
  df_mu_by_group <- split(df_mu, df_mu$Group)

  # === Build group_data ===
  groups <- unique(data$Group)
  group_data <- dplyr::tibble(
    group = groups,
    muk = purrr::map(groups, function(g) {
      subset <- df_mu_by_group[[g]]
      stats::setNames(subset$muk_vector, subset$ID)
    }),
    # n_ids x D matrix, rownames = ID (D = 1 for the legacy scalar-Input
    # case). input_matrix_by_id() itself errors clearly on a duplicated or
    # ragged Input_ID design within this group.
    id_to_input = purrr::map(groups, function(g) {
      input_matrix_by_id(id_input_by_group[[g]], key_cols = "ID")$matrix
    }),
    n_ids = purrr::map_int(groups, function(g) nrow(df_mu_by_group[[g]])),
    # n_obs = observations per ID (= nb_sample); same N used in muk numerator.
    # All IDs within a group must share this count: the group-level posterior
    # scale (n_obs + lambda_0) is a single divisor applied to every ID, so a
    # silently-averaged n_obs would bias the covariance for groups with
    # unequal per-ID observation counts.
    n_obs = purrr::map_int(groups, function(g) {
      lens <- df_mu_by_group[[g]]$lengths
      if (length(unique(lens)) > 1) {
        stop(paste0(
          "Group '", g, "' has IDs with different numbers of observations (",
          paste(sort(unique(lens)), collapse = ", "), "); multi_posterior_mean() ",
          "requires every ID within a group to have the same number of observations."
        ))
      }
      as.integer(lens[1])
    })
  )

  # === Check muk / sigmak dimension consistency: each ID must map to a
  # DISTINCT Input position within its group (no two IDs sharing the same
  # D-dimensional coordinates) ===
  muk_lengths <- vapply(group_data$muk, length, integer(1))
  distinct_positions <- vapply(group_data$id_to_input, function(m) {
    length(unique(row_input_key(m)))
  }, integer(1))
  mismatch <- which(muk_lengths != distinct_positions)
  if (length(mismatch) > 0) {
    bad <- group_data$group[mismatch]
    stop(paste0(
      "In group(s) [", paste(bad, collapse = ", "), "]: ",
      "number of unique IDs (", paste(muk_lengths[mismatch], collapse = ", "), ") ",
      "does not match number of unique Input positions (",
      paste(distinct_positions[mismatch], collapse = ", "), "). ",
      "Each ID must map to a distinct Input position within its group."
    ))
  }

  # === Compute kernel matrices ===
  # Kernel matrices are cached by their (sorted) set of Input positions: two
  # groups with the same set of positions share the exact same matrix object
  # (no duplication), regardless of which IDs use those positions in each
  # group. The matrix's dimnames are the positions themselves (via
  # row_input_key(), not IDs), so it can be safely re-indexed by ID per group
  # via get_sigmak() without ever mixing up which row/col belongs to which ID.
  cache <- new.env(hash = TRUE, parent = emptyenv())

  groups_list <- lapply(seq_len(nrow(group_data)), function(i) {
    id_to_input <- group_data$id_to_input[[i]]
    row_keys    <- row_input_key(id_to_input)
    ord         <- order(row_keys)
    keep        <- !duplicated(row_keys[ord])
    unique_rows <- id_to_input[ord[keep], , drop = FALSE]
    unique_keys <- row_keys[ord[keep]]
    # A canonical, deterministic cache key built from the exact same strings
    # used for the matrix's own dimnames just below -- unlike the previous
    # toString(sort(vec)) key, the two can never disagree.
    vec_hash <- paste(unique_keys, collapse = ";;")

    if (!exists(vec_hash, envir = cache, inherits = FALSE)) {
      kern_mat <- keRnel::evaluate(kern, unique_rows, unique_rows)
      dimnames(kern_mat) <- list(unique_keys, unique_keys)
      assign(vec_hash, kern_mat, envir = cache)
    }
    list(
      muk         = group_data$muk[[i]],
      id_to_input = id_to_input,
      kernel_key  = vec_hash,
      scale       = group_data$n_obs[i] + lambda_0,
      obs_noise   = obs_noise
    )
  })
  names(groups_list) <- group_data$group

  structure(
    list(kernels = as.list(cache), groups = groups_list),
    class = "bayesomics_posterior"
  )
}

#' Print a BayesOmics Posterior Object
#'
#' @description
#' Pretty-prints the result of \code{\link{multi_posterior_mean}}: for each
#' group, the posterior mean vector (\code{muk}) and the reconstructed
#' (ID-aligned) posterior covariance matrix, obtained via the internal
#' \code{get_sigmak()} helper.
#'
#' @param x A list returned by \code{\link{multi_posterior_mean}}.
#' @param digits Number of significant digits used when rounding the
#'   displayed mean vector and covariance matrix. Defaults to \code{3}.
#' @param ... Unused, included for S3 consistency.
#'
#' @return \code{x}, invisibly.
#' @export
print.bayesomics_posterior <- function(x, digits = 3, ...) {
  cat(sprintf(
    "<BayesOmics posterior> %d group(s), %d cached kernel matrix/matrices\n",
    length(x$groups), length(x$kernels)
  ))
  for (g in names(x$groups)) {
    entry <- x$groups[[g]]
    sigmak <- get_sigmak(entry, x$kernels)
    cat(sprintf("\n-- Group %s (%d IDs) --\n", g, length(entry$muk)))
    cat("Posterior mean:\n")
    print(round(entry$muk, digits))
    cat("Posterior covariance:\n")
    print(round(sigmak, digits))
  }
  invisible(x)
}

#' @title Sample from a Normal multivariate distribution
#'
#' @description
#' Sample n elements from the posterior distribution of each group, and
#' reshape the result directly into the long-format data frame expected by
#' \code{\link{plot_distrib}}.
#'
#' @importFrom mvtnorm rmvnorm
#' @param results A list returned by \code{\link{multi_posterior_mean}}, with
#'   elements \code{kernels} and \code{groups}.
#' @param n A number indicating the number of samples
#'
#' @return A data frame with columns \code{ID}, \code{Group}, and
#'   \code{Sample} (one row per draw, per ID, per group).
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 100)
#' head(samples)
sample_posterior <- function(results, n) {
  if (!is.list(results) || !all(c("kernels", "groups") %in% names(results)) || length(results$groups) == 0) {
    stop("'results' must be the list returned by multi_posterior_mean() (with 'kernels' and 'groups').")
  }
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  for (group in names(results$groups)) {
    g <- results$groups[[group]]
    if (!all(c("muk", "id_to_input", "kernel_key", "scale") %in% names(g))) {
      stop(paste("Group", group, "is missing required elements (muk, id_to_input, kernel_key, scale)."))
    }
  }

  samples_list <- lapply(results$groups, function(g) {
    sigmak <- get_sigmak(g, results$kernels)
    mat <- mvtnorm::rmvnorm(n = n, mean = g$muk, sigma = sigmak)
    colnames(mat) <- names(g$muk)
    mat
  })
  names(samples_list) <- names(results$groups)

  # Melt each group's (n_draws x n_ids) matrix to long format in one shot
  # (rep()/as.vector() are O(n_draws * n_ids)) instead of rbind-ing one tiny
  # data frame per ID, which is O(n_ids) rbind calls each re-copying an
  # ever-growing data frame (O(n_ids^2) overall).
  do.call(rbind, lapply(names(samples_list), function(g) {
    mat <- samples_list[[g]]
    data.frame(
      ID     = rep(colnames(mat), each = n),
      Group  = g,
      Sample = as.vector(mat),
      stringsAsFactors = FALSE
    )
  }))
}
