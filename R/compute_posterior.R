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
#'   n_ids x D matrix with `rownames = ID` (the shape `posterior_mean()`
#'   itself now always produces, D >= 1). Both are turned into one
#'   `row_input_key()` string per ID, so re-indexing the shared kernel matrix
#'   works identically either way.
#' @details `group_entry$obs_noise` (added by `posterior_mean()`, default
#'   0) is a fixed, known observation-noise variance to add back to the
#'   kernel-derived covariance before dividing by `scale`. It exists because
#'   HP-fitting workflows built around `fit_kernel(..., prior_cov = )` (see
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
#' @param data A data frame containing the data to be analyzed. Must include
#'   columns 'Group' and 'Output'. 'ID' and 'Input' are optional, but must be
#'   supplied together or omitted together: if BOTH are missing,
#'   \code{posterior_mean()} runs in **univariate mode** (a warning is issued)
#'   -- each group is treated as a single feature (a dummy constant
#'   'ID'/'Input' is added internally), with no cross-feature correlation
#'   structure, exactly the \code{nb_id = 1} case validated in
#'   \code{dev/univariate/NOTES_univariate.md}. Supplying only one of the two
#'   is ambiguous and errors instead. Supply both 'ID' and 'Input' (plus
#'   'Input_ID' for a multi-dimensional Input) for the general multi-feature
#'   case.
#' @param kern A kernel object (from the keRnel package) used to compute
#'   pairwise covariances, or a named list of kernel objects (one entry per
#'   group in \code{data}, e.g. from \code{\link{fit_kernel}(..., group_col =
#'   "Group", pooled = FALSE)}) to give each group its own, independently
#'   fitted kernel -- every group then gets its own \code{kernel_key}, so
#'   metrics with \code{requires_shared_kernel() == TRUE} (e.g.
#'   \code{\link{ovl_metric}}) cannot compare them (same caveat as
#'   \code{pooled = FALSE} below). Defaults to \code{NULL}, in which case a diagonal
#'   (\code{keRnel::white_noise_kernel()}) kernel is built automatically from a
#'   closed-form residual-variance estimate -- the exact REML minimizer
#'   \code{fit_kernel(..., group_col = "Group")} would converge to numerically,
#'   computed directly instead (see \code{resolve_closed_form_kernel()} in
#'   \code{R/optim_kernel.R} and \code{dev/univariate/NOTES_univariate.md} for
#'   the derivation). This is the natural default in univariate mode, but also
#'   works with a real 'ID'/'Input' design (every feature is then treated as
#'   independent -- no spatial structure -- unlike a real kernel fit via
#'   \code{fit_kernel()}, which must be supplied explicitly via \code{kern} for
#'   that).
#' @param pooled Only used when \code{kern = NULL}. If \code{TRUE} (default), a
#'   single noise variance is estimated and shared by every group (more
#'   residual degrees of freedom, assumes homogeneous noise across groups). If
#'   \code{FALSE}, each group gets its own independently-estimated variance
#'   (heteroscedastic); groups then never share a \code{kernel_key}, so
#'   metrics with \code{requires_shared_kernel() == TRUE} (e.g.
#'   \code{\link{ovl_metric}}) cannot compare them. Ignored when \code{kern} is
#'   supplied directly.
#' @param df_warn Only used when \code{kern = NULL}. A \code{warning()} is
#'   issued whenever the residual degrees of freedom backing a closed-form
#'   variance estimate fall below this (the estimate runs but is noisy); an
#'   \code{error()} is always raised at 0 or fewer degrees of freedom
#'   (the estimate is undefined). Defaults to \code{8} (see
#'   \code{dev/optim_exploration/14_group_and_id_demean/README.md}, point 3:
#'   empirically, the fitted variance's coefficient of variation drops below
#'   0.5 around 8 residual degrees of freedom).
#' @param mu_0 Prior mean parameter.
#' @param lambda_0 Prior precision parameter.
#' @param obs_noise Fixed, known observation-noise variance to add back into
#'   the posterior covariance (see `get_sigmak()`'s `@details`). Needed
#'   whenever `kern`'s hyperparameters were fit with `fit_kernel(..., prior_cov
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
#' posterior <- posterior_mean(data, kern)
#' posterior$groups[["1"]]$muk
#'
#' # Univariate mode (single CpG/feature): no 'ID'/'Input' columns, no kernel
#' # object needed -- a warning is issued and the noise variance is estimated
#' # in closed form (pooled across every group by default):
#' cpg <- data.frame(
#'   Group  = rep(c("A", "B"), each = 5),
#'   Sample = rep(1:5, 2),
#'   Output = c(rnorm(5, 0, 1), rnorm(5, 3, 1))
#' )
#' uni_posterior <- posterior_mean(cpg)
posterior_mean <- function(data, kern = NULL, mu_0 = 1, lambda_0 = 1, obs_noise = 0,
                            pooled = TRUE, df_warn = 8) {
  # FIXME: mu_0 defaults to 1, which is unusual for a Gaussian prior mean
  # (0 would be the conventional uninformative default). Verify this is
  # intentional (e.g. tied to a specific use case) and document the
  # rationale in @details, or change the default to 0.
  # === Initial checks ===
  if (!all(c("Group", "Output") %in% names(data))) {
    stop(paste0(
      "The following columns are missing: ",
      paste(setdiff(c("Group", "Output"), names(data)), collapse = ", ")
    ))
  }

  # === Univariate mode: no 'ID'/'Input' -> one feature per group ===
  # Shared with fit_kernel()'s kern = NULL branch (inject_univariate_dummy_cols(),
  # see R/optim_kernel.R) so both switch into univariate mode under the exact
  # same condition, with the exact same warning/error text.
  data <- inject_univariate_dummy_cols(data, "posterior_mean")

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
  # A named list of kernel objects (one per group, e.g. from
  # fit_kernel(..., group_col = "Group", pooled = FALSE)) is accepted here
  # too, in addition to a single shared kernel object -- mirroring the shape
  # resolve_closed_form_kernel(pooled = FALSE) already produces internally
  # for the kern = NULL path. kern_is_shared below then dispatches per group.
  kern_is_list <- is.list(kern) && !inherits(kern, "kernel") && length(kern) > 0 &&
    all(vapply(kern, inherits, logical(1), "kernel"))
  if (!is.null(kern) && !inherits(kern, "kernel") && !kern_is_list) {
    stop("The 'kern' argument must be a valid kernel object from the keRnel package, ",
         "or a named list of kernel objects (one per group).")
  }
  if (kern_is_list) {
    missing_groups <- setdiff(unique(as.character(data$Group)), names(kern))
    if (length(missing_groups) > 0) {
      stop("'kern' is a named list but is missing an entry for group(s): ",
           paste(missing_groups, collapse = ", "))
    }
  }

  # === Convert Group to character if necessary ===
  if (!is.character(data$Group)) {
    data$Group <- as.character(data$Group)
  }

  # === Closed-form kernel when none is supplied ===
  # Must run after the Group -> character conversion above so the per-group
  # kernel list resolve_closed_form_kernel() returns (pooled = FALSE) is keyed
  # by the exact same group labels group_data$group is built from below.
  if (is.null(kern)) {
    kern <- resolve_closed_form_kernel(data, pooled = pooled, df_warn = df_warn)
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
  #
  # This n_dim division cannot, by itself, distinguish clean data from every
  # row being uniformly duplicated (e.g. 2x): n_distinct(Input_ID) and n()
  # both double, so their ratio (and the later "is lengths near-integer?"
  # check below) is unchanged either way -- a duplicated dataset silently
  # returns a wrong muk with no error (found by dev/multi_agent_audit,
  # reproduced: muk = 12 instead of 10 under a uniform 2x row duplication).
  # When a 'Sample' column identifies individual replicates, guard against
  # this directly: within a single (Group, ID, Sample) replicate, each
  # Input_ID must appear exactly once (n_distinct(Input_ID) == n()) -- a
  # uniform duplication violates this even though the (Group, ID)-level
  # ratio stays a clean integer. No 'Sample' column means there is no way to
  # identify which rows form one replicate, so this check is skipped (same
  # convention as sum_logGaussian()'s is_replicated detection elsewhere).
  if ("Sample" %in% names(data)) {
    replicate_check <- data %>%
      dplyr::group_by(.data$Group, .data$ID, .data$Sample) %>%
      dplyr::summarise(
        n_rows            = dplyr::n(),
        n_distinct_input  = dplyr::n_distinct(.data$Input_ID),
        .groups = "drop"
      )
    if (any(replicate_check$n_rows != replicate_check$n_distinct_input)) {
      stop(paste0(
        "Duplicated (or missing) Input_ID row(s) within at least one ",
        "(Group, ID, Sample) replicate: each replicate must have exactly ",
        "one row per Input_ID dimension, with no repeats."
      ))
    }
  }
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
          paste(sort(unique(lens)), collapse = ", "), "); posterior_mean() ",
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
  # `kern` is either a single kernel object shared by every group (the normal
  # case: user-supplied, or the closed-form fit when pooled = TRUE), or a
  # named list of per-group kernel objects (closed-form fit, pooled = FALSE,
  # from resolve_closed_form_kernel() in R/optim_kernel.R) -- resolved once
  # here rather than per iteration below.
  kern_is_shared <- inherits(kern, "kernel")

  groups_list <- lapply(seq_len(nrow(group_data)), function(i) {
    id_to_input <- group_data$id_to_input[[i]]
    row_keys    <- row_input_key(id_to_input)
    ord         <- order(row_keys)
    keep        <- !duplicated(row_keys[ord])
    unique_rows <- id_to_input[ord[keep], , drop = FALSE]
    unique_keys <- row_keys[ord[keep]]

    group_name <- group_data$group[i]
    kern_i     <- if (kern_is_shared) kern else kern[[group_name]]

    # A canonical, deterministic cache key built from the exact same strings
    # used for the matrix's own dimnames just below -- unlike the previous
    # toString(sort(vec)) key, the two can never disagree. Hashed (not used
    # raw) because assign()/exists() use this as an R *symbol* name, which
    # has a hard 10,000-byte limit enforced by R itself -- the raw
    # collapsed-keys string blows past that well before nb_id reaches a few
    # thousand (found running E14_scale_server.R at nb_id=2000: "variable
    # names are limited to 10000 bytes"). rlang::hash() (already a
    # dependency) gives a fixed-length 32-char digest regardless of input
    # size, deterministic within a session -- exactly what this
    # single-call-local cache needs (no cross-session persistence implied
    # or required, see this function's own cache-scope documentation).
    # When kern varies per group (pooled = FALSE), the same Input positions
    # must NOT be shared across groups' cached matrices -- they come from
    # different fitted variances -- so the group name is folded into the key,
    # guaranteeing every group gets its own cache entry (and therefore its own
    # kernel_key, which is exactly what requires_shared_kernel() metrics like
    # ovl_metric() need to correctly refuse to compare them, see
    # compute_group_diff()).
    cache_key <- if (kern_is_shared) {
      paste(unique_keys, collapse = ";;")
    } else {
      paste0(group_name, "", paste(unique_keys, collapse = ";;"))
    }
    vec_hash <- rlang::hash(cache_key)

    if (!exists(vec_hash, envir = cache, inherits = FALSE)) {
      kern_mat <- keRnel::evaluate(kern_i, unique_rows, unique_rows)
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
#' Pretty-prints the result of \code{\link{posterior_mean}}: for each
#' group, the posterior mean vector (\code{muk}) and the reconstructed
#' (ID-aligned) posterior covariance matrix, obtained via the internal
#' \code{get_sigmak()} helper.
#'
#' @param x A list returned by \code{\link{posterior_mean}}.
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
#' @param results A list returned by \code{\link{posterior_mean}}, with
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
#' posterior <- posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 100)
#' head(samples)
sample_posterior <- function(results, n) {
  if (!is.list(results) || !all(c("kernels", "groups") %in% names(results)) || length(results$groups) == 0) {
    stop("'results' must be the list returned by posterior_mean() (with 'kernels' and 'groups').")
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

## ===========================================================================
## Block wrapper -- split one BayesOmics problem into several independent,
## smaller ones (one per block, or per (Group, block)), each solved by the
## UNCHANGED posterior_mean()/fit_kernel() pipeline above, then consulted as if
## the result were a single block-diagonal covariance. See
## dev/block_wrapper_demo/README.md for the full design discussion; this is
## the "Option B" integration into R/ decided there.
## ===========================================================================

#' Resolve a Partition Strategy into a Block Assignment
#'
#' @description
#' Generic dispatching on the S3 class of \code{strategy} (a
#' \code{partition_strategy} object, e.g. from \code{\link{partition_by_id}}).
#' Not normally called directly -- \code{\link{fit_block_posterior}} calls it
#' internally -- but exported so custom partition strategies can be tested in
#' isolation, or implemented by defining a new \code{resolve_partition.*}
#' method.
#'
#' @param strategy A \code{partition_strategy} object, or \code{NULL} (one
#'   single dense block, see \code{\link{partition_dense}}).
#' @param ids A character vector of all IDs to assign.
#' @param db The full data frame being partitioned (only needed by strategies
#'   that inspect \code{Input}/annotation columns, e.g.
#'   \code{\link{partition_by_range}}).
#' @return A named character vector, \code{block_of_id[id]}, giving each ID's
#'   block label.
#' @export
resolve_partition <- function(strategy, ids, db = NULL) UseMethod("resolve_partition")

#' @export
resolve_partition.default <- function(strategy, ids, db = NULL) {
  if (is.null(strategy)) return(stats::setNames(rep("ALL", length(ids)), ids))
  if (is.character(strategy) && !is.null(names(strategy))) {
    missing_ids <- setdiff(ids, names(strategy))
    if (length(missing_ids) > 0) {
      stop("resolve_partition(): missing block assignment for ID(s): ",
           paste(missing_ids, collapse = ", "))
    }
    return(as.character(strategy[ids]))
  }
  stop("resolve_partition(): unrecognized partition specification (class: ",
       paste(class(strategy), collapse = "/"), ").")
}

## Natural sort of block labels: compares the final numeric suffix as a
## NUMBER rather than character-by-character -- plain sort() would place
## "blk10" before "blk2". Never affects correctness (every consumer below is
## invariant to block order) but avoids a trap for future code that displays
## or groups adjacent blocks in this order.
#' @noRd
natural_sort <- function(x) {
  suffix_num <- suppressWarnings(as.numeric(sub("^.*?([0-9]+)$", "\\1", x)))
  prefix <- sub("[0-9]+$", "", x)
  x[order(prefix, ifelse(is.na(suffix_num), Inf, suffix_num), x)]
}

#' Explicit Manual ID-to-Block Assignment
#'
#' @param assignment A named character (or coercible) vector, \code{assignment[id]}
#'   giving that ID's block label. Must cover every ID being partitioned.
#' @return A \code{partition_strategy} object usable as
#'   \code{\link{fit_block_posterior}}'s \code{partition} argument.
#' @export
partition_by_id <- function(assignment) {
  structure(list(assignment = assignment), class = c("by_id_partition", "partition_strategy"))
}

#' @export
resolve_partition.by_id_partition <- function(strategy, ids, db = NULL) {
  assignment <- strategy$assignment
  missing_ids <- setdiff(ids, names(assignment))
  if (length(missing_ids) > 0) {
    stop("partition_by_id(): missing assignment for ID(s): ", paste(missing_ids, collapse = ", "))
  }
  if (any(is.na(assignment[ids]))) {
    stop("partition_by_id(): NA block label(s) not allowed.")
  }
  stats::setNames(as.character(assignment[ids]), ids)
}

#' Partition IDs by Slicing Input into Ranges
#'
#' @description
#' Cuts \code{Input} values into intervals (via \code{\link[base]{cut}}); IDs
#' whose Input falls in the same interval go into the same block.
#'
#' @param breaks Explicit interval breakpoints, passed to \code{cut()}. If
#'   \code{NULL} (default), \code{n_bins} equal-width bins spanning the
#'   observed range of \code{Input} are used instead.
#' @param n_bins Number of equal-width bins when \code{breaks} is \code{NULL}.
#'   Defaults to \code{3}.
#' @param input_col Name of the Input column to slice. Defaults to
#'   \code{"Input"}.
#' @param input_id When \code{db} has a multi-dimensional Input (an
#'   \code{Input_ID} column with more than one distinct value), which axis to
#'   slice on. Required in that case -- without it, guessing an axis could
#'   silently cut on the wrong one, so \code{\link{resolve_partition}} errors
#'   instead.
#' @return A \code{partition_strategy} object usable as
#'   \code{\link{fit_block_posterior}}'s \code{partition} argument.
#' @export
partition_by_range <- function(breaks = NULL, n_bins = NULL, input_col = "Input", input_id = NULL) {
  structure(list(breaks = breaks, n_bins = n_bins, input_col = input_col, input_id = input_id),
            class = c("by_range_partition", "partition_strategy"))
}

#' @export
resolve_partition.by_range_partition <- function(strategy, ids, db) {
  if (is.null(db)) stop("partition_by_range(): requires 'db' (needs the Input column).")
  db_axis <- db
  if (!is.null(strategy$input_id)) {
    if (!"Input_ID" %in% names(db)) {
      stop("partition_by_range(): 'input_id' was specified but 'db' has no 'Input_ID' column.")
    }
    db_axis <- db[db$Input_ID == strategy$input_id, , drop = FALSE]
  } else if ("Input_ID" %in% names(db) && length(unique(db$Input_ID)) > 1) {
    stop("partition_by_range(): 'db' has more than one Input_ID axis -- specify which one via ",
         "'input_id' (e.g. partition_by_range(n_bins = 3, input_id = 1)); silently guessing an ",
         "axis would risk cutting on the wrong one.")
  }
  pos <- vapply(ids, function(id) db_axis[[strategy$input_col]][match(id, db_axis$ID)], numeric(1))
  breaks <- strategy$breaks
  if (is.null(breaks)) {
    n_bins <- if (is.null(strategy$n_bins)) 3 else strategy$n_bins
    breaks <- seq(min(pos), max(pos), length.out = n_bins + 1)
  }
  bin <- cut(pos, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  if (any(is.na(bin))) stop("partition_by_range(): some Input value(s) fall outside 'breaks'.")
  stats::setNames(paste0("range_", bin), ids)
}

#' Partition IDs from an External Annotation
#'
#' @param mapping Either a data frame with a \code{by} column and a
#'   \code{Block} column, or a named vector (ID -> block label).
#' @param by Name of the ID column in \code{mapping} when it is a data frame.
#'   Defaults to \code{"ID"}.
#' @return A \code{partition_strategy} object usable as
#'   \code{\link{fit_block_posterior}}'s \code{partition} argument.
#' @export
partition_from_annotation <- function(mapping, by = "ID") {
  structure(list(mapping = mapping, by = by), class = c("annotation_partition", "partition_strategy"))
}

#' @export
resolve_partition.annotation_partition <- function(strategy, ids, db = NULL) {
  mapping <- strategy$mapping
  if (is.data.frame(mapping)) {
    idx <- match(ids, mapping[[strategy$by]])
    if (any(is.na(idx))) {
      stop("partition_from_annotation(): no annotation for ID(s): ",
           paste(ids[is.na(idx)], collapse = ", "))
    }
    return(stats::setNames(as.character(mapping$Block[idx]), ids))
  }
  missing_ids <- setdiff(ids, names(mapping))
  if (length(missing_ids) > 0) {
    stop("partition_from_annotation(): missing mapping for ID(s): ", paste(missing_ids, collapse = ", "))
  }
  stats::setNames(as.character(mapping[ids]), ids)
}

#' Partition IDs by Density-Based Clustering (Reserved, Not Implemented)
#'
#' @description
#' Placeholder for an automatic, density-based partition strategy.
#' Constructing this object is allowed (so it can be passed around and
#' inspected), but \code{\link{resolve_partition}} always errors on it --
#' auto-detection is not implemented yet.
#'
#' @param eps Reserved.
#' @param minPts Reserved.
#' @param ... Reserved.
#' @return A \code{partition_strategy} object that always errors when resolved.
#' @export
partition_dbscan <- function(eps, minPts, ...) {
  structure(list(eps = eps, minPts = minPts, ...), class = c("dbscan_partition", "partition_strategy"))
}

#' @export
resolve_partition.dbscan_partition <- function(strategy, ids, db = NULL) {
  stop("partition_dbscan(): auto-detection is not implemented yet.")
}

#' Dense Partition (a Single Block)
#'
#' @description
#' The "no partitioning at all" extreme: every ID goes into one block, so
#' \code{\link{fit_block_posterior}} degrades exactly to a single call to
#' \code{\link{posterior_mean}} on the whole dataset. Equivalent to passing
#' \code{partition = NULL}.
#'
#' @return \code{NULL} (resolved by \code{\link{resolve_partition}()}'s default method).
#' @export
partition_dense <- function() NULL

#' Diagonal Partition (One Block per ID)
#'
#' @description
#' The opposite extreme from \code{\link{partition_dense}}: every ID is its
#' own singleton block, so every block is a \code{nb_id = 1} problem.
#'
#' @param ids Character vector of IDs to assign one-per-block.
#' @return A \code{partition_strategy} object usable as
#'   \code{\link{fit_block_posterior}}'s \code{partition} argument.
#' @export
partition_diagonal <- function(ids) partition_by_id(stats::setNames(ids, ids))

#' Fit Independent Posteriors per Block
#'
#' @description
#' Splits \code{data} according to \code{partition}, then calls the unchanged
#' \code{\link{posterior_mean}}/\code{\link{fit_kernel}} pipeline independently
#' per sub-problem (per block, or per \code{(Group, block)} when
#' \code{pooled = FALSE}), collecting the results into a list consultable via
#' \code{\link{get_sigmak_block}}/\code{\link{get_muk_block}}. No result
#' fusion happens here -- each sub-problem stays a genuine
#' \code{bayesomics_posterior} object, so no numeric primitive is
#' reimplemented.
#'
#' Passing \code{partition = \link{partition_dense}()} (or, equivalently,
#' \code{NULL}) collapses this to a single sub-problem covering the whole
#' dataset -- with the same \code{pooled} default as \code{posterior_mean()}
#' (\code{TRUE}), this reproduces \code{posterior_mean(data)} exactly, so
#' calling \code{fit_block_posterior()} with nothing else specified stays
#' agnostic to the non-partitioned baseline.
#'
#' @param data A data frame, same requirements as \code{\link{posterior_mean}}.
#' @param partition A \code{partition_strategy} object (e.g.
#'   \code{\link{partition_by_id}}, \code{\link{partition_by_range}},
#'   \code{\link{partition_from_annotation}}, \code{\link{partition_diagonal}}),
#'   or \code{NULL} / \code{\link{partition_dense}()} for a single block.
#' @param kern \code{NULL} (closed form, as in \code{posterior_mean()}), or an
#'   UNFITTED kernel object (template) -- in that case \code{fit_kernel()} is
#'   called once per sub-problem before building its posterior (mirroring how
#'   \code{posterior_mean()} never fits \code{kern} itself when it is
#'   supplied, see \code{R/optim_kernel.R}).
#' @param pooled A single parameter threaded to both branches: when
#'   \code{kern = NULL}, passed straight through to
#'   \code{posterior_mean(kern = NULL, pooled = pooled)}; when \code{kern} is
#'   an unfitted template, \code{TRUE} fits one \code{fit_kernel(group_col =
#'   "Group")} shared across every Group within a block, \code{FALSE} fits an
#'   independent \code{fit_kernel()} per \code{(Group, block)}. Defaults to
#'   \code{TRUE}, matching \code{posterior_mean()}'s own default. See also
#'   \code{\link{fit_kernel}}'s own \code{pooled} argument for the same
#'   pooled/non-pooled choice without the block-partitioning machinery.
#' @param mu_0,lambda_0,obs_noise,df_warn Forwarded to \code{posterior_mean()}
#'   on every sub-problem; see its documentation.
#' @param prior_mean,prior_cov,pen_diag Forwarded to \code{fit_kernel()} when
#'   \code{kern} is an unfitted template; see its documentation.
#' @return A \code{block_posterior_fit} object: a list with
#'   \code{block_results} (named list of per-sub-problem
#'   \code{bayesomics_posterior} objects), \code{block_index} (named list per
#'   Group, mapping block label -> \code{block_results} key), \code{block_of_id}
#'   (the resolved ID -> block assignment), \code{blocks} (block labels, in
#'   natural-sort order), and \code{pooled}.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 12, nb_group = 2, nb_sample = 5)
#' ids <- unique(data$ID)
#' partition <- partition_by_id(stats::setNames(rep(c("A", "B"), length.out = length(ids)), ids))
#' fit <- fit_block_posterior(data, partition, kern = NULL, pooled = TRUE)
#' fit$blocks
fit_block_posterior <- function(data, partition, kern = NULL, pooled = TRUE,
                                 mu_0 = 1, lambda_0 = 1, obs_noise = 0, df_warn = 8,
                                 prior_mean = 0, prior_cov = 1, pen_diag = 1e-6) {
  ## posterior_mean() forces Group to character but never ID -- setdiff()/
  ## match()/paste() in resolve_partition() implicitly assume textual IDs; an
  ## ID stored as a factor/integer could otherwise behave unpredictably
  ## (silent coercions). Coerced HERE, before anything else, so
  ## resolve_partition() and the split below all see the same textual IDs.
  data$ID <- as.character(data$ID)
  ## Group must be coerced HERE too, for a different reason: posterior_mean()
  ## coerces it internally, so every sub-problem's res$groups is always keyed
  ## by CHARACTER group labels -- but the bookkeeping below (block_index,
  ## `for (g in groups_all) ... res$groups[[g]]`) reads data$Group directly.
  ## With an un-coerced integer/numeric Group (simu_db()'s default), `[[`
  ## silently switches to POSITIONAL indexing instead of name lookup, which
  ## only "works" by accident when group values happen to equal 1, 2, ... in
  ## their natural order, and throws "subscript out of bounds" as soon as a
  ## group's sub-problem has fewer entries than that positional index (caught
  ## by tests/testthat/test-compute_posterior.R's pooled = FALSE regression
  ## test, using simu_db()'s integer Group as-is).
  data$Group <- as.character(data$Group)
  ids_all <- unique(data$ID)
  block_of_id <- resolve_partition(partition, ids_all, data)
  blocks <- natural_sort(unique(block_of_id))
  groups_all <- unique(data$Group)

  block_results <- list()
  block_index <- stats::setNames(vector("list", length(groups_all)), groups_all)
  for (g in groups_all) block_index[[g]] <- stats::setNames(rep(NA_character_, length(blocks)), blocks)

  fit_one <- function(sub, group_col = NULL) {
    if (is.null(kern)) {
      posterior_mean(sub, kern = NULL, mu_0 = mu_0, lambda_0 = lambda_0,
                      obs_noise = obs_noise, pooled = pooled, df_warn = df_warn)
    } else {
      opt <- fit_kernel(kern, sub, prior_mean = if (is.null(group_col)) prior_mean else NULL,
                        prior_cov = prior_cov, pen_diag = pen_diag, verbose = TRUE,
                        group_col = group_col)
      res <- posterior_mean(sub, kern = opt$kern, mu_0 = mu_0, lambda_0 = lambda_0, obs_noise = obs_noise)
      ## posterior_mean()$kernels stores the CACHED kernel matrix, not the
      ## fitted kernel object -- kept as a separate attribute so callers can
      ## inspect the per-block fitted hyperparameters
      ## (keRnel::get_trainable_params(attr(res, "fitted_kern"))).
      attr(res, "fitted_kern") <- opt$kern
      res
    }
  }

  for (b in blocks) {
    ids_b <- names(block_of_id)[block_of_id == b]
    if (pooled) {
      sub <- data[data$ID %in% ids_b, , drop = FALSE]
      key <- b
      block_results[[key]] <- fit_one(sub, group_col = if (is.null(kern)) NULL else "Group")
      for (g in groups_all) if (g %in% unique(sub$Group)) block_index[[g]][b] <- key
    } else {
      for (g in groups_all) {
        sub <- data[data$ID %in% ids_b & data$Group == g, , drop = FALSE]
        if (nrow(sub) == 0) next
        key <- paste(b, g, sep = "__")
        block_results[[key]] <- fit_one(sub, group_col = NULL)
        block_index[[g]][b] <- key
      }
    }
  }

  structure(
    list(block_results = block_results, block_index = block_index,
         block_of_id = block_of_id, blocks = blocks, pooled = pooled),
    class = "block_posterior_fit"
  )
}

#' Reconstruct Per-Block Posterior Covariances for a Group
#'
#' @description
#' Calls the real, unchanged (internal) \code{get_sigmak()} on each of a group's
#' sub-problem results, collecting one (ID-aligned) covariance matrix per
#' block.
#'
#' @param fit A \code{block_posterior_fit} object from
#'   \code{\link{fit_block_posterior}}.
#' @param group Which group's blocks to reconstruct.
#' @return A named list of covariance matrices, one per block that group
#'   participates in.
#' @export
get_sigmak_block <- function(fit, group) {
  Sigma_blocks <- list()
  for (b in fit$blocks) {
    key <- fit$block_index[[group]][b]
    if (is.na(key)) next
    res <- fit$block_results[[key]]
    if (!(group %in% names(res$groups))) next
    Sigma_blocks[[b]] <- get_sigmak(res$groups[[group]], res$kernels)
  }
  Sigma_blocks
}

#' Concatenate Per-Block Posterior Means for a Group
#'
#' @param fit A \code{block_posterior_fit} object from
#'   \code{\link{fit_block_posterior}}.
#' @param group Which group's blocks to concatenate.
#' @return A single named numeric vector (posterior mean per ID), concatenated
#'   across every block that group participates in.
#' @export
get_muk_block <- function(fit, group) {
  mus <- list()
  for (b in fit$blocks) {
    key <- fit$block_index[[group]][b]
    if (is.na(key)) next
    res <- fit$block_results[[key]]
    if (!(group %in% names(res$groups))) next
    mus[[b]] <- res$groups[[group]]$muk
  }
  ## do.call(c, unname(mus)) rather than unlist(mus): unlist() would prefix
  ## each ID name with its block name (named list of named vectors), which
  ## would break ID matching against Sigma_blocks' rownames.
  do.call(c, unname(mus))
}

#' Assemble a Dense Block-Diagonal Covariance Matrix
#'
#' @description
#' Assembles a list of per-block covariance matrices (as returned by
#' \code{\link{get_sigmak_block}}) into one dense block-diagonal matrix --
#' e.g. to feed into the dense \code{\link{evaluate_metric}} family (such as
#' \code{\link{wasserstein_metric}}) when a metric has no block-decoupled
#' shortcut of its own. Off-block-diagonal entries are exactly \code{0}, not
#' merely small.
#'
#' @param sigma_blocks A named list of covariance matrices, e.g. from
#'   \code{\link{get_sigmak_block}}.
#' @return A dense \code{p x p} matrix (\code{p} = total IDs across blocks),
#'   with \code{dimnames} set to the IDs in block order.
#' @export
as_block_diag_matrix <- function(sigma_blocks) {
  sizes <- vapply(sigma_blocks, nrow, integer(1))
  p <- sum(sizes)
  M <- matrix(0, p, p)
  offset <- 0
  all_ids <- character(0)
  for (Sb in sigma_blocks) {
    n_b <- nrow(Sb)
    M[(offset + 1):(offset + n_b), (offset + 1):(offset + n_b)] <- Sb
    all_ids <- c(all_ids, rownames(Sb))
    offset <- offset + n_b
  }
  dimnames(M) <- list(all_ids, all_ids)
  M
}
