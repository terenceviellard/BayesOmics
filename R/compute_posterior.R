#' @noRd
# TODO: format_input_key() and get_sigmak() are @noRd (internal) but are
# called directly from tests/testthat/ — their signature is therefore a de
# facto stable contract without being part of the public API. Either export
# them officially or accept that test breakage on refactor is expected.
format_input_key <- function(x) sprintf("%.17g", x)

#' @noRd
get_sigmak <- function(group_entry, kernels) {
  ids <- names(group_entry$muk)
  input_keys <- format_input_key(group_entry$id_to_input[ids])
  kern_mat <- kernels[[group_entry$kernel_key]]
  sigmak <- kern_mat[input_keys, input_keys, drop = FALSE] / group_entry$scale
  dimnames(sigmak) <- list(ids, ids)
  sigmak
}

#' Compute Posterior Means for Different Groups
#'
#' This function computes the posterior means for different groups within a dataset.
#' It expects the input data frame to contain a column named 'Group' that identifies
#' the groups, as well as columns 'ID', 'Output', and 'Input'.
#'
#' The posterior distribution is given by:
#'
#' \eqn{p(\mathbf{\mu} \mid y_1, \dots, y_N, \Sigma_{\hat{\theta}}) = \mathcal{N}\left(\mathbf{\mu}; \ \dfrac{\lambda_0 \mu_0 + \sum_{n=1}^{N} y_n}{N + \lambda_0}, \dfrac{1}{N + \lambda_0} \Sigma_{\hat{\theta}}\right)}
#'
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @param data A data frame containing the data to be analyzed. Must include columns 'Group', 'ID', 'Output', and 'Input'.
#' @param kern A kernel function or object used to compute pairwise kernels.
#' @param mu_0 Prior mean parameter.
#' @param lambda_0 Prior precision parameter.
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{kernels}}{A named list of kernel/correlation matrices, one per
#'       distinct set of Input values found across groups. Each matrix has
#'       \code{dimnames} set to the (string-formatted) Input values it was built
#'       from, and is shared by reference across every group with that same set
#'       of Input values (no duplication).}
#'     \item{\code{groups}}{A named list (one entry per group) with: \code{muk},
#'       a named vector of posterior means keyed by ID; \code{id_to_input}, the
#'       ID -> Input mapping for that group; \code{kernel_key}, which entry of
#'       \code{kernels} to use; and \code{scale}, the divisor (\code{n_obs + lambda_0})
#'       applied to that kernel matrix to get the posterior covariance. The
#'       (internal) \code{get_sigmak()} helper reconstructs the actual
#'       (ID-aligned) posterior covariance matrix for a group.}
#'   }
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' posterior$groups[["1"]]$muk
multi_posterior_mean <- function(data, kern, mu_0 = 1, lambda_0 = 1) {
  # FIXME: mu_0 defaults to 1, which is unusual for a Gaussian prior mean
  # (0 would be the conventional uninformative default). Verify this is
  # intentional (e.g. tied to a specific use case) and document the
  # rationale in @details, or change the default to 0.
  # === Initial checks ===
  required_cols <- c("Group", "ID", "Output", "Input")
  if (!all(required_cols %in% names(data))) {
    stop(paste0("The following columns are missing: ", paste(setdiff(required_cols, names(data)), collapse = ", ")))
  }
  if (!is.numeric(data$Input) || !is.numeric(data$Output) || !is.numeric(mu_0)) {
    stop("The 'Input' and 'Output' columns must be numeric, and mu_0 must be numeric.")
  }
  # TODO: no check that data$Input / data$Output are finite (no NaN/Inf).
  # NaN/Inf currently propagate silently into dplyr::summarise() (sum_outputs,
  # unique_input) and surface as a cryptic failure much later (e.g. inside
  # chol()/pairwise_kernel()) rather than a clear error at the input boundary.
  # Consider: stop() if any(!is.finite(data$Input)) || any(!is.finite(data$Output)).
  if (lambda_0 <= 0 || !is.numeric(lambda_0)) {
    stop("lambda_0 must be numeric and strictly positive.")
  }
  if (any(is.na(data$Group)) || any(is.na(data$ID))) {
    stop("The 'Group' and 'ID' columns must not contain NA values.")
  }
  if (!is(kern, "AbstractKernel")) {
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
  df_mu <- data %>%
    dplyr::group_by(.data$Group, .data$ID) %>%
    dplyr::summarise(
      sum_outputs = sum(.data$Output, na.rm = TRUE),
      lengths = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(muk_vector = (lambda_0 * mu_0 + .data$sum_outputs) / (lengths + lambda_0))

  # === Calculate unique_inputs (list of sorted vectors) ===
  unique_inputs <- data %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::summarise(unique_input = list(sort(unique(.data$Input))), .groups = "drop")

  # === Per-(Group, ID) Input lookup (Input is unique per ID within a group,
  # already enforced below) so muk and the kernel matrix can be re-aligned by
  # ID rather than by position ===
  id_input <- data %>%
    dplyr::distinct(.data$Group, .data$ID, .data$Input)

  # Catch an ID mapping to several distinct Input values directly (rather than
  # relying on the muk/vec_name length-mismatch check below, which can miss
  # this if the extra Input value happens to coincide with another ID's).
  id_input_counts <- id_input %>%
    dplyr::count(.data$Group, .data$ID, name = "n_inputs")
  bad_ids <- id_input_counts[id_input_counts$n_inputs > 1, ]
  if (nrow(bad_ids) > 0) {
    stop(paste0(
      "Each ID must map to a single Input value within its group, but ",
      "Group '", bad_ids$Group[1], "', ID '", bad_ids$ID[1], "' has ",
      bad_ids$n_inputs[1], " distinct Input values."
    ))
  }

  id_input_by_group <- split(id_input, id_input$Group)

  # === Pre-compute subsets for group_data ===
  df_mu_by_group <- split(df_mu, df_mu$Group)
  unique_inputs_by_group <- stats::setNames(unique_inputs$unique_input, unique_inputs$Group)

  # === Build group_data ===
  groups <- unique(data$Group)
  group_data <- dplyr::tibble(
    group = groups,
    muk = purrr::map(groups, function(g) {
      subset <- df_mu_by_group[[g]]
      stats::setNames(subset$muk_vector, subset$ID)
    }),
    vec_name = purrr::map(groups, function(g) unique_inputs_by_group[[g]]),
    id_to_input = purrr::map(groups, function(g) {
      sub <- id_input_by_group[[g]]
      stats::setNames(sub$Input, sub$ID)
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

  # === Check muk / sigmak dimension consistency ===
  muk_lengths <- vapply(group_data$muk,      length, integer(1))
  vec_lengths <- vapply(group_data$vec_name, length, integer(1))
  mismatch    <- which(muk_lengths != vec_lengths)
  if (length(mismatch) > 0) {
    bad <- group_data$group[mismatch]
    stop(paste0(
      "In group(s) [", paste(bad, collapse = ", "), "]: ",
      "number of unique IDs (", paste(muk_lengths[mismatch], collapse = ", "), ") ",
      "does not match number of unique Input values (",
      paste(vec_lengths[mismatch], collapse = ", "), "). ",
      "Each ID must map to a distinct Input value within its group."
    ))
  }

  # === Compute kernel matrices ===
  # Kernel matrices are cached by their (sorted) set of Input values: two
  # groups with the same set of Input values share the exact same matrix
  # object (no duplication), regardless of which IDs use those values in each
  # group. The matrix's dimnames are the Input values themselves (not IDs, and
  # not positions), so it can be safely re-indexed by ID per group via
  # get_sigmak() without ever mixing up which row/col belongs to which ID.
  cache <- new.env(hash = TRUE, parent = emptyenv())

  groups_list <- lapply(seq_len(nrow(group_data)), function(i) {
    vec      <- group_data$vec_name[[i]]
    # TODO: vec_hash is a plain toString(sort(vec)) of the raw doubles, not a
    # proper hash (e.g. digest::digest()) and not normalized like
    # format_input_key()'s "%.17g" formatting used for dimnames just below.
    # Two numerically-equal-but-differently-represented vectors could in
    # theory collide or fail to collide inconsistently with the dimnames key.
    # No test currently exercises this edge case.
    vec_hash <- toString(sort(vec))
    if (!exists(vec_hash, envir = cache, inherits = FALSE)) {
      kern_mat <- keRnel::pairwise_kernel(kern, as.matrix(vec), as.matrix(vec))
      dimnames(kern_mat) <- list(format_input_key(vec), format_input_key(vec))
      assign(vec_hash, kern_mat, envir = cache)
    }
    list(
      muk         = group_data$muk[[i]],
      id_to_input = group_data$id_to_input[[i]],
      kernel_key  = vec_hash,
      scale       = group_data$n_obs[i] + lambda_0
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
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
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
