#' @noRd
#'
#' @details Shared diagonal-jitter loop used by both \code{chol_inv_jitter()}
#'   (returns the inverse) and \code{chol_inv_jitter_diag()} in simu_db.R
#'   (returns the jittered matrix itself). \code{transform} is applied to the
#'   jittered matrix on each attempt and must throw an error (or return NULL)
#'   if the matrix is not positive-definite. Handles \code{pen_diag = 0} by
#'   bumping to \code{1e-6} only once jitter is actually needed, instead of
#'   multiplying zero by ten forever.
jitter_until_pd <- function(mat, pen_diag, transform, max_tries = 20,
                             warn_ratio = 100, label = "jitter_until_pd") {
  initial_pen_diag <- if (pen_diag == 0) 1e-6 else pen_diag
  current <- pen_diag
  for (i in seq_len(max_tries)) {
    mat_j <- mat
    diag(mat_j) <- diag(mat_j) + current
    result <- tryCatch(transform(mat_j), error = function(e) NULL)
    if (!is.null(result)) {
      if (current > warn_ratio * initial_pen_diag) {
        warning(
          label, ": required a jitter of ", signif(current, 3),
          " (", signif(current / initial_pen_diag, 3), "x the initial pen_diag = ",
          signif(initial_pen_diag, 3), ") to reach positive-definiteness; ",
          "this may indicate ill-conditioned kernel hyperparameters."
        )
      }
      return(result)
    }
    last_attempted <- current
    current <- if (current == 0) 1e-6 else 10 * current
  }
  stop(label, ": matrix could not be made positive-definite after ",
       max_tries, " jitter attempts (last attempted pen_diag = ", last_attempted, ").")
}

#' @noRd
#'
#' @details Memoizes the result of `compute()` (expected to be expensive, e.g.
#'   building and inverting a kernel covariance matrix) keyed by exact equality
#'   of `free` (the free/unconstrained hyperparameter vector) against the last
#'   call. `cache` must be a single environment reused across calls (e.g. one
#'   created per `optim_hp()` run); if `cache` is `NULL`, `compute()` is always
#'   called fresh (no memoization).
cached_cov_inv <- function(cache, free, compute) {
  if (is.null(cache)) {
    return(compute())
  }
  if (!is.null(cache$free) && identical(cache$free, free)) {
    return(cache$inv)
  }
  inv <- compute()
  cache$free <- free
  cache$inv  <- inv
  inv
}

#' @noRd
chol_inv_jitter <- function(mat, pen_diag, max_tries = 20, warn_ratio = 100) {
  jitter_until_pd(mat, pen_diag, function(m) chol2inv(chol(m)),
                  max_tries, warn_ratio, label = "chol_inv_jitter")
}

#' @noRd
#'
#' @details Centers `output` on its own `group`-specific empirical mean --
#'   used by `optim_hp()`'s `group_col` argument to pool replicates from
#'   multiple groups without an unmodeled between-group mean shift
#'   contaminating the fitted kernel hyperparameters (see
#'   `dev/optim_exploration/NOTES_math.md`, sec. "Confusion du decalage de
#'   moyenne entre groupes"). `tapply()` returns a 1-d array; naive
#'   arithmetic on its result (`output - grp_means[group]`) silently keeps a
#'   stray `dim` attribute that breaks `dmnorm()`'s `is.vector(x)` check
#'   downstream, producing a cryptic "missing value where TRUE/FALSE needed"
#'   error with no hint of the real cause -- guarded here with
#'   `unname(as.numeric(...))` and asserted via `stopifnot()`.
demean_by_group <- function(output, group) {
  grp_means <- tapply(output, group, mean)
  demeaned  <- unname(as.numeric(output - grp_means[group]))
  stopifnot(is.null(dim(demeaned)))
  list(output = demeaned, n_groups = length(unique(group)))
}

#' @noRd
#'
#' @details Appends one row to `trace_log$rows` (a plain list, grown by
#'   index) recording a single objective/gradient evaluation. No-ops when
#'   `trace_log` is `NULL`, so callers can pass it unconditionally without
#'   branching, and `track_trace = FALSE` costs nothing beyond the `NULL`
#'   check. `value` is the NLL for `type = "fn"` rows or the gradient norm
#'   for `type = "gr"` rows; `hp` is recorded as a named list so rows can be
#'   reassembled into a data.frame with one column per hyperparameter.
record_trace <- function(trace_log, type, hp, value, t0) {
  if (is.null(trace_log)) {
    return(invisible(NULL))
  }
  i <- length(trace_log$rows) + 1L
  trace_log$rows[[i]] <- c(
    list(eval_type = type, eval_index = i),
    as.list(hp),
    list(value = value, elapsed_sec = as.numeric(proc.time()["elapsed"] - t0))
  )
  invisible(NULL)
}

#' @noRd
#'
#' @details Reassembles `trace_log$rows` (set by `record_trace()`) into a
#'   single data.frame, one row per evaluation. Returns `NULL` when
#'   `trace_log` is `NULL` or no evaluations were recorded.
build_trace_df <- function(trace_log) {
  if (is.null(trace_log) || length(trace_log$rows) == 0) {
    return(NULL)
  }
  do.call(rbind, lapply(trace_log$rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
}


#' @noRd
#'
#' @details Accepts `prior_cov` as either a single numeric value (applied as a
#'   uniform diagonal of size `n`) or a square matrix of size `n x n`. Any other
#'   shape produces an informative error naming the offending class/dimensions,
#'   instead of the cryptic failure that would occur from indexing `dim()` on a
#'   non-matrix object.
resolve_prior_cov <- function(prior_cov, n) {
  if (length(prior_cov) == 1 && is.numeric(prior_cov)) {
    return(diag(n) * prior_cov)
  }
  if (!is.matrix(prior_cov)) {
    stop(
      "prior_cov must be a single numeric value (used as a uniform diagonal) ",
      "or a square matrix of size ", n, "x", n, ", but got an object of class '",
      class(prior_cov)[1], "' with length ", length(prior_cov), "."
    )
  }
  if (nrow(prior_cov) != n || ncol(prior_cov) != n) {
    stop(
      "prior_cov must be a square matrix of size ", n, "x", n,
      ", but got a ", nrow(prior_cov), "x", ncol(prior_cov), " matrix."
    )
  }
  prior_cov
}

#' @noRd
# TODO: dmnorm()'s mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE) below
# duplicates mu into a full n x p matrix. Currently harmless because the only
# call site (sum_logGaussian(), passing db$Output as a single vector) always
# has n = 1, but this is an un-enforced assumption -- if dmnorm() is ever
# called with n > 1 this becomes an O(n*p) allocation with no guard or test.
dmnorm <- function(x, mu, inv_Sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  n <- nrow(x)
  p <- ncol(x)

  if (is.vector(mu)) {
    if (length(mu) != p) {
      stop("The length of mu must match the number of columns in x.")
    }
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  } else if (is.matrix(mu)) {
    if (ncol(mu) != p || nrow(mu) != n) {
      stop("The dimensions of mu must match the dimensions of x.")
    }
  } else {
    stop("mu must be a vector or a matrix.")
  }

  z <- x - mu


  if (ncol(z) != nrow(inv_Sigma)) {
    if (nrow(z) == nrow(inv_Sigma)) {
      # Transpose z if necessary
      z <- t(z)
    } else {
      stop("The number of columns in z must match the number of rows in inv_Sigma.")
    }
  }

  logdetS <- tryCatch(
    -determinant(inv_Sigma, logarithm = TRUE)$modulus,
    error = function(e) {
      stop(
        "dmnorm: failed to compute the determinant of 'inv_Sigma' (", e$message, "). ",
        "This usually means 'inv_Sigma' is not a valid precision matrix; check the ",
        "kernel covariance and the 'pen_diag' jitter used to build it.",
        call. = FALSE
      )
    }
  )
  attributes(logdetS) <- NULL

  ssq <- sum((z %*% inv_Sigma) * z)

  loglik <- as.vector(-(p * log(2 * pi) + logdetS + ssq) / 2)

  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}

#' @noRd
#'
#' @details Recovers one Output value per (grouping key) from a `db` whose
#'   Output is repeated identically across the D `Input_ID` rows of one
#'   observation, aligning the result to `row_order` (a character vector of
#'   `ID` values, e.g. `rownames(input_mat)`). Errors if any ID in `row_order`
#'   is missing from `sub`, or if `sub` doesn't have exactly one Output value
#'   per ID (a contract violation -- Output must not vary across the
#'   Input_ID rows of a single observation).
aligned_output <- function(sub, row_order) {
  out <- dplyr::distinct(sub, .data$ID, .data$Output)
  if (nrow(out) != length(row_order) || !setequal(out$ID, row_order)) {
    stop(
      "Each observation must have exactly one Output value per ID, matching ",
      "the ID set of the Input positions -- check for a ragged design or an ",
      "Output value that varies across the Input_ID rows of one observation."
    )
  }
  out$Output[match(row_order, out$ID)]
}

#' @noRd
#'
#' @details If `cache` is a (single, fresh-per-optimization-run) environment,
#'   the kernel covariance and its inverse are memoized by `free`: since
#'   `optim_hp()`'s objective and gradient are evaluated at the same
#'   `free` within one L-BFGS-B iteration, this avoids recomputing
#'   `evaluate()` (an O(n^3) operation) twice for the same point.
#'   `cache = NULL` (the default) disables memoization and recomputes as before.
#'
#'   `db` must already contain `Input_ID`/`Input` columns (see
#'   `normalize_input_cols()`, called by `optim_hp()` before this function is
#'   ever reached). When `db` has exactly one distinct `Input_ID` (the legacy
#'   scalar-Input case, or an explicitly single-dimensional design), the
#'   original scalar-position logic runs unchanged: if `db` also contains a
#'   `Sample` column and duplicate `Input` values (i.e. multiple replicates
#'   share the same input positions), the function uses the correct
#'   replicated likelihood: it builds a pxp kernel matrix on the *unique*
#'   input positions and sums independent log N(y_s; mu, K_p) terms over
#'   each sample s. When `db` has more than one distinct `Input_ID` (a
#'   genuinely multi-dimensional Input), the same replicated/non-replicated
#'   logic is generalized over a p x D matrix of unique ID positions,
#'   requiring an `ID` column to group each observation's Input_ID rows.
#'
#'   `n_groups`, when non-`NULL`, adds the closed-form REML correction for
#'   `n_groups` group-specific means already profiled out of `db$Output`
#'   upstream (by `optim_hp()`'s `group_col` demeaning step -- see there):
#'   `NLL_REML(theta) = NLL(theta) - (n_groups/2) * log|Sigma'_theta|`. This
#'   is the exact restricted-likelihood correction for `n_groups` fixed,
#'   distinct group means estimated by their own empirical mean (proved in
#'   `dev/optim_exploration/NOTES_math.md`, sec. "Confusion du decalage de
#'   moyenne entre groupes": since Sigma'_theta weighs every replicate
#'   identically, the GLS group-mean estimator collapses exactly to the
#'   plain empirical mean, so the only missing piece versus the true
#'   restricted likelihood is this log-determinant term, which compensates
#'   for the `n_groups` degrees of freedom spent estimating those means).
#'   Requires the `is_replicated` branch (a balanced Group x Sample design);
#'   errors if `n_groups` is supplied for a non-replicated `db`, since the
#'   correction's derivation assumes a well-defined shared p x p Sigma'_theta
#'   reused identically across every replicate.
sum_logGaussian <- function(free, db, prior_mean, kern, prior_cov, pen_diag, cache = NULL,
                            trace_log = NULL, t0 = NULL, n_groups = NULL) {
  kern <- keRnel::set_free_params(kern, free)
  db   <- normalize_input_cols(db)
  n_dim <- length(unique(db$Input_ID))

  if (n_dim == 1) {
    # -- D = 1: legacy scalar-Input logic, unchanged --------------------------
    unique_inputs <- sort(unique(db$Input))
    p             <- length(unique_inputs)
    is_replicated <- FALSE
    if ("Sample" %in% names(db) && nrow(db) > p) {
      grp_sample_id <- if ("Group" %in% names(db))
        interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
      sample_counts <- table(grp_sample_id)
      is_replicated <- all(sample_counts == p)
    }

    if (is_replicated) {
      input_mat <- as.matrix(unique_inputs)

      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, p)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, p) else prior_mean

      grp_sample_id <- if ("Group" %in% names(db))
        interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
      replicates <- levels(grp_sample_id)
      log_lik <- 0
      for (gs in replicates) {
        sub <- db[grp_sample_id == gs, ]
        y_s <- sub$Output[order(sub$Input)]
        log_lik <- log_lik + dmnorm(y_s, mean_vec, inv, log = TRUE)
      }
      neg_sum_log_likelihood <- -log_lik

      if (!is.null(n_groups)) {
        log_det_sigma <- -as.numeric(determinant(inv, logarithm = TRUE)$modulus)
        neg_sum_log_likelihood <- neg_sum_log_likelihood - (n_groups / 2) * log_det_sigma
      }

    } else {
      if (!is.null(n_groups)) {
        stop(
          "n_groups (REML group correction) requires a balanced replicated ",
          "design (the same features observed for every Group x Sample ",
          "combination); this dataset does not qualify (ragged or ",
          "single-replicate design)."
        )
      }
      input <- as.matrix(db$Input)
      n     <- nrow(input)

      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input, input)
        cov <- cov + resolve_prior_cov(prior_cov, n)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      if (length(prior_mean) == 1) {
        mean_vec <- rep(prior_mean, n)
      } else {
        mean_vec <- prior_mean
        if (length(mean_vec) != n) {
          stop("The length of 'prior_mean' must match the number of rows in db$Input")
        }
      }

      neg_sum_log_likelihood <- -sum(dmnorm(db$Output, mean_vec, inv, log = TRUE))
    }

  } else {
    # -- D > 1: multi-dimensional Input, generalized over the same logic -----
    if (!("ID" %in% names(db))) {
      stop(
        "Multi-dimensional Input (more than one distinct Input_ID) requires ",
        "an 'ID' column to group each observation's Input_ID rows."
      )
    }
    positions <- input_matrix_by_id(
      dplyr::distinct(db, .data$ID, .data$Input_ID, .data$Input),
      key_cols = "ID"
    )
    input_mat <- positions$matrix
    p         <- nrow(input_mat)

    grp_sample_id <- if ("Sample" %in% names(db)) {
      if ("Group" %in% names(db)) interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
    } else {
      NULL
    }
    replicates    <- if (!is.null(grp_sample_id)) levels(factor(grp_sample_id)) else NULL
    obs_per_rep   <- if (!is.null(replicates)) {
      vapply(replicates, function(gs) dplyr::n_distinct(db$ID[grp_sample_id == gs]), integer(1))
    } else {
      integer(0)
    }
    is_replicated <- length(replicates) > 1 && all(obs_per_rep == p)

    if (is_replicated) {
      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, p)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, p) else prior_mean

      log_lik <- 0
      for (gs in replicates) {
        sub <- db[grp_sample_id == gs, ]
        y_s <- aligned_output(sub, rownames(input_mat))
        log_lik <- log_lik + dmnorm(y_s, mean_vec, inv, log = TRUE)
      }
      neg_sum_log_likelihood <- -log_lik

      if (!is.null(n_groups)) {
        log_det_sigma <- -as.numeric(determinant(inv, logarithm = TRUE)$modulus)
        neg_sum_log_likelihood <- neg_sum_log_likelihood - (n_groups / 2) * log_det_sigma
      }

    } else {
      if (!is.null(n_groups)) {
        stop(
          "n_groups (REML group correction) requires a balanced replicated ",
          "design (the same features observed for every Group x Sample ",
          "combination); this dataset does not qualify (ragged or ",
          "single-replicate design)."
        )
      }
      n <- p
      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, n)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      if (length(prior_mean) == 1) {
        mean_vec <- rep(prior_mean, n)
      } else {
        mean_vec <- prior_mean
        if (length(mean_vec) != n) {
          stop("The length of 'prior_mean' must match the number of rows in db$Input")
        }
      }

      output <- aligned_output(db, rownames(input_mat))
      neg_sum_log_likelihood <- -sum(dmnorm(output, mean_vec, inv, log = TRUE))
    }
  }

  # get_free_params() is deliberately unnamed/positional; label a copy with
  # the corresponding natural-space hyperparameter names (same length, same
  # tree-walk order) purely for a readable trace log -- `free` itself stays
  # untouched, since cached_cov_inv() keys off its exact (unnamed) identity.
  named_free <- free
  names(named_free) <- names(keRnel::get_trainable_params(kern))
  record_trace(trace_log, "fn", named_free, neg_sum_log_likelihood, t0)
  return(neg_sum_log_likelihood)
}

#' @noRd
#'
#' @details Analytic gradient counterpart of `sum_logGaussian()`, in the same
#'   free/unconstrained space: a single `keRnel::kernel_grad_free()` call
#'   returns the covariance derivative with respect to every non-frozen
#'   hyperparameter at once (positional, in `keRnel::get_free_params()`'s
#'   tree-walk order) -- replacing the previous per-hyperparameter
#'   `keRnel::kernel_deriv()` loop entirely. Working positionally in free
#'   space (rather than by natural-space hyperparameter name) is what lets
#'   this stay correct even when `kern` has two sibling sub-kernels sharing a
#'   bare hyperparameter name (e.g. `se_kernel(length_scale=1) +
#'   se_kernel(length_scale=3)`), which `keRnel::kupdate()` cannot address
#'   independently -- see `?keRnel::kernel_grad_free`.
#'
#'   `n_groups`: gradient counterpart of `sum_logGaussian()`'s REML
#'   correction. `d/dtheta_j[-(n_groups/2) * log|Sigma'_theta|] =
#'   -(n_groups/2) * tr(Sigma'^-1 %*% dK/dtheta_j)`, folded into
#'   `common_term` as `+ n_groups * inv` before the existing
#'   `-0.5 * sum(common_term * t(dK))` contraction (so `common_term`'s
#'   sign convention there already accounts for the two minus signs
#'   cancelling) -- see `dev/optim_exploration/11_reml_group_gradient_check.R`
#'   for the finite-difference verification of this term specifically.
gr_sum_logGaussian <- function(free, db, prior_mean, kern, prior_cov, pen_diag, cache = NULL,
                               trace_log = NULL, t0 = NULL, n_groups = NULL) {
  kern  <- keRnel::set_free_params(kern, free)
  db    <- normalize_input_cols(db)
  n_dim <- length(unique(db$Input_ID))

  if (n_dim == 1) {
    # -- D = 1: legacy scalar-Input logic ------------------------------------
    unique_inputs <- sort(unique(db$Input))
    p             <- length(unique_inputs)
    is_replicated <- FALSE
    if ("Sample" %in% names(db) && nrow(db) > p) {
      grp_sample_id <- if ("Group" %in% names(db))
        interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
      sample_counts <- table(grp_sample_id)
      is_replicated <- all(sample_counts == p)
    }

    if (is_replicated) {
      input_mat <- as.matrix(unique_inputs)

      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, p)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, p) else prior_mean
      grp_sample_id <- if ("Group" %in% names(db))
        interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
      replicates <- levels(grp_sample_id)
      n_s        <- length(replicates)

      outer_sum <- matrix(0, p, p)
      for (gs in replicates) {
        sub  <- db[grp_sample_id == gs, ]
        y_s  <- sub$Output[order(sub$Input)]
        r_s  <- inv %*% (y_s - mean_vec)
        outer_sum <- outer_sum + r_s %*% t(r_s)
      }
      common_term <- outer_sum - n_s * inv
      if (!is.null(n_groups)) {
        common_term <- common_term + n_groups * inv
      }

      grads <- keRnel::kernel_grad_free(kern, input_mat, input_mat)
      grad  <- vapply(grads, function(dK) -0.5 * sum(common_term * t(dK)), numeric(1))

    } else {
      if (!is.null(n_groups)) {
        stop(
          "n_groups (REML group correction) requires a balanced replicated ",
          "design (the same features observed for every Group x Sample ",
          "combination); this dataset does not qualify (ragged or ",
          "single-replicate design)."
        )
      }
      output <- db$Output
      input  <- as.matrix(db$Input)
      n      <- nrow(input)

      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input, input)
        cov <- cov + resolve_prior_cov(prior_cov, n)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, n) else prior_mean

      prod_inv    <- inv %*% (output - mean_vec)
      common_term <- prod_inv %*% t(prod_inv) - inv

      grads <- keRnel::kernel_grad_free(kern, input, input)
      grad  <- vapply(grads, function(dK) -0.5 * sum(common_term * t(dK)), numeric(1))
    }

  } else {
    # -- D > 1: multi-dimensional Input ---------------------------------------
    if (!("ID" %in% names(db))) {
      stop(
        "Multi-dimensional Input (more than one distinct Input_ID) requires ",
        "an 'ID' column to group each observation's Input_ID rows."
      )
    }
    positions <- input_matrix_by_id(
      dplyr::distinct(db, .data$ID, .data$Input_ID, .data$Input),
      key_cols = "ID"
    )
    input_mat <- positions$matrix
    p         <- nrow(input_mat)

    grp_sample_id <- if ("Sample" %in% names(db)) {
      if ("Group" %in% names(db)) interaction(db$Group, db$Sample, drop = TRUE) else db$Sample
    } else {
      NULL
    }
    replicates  <- if (!is.null(grp_sample_id)) levels(factor(grp_sample_id)) else NULL
    obs_per_rep <- if (!is.null(replicates)) {
      vapply(replicates, function(gs) dplyr::n_distinct(db$ID[grp_sample_id == gs]), integer(1))
    } else {
      integer(0)
    }
    is_replicated <- length(replicates) > 1 && all(obs_per_rep == p)

    if (is_replicated) {
      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, p)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, p) else prior_mean
      n_s      <- length(replicates)

      outer_sum <- matrix(0, p, p)
      for (gs in replicates) {
        sub <- db[grp_sample_id == gs, ]
        y_s <- aligned_output(sub, rownames(input_mat))
        r_s <- inv %*% (y_s - mean_vec)
        outer_sum <- outer_sum + r_s %*% t(r_s)
      }
      common_term <- outer_sum - n_s * inv
      if (!is.null(n_groups)) {
        common_term <- common_term + n_groups * inv
      }

      grads <- keRnel::kernel_grad_free(kern, input_mat, input_mat)
      grad  <- vapply(grads, function(dK) -0.5 * sum(common_term * t(dK)), numeric(1))

    } else {
      if (!is.null(n_groups)) {
        stop(
          "n_groups (REML group correction) requires a balanced replicated ",
          "design (the same features observed for every Group x Sample ",
          "combination); this dataset does not qualify (ragged or ",
          "single-replicate design)."
        )
      }
      n <- p
      inv <- cached_cov_inv(cache, free, function() {
        cov <- keRnel::evaluate(kern, input_mat, input_mat)
        cov <- cov + resolve_prior_cov(prior_cov, n)
        chol_inv_jitter(cov, pen_diag = pen_diag)
      })

      mean_vec <- if (length(prior_mean) == 1) rep(prior_mean, n) else prior_mean
      output   <- aligned_output(db, rownames(input_mat))

      prod_inv    <- inv %*% (output - mean_vec)
      common_term <- prod_inv %*% t(prod_inv) - inv

      grads <- keRnel::kernel_grad_free(kern, input_mat, input_mat)
      grad  <- vapply(grads, function(dK) -0.5 * sum(common_term * t(dK)), numeric(1))
    }
  }

  named_free <- free
  names(named_free) <- names(keRnel::get_trainable_params(kern))
  record_trace(trace_log, "gr", named_free, sqrt(sum(grad^2)), t0)
  return(grad)
}



#' Optimize Hyperparameters for a kernel (with additive kernel support)
#'
#' @param kern A kernel object inheriting from keRnel's `kernel` class. Its
#'   construction values (e.g. \code{se_kernel(length_scale = 2)}) are the
#'   optimization's starting point -- there is no separate `hp` argument.
#'   Optimizing in \code{keRnel::get_free_params(kern)}'s free/unconstrained
#'   space (rather than by natural-space hyperparameter name) is what lets
#'   two sibling sub-kernels sharing a bare hyperparameter name (e.g.
#'   \code{se_kernel(length_scale=1) + se_kernel(length_scale=3)}) be given
#'   genuinely different starting values and be fit independently -- a case
#'   \code{keRnel::kupdate()} cannot express, since it updates every
#'   occurrence of a bare name at once.
#' @param db The dataset used for optimization. Must contain columns `Input` and
#'   `Output` (plus `Input_ID` for a multi-dimensional Input; see
#'   \code{\link{normalize_input_cols}}). If it also contains a `Sample` column
#'   and the same `Input` positions are repeated across samples (i.e.
#'   replicated measurements), the function automatically uses the correct
#'   replicated likelihood: it builds a pxp kernel matrix on the *unique*
#'   input positions and sums independent log N(y_s; mu, K_p) terms over each
#'   sample s. Without this correction, same-position observations in
#'   different samples would be treated as perfectly correlated, biasing all
#'   hyperparameter estimates.
#' @param prior_mean Prior mean: either a scalar or a vector of length equal to
#'   the number of unique `Input` positions in `db`. Required unless
#'   `group_col` is supplied (in which case it is ignored -- see `group_col`).
#' @param prior_cov Prior covariance: either a scalar (diagonal value) or a
#'   square matrix of size equal to the number of unique `Input` positions in `db`.
#' @param pen_diag Jitter added to the diagonal for numerical stability. Defaults to `1e-6`.
#' @param verbose If `FALSE` (default), returns the optimized parameter vector;
#'   if `TRUE`, returns the full `optim` result list (with the fitted kernel
#'   object attached as `result$kern`).
#' @param max_iter Maximum number of L-BFGS-B iterations. Defaults to `1000`.
#' @param factr L-BFGS-B relative convergence tolerance, passed straight through
#'   to `stats::optim()`'s `control$factr`. Smaller values demand tighter
#'   convergence (more iterations); defaults to `1e7` (the value previously
#'   hardcoded).
#' @param pgtol L-BFGS-B projected-gradient convergence tolerance, passed
#'   straight through to `stats::optim()`'s `control$pgtol`. Defaults to `0`
#'   (R's own `optim()` default, previously left unset).
#' @param track_trace If `TRUE`, records every objective/gradient evaluation
#'   during the optimization (hyperparameter values, NLL or gradient norm,
#'   elapsed time) and attaches it as a `trace` data.frame. Defaults to
#'   `FALSE`, in which case nothing is recorded and the return value is
#'   identical to before this parameter existed.
#' @param group_col Name of a column in `db` (e.g. `"Group"`) identifying
#'   which experimental group each observation belongs to. When supplied,
#'   pools replicates from every group into a single fit -- centering each
#'   group on its own empirical mean first (`prior_mean` is then ignored,
#'   forced to `0`) -- and adds the exact closed-form REML correction for
#'   the degrees of freedom spent estimating those group means:
#'   \eqn{\text{NLL}_\text{REML}(\theta) = \text{NLL}_\text{demeaned}(\theta) -
#'   (G/2)\log|\Sigma'_\theta|}, \eqn{G} = number of distinct groups. Without
#'   this, pooling multiple groups under one shared `prior_mean` lets an
#'   unmodeled between-group mean shift bias the fitted kernel variance
#'   (severely, and roughly quadratically in the shift) -- see
#'   `dev/optim_exploration/NOTES_math.md` and
#'   `dev/optim_exploration/10_group_mean_shift_confound.R` for the full
#'   derivation and empirical quantification. Requires a balanced Group x
#'   Sample design (the same features observed for every group/sample
#'   combination); errors otherwise. Defaults to `NULL` (no group handling,
#'   behavior identical to before this parameter existed).
#'
#' @return If `verbose` is `FALSE`, a named vector of optimized hyperparameters
#'   in natural (constrained) space (via `keRnel::get_trainable_params()`),
#'   with the optimizer's `convergence` code and final objective `value`
#'   attached as attributes (`attr(result, "convergence")`, `attr(result,
#'   "value")`) so convergence can be checked without re-running with `verbose
#'   = TRUE`; a `convergence` of `0` means success. If `track_trace = TRUE`, a
#'   `trace` data.frame (one row per objective/gradient evaluation) is also
#'   attached as `attr(result, "trace")`. Otherwise (`verbose = TRUE`) the full
#'   list returned by [stats::optim()], with `result$kern` (the fitted kernel
#'   object) and `result$trace` (when `track_trace = TRUE`) added.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 1, nb_sample = 5)
#' kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
#' optim_hp(kern, data, prior_mean = 0, prior_cov = 1)
#'
#' # Pooling replicates from multiple groups into a single fit, without an
#' # unmodeled group mean shift biasing the fitted variance:
#' data2 <- simu_db(nb_id = 8, nb_group = 3, nb_sample = 5, diff_group = 4)
#' optim_hp(kern, data2, prior_cov = 1, group_col = "Group")
optim_hp <- function(kern, db, prior_mean = NULL, prior_cov,
                      pen_diag = 1e-6, verbose = FALSE,
                      max_iter = 1000,
                      factr = 1e7, pgtol = 0, track_trace = FALSE,
                      group_col = NULL) {
  if (!inherits(kern, "kernel")) {
    stop("'kern' must be a valid kernel object from the keRnel package.")
  }
  if (!all(c("Input", "Output") %in% names(db))) {
    stop("'db' must contain columns 'Input' and 'Output'.")
  }
  db <- normalize_input_cols(db)
  if (any(!is.finite(db$Input)) || any(!is.finite(db$Output))) {
    stop("'db$Input' and 'db$Output' must not contain NaN, Inf, or NA values.")
  }

  n_groups <- NULL
  if (!is.null(group_col)) {
    if (!group_col %in% names(db)) {
      stop("'group_col' (\"", group_col, "\") is not a column of 'db'.")
    }
    if (!is.null(prior_mean) && !isTRUE(all.equal(unname(prior_mean), 0))) {
      warning(
        "'prior_mean' is ignored when 'group_col' is supplied: each group is ",
        "centered on its own empirical mean first, forcing an effective ",
        "prior_mean of 0."
      )
    }
    demeaned   <- demean_by_group(db$Output, db[[group_col]])
    db$Output  <- demeaned$output
    n_groups   <- demeaned$n_groups
    prior_mean <- 0
  } else if (is.null(prior_mean)) {
    stop("'prior_mean' is required unless 'group_col' is supplied.")
  }

  if (length(prior_mean) != 1) {
    n_positions <- if (length(unique(db$Input_ID)) == 1) {
      length(unique(db$Input))
    } else {
      nrow(input_matrix_by_id(dplyr::distinct(db, .data$ID, .data$Input_ID, .data$Input), "ID")$matrix)
    }
    if (length(prior_mean) != n_positions) {
      stop("'prior_mean' must be a scalar or a vector of length equal to the number of unique Input positions in 'db'.")
    }
  }

  # Shared across objective/gradient for this optimization run only: avoids
  # recomputing the O(n^3) kernel covariance + inverse when L-BFGS-B evaluates
  # both at the same free-space point (the common case).
  cov_cache <- new.env(parent = emptyenv())

  # trace_log stays NULL (and record_trace() no-ops) when track_trace = FALSE,
  # so the default path costs nothing beyond this one check per evaluation.
  trace_log <- if (track_trace) list2env(list(rows = list()), parent = emptyenv()) else NULL
  t0 <- proc.time()["elapsed"]

  objective <- function(free) {
    sum_logGaussian(free, db, prior_mean, kern, prior_cov, pen_diag, cache = cov_cache,
                     trace_log = trace_log, t0 = t0, n_groups = n_groups)
  }

  gradient <- function(free) {
    gr_sum_logGaussian(free, db, prior_mean, kern, prior_cov, pen_diag, cache = cov_cache,
                        trace_log = trace_log, t0 = t0, n_groups = n_groups)
  }

  free0 <- keRnel::get_free_params(kern)

  result <- stats::optim(
    par = free0,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B",
    control = list(factr = factr, pgtol = pgtol, maxit = max_iter)
  )

  if (result$convergence != 0) {
    warning(
      "L-BFGS-B did not converge (code ", result$convergence, "). ",
      "Results may be unreliable. Consider increasing max_iter or adjusting pen_diag."
    )
  }

  kern_fit <- keRnel::set_free_params(kern, result$par)
  trace_df <- build_trace_df(trace_log)

  if (!verbose) {
    par <- keRnel::get_trainable_params(kern_fit)
    attr(par, "convergence") <- result$convergence
    attr(par, "value") <- result$value
    if (track_trace) attr(par, "trace") <- trace_df
    return(par)
  } else {
    result$kern <- kern_fit
    if (track_trace) result$trace <- trace_df
    return(result)
  }
}
