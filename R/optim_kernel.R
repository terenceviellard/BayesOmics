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
#'   of `hp` against the last call. `cache` must be a single environment reused
#'   across calls (e.g. one created per `optim_hp()` run); if `cache`
#'   is `NULL`, `compute()` is always called fresh (no memoization).
cached_cov_inv <- function(cache, hp, compute) {
  if (is.null(cache)) {
    return(compute())
  }
  if (!is.null(cache$hp) && identical(cache$hp, hp)) {
    return(cache$inv)
  }
  inv <- compute()
  cache$hp  <- hp
  cache$inv <- inv
  inv
}

#' @noRd
chol_inv_jitter <- function(mat, pen_diag, max_tries = 20, warn_ratio = 100) {
  jitter_until_pd(mat, pen_diag, function(m) chol2inv(chol(m)),
                  max_tries, warn_ratio, label = "chol_inv_jitter")
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
#' @details If `cache` is a (single, fresh-per-optimization-run) environment,
#'   the kernel covariance and its inverse are memoized by `hp`: since
#'   `optim_hp()`'s objective and gradient are evaluated at the same
#'   `hp` within one L-BFGS-B iteration, this avoids recomputing
#'   `pairwise_kernel()` (an O(n^3) operation) twice for the same point.
#'   `cache = NULL` (the default) disables memoization and recomputes as before.
sum_logGaussian <- function(hp, db, prior_mean, kern, prior_cov, pen_diag, cache = NULL,
                            trace_log = NULL, t0 = NULL) {
  kern <- keRnel::set_hyperparameters(kern, hp)
  input <- db$Input
  input <- as.matrix(input)
  n <- nrow(input)

  inv <- cached_cov_inv(cache, hp, function() {
    cov <- keRnel::pairwise_kernel(kern, input, input)
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

  log_likelihoods <- dmnorm(db$Output, mean_vec, inv, log = TRUE)
  neg_sum_log_likelihood <- -sum(log_likelihoods)

  record_trace(trace_log, "fn", hp, neg_sum_log_likelihood, t0)

  return(neg_sum_log_likelihood)
}

#' @noRd
#'
#' @importFrom methods is
#' @importFrom keRnel gt_HPs
gr_sum_logGaussian <- function(hp, db, prior_mean, kern, prior_cov, pen_diag, cache = NULL,
                               trace_log = NULL, t0 = NULL) {
  kern <- keRnel::set_hyperparameters(kern, hp)
  list_hp <- keRnel::get_hyperparameter_names(kern)
  output <- db$Output
  input <- db$Input
  input <- as.matrix(input)
  n <- nrow(input)

  inv <- cached_cov_inv(cache, hp, function() {
    cov <- keRnel::pairwise_kernel(kern, input, input)
    cov <- cov + resolve_prior_cov(prior_cov, n)
    chol_inv_jitter(cov, pen_diag = pen_diag)
  })

  if (length(prior_mean) == 1) {
    mean_vec <- rep(prior_mean, n)
  } else {
    mean_vec <- prior_mean
  }

  prod_inv <- inv %*% (output - mean_vec)
  common_term <- prod_inv %*% t(prod_inv) - inv

  grad <- numeric(length(list_hp))
  names(grad) <- list_hp

  for (i in seq_along(list_hp)) {
    hp_name <- list_hp[i]
    kern_deriv <- keRnel::kernel_deriv(kern, input, input, hp_name)
    grad[i] <- -0.5 * sum(common_term * t(kern_deriv))
  }

  record_trace(trace_log, "gr", hp, sqrt(sum(grad^2)), t0)

  return(grad)
}



#' Optimize Hyperparameters for a kernel (with additive kernel support)
#'
#' @param hp A vector of initial hyperparameters to be optimized. Length must
#'   match the number of hyperparameters of `kern`.
#' @param db The dataset used for optimization. Must contain columns `Input` and `Output`.
#' @param prior_mean Prior mean: either a scalar or a vector of length `nrow(db)`.
#' @param kern A kernel object inheriting from `AbstractKernel` (keRnel package).
#' @param prior_cov Prior covariance: either a scalar (diagonal value) or a
#'   square matrix of size `nrow(db)`.
#' @param pen_diag Jitter added to the diagonal for numerical stability. Defaults to `1e-6`.
#' @param verbose If `FALSE` (default), returns the optimized parameter vector;
#'   if `TRUE`, returns the full `optim` result list.
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
#'
#' @return If `verbose` is `FALSE`, a named vector of optimized hyperparameters,
#'   with the optimizer's `convergence` code and final objective `value` attached
#'   as attributes (`attr(result, "convergence")`, `attr(result, "value")`) so
#'   convergence can be checked without re-running with `verbose = TRUE`; a
#'   `convergence` of `0` means success. If `track_trace = TRUE`, a `trace`
#'   data.frame (one row per objective/gradient evaluation) is also attached
#'   as `attr(result, "trace")`. Otherwise (`verbose = TRUE`) the full list
#'   returned by [stats::optim()], with `result$trace` added when
#'   `track_trace = TRUE`.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 1, nb_sample = 5)
#' kern <- methods::new("SEKernel")
#' hp <- unlist(keRnel::gt_HPs(kern))
#' optim_hp(hp, data, prior_mean = 0, kern = kern, prior_cov = 1)
optim_hp <- function(hp, db, prior_mean, kern, prior_cov,
                               pen_diag = 1e-6, verbose = FALSE,
                               max_iter = 1000,
                               factr = 1e7, pgtol = 0, track_trace = FALSE) {
  if (!methods::is(kern, "AbstractKernel")) {
    stop("'kern' must be a valid kernel object from the keRnel package.")
  }
  if (!all(c("Input", "Output") %in% names(db))) {
    stop("'db' must contain columns 'Input' and 'Output'.")
  }
  # TODO: no check that db$Input / db$Output are finite -- NaN/Inf would
  # propagate into pairwise_kernel()/chol() and fail there with an opaque
  # error instead of a clear message at the function boundary.
  n_obs <- nrow(as.matrix(db$Input))
  if (length(prior_mean) != 1 && length(prior_mean) != n_obs) {
    stop("'prior_mean' must be a scalar or a vector of length nrow(db).")
  }

  # Shared across objective/gradient for this optimization run only: avoids
  # recomputing the O(n^3) kernel covariance + inverse when L-BFGS-B evaluates
  # both at the same hp (the common case).
  cov_cache <- new.env(parent = emptyenv())

  # trace_log stays NULL (and record_trace() no-ops) when track_trace = FALSE,
  # so the default path costs nothing beyond this one check per evaluation.
  trace_log <- if (track_trace) list2env(list(rows = list()), parent = emptyenv()) else NULL
  t0 <- proc.time()["elapsed"]

  objective <- function(hp) {
    sum_logGaussian(hp, db, prior_mean, kern, prior_cov, pen_diag, cache = cov_cache,
                     trace_log = trace_log, t0 = t0)
  }

  gradient <- function(hp) {
    gr_sum_logGaussian(hp, db, prior_mean, kern, prior_cov, pen_diag, cache = cov_cache,
                        trace_log = trace_log, t0 = t0)
  }

  lower <- rep(1e-6, length(hp))
  names(lower) <- names(hp)
  upper <- rep(1e6, length(hp))
  names(upper) <- names(hp)

  result <- stats::optim(
    par = hp,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(factr = factr, pgtol = pgtol, maxit = max_iter)
  )

  if (result$convergence != 0) {
    warning(
      "L-BFGS-B did not converge (code ", result$convergence, "). ",
      "Results may be unreliable. Consider increasing max_iter or adjusting pen_diag."
    )
  }

  trace_df <- build_trace_df(trace_log)

  if (!verbose) {
    par <- result$par
    attr(par, "convergence") <- result$convergence
    attr(par, "value") <- result$value
    if (track_trace) attr(par, "trace") <- trace_df
    return(par)
  } else {
    if (track_trace) result$trace <- trace_df
    return(result)
  }
}

