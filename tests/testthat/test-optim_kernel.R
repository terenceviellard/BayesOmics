# ── chol_inv_jitter ───────────────────────────────────────────────────────────

test_that("chol_inv_jitter returns inverse for PD matrix", {
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  inv <- BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6)
  expect_true(is.matrix(inv))
  expect_equal(dim(inv), c(2, 2))
  # A %*% A^-1 ≈ I
  expect_equal(mat %*% inv, diag(2), tolerance = 1e-6)
})

test_that("chol_inv_jitter handles near-singular matrix via jitter", {
  mat <- matrix(c(1, 1, 1, 1), nrow = 2)  # singular
  expect_no_error(BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6))
})

test_that("chol_inv_jitter result is symmetric", {
  mat <- matrix(c(4, 1, 1, 3), nrow = 2)
  inv <- BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6)
  expect_equal(inv, t(inv), tolerance = 1e-10)
})

test_that("chol_inv_jitter works for 1x1 matrix", {
  mat <- matrix(4.0, 1, 1)
  inv <- BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6)
  expect_equal(inv[1, 1], 1 / 4, tolerance = 1e-6)
})

test_that("chol_inv_jitter warns when a large jitter (relative to pen_diag) was needed", {
  mat <- diag(c(-50, 5, 5))  # one large negative eigenvalue forces many jitter doublings
  expect_warning(
    BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6),
    "ill-conditioned"
  )
})

test_that("chol_inv_jitter does not warn when no/little jitter was needed", {
  mat <- diag(c(4, 3, 2))  # already positive-definite
  expect_warning(BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6), NA)
})

test_that("chol_inv_jitter stops with an informative error instead of recursing forever", {
  # A matrix containing NaN can never become positive-definite via diagonal jitter;
  # this must not crash with a C stack overflow.
  mat <- matrix(NaN, 2, 2)
  expect_error(
    BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6, max_tries = 5),
    "could not be made positive-definite"
  )
})

# ── dmnorm ────────────────────────────────────────────────────────────────────

test_that("dmnorm returns a scalar", {
  x   <- matrix(c(0, 0), nrow = 1)
  mu  <- c(0, 0)
  inv <- diag(2)
  val <- BayesOmics:::dmnorm(x, mu, inv)
  expect_true(is.numeric(val))
  expect_length(val, 1)
})

test_that("dmnorm returns value in (0, 1] for density", {
  x   <- matrix(c(0, 0), nrow = 1)
  mu  <- c(0, 0)
  inv <- diag(2)
  val <- BayesOmics:::dmnorm(x, mu, inv)
  expect_gt(val, 0)
  expect_lte(val, 1)
})

test_that("dmnorm log=TRUE returns log of log=FALSE", {
  x   <- matrix(c(1, -1), nrow = 1)
  mu  <- c(0, 0)
  inv <- diag(2) * 2
  v   <- BayesOmics:::dmnorm(x, mu, inv, log = FALSE)
  lv  <- BayesOmics:::dmnorm(x, mu, inv, log = TRUE)
  expect_equal(lv, log(v), tolerance = 1e-10)
})

test_that("dmnorm is maximized at the mean", {
  mu  <- c(1, 2)
  inv <- diag(2)
  at_mean <- BayesOmics:::dmnorm(matrix(mu, 1), mu, inv, log = TRUE)
  off     <- BayesOmics:::dmnorm(matrix(mu + 1, 1), mu, inv, log = TRUE)
  expect_gt(at_mean, off)
})

test_that("dmnorm accepts vector x (coerced to 1-row matrix)", {
  x   <- c(0, 0)
  mu  <- c(0, 0)
  inv <- diag(2)
  expect_no_error(BayesOmics:::dmnorm(x, mu, inv))
})

test_that("dmnorm errors when mu length mismatches x ncol", {
  x   <- matrix(c(1, 2, 3), nrow = 1)
  mu  <- c(1, 2)
  inv <- diag(3)
  expect_error(BayesOmics:::dmnorm(x, mu, inv))
})

test_that("dmnorm accepts mu as a matrix matching x's dimensions", {
  x   <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  mu  <- matrix(0, nrow = 2, ncol = 2)
  inv <- diag(2)
  val_vec <- BayesOmics:::dmnorm(x, c(0, 0), inv)  # same mu via vector recycled per row
  val_mat <- BayesOmics:::dmnorm(x, mu, inv)
  expect_equal(val_mat, val_vec)
})

test_that("dmnorm errors when mu matrix dimensions do not match x", {
  x   <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  mu  <- matrix(0, nrow = 3, ncol = 2)
  inv <- diag(2)
  expect_error(BayesOmics:::dmnorm(x, mu, inv))
})

test_that("dmnorm transposes z when ncol(z) doesn't match inv_Sigma but nrow does", {
  inv_Sigma <- diag(2)
  x  <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)  # n = 2, p = 3
  mu <- matrix(0, nrow = 2, ncol = 3)
  expect_no_error(BayesOmics:::dmnorm(x, mu, inv_Sigma))
})

test_that("dmnorm raises an informative error instead of silently propagating a try-error when determinant() fails", {
  # determinant() errors on a non-square input; previously this was caught by a
  # silent try() whose try-error attributes were stripped, letting NA/garbage
  # flow into the log-likelihood instead of failing loudly.
  x   <- matrix(c(0, 0), nrow = 1)
  mu  <- c(0, 0)
  bad_inv_Sigma <- matrix(1, nrow = 2, ncol = 3)  # not square -> determinant() errors
  expect_error(
    BayesOmics:::dmnorm(x, mu, bad_inv_Sigma),
    "failed to compute the determinant"
  )
})

# ── resolve_prior_cov ─────────────────────────────────────────────────────────

test_that("resolve_prior_cov turns a numeric scalar into a uniform diagonal matrix", {
  out <- BayesOmics:::resolve_prior_cov(2.5, 3)
  expect_equal(out, diag(3) * 2.5)
})

test_that("resolve_prior_cov accepts a matrix of the right size unchanged", {
  mat <- diag(4) * 3
  out <- BayesOmics:::resolve_prior_cov(mat, 4)
  expect_identical(out, mat)
})

test_that("resolve_prior_cov errors with an informative message for a non-matrix, non-scalar input", {
  expect_error(
    BayesOmics:::resolve_prior_cov(c(1, 2, 3), 3),
    "prior_cov must be a single numeric value.*or a square matrix"
  )
})

test_that("resolve_prior_cov errors with an informative message for a wrong-size matrix", {
  expect_error(
    BayesOmics:::resolve_prior_cov(diag(2), 5),
    "prior_cov must be a square matrix of size 5x5"
  )
})

# ── cached_cov_inv ─────────────────────────────────────────────────────────────

test_that("cached_cov_inv calls compute() fresh every time when cache is NULL", {
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }
  v1 <- BayesOmics:::cached_cov_inv(NULL, c(1, 2), compute)
  v2 <- BayesOmics:::cached_cov_inv(NULL, c(1, 2), compute)
  expect_equal(c(v1, v2), c(1, 2))
})

test_that("cached_cov_inv reuses the cached result for an identical hp", {
  cache <- new.env()
  calls <- 0
  compute <- function() { calls <<- calls + 1; "result" }
  BayesOmics:::cached_cov_inv(cache, c(1, 2), compute)
  BayesOmics:::cached_cov_inv(cache, c(1, 2), compute)
  expect_equal(calls, 1)
})

test_that("cached_cov_inv recomputes when hp changes", {
  cache <- new.env()
  calls <- 0
  compute <- function() { calls <<- calls + 1; calls }
  BayesOmics:::cached_cov_inv(cache, c(1, 2), compute)
  BayesOmics:::cached_cov_inv(cache, c(3, 4), compute)
  expect_equal(calls, 2)
})

# ── sum_logGaussian ───────────────────────────────────────────────────────────

test_that("sum_logGaussian returns a finite scalar", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- keRnel::gt_HPs(kern)
  val  <- BayesOmics:::sum_logGaussian(unlist(hp), data, 0, kern, 1, 1e-6)
  expect_true(is.numeric(val))
  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("sum_logGaussian decreases when kernel fits data better", {
  set.seed(1)
  data <- make_data(nb_id = 20)
  kern <- make_kernel()
  hp0 <- unlist(keRnel::gt_HPs(kern))

  v_init <- BayesOmics:::sum_logGaussian(hp0, data, 0, kern, 1, 1e-6)
  hp_opt <- optim_hp(hp0, data, 0, kern, 1)
  v_opt  <- BayesOmics:::sum_logGaussian(hp_opt, data, 0, kern, 1, 1e-6)
  expect_lte(v_opt, v_init)
})

test_that("sum_logGaussian is larger (more negative) for mismatched prior_mean", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  val_right  <- BayesOmics:::sum_logGaussian(hp, data, mean(data$Output), kern, 1, 1e-6)
  val_wrong  <- BayesOmics:::sum_logGaussian(hp, data, mean(data$Output) + 1000, kern, 1, 1e-6)
  expect_gt(val_wrong, val_right)
})

test_that("sum_logGaussian gives the same result with and without a cache", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  no_cache <- BayesOmics:::sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  cache <- new.env()
  with_cache <- BayesOmics:::sum_logGaussian(hp, data, 0, kern, 1, 1e-6, cache = cache)
  expect_equal(with_cache, no_cache)
})

test_that("gr_sum_logGaussian gives the same result with and without a cache", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  no_cache <- BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  cache <- new.env()
  with_cache <- BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6, cache = cache)
  expect_equal(with_cache, no_cache)
})

test_that("a shared cache makes sum_logGaussian then gr_sum_logGaussian reuse the same cov inverse", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  cache <- new.env()
  BayesOmics:::sum_logGaussian(hp, data, 0, kern, 1, 1e-6, cache = cache)
  expect_identical(cache$hp, hp)
  cached_inv <- cache$inv
  BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6, cache = cache)
  expect_identical(cache$inv, cached_inv)
})

test_that("sum_logGaussian accepts a full prior_cov matrix (not just a scalar)", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  n    <- nrow(data)
  prior_cov_mat <- diag(n) * 1.5
  val_scalar <- BayesOmics:::sum_logGaussian(hp, data, 0, kern, 1.5, 1e-6)
  val_matrix <- BayesOmics:::sum_logGaussian(hp, data, 0, kern, prior_cov_mat, 1e-6)
  expect_equal(val_matrix, val_scalar, tolerance = 1e-8)
})

test_that("sum_logGaussian errors when prior_cov matrix has the wrong size", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  bad_prior_cov <- diag(2)  # wrong size vs nrow(data)
  expect_error(BayesOmics:::sum_logGaussian(hp, data, 0, kern, bad_prior_cov, 1e-6),
               "prior_cov must be a square matrix")
})

# ── gr_sum_logGaussian ────────────────────────────────────────────────────────

test_that("gr_sum_logGaussian returns numeric vector matching hp length", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  g    <- BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  expect_true(is.numeric(g))
  expect_length(g, length(hp))
})

test_that("gr_sum_logGaussian gradient is finite at reasonable hp", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  g    <- BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  expect_true(all(is.finite(g)))
})

test_that("gr_sum_logGaussian numerical gradient matches analytic gradient", {
  set.seed(2)
  data <- make_data(nb_id = 8)
  kern <- make_kernel(hp = c(2.0, 1.5))
  hp   <- unlist(keRnel::gt_HPs(kern))
  eps  <- 1e-5
  analytic <- BayesOmics:::gr_sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(hp), function(i) {
    hp_up <- hp; hp_up[i] <- hp[i] + eps
    hp_dn <- hp; hp_dn[i] <- hp[i] - eps
    (BayesOmics:::sum_logGaussian(hp_up, data, 0, kern, 1, 1e-6) -
     BayesOmics:::sum_logGaussian(hp_dn, data, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

# ── optim_hp ────────────────────────────────────────────────────────

test_that("optim_hp returns named numeric vector by default", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1)
  expect_true(is.numeric(res))
  expect_length(res, length(hp))
})

test_that("optim_hp returns full optim list when verbose=TRUE", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1, verbose = TRUE)
  expect_true(is.list(res))
  expect_true("par" %in% names(res))
  expect_true("convergence" %in% names(res))
})

test_that("optim_hp result keeps hp positive (lower = 1e-6)", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1)
  expect_true(all(res >= 1e-7))  # tolerance for floating-point at boundary
})

test_that("optim_hp errors on wrong kern type", {
  data <- make_data()
  expect_error(optim_hp(c(1, 1), data, 0, list(a = 1), 1),
               "keRnel")
})

test_that("optim_hp errors when db lacks Input/Output", {
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  expect_error(optim_hp(hp, data.frame(x = 1:5), 0, kern, 1),
               "Input.*Output|Output.*Input")
})

test_that("optim_hp errors when prior_mean has wrong length", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  expect_error(optim_hp(hp, data, c(1, 2, 3), kern, 1),
               "prior_mean")
})

test_that("track_trace = FALSE (default) attaches no trace attribute", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1)
  expect_null(attr(res, "trace"))
})

test_that("track_trace = TRUE attaches a trace data.frame with fn and gr rows", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1, track_trace = TRUE)
  tr   <- attr(res, "trace")
  expect_s3_class(tr, "data.frame")
  expect_true(all(c("eval_type", "eval_index", "value", "elapsed_sec", names(hp)) %in% names(tr)))
  expect_true("fn" %in% tr$eval_type)
  expect_true("gr" %in% tr$eval_type)
})

test_that("track_trace = TRUE does not change the optimized hyperparameters", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res_plain <- optim_hp(hp, data, 0, kern, 1)
  res_trace <- optim_hp(hp, data, 0, kern, 1, track_trace = TRUE)
  expect_equal(as.numeric(res_trace), as.numeric(res_plain))
})

test_that("track_trace = TRUE also populates result$trace when verbose = TRUE", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1, track_trace = TRUE, verbose = TRUE)
  expect_s3_class(res$trace, "data.frame")
})

test_that("factr and pgtol are accepted and influence convergence speed", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res_default <- optim_hp(hp, data, 0, kern, 1, verbose = TRUE)
  res_loose   <- optim_hp(hp, data, 0, kern, 1, factr = 1e12, pgtol = 1e-2, verbose = TRUE)
  expect_true(res_loose$counts[1] <= res_default$counts[1])
})


test_that("optim_hp attaches convergence and value as attributes by default", {
  data <- make_data()
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  res  <- optim_hp(hp, data, 0, kern, 1)
  full <- optim_hp(hp, data, 0, kern, 1, verbose = TRUE)
  expect_equal(attr(res, "convergence"), full$convergence)
  expect_equal(attr(res, "value"), full$value)
})

test_that("optim_hp warns when max_iter is too low to converge", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel(hp = c(0.01, 0.01))
  hp   <- unlist(keRnel::gt_HPs(kern))
  expect_warning(
    optim_hp(hp, data, 0, kern, 1, max_iter = 1),
    "did not converge"
  )
})

test_that("optim_hp respects max_iter (stops early instead of fully converging)", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel(hp = c(0.01, 0.01))
  hp   <- unlist(keRnel::gt_HPs(kern))
  res_capped <- suppressWarnings(
    optim_hp(hp, data, 0, kern, 1, max_iter = 1, verbose = TRUE)
  )
  res_full <- optim_hp(hp, data, 0, kern, 1, verbose = TRUE)
  expect_equal(res_capped$convergence, 1)  # 1 = hit the iteration limit (see ?optim)
  expect_true(res_capped$value >= res_full$value)
})

test_that("optim_hp actually reduces objective vs initial hp", {
  set.seed(5)
  data <- make_data(nb_id = 15)
  kern <- make_kernel()
  hp   <- unlist(keRnel::gt_HPs(kern))
  v0   <- BayesOmics:::sum_logGaussian(hp, data, 0, kern, 1, 1e-6)
  hp2  <- optim_hp(hp, data, 0, kern, 1)
  v1   <- BayesOmics:::sum_logGaussian(hp2, data, 0, kern, 1, 1e-6)
  expect_lte(v1, v0)
})
