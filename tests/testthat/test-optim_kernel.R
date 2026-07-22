# -- chol_inv_jitter -----------------------------------------------------------

test_that("chol_inv_jitter returns inverse for PD matrix", {
  mat <- matrix(c(4, 2, 2, 3), nrow = 2)
  inv <- BayesOmics:::chol_inv_jitter(mat, pen_diag = 1e-6)
  expect_true(is.matrix(inv))
  expect_equal(dim(inv), c(2, 2))
  # A %*% A^-1 ~ I
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

# -- dmnorm --------------------------------------------------------------------

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

# -- resolve_prior_cov ---------------------------------------------------------

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

# -- cached_cov_inv -------------------------------------------------------------

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

# -- sum_logGaussian -----------------------------------------------------------

test_that("sum_logGaussian returns a finite scalar", {
  data <- make_data()
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  val  <- BayesOmics:::sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  expect_true(is.numeric(val))
  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("sum_logGaussian decreases when kernel fits data better", {
  set.seed(1)
  data <- make_data(nb_id = 20)
  kern <- make_kernel()
  free0 <- keRnel::get_free_params(kern)

  v_init  <- BayesOmics:::sum_logGaussian(free0, data, 0, kern, 1, 1e-6)
  hp_opt  <- optim_hp(kern, data, 0, 1)
  kern_opt <- do.call(keRnel::kupdate, c(list(kern), as.list(hp_opt)))
  free_opt <- keRnel::get_free_params(kern_opt)
  v_opt   <- BayesOmics:::sum_logGaussian(free_opt, data, 0, kern, 1, 1e-6)
  expect_lte(v_opt, v_init)
})

test_that("sum_logGaussian is larger (more negative) for mismatched prior_mean", {
  data <- make_data()
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  val_right  <- BayesOmics:::sum_logGaussian(free, data, mean(data$Output), kern, 1, 1e-6)
  val_wrong  <- BayesOmics:::sum_logGaussian(free, data, mean(data$Output) + 1000, kern, 1, 1e-6)
  expect_gt(val_wrong, val_right)
})

test_that("sum_logGaussian gives the same result with and without a cache", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  no_cache <- BayesOmics:::sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  cache <- new.env()
  with_cache <- BayesOmics:::sum_logGaussian(free, data, 0, kern, 1, 1e-6, cache = cache)
  expect_equal(with_cache, no_cache)
})

test_that("gr_sum_logGaussian gives the same result with and without a cache", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  no_cache <- BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  cache <- new.env()
  with_cache <- BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6, cache = cache)
  expect_equal(with_cache, no_cache)
})

test_that("a shared cache makes sum_logGaussian then gr_sum_logGaussian reuse the same cov inverse", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  cache <- new.env()
  BayesOmics:::sum_logGaussian(free, data, 0, kern, 1, 1e-6, cache = cache)
  expect_identical(cache$free, free)
  cached_inv <- cache$inv
  BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6, cache = cache)
  expect_identical(cache$inv, cached_inv)
})

test_that("sum_logGaussian accepts a full prior_cov matrix (not just a scalar)", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  p    <- length(unique(data$Input))  # kernel matrix is p x p (unique inputs)
  prior_cov_mat <- diag(p) * 1.5
  val_scalar <- BayesOmics:::sum_logGaussian(free, data, 0, kern, 1.5, 1e-6)
  val_matrix <- BayesOmics:::sum_logGaussian(free, data, 0, kern, prior_cov_mat, 1e-6)
  expect_equal(val_matrix, val_scalar, tolerance = 1e-8)
})

test_that("sum_logGaussian errors when prior_cov matrix has the wrong size", {
  data <- make_data(nb_id = 6)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  bad_prior_cov <- diag(2)  # wrong size vs nrow(data)
  expect_error(BayesOmics:::sum_logGaussian(free, data, 0, kern, bad_prior_cov, 1e-6),
               "prior_cov must be a square matrix")
})

# -- gr_sum_logGaussian --------------------------------------------------------

test_that("gr_sum_logGaussian returns numeric vector matching the number of free params", {
  data <- make_data()
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  g    <- BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  expect_true(is.numeric(g))
  expect_length(g, length(free))
})

test_that("gr_sum_logGaussian gradient is finite at reasonable hp", {
  data <- make_data()
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  g    <- BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  expect_true(all(is.finite(g)))
})

test_that("gr_sum_logGaussian numerical gradient matches analytic gradient", {
  set.seed(2)
  data <- make_data(nb_id = 8)
  kern <- make_kernel(hp = c(2.0, 1.5))
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic <- BayesOmics:::gr_sum_logGaussian(free, data, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, data, 0, kern, 1, 1e-6) -
     BayesOmics:::sum_logGaussian(free_dn, data, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

# -- optim_hp --------------------------------------------------------

test_that("optim_hp returns named numeric vector by default", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1)
  expect_true(is.numeric(res))
  expect_length(res, length(keRnel::get_trainable_params(kern)))
})

test_that("optim_hp returns full optim list when verbose=TRUE", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1, verbose = TRUE)
  expect_true(is.list(res))
  expect_true("par" %in% names(res))
  expect_true("convergence" %in% names(res))
  expect_true(inherits(res$kern, "kernel"))
})

test_that("optim_hp result keeps hyperparameters positive", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1)
  expect_true(all(res > 0))
})

test_that("optim_hp errors on wrong kern type", {
  data <- make_data()
  expect_error(optim_hp(list(a = 1), data, 0, 1),
               "keRnel")
})

test_that("optim_hp errors when db lacks Input/Output", {
  kern <- make_kernel()
  expect_error(optim_hp(kern, data.frame(x = 1:5), 0, 1),
               "Input.*Output|Output.*Input")
})

test_that("optim_hp errors when db$Input contains NaN", {
  data <- make_data()
  kern <- make_kernel()
  data$Input[1] <- NaN
  expect_error(optim_hp(kern, data, 0, 1), "NaN|Inf|NA")
})

test_that("optim_hp errors when db$Output contains Inf", {
  data <- make_data()
  kern <- make_kernel()
  data$Output[1] <- Inf
  expect_error(optim_hp(kern, data, 0, 1), "NaN|Inf|NA")
})

test_that("optim_hp errors when db$Output contains NA", {
  data <- make_data()
  kern <- make_kernel()
  data$Output[1] <- NA_real_
  expect_error(optim_hp(kern, data, 0, 1), "NaN|Inf|NA")
})

test_that("optim_hp errors when prior_mean has wrong length", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(optim_hp(kern, data, c(1, 2, 3), 1),
               "prior_mean")
})

test_that("track_trace = FALSE (default) attaches no trace attribute", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1)
  expect_null(attr(res, "trace"))
})

test_that("track_trace = TRUE attaches a trace data.frame with fn and gr rows", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1, track_trace = TRUE)
  tr   <- attr(res, "trace")
  expect_s3_class(tr, "data.frame")
  expect_true(all(c("eval_type", "eval_index", "value", "elapsed_sec") %in% names(tr)))
  expect_true("fn" %in% tr$eval_type)
  expect_true("gr" %in% tr$eval_type)
})

test_that("track_trace = TRUE does not change the optimized hyperparameters", {
  data <- make_data()
  kern <- make_kernel()
  res_plain <- optim_hp(kern, data, 0, 1)
  res_trace <- optim_hp(kern, data, 0, 1, track_trace = TRUE)
  expect_equal(as.numeric(res_trace), as.numeric(res_plain))
})

test_that("track_trace = TRUE also populates result$trace when verbose = TRUE", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1, track_trace = TRUE, verbose = TRUE)
  expect_s3_class(res$trace, "data.frame")
})

test_that("factr and pgtol are accepted and influence convergence speed", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel()
  res_default <- optim_hp(kern, data, 0, 1, verbose = TRUE)
  res_loose   <- optim_hp(kern, data, 0, 1, factr = 1e12, pgtol = 1e-2, verbose = TRUE)
  expect_true(res_loose$counts[1] <= res_default$counts[1])
})


test_that("optim_hp attaches convergence and value as attributes by default", {
  data <- make_data()
  kern <- make_kernel()
  res  <- optim_hp(kern, data, 0, 1)
  full <- optim_hp(kern, data, 0, 1, verbose = TRUE)
  expect_equal(attr(res, "convergence"), full$convergence)
  expect_equal(attr(res, "value"), full$value)
})

test_that("optim_hp warns when max_iter is too low to converge", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel(hp = c(0.01, 0.01))
  expect_warning(
    optim_hp(kern, data, 0, 1, max_iter = 1),
    "did not converge"
  )
})

test_that("optim_hp respects max_iter (stops early instead of fully converging)", {
  data <- make_data(nb_id = 20)
  kern <- make_kernel(hp = c(0.01, 0.01))
  res_capped <- suppressWarnings(
    optim_hp(kern, data, 0, 1, max_iter = 1, verbose = TRUE)
  )
  res_full <- optim_hp(kern, data, 0, 1, verbose = TRUE)
  expect_equal(res_capped$convergence, 1)  # 1 = hit the iteration limit (see ?optim)
  expect_true(res_capped$value >= res_full$value)
})

test_that("optim_hp actually reduces objective vs initial hp", {
  set.seed(5)
  data <- make_data(nb_id = 15)
  kern <- make_kernel()
  free0 <- keRnel::get_free_params(kern)
  v0   <- BayesOmics:::sum_logGaussian(free0, data, 0, kern, 1, 1e-6)
  hp2  <- optim_hp(kern, data, 0, 1)
  kern2 <- do.call(keRnel::kupdate, c(list(kern), as.list(hp2)))
  v1   <- BayesOmics:::sum_logGaussian(keRnel::get_free_params(kern2), data, 0, kern, 1, 1e-6)
  expect_lte(v1, v0)
})

# -- sum_logGaussian / gr_sum_logGaussian: cas replique (nb_sample > 1) -------

test_that("sum_logGaussian triggers the replicated branch when nb_sample > 1", {
  data_rep <- make_data(nb_id = 5, nb_sample = 3)
  data_rep <- data_rep[data_rep$Group == "1", ]
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  expect_true(nrow(data_rep) > length(unique(data_rep$Input)))
  expect_true("Sample" %in% names(data_rep))
  val <- BayesOmics:::sum_logGaussian(free, data_rep, 0, kern, 1, 1e-6)
  expect_true(is.finite(val))
})

test_that("gr_sum_logGaussian analytic gradient matches numerical gradient (replicated case)", {
  set.seed(3)
  data_rep <- make_data(nb_id = 8, nb_sample = 4)
  data_rep <- data_rep[data_rep$Group == "1", ]
  kern <- make_kernel(hp = c(1.5, 2.0))
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic <- BayesOmics:::gr_sum_logGaussian(free, data_rep, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, data_rep, 0, kern, 1, 1e-6) -
       BayesOmics:::sum_logGaussian(free_dn, data_rep, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

test_that("sum_logGaussian replicated with 1 sample gives finite result (non-replicated branch)", {
  data_1 <- make_data(nb_id = 5, nb_sample = 1)
  data_1 <- data_1[data_1$Group == "1", ]
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  expect_equal(nrow(data_1), length(unique(data_1$Input)))
  val <- BayesOmics:::sum_logGaussian(free, data_1, 0, kern, 1, 1e-6)
  expect_true(is.finite(val))
})

test_that("optim_hp recovers HPs for replicated data generated with known kernel", {
  skip_on_cran()
  set.seed(11)
  kern_true <- keRnel::variance_kernel(variance = 3) * keRnel::se_kernel(length_scale = 5)
  data_rep  <- simu_db_kernel(kernel = kern_true, nb_id = 10, nb_group = 1,
                               nb_sample = 8, diff_group = 0, var_sample = 1e-4)
  data_rep  <- data_rep[data_rep$Group == "1", ]
  kern_fit  <- make_kernel()
  hp_opt    <- optim_hp(kern_fit, data_rep, mean(data_rep$Output), 1e-4)
  expect_equal(unname(hp_opt["variance"]),     3, tolerance = 0.8)
  expect_equal(unname(hp_opt["length_scale"]), 5, tolerance = 2.0)
})

# -- Multi-groupe : cas limites (benchmark 07) ----------------------------------

test_that("sum_logGaussian retourne un scalaire fini pour donnees multi-groupes", {
  set.seed(20)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_mg   <- simu_db_kernel(nb_id = 6, nb_group = 3, nb_sample = 4, diff_group = 0,
                               var_sample = 1, kernel = kern_true)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  val  <- BayesOmics:::sum_logGaussian(free, data_mg, 0, kern, 1, 1e-6)
  expect_length(val, 1)
  expect_true(is.finite(val))
})

test_that("gr_sum_logGaussian gradient analytique = gradient numerique (multi-groupe, equilibre)", {
  set.seed(21)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_mg   <- simu_db_kernel(nb_id = 6, nb_group = 3, nb_sample = 3, diff_group = 0,
                               var_sample = 1, kernel = kern_true)
  kern <- make_kernel(hp = c(1.5, 2.0))
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic     <- BayesOmics:::gr_sum_logGaussian(free, data_mg, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, data_mg, 0, kern, 1, 1e-6) -
     BayesOmics:::sum_logGaussian(free_dn, data_mg, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

test_that("gr_sum_logGaussian gradient analytique = gradient numerique (design desequilibre)", {
  set.seed(22)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_unbal <- simu_db_kernel(nb_id = 6, nb_group = 3, nb_sample = c(2, 5, 8),
                                diff_group = 0, var_sample = 1, kernel = kern_true)
  kern <- make_kernel(hp = c(1.5, 2.0))
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic     <- BayesOmics:::gr_sum_logGaussian(free, data_unbal, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, data_unbal, 0, kern, 1, 1e-6) -
     BayesOmics:::sum_logGaussian(free_dn, data_unbal, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

test_that("sum_logGaussian branche repliquee declenche pour nb_group > 1, nb_sample = 1", {
  set.seed(23)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_g5s1 <- simu_db_kernel(nb_id = 5, nb_group = 5, nb_sample = 1, diff_group = 0,
                               var_sample = 1, kernel = kern_true)
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  # nrow > p_unique -> branche repliquee avec 5 replicats (un par groupe)
  expect_true(nrow(data_g5s1) > length(unique(data_g5s1$Input)))
  val <- BayesOmics:::sum_logGaussian(free, data_g5s1, 0, kern, 1, 1e-6)
  expect_true(is.finite(val))
})

test_that("optim_hp converge sur donnees multi-groupes (diff_group = 0)", {
  set.seed(24)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_mg   <- simu_db_kernel(nb_id = 8, nb_group = 3, nb_sample = 4, diff_group = 0,
                               var_sample = 1, kernel = kern_true)
  kern <- make_kernel()
  res  <- optim_hp(kern, data_mg, prior_mean = 0, prior_cov = 1, verbose = TRUE)
  expect_equal(res$convergence, 0)
  expect_true(all(is.finite(res$par)))
})

test_that("doubler les groupes identiques (diff_group=0) double la NLL a hp fixes", {
  # Propriete de la vraisemblance repliquee : si on duplique tous les groupes
  # (meme donnees, double nb_group), la NLL doit etre ~ 2x la NLL d'origine.
  set.seed(25)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data1 <- simu_db_kernel(nb_id = 5, nb_group = 2, nb_sample = 3, diff_group = 0,
                           var_sample = 1, kernel = kern_true)
  # Duplique les groupes : 4 groupes avec les memes donnees que les 2 d'origine
  data2 <- rbind(data1, transform(data1, Group = paste0(data1$Group, "_b")))
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  nll1 <- BayesOmics:::sum_logGaussian(free, data1, 0, kern, 1, 1e-6)
  nll2 <- BayesOmics:::sum_logGaussian(free, data2, 0, kern, 1, 1e-6)
  expect_equal(nll2, 2 * nll1, tolerance = 1e-6)
})

# -- Edge cases from optim_exploration/ -----------------------------------

test_that("chol_inv_jitter with pen_diag = 0 bumps to 1e-6 only once jitter is needed", {
  # jitter_until_pd() must not multiply 0 by 10 forever: a singular (rank-1)
  # matrix should still be rescued via the documented "bump to 1e-6 once
  # jitter is actually needed" special case.
  singular_mat <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_no_error(BayesOmics:::chol_inv_jitter(singular_mat, pen_diag = 0))
})

test_that("optim_hp recovers hyperparameters on every well-behaved kernel family", {
  # Smoke test only (runs without error, gives finite positive output) --
  # not a tight RMSE tolerance, see optim_exploration/07_kernel_family_zoo.R
  # for the full Monte-Carlo recovery study. All families listed here are
  # supported by keRnel::kernel_grad_free() (including linear_kernel-family
  # and periodic_kernel, previously excluded due to a keRnel bug that has
  # since been fixed).
  zoo <- list(
    SE                = keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 1.5),
    Matern12          = keRnel::variance_kernel(variance = 2) * keRnel::matern12_kernel(length_scale = 1.5),
    Matern32          = keRnel::variance_kernel(variance = 2) * keRnel::matern32_kernel(length_scale = 1.5),
    Matern52          = keRnel::variance_kernel(variance = 2) * keRnel::matern52_kernel(length_scale = 1.5),
    RationalQuadratic = keRnel::variance_kernel(variance = 2) *
      keRnel::rational_quadratic_kernel(length_scale = 1.5, alpha = 1),
    Constant          = keRnel::constant_kernel(value = 2),
    Noise             = keRnel::white_noise_kernel(noise = 1),
    Linear            = keRnel::affine_kernel(slope_var = 1, offset = 2),
    Periodic          = keRnel::variance_kernel(variance = 2) *
      keRnel::periodic_kernel(length_scale = 1.5, period = 4)
  )
  for (kname in names(zoo)) {
    set.seed(50)
    kern <- zoo[[kname]]
    data <- simu_db_kernel(nb_id = 8, nb_group = 1, nb_sample = 5, diff_group = 0,
                            var_sample = 0.5, kernel = kern, range_input = c(0, 10))
    data <- data[data$Group == "1", ]
    res <- optim_hp(kern, data, mean(data$Output), prior_cov = 0.5, pen_diag = 1e-3)
    expect_true(all(is.finite(res)), info = kname)
  }
})

test_that("sum_logGaussian falls back to the non-replicated path when replicate counts are unequal", {
  # is_replicated is gated on every (Group, Sample) combination containing
  # exactly one row per unique Input position; a ragged design (one
  # replicate missing an id) must NOT be misinterpreted as replicated.
  set.seed(52)
  kern_true <- keRnel::variance_kernel(variance = 4) * keRnel::se_kernel(length_scale = 15)
  data <- simu_db_kernel(nb_id = 8, nb_group = 1, nb_sample = 4, diff_group = 0,
                         var_sample = 0.5, kernel = kern_true, range_input = c(0, 50))
  data <- data[data$Group == "1", ]
  drop_idx <- which(data$Sample == data$Sample[1])[1]
  data_unequal <- data[-drop_idx, ]

  p_uniq <- length(unique(data_unequal$Input))
  sample_counts <- table(data_unequal$Sample)
  expect_false(all(sample_counts == p_uniq))

  free <- keRnel::get_free_params(kern_true)
  val <- BayesOmics:::sum_logGaussian(free, data_unequal,
                                       mean(data_unequal$Output), kern_true, 0.5, 1e-6)
  expect_true(is.finite(val))
})

# -- Multi-dimensional Input (D > 1) --------------------------------------

test_that("sum_logGaussian handles a genuinely multi-dimensional Input (D = 2)", {
  set.seed(60)
  ids <- paste0("ID_", 1:6)
  base <- data.frame(ID = ids, dim1 = seq(0, 10, length.out = 6), dim2 = seq(5, 0, length.out = 6))
  long <- do.call(rbind, lapply(seq_len(3), function(s) {
    y <- rnorm(6, mean = base$dim1)
    rbind(
      data.frame(ID = ids, Sample = s, Input_ID = 1, Input = base$dim1, Output = y),
      data.frame(ID = ids, Sample = s, Input_ID = 2, Input = base$dim2, Output = y)
    )
  }))
  kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 2)
  free <- keRnel::get_free_params(kern)
  val <- BayesOmics:::sum_logGaussian(free, long, 0, kern, 1, 1e-6)
  expect_true(is.finite(val))
})

test_that("sum_logGaussian errors on multi-dimensional Input without an ID column", {
  db <- data.frame(Input_ID = rep(1:2, 3), Input = runif(6), Output = rnorm(6))
  kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
  free <- keRnel::get_free_params(kern)
  expect_error(BayesOmics:::sum_logGaussian(free, db, 0, kern, 1, 1e-6), "ID")
})

test_that("gr_sum_logGaussian analytic gradient matches numerical gradient (D = 2)", {
  set.seed(61)
  ids <- paste0("ID_", 1:6)
  base <- data.frame(ID = ids, dim1 = seq(0, 10, length.out = 6), dim2 = seq(5, 0, length.out = 6))
  long <- do.call(rbind, lapply(seq_len(3), function(s) {
    y <- rnorm(6, mean = base$dim1)
    rbind(
      data.frame(ID = ids, Sample = s, Input_ID = 1, Input = base$dim1, Output = y),
      data.frame(ID = ids, Sample = s, Input_ID = 2, Input = base$dim2, Output = y)
    )
  }))
  kern <- keRnel::variance_kernel(variance = 1.5) * keRnel::se_kernel(length_scale = 2.5)
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic <- BayesOmics:::gr_sum_logGaussian(free, long, 0, kern, 1, 1e-6)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, long, 0, kern, 1, 1e-6) -
     BayesOmics:::sum_logGaussian(free_dn, long, 0, kern, 1, 1e-6)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)
})

test_that("optim_hp fits a genuinely multi-dimensional (D = 2) isotropic kernel", {
  skip_on_cran()
  set.seed(62)
  ids <- paste0("ID_", 1:10)
  base <- data.frame(ID = ids, dim1 = runif(10, 0, 10), dim2 = runif(10, 0, 10))
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 3)
  input_mat <- as.matrix(base[, c("dim1", "dim2")])
  cov <- keRnel::evaluate(kern_true, input_mat, input_mat) + diag(1e-4, 10)
  L   <- t(chol(cov))
  long <- do.call(rbind, lapply(seq_len(6), function(s) {
    y <- as.vector(L %*% rnorm(10))
    rbind(
      data.frame(ID = ids, Sample = s, Input_ID = 1, Input = base$dim1, Output = y),
      data.frame(ID = ids, Sample = s, Input_ID = 2, Input = base$dim2, Output = y)
    )
  }))
  kern_fit <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
  res <- optim_hp(kern_fit, long, mean(long$Output), 1e-4)
  expect_true(all(is.finite(res)))
  expect_equal(unname(res["variance"]), 2, tolerance = 1.5)
  expect_equal(unname(res["length_scale"]), 3, tolerance = 2.0)
})

test_that("optim_hp fits two sibling sub-kernels sharing a bare hyperparameter name independently", {
  # The direct validation of the free-space-gradient design: kupdate() cannot
  # give se_kernel(length_scale=1) + se_kernel(length_scale=3) two different
  # starting values that stay independent during optimization (it updates
  # every occurrence of a bare name at once), but the free-space path
  # (get_free_params()/set_free_params()/kernel_grad_free()) is positional
  # and therefore unaffected by the name collision.
  skip_on_cran()
  set.seed(63)
  data <- make_data(nb_id = 25, nb_sample = 6)
  data <- data[data$Group == "1", ]
  kern <- keRnel::se_kernel(length_scale = 1) + keRnel::se_kernel(length_scale = 8)
  res  <- optim_hp(kern, data, mean(data$Output), 1)
  expect_length(res, 2)
  expect_true(all(is.finite(res)))
  expect_gt(abs(res[1] - res[2]), 1e-6)
})

# -- Group pooling / REML correction (group_col, "Solution S2") ----------------

test_that("demean_by_group centers each group on its own mean, for any G", {
  output <- c(1, 2, 3, 10, 12, 14, 100, 106)
  group  <- c("a", "a", "a", "b", "b", "b", "c", "c")
  res <- BayesOmics:::demean_by_group(output, group)
  expect_equal(res$n_groups, 3)
  expect_equal(as.numeric(tapply(res$output, group, mean)), c(0, 0, 0), tolerance = 1e-10)
  # No dim attribute leaks through (the tapply-array footgun documented in
  # dev/optim_exploration/NOTES_math.md and 10_group_mean_shift_confound.R).
  expect_null(dim(res$output))
  expect_true(is.vector(res$output))
})

test_that("optim_hp requires prior_mean unless group_col is supplied", {
  data <- make_data(nb_id = 5, nb_group = 2, nb_sample = 3)
  kern <- make_kernel()
  expect_error(optim_hp(kern, data, prior_cov = 1), "prior_mean.*required")
})

test_that("optim_hp errors when group_col is not a column of db", {
  data <- make_data(nb_id = 5, nb_group = 2, nb_sample = 3)
  kern <- make_kernel()
  expect_error(
    optim_hp(kern, data, prior_cov = 1, group_col = "NotAColumn"),
    "not a column"
  )
})

test_that("optim_hp warns when a nonzero prior_mean is supplied alongside group_col", {
  skip_on_cran()
  set.seed(70)
  data <- make_data(nb_id = 5, nb_group = 2, nb_sample = 3)
  kern <- make_kernel()
  expect_warning(
    optim_hp(kern, data, prior_mean = 5, prior_cov = 1, group_col = "Group"),
    "ignored"
  )
  # prior_mean = 0 alongside group_col is the documented no-op case: no warning.
  expect_warning(
    optim_hp(kern, data, prior_mean = 0, prior_cov = 1, group_col = "Group"),
    NA
  )
})

test_that("optim_hp(group_col=...) internally demeans and applies n_groups (matches the manual equivalent)", {
  skip_on_cran()
  set.seed(71)
  data <- make_data(nb_id = 6, nb_group = 3, nb_sample = 4)
  kern <- make_kernel(hp = c(1.5, 2.0))

  res_grouped <- optim_hp(kern, data, prior_cov = 1, group_col = "Group", verbose = TRUE)

  # Re-evaluate the NLL that optim_hp() must have minimized, by reproducing
  # its two documented steps by hand: demean_by_group() then n_groups = G.
  demeaned    <- BayesOmics:::demean_by_group(data$Output, data$Group)
  data_manual <- data
  data_manual$Output <- demeaned$output
  free_at_opt <- keRnel::get_free_params(res_grouped$kern)
  nll_manual  <- BayesOmics:::sum_logGaussian(free_at_opt, data_manual, 0, kern, 1, 1e-6,
                                              n_groups = demeaned$n_groups)
  expect_equal(nll_manual, res_grouped$value, tolerance = 1e-8)
  expect_equal(demeaned$n_groups, 3)
})

test_that("sum_logGaussian errors when n_groups is supplied for a non-replicated design", {
  set.seed(72)
  # A single group, single sample -- the genuinely non-replicated branch (see
  # "sum_logGaussian replicated with 1 sample gives finite result" above:
  # nb_group > 1 with nb_sample = 1 still hits the *replicated* branch, since
  # each group's own sample already counts as one of >1 Group x Sample
  # replicates -- only nb_group = 1, nb_sample = 1 is truly non-replicated).
  data_1 <- make_data(nb_id = 5, nb_group = 2, nb_sample = 1)
  data_1 <- data_1[data_1$Group == "1", ]
  kern <- make_kernel()
  free <- keRnel::get_free_params(kern)
  expect_error(
    BayesOmics:::sum_logGaussian(free, data_1, 0, kern, 1, 1e-6, n_groups = 2),
    "balanced replicated design"
  )
  expect_error(
    BayesOmics:::gr_sum_logGaussian(free, data_1, 0, kern, 1, 1e-6, n_groups = 2),
    "balanced replicated design"
  )
})

test_that("gr_sum_logGaussian REML correction: analytic gradient matches numerical gradient", {
  set.seed(73)
  kern_true <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 2)
  data_mg   <- simu_db_kernel(nb_id = 6, nb_group = 4, nb_sample = 3, diff_group = 5,
                               var_sample = 1, kernel = kern_true)
  demeaned  <- BayesOmics:::demean_by_group(data_mg$Output, data_mg$Group)
  data_mg$Output <- demeaned$output

  kern <- make_kernel(hp = c(1.5, 2.0))
  free <- keRnel::get_free_params(kern)
  eps  <- 1e-5
  analytic <- BayesOmics:::gr_sum_logGaussian(free, data_mg, 0, kern, 1, 1e-6, n_groups = 4)
  numeric_grad <- vapply(seq_along(free), function(i) {
    free_up <- free; free_up[i] <- free[i] + eps
    free_dn <- free; free_dn[i] <- free[i] - eps
    (BayesOmics:::sum_logGaussian(free_up, data_mg, 0, kern, 1, 1e-6, n_groups = 4) -
     BayesOmics:::sum_logGaussian(free_dn, data_mg, 0, kern, 1, 1e-6, n_groups = 4)) / (2 * eps)
  }, numeric(1))
  expect_equal(unname(analytic), numeric_grad, tolerance = 1e-3)

  # The correction should have a real (non-negligible) effect, not be a
  # rounding-level no-op -- see dev/optim_exploration/11_reml_group_gradient_check.R
  # for the full cross-kernel/cross-G sweep of this same comparison.
  analytic_no_reml <- BayesOmics:::gr_sum_logGaussian(free, data_mg, 0, kern, 1, 1e-6, n_groups = NULL)
  expect_gt(max(abs(analytic - analytic_no_reml)), 1e-2)
})

test_that("optim_hp(group_col=...) recovers kernel variance closer to truth than naive pooling under a group mean shift", {
  skip_on_cran()
  set.seed(74)
  kern_true <- keRnel::variance_kernel(variance = 3) * keRnel::se_kernel(length_scale = 20)
  data <- simu_db_kernel(kernel = kern_true, nb_id = 15, nb_group = 3, nb_sample = 8,
                          diff_group = 10, var_sample = 1, mu_0 = 40,
                          range_input = c(0, 150), integer_input = TRUE)

  kern_fit <- make_kernel()
  hp_naive   <- optim_hp(kern_fit, data, prior_mean = mean(data$Output), prior_cov = 1e-6)
  hp_grouped <- optim_hp(kern_fit, data, prior_cov = 1e-6, group_col = "Group")

  expect_true(abs(unname(hp_grouped["variance"]) - 3) < abs(unname(hp_naive["variance"]) - 3))
})
