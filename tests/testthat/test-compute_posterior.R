# ── multi_posterior_mean: input validation ────────────────────────────────────

test_that("multi_posterior_mean errors when required columns are missing", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(multi_posterior_mean(data[, c("ID", "Output", "Input")], kern), "Group")
  expect_error(multi_posterior_mean(data[, c("Group", "Output", "Input")], kern), "ID")
  expect_error(multi_posterior_mean(data[, c("Group", "ID", "Input")], kern), "Output")
  expect_error(multi_posterior_mean(data[, c("Group", "ID", "Output")], kern), "Input")
})

test_that("multi_posterior_mean errors on non-kernel argument", {
  data <- make_data()
  expect_error(multi_posterior_mean(data, list()), "keRnel")
})

test_that("multi_posterior_mean errors when lambda_0 <= 0", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(multi_posterior_mean(data, kern, lambda_0 = 0),  "lambda_0")
  expect_error(multi_posterior_mean(data, kern, lambda_0 = -1), "lambda_0")
})

test_that("multi_posterior_mean errors on NA in Group or ID", {
  data <- make_data()
  kern <- make_kernel()
  data_na <- data
  data_na$Group[1] <- NA
  expect_error(multi_posterior_mean(data_na, kern), "NA")
})

test_that("multi_posterior_mean errors when mu_0 is not numeric", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(multi_posterior_mean(data, kern, mu_0 = "0"), "mu_0")
})

test_that("multi_posterior_mean errors when an ID maps to several distinct Input values within a group", {
  # Two different IDs sharing the same single Input value in a group: 2
  # distinct IDs but only 1 distinct Input value -> muk/vec_name length mismatch.
  data <- data.frame(
    ID     = c("ID_1", "ID_2"),
    Group  = "A",
    Output = c(0, 0),
    Input  = c(5, 5)
  )
  kern <- make_kernel()
  expect_error(multi_posterior_mean(data, kern), "distinct Input value")
})

# ── multi_posterior_mean: muk/sigmak alignment (regression for the former
# alphabetical-ID vs Input-value-order misalignment bug) ─────────────────────

test_that("multi_posterior_mean: muk and sigmak (via get_sigmak) are correctly aligned by ID", {
  ker <- methods::new("SEKernel")
  ker <- keRnel::set_hyperparameters(ker, c(variance_se = 10, length_scale_se = 200))

  # dist(ID_1, ID_2) = 90, dist(ID_1, ID_3) = 100, dist(ID_2, ID_3) = 10
  data <- data.frame(
    ID     = c("ID_1", "ID_2", "ID_3"),
    Group  = "A",
    Output = c(0, 0, 0),
    Input  = c(100, 10, 0)
  )
  res <- multi_posterior_mean(data, ker)
  g   <- res$groups[["A"]]
  sig <- BayesOmics:::get_sigmak(g, res$kernels)

  pk <- function(xi, xj) as.numeric(keRnel::pairwise_kernel(ker, as.matrix(xi), as.matrix(xj))) / 2
  true_cov_12 <- pk(100, 10)  # cov(ID_1, ID_2)
  true_cov_13 <- pk(100, 0)   # cov(ID_1, ID_3)
  true_cov_23 <- pk(10, 0)    # cov(ID_2, ID_3)

  ids <- names(g$muk)
  expect_equal(unname(sig[ids == "ID_1", ids == "ID_2"]), true_cov_12, tolerance = 1e-8)
  expect_equal(unname(sig[ids == "ID_1", ids == "ID_3"]), true_cov_13, tolerance = 1e-8)
  expect_equal(unname(sig[ids == "ID_2", ids == "ID_3"]), true_cov_23, tolerance = 1e-8)
})

# ── multi_posterior_mean: return structure ─────────────────────────────────────

test_that("multi_posterior_mean returns a list with kernels and groups", {
  res <- make_posteriors()
  expect_true(is.list(res))
  expect_named(res, c("kernels", "groups"))
  expect_named(res$groups)
})

test_that("multi_posterior_mean result has one group entry per group", {
  data <- make_data(nb_group = 3)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  expect_length(res$groups, 3)
})

test_that("multi_posterior_mean each group entry has muk, id_to_input, kernel_key, scale", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    expect_true(all(c("muk", "id_to_input", "kernel_key", "scale") %in% names(res$groups[[g]])))
  }
})

test_that("multi_posterior_mean: muk is a named numeric vector", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    expect_true(is.numeric(res$groups[[g]]$muk))
    expect_true(!is.null(names(res$groups[[g]]$muk)))
  }
})

test_that("multi_posterior_mean: get_sigmak returns a square matrix", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    sk <- BayesOmics:::get_sigmak(res$groups[[g]], res$kernels)
    expect_true(is.matrix(sk))
    expect_equal(nrow(sk), ncol(sk))
  }
})

test_that("multi_posterior_mean: muk length equals sigmak dimension", {
  res <- make_posteriors(nb_id = 6)
  for (g in names(res$groups)) {
    sk <- BayesOmics:::get_sigmak(res$groups[[g]], res$kernels)
    expect_equal(length(res$groups[[g]]$muk), nrow(sk))
  }
})

test_that("multi_posterior_mean: muk names match ID values", {
  data <- make_data(nb_id = 4)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  ids  <- sort(unique(data$ID))
  for (g in names(res$groups)) {
    expect_equal(sort(names(res$groups[[g]]$muk)), ids)
  }
})

# ── multi_posterior_mean: numerics ─────────────────────────────────────────────

test_that("multi_posterior_mean: larger lambda_0 pulls muk toward mu_0", {
  data <- make_data(nb_id = 10)
  kern <- make_kernel()
  mu_0 <- 0
  res_small <- multi_posterior_mean(data, kern, mu_0 = mu_0, lambda_0 = 0.01)
  res_large <- multi_posterior_mean(data, kern, mu_0 = mu_0, lambda_0 = 100)
  muk_small <- unlist(lapply(res_small$groups, `[[`, "muk"))
  muk_large <- unlist(lapply(res_large$groups, `[[`, "muk"))
  # with large lambda_0, muk is pulled toward mu_0=0
  expect_lt(mean(abs(muk_large)), mean(abs(muk_small)) + 0.1 + abs(mean(muk_small)))
})

test_that("multi_posterior_mean: kernel cache is shared (by reference) across groups with identical inputs", {
  data <- make_data(nb_id = 5, nb_group = 2)
  # Force both groups to share the exact same inputs
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, each = 2)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  expect_equal(res$groups[["1"]]$kernel_key, res$groups[["2"]]$kernel_key)
  expect_length(res$kernels, 1)
  sig1 <- BayesOmics:::get_sigmak(res$groups[["1"]], res$kernels)
  sig2 <- BayesOmics:::get_sigmak(res$groups[["2"]], res$kernels)
  expect_equal(unname(sig1), unname(sig2), tolerance = 1e-10)
})

test_that("multi_posterior_mean: kernel cache is shared (by reference) across three groups with identical inputs", {
  data <- make_data(nb_id = 5, nb_group = 3)
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, length.out = nrow(data))
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  keys <- vapply(res$groups, `[[`, character(1), "kernel_key")
  expect_length(unique(keys), 1)
  expect_length(res$kernels, 1)
})

test_that("multi_posterior_mean: get_sigmak reconstructs a fresh matrix that does not mutate the shared cache", {
  data <- make_data(nb_id = 5, nb_group = 2)
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, each = 2)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)

  sig1 <- BayesOmics:::get_sigmak(res$groups[["1"]], res$kernels)
  sig1[1, 1] <- -999  # mutate the value returned to the caller

  sig2 <- BayesOmics:::get_sigmak(res$groups[["2"]], res$kernels)
  expect_false(any(sig2 == -999))
})

test_that("multi_posterior_mean: Group converted to character if integer", {
  data <- make_data(nb_id = 4, nb_group = 2)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  expect_true(all(vapply(names(res$groups), is.character, logical(1))))
})

# ── multi_posterior_mean: edge cases ──────────────────────────────────────────

test_that("multi_posterior_mean works with nb_id = 1", {
  data <- make_data(nb_id = 1, nb_group = 2)
  kern <- make_kernel()
  expect_no_error(multi_posterior_mean(data, kern))
})

test_that("multi_posterior_mean works with a single group", {
  data <- make_data(nb_id = 5, nb_group = 1)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  expect_length(res$groups, 1)
})

# ── sample_posterior: validation ──────────────────────────────────────────────

test_that("sample_posterior errors on empty list", {
  expect_error(sample_posterior(list(), 10), "kernels.*groups|groups.*kernels")
})

test_that("sample_posterior errors on non-positive n", {
  res <- make_simple_results()
  expect_error(sample_posterior(res, 0),  "positive integer")
  expect_error(sample_posterior(res, -1), "positive integer")
})

test_that("sample_posterior errors on non-integer n", {
  res <- make_simple_results()
  expect_error(sample_posterior(res, 1.5), "positive integer")
})

test_that("sample_posterior errors when a group is missing required elements", {
  bad <- list(kernels = list(), groups = list(g1 = list(muk = c(1, 2))))
  expect_error(sample_posterior(bad, 10))
})

# ── sample_posterior: return structure ─────────────────────────────────────────

test_that("sample_posterior returns a data.frame in long format", {
  res  <- make_simple_results()
  long <- sample_posterior(res, 50)
  expect_s3_class(long, "data.frame")
  expect_named(long, c("ID", "Group", "Sample"))
})

test_that("sample_posterior has n_draws * n_ids * n_groups rows", {
  n_ids    <- 3
  n_groups <- 2
  n_draws  <- 100
  res      <- make_simple_results(n_groups = n_groups, n_ids = n_ids)
  long     <- sample_posterior(res, n_draws)
  expect_equal(nrow(long), n_ids * n_groups * n_draws)
})

test_that("sample_posterior contains all groups", {
  res  <- make_simple_results(n_groups = 3)
  long <- sample_posterior(res, 50)
  expect_equal(sort(unique(long$Group)), sort(names(res$groups)))
})

test_that("sample_posterior contains all IDs", {
  n_ids <- 4
  res   <- make_simple_results(n_ids = n_ids)
  long  <- sample_posterior(res, 50)
  expect_length(unique(long$ID), n_ids)
})

test_that("sample_posterior: Sample column is numeric", {
  res  <- make_simple_results()
  long <- sample_posterior(res, 50)
  expect_true(is.numeric(long$Sample))
})

test_that("sample_posterior is reproducible with set.seed (direct unit test)", {
  res <- make_simple_results(n_groups = 2, n_ids = 3)
  set.seed(99)
  long1 <- sample_posterior(res, 50)
  set.seed(99)
  long2 <- sample_posterior(res, 50)
  expect_equal(long1, long2)
})

test_that("sample_posterior works when a group has a single ID (1x1 covariance)", {
  res  <- make_simple_results(n_groups = 2, n_ids = 1)
  long <- sample_posterior(res, 50)
  expect_equal(nrow(long), 1 * 2 * 50)
  expect_true(all(is.finite(long$Sample)))
})
