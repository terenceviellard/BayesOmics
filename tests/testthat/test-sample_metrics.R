# -- compute_group_diff_samples: input validation ------------------------------

test_that("compute_group_diff_samples errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(compute_group_diff_samples(bad, energy_metric()), "following columns are missing")
})

test_that("compute_group_diff_samples errors when 'metric' is not a sample_metric", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
  ))
  expect_error(compute_group_diff_samples(sd, mahalanobis_metric()), "sample_metric")
})

test_that("compute_group_diff_samples errors when two groups don't share the same set of IDs", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_3")))
  ))
  expect_error(compute_group_diff_samples(sd, energy_metric()), "do not share the same set of IDs")
})

# -- compute_group_diff_samples: return structure -------------------------------

test_that("compute_group_diff_samples returns a square, symmetric matrix with group dimnames", {
  set.seed(1)
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = c("ID_1", "ID_2"), n = 100)
  mat <- compute_group_diff_samples(sd, energy_metric())
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(3, 3))
  expect_equal(rownames(mat), c("G1", "G2", "G3"))
  expect_equal(mat, t(mat))
})

test_that("compute_group_diff_samples with a single group returns a 1x1 zero matrix", {
  sd <- make_sample_distrib(groups = c("G1"), ids = c("ID_1"), n = 50)
  mat <- compute_group_diff_samples(sd, energy_metric())
  expect_equal(dim(mat), c(1, 1))
  expect_equal(mat[1, 1], 0)
})

test_that("compute_group_diff_samples fixes the diagonal at exactly 0 for all sample_metrics", {
  set.seed(1)
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1", "ID_2"), n = 80)
  for (m in list(energy_metric(), mmd_metric(), sliced_wasserstein_metric(n_projections = 10))) {
    mat <- compute_group_diff_samples(sd, m)
    expect_equal(unname(diag(mat)), c(0, 0), info = class(m)[1])
  }
})

# -- energy_metric / mmd_metric: separated groups score higher than close ones --

test_that("energy_metric is larger for well-separated groups than for nearly identical ones", {
  set.seed(1)
  sd_close <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1"), n = 300)
  set.seed(1)
  n <- 300
  sd_far <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n, mean = 0), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n, mean = 50), stringsAsFactors = FALSE)
  ))
  e_close <- compute_group_diff_samples(sd_close, energy_metric())["G1", "G2"]
  e_far   <- compute_group_diff_samples(sd_far, energy_metric())["G1", "G2"]
  expect_gt(e_far, e_close)
})

test_that("energy_metric is near 0 for two samples drawn from the same distribution", {
  set.seed(7)
  n <- 500
  sd <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n), stringsAsFactors = FALSE)
  ))
  e <- compute_group_diff_samples(sd, energy_metric())["G1", "G2"]
  expect_true(abs(e) < 0.5)
})

test_that("mmd_metric is larger for well-separated groups than for identical-distribution groups", {
  set.seed(3)
  n <- 300
  sd_same <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n), stringsAsFactors = FALSE)
  ))
  sd_far <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n, mean = 0), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n, mean = 20), stringsAsFactors = FALSE)
  ))
  mmd_same <- compute_group_diff_samples(sd_same, mmd_metric())["G1", "G2"]
  mmd_far  <- compute_group_diff_samples(sd_far, mmd_metric())["G1", "G2"]
  expect_gt(mmd_far, mmd_same)
})

test_that("mmd_metric accepts an explicit bandwidth", {
  set.seed(1)
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1"), n = 100)
  expect_no_error(mat <- compute_group_diff_samples(sd, mmd_metric(bandwidth = 2)))
  expect_true(is.finite(mat["G1", "G2"]))
})

# -- sliced_wasserstein_metric: separation-sensitive, reproducible with a seed --

test_that("sliced_wasserstein_metric is larger for well-separated groups than for close ones", {
  set.seed(1)
  n <- 300
  sd_close <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n, mean = 0), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n, mean = 0.1), stringsAsFactors = FALSE)
  ))
  sd_far <- do.call(rbind, list(
    data.frame(ID = "ID_1", Group = "G1", Sample = stats::rnorm(n, mean = 0), stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = stats::rnorm(n, mean = 20), stringsAsFactors = FALSE)
  ))
  set.seed(99)
  sw_close <- compute_group_diff_samples(sd_close, sliced_wasserstein_metric(n_projections = 30))["G1", "G2"]
  set.seed(99)
  sw_far <- compute_group_diff_samples(sd_far, sliced_wasserstein_metric(n_projections = 30))["G1", "G2"]
  expect_gt(sw_far, sw_close)
})

test_that("sliced_wasserstein_metric is reproducible given the same random seed", {
  set.seed(1)
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1", "ID_2"), n = 100)
  set.seed(123)
  mat1 <- compute_group_diff_samples(sd, sliced_wasserstein_metric(n_projections = 10))
  set.seed(123)
  mat2 <- compute_group_diff_samples(sd, sliced_wasserstein_metric(n_projections = 10))
  expect_equal(mat1, mat2)
})

# -- 1D sanity: sliced Wasserstein with 1 shared ID matches direct 1D W1 -------

test_that("sliced_wasserstein_metric with a single shared ID matches the direct 1D Wasserstein distance", {
  set.seed(5)
  n <- 200
  x <- stats::rnorm(n, mean = 0)
  y <- stats::rnorm(n, mean = 3)
  sd <- rbind(
    data.frame(ID = "ID_1", Group = "G1", Sample = x, stringsAsFactors = FALSE),
    data.frame(ID = "ID_1", Group = "G2", Sample = y, stringsAsFactors = FALSE)
  )
  sw <- compute_group_diff_samples(sd, sliced_wasserstein_metric(n_projections = 5))["G1", "G2"]

  qx <- stats::quantile(x, probs = seq(0, 1, length.out = n + 1), type = 7, names = FALSE)
  qy <- stats::quantile(y, probs = seq(0, 1, length.out = n + 1), type = 7, names = FALSE)
  w1_direct <- mean(abs(qx - qy))
  # With a single shared ID, every random projection direction is +/-1, so
  # every one of the n_projections terms reduces to the same 1D W1 distance.
  expect_equal(sw, w1_direct, tolerance = 1e-8)
})
