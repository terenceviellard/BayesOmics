# ── Structure ─────────────────────────────────────────────────────────────────

test_that("simu_db returns a data.frame", {
  expect_s3_class(simu_db(), "data.frame")
})

test_that("simu_db has exactly the required columns in order", {
  expect_named(simu_db(), c("ID", "Group", "Sample", "Input", "Output"))
})

test_that("simu_db row count = nb_id * nb_group * nb_sample", {
  expect_equal(nrow(simu_db(nb_id = 4, nb_group = 3, nb_sample = 2)), 4 * 3 * 2)
  expect_equal(nrow(simu_db(nb_id = 1, nb_group = 1, nb_sample = 1)), 1)
  expect_equal(nrow(simu_db(nb_id = 10, nb_group = 5, nb_sample = 3)), 150)
})

test_that("simu_db has correct number of unique IDs", {
  data <- simu_db(nb_id = 7, nb_group = 2, nb_sample = 3)
  expect_equal(length(unique(data$ID)), 7)
})

test_that("simu_db has correct number of unique Groups", {
  data <- simu_db(nb_id = 4, nb_group = 5, nb_sample = 2)
  expect_equal(length(unique(data$Group)), 5)
})

test_that("simu_db has correct Sample values per (ID, Group)", {
  data <- simu_db(nb_id = 3, nb_group = 2, nb_sample = 4)
  for (id in unique(data$ID)) {
    for (g in unique(data$Group)) {
      s <- data$Sample[data$ID == id & data$Group == g]
      expect_equal(sort(s), seq_len(4))
    }
  }
})

# ── Types ──────────────────────────────────────────────────────────────────────

test_that("simu_db: Input and Output are numeric", {
  data <- simu_db()
  expect_true(is.numeric(data$Input))
  expect_true(is.numeric(data$Output))
})

test_that("simu_db: no NA values", {
  expect_false(anyNA(simu_db()))
})

# ── Input ─────────────────────────────────────────────────────────────────────

test_that("simu_db: Input is constant within each (ID, Group) pair across samples", {
  data <- simu_db(nb_id = 4, nb_group = 2, nb_sample = 5)
  for (id in unique(data$ID)) {
    for (g in unique(data$Group)) {
      inputs <- data$Input[data$ID == id & data$Group == g]
      expect_equal(length(unique(inputs)), 1,
                   info = paste("ID =", id, ", Group =", g))
    }
  }
})

test_that("simu_db: Input stays within range_input", {
  data <- simu_db(nb_id = 20, nb_group = 3, nb_sample = 2,
                  range_input = c(10, 20))
  expect_true(all(data$Input >= 10 & data$Input <= 20))
})

test_that("simu_db: each (ID, Group) pair has a distinct Input", {
  data <- simu_db(nb_id = 10, nb_group = 2, nb_sample = 1)
  for (g in unique(data$Group)) {
    sub  <- data[data$Group == g, ]
    n_ids    <- length(unique(sub$ID))
    n_inputs <- length(unique(sub$Input))
    expect_equal(n_inputs, n_ids)
  }
})

test_that("simu_db: Input is shared across groups for the same ID", {
  # So that two groups built from the same simu_db() call end up with the
  # same kernel_key in multi_posterior_mean(), enabling calculate_group_overlaps().
  data <- simu_db(nb_id = 6, nb_group = 3, nb_sample = 2)
  for (id in unique(data$ID)) {
    inputs <- data$Input[data$ID == id]
    expect_equal(length(unique(inputs)), 1, info = paste("ID =", id))
  }
})

# ── Output ────────────────────────────────────────────────────────────────────

test_that("simu_db: group effect is diff_group per group level", {
  set.seed(1)
  data <- simu_db(nb_id = 200, nb_group = 2, nb_sample = 1,
                  diff_group = 10, var_sample = 0)
  mean_g1 <- mean(data$Output[data$Group == 1])
  mean_g2 <- mean(data$Output[data$Group == 2])
  expect_equal(mean_g2 - mean_g1, 10, tolerance = 0.5)
})

test_that("simu_db: diff_group=0 removes systematic group difference", {
  set.seed(1)
  data <- simu_db(nb_id = 500, nb_group = 2, nb_sample = 1,
                  diff_group = 0, var_sample = 0)
  mean_g1 <- mean(data$Output[data$Group == 1])
  mean_g2 <- mean(data$Output[data$Group == 2])
  expect_equal(mean_g2 - mean_g1, 0, tolerance = 0.1)
})

test_that("simu_db: var_sample=0 produces no noise (deterministic given seed)", {
  set.seed(7)
  d1 <- simu_db(nb_id = 5, nb_group = 2, nb_sample = 3, var_sample = 0)
  set.seed(7)
  d2 <- simu_db(nb_id = 5, nb_group = 2, nb_sample = 3, var_sample = 0)
  expect_equal(d1$Output, d2$Output)
})

test_that("simu_db: negative diff_group decreases output across groups", {
  set.seed(1)
  data <- simu_db(nb_id = 200, nb_group = 2, nb_sample = 1,
                  diff_group = -10, var_sample = 0)
  mean_g1 <- mean(data$Output[data$Group == 1])
  mean_g2 <- mean(data$Output[data$Group == 2])
  expect_equal(mean_g2 - mean_g1, -10, tolerance = 0.5)
})

test_that("simu_db: larger var_sample increases output variance", {
  set.seed(1)
  data_low  <- simu_db(nb_id = 100, var_sample = 0.01)
  data_high <- simu_db(nb_id = 100, var_sample = 10)
  expect_gt(stats::var(data_high$Output), stats::var(data_low$Output))
})

# ── Reproducibility ────────────────────────────────────────────────────────────

test_that("simu_db is reproducible with set.seed", {
  set.seed(123)
  d1 <- simu_db(nb_id = 10, nb_group = 3, nb_sample = 4)
  set.seed(123)
  d2 <- simu_db(nb_id = 10, nb_group = 3, nb_sample = 4)
  expect_equal(d1, d2)
})

# ── Edge cases ─────────────────────────────────────────────────────────────────

test_that("simu_db works with nb_id = nb_group = nb_sample = 1", {
  data <- simu_db(nb_id = 1, nb_group = 1, nb_sample = 1)
  expect_equal(nrow(data), 1)
  expect_named(data, c("ID", "Group", "Sample", "Input", "Output"))
})

test_that("simu_db works with large nb_id", {
  expect_no_error(simu_db(nb_id = 500, nb_group = 3, nb_sample = 1))
})

# ── simu_db: input validation (mirrors simu_db_kernel) ──────────────────────

test_that("simu_db validates nb_id, nb_group, nb_sample", {
  expect_error(simu_db(nb_id = 0), "nb_id")
  expect_error(simu_db(nb_id = 2.5), "nb_id")
  expect_error(simu_db(nb_group = 0), "nb_group")
  expect_error(simu_db(nb_sample = -1), "nb_sample")
})

test_that("simu_db validates range_input and range_output", {
  expect_error(simu_db(range_output = c(10, 1)), "range_output")
  expect_error(simu_db(range_input = c(1, 2, 3)), "range_input")
})

test_that("simu_db validates diff_group", {
  expect_error(simu_db(diff_group = c(1, 2)), "diff_group")
})

test_that("simu_db allows var_sample = 0 but rejects negative var_sample", {
  expect_no_error(simu_db(var_sample = 0))
  expect_error(simu_db(var_sample = -1), "var_sample")
})

# ── simu_db_kernel: structure ───────────────────────────────────────────────

test_that("simu_db_kernel returns a data.frame with the expected columns and row count", {
  data <- simu_db_kernel(nb_id = 4, nb_group = 3, nb_sample = 2, kernel = make_kernel())
  expect_s3_class(data, "data.frame")
  expect_named(data, c("ID", "Group", "Sample", "Input", "Output"))
  expect_equal(nrow(data), 4 * 3 * 2)
})

test_that("simu_db_kernel: no NA values, even with nb_sample > 1 (repeated Input)", {
  data <- simu_db_kernel(nb_id = 5, nb_group = 2, nb_sample = 5, kernel = make_kernel())
  expect_false(anyNA(data))
})

test_that("simu_db_kernel: Input is constant within each (ID, Group) pair across samples", {
  data <- simu_db_kernel(nb_id = 4, nb_group = 2, nb_sample = 5, kernel = make_kernel())
  for (id in unique(data$ID)) {
    for (g in unique(data$Group)) {
      inputs <- data$Input[data$ID == id & data$Group == g]
      expect_equal(length(unique(inputs)), 1,
                   info = paste("ID =", id, ", Group =", g))
    }
  }
})

test_that("simu_db_kernel: Input is shared across groups for the same ID", {
  data <- simu_db_kernel(nb_id = 5, nb_group = 3, nb_sample = 2, kernel = make_kernel())
  for (id in unique(data$ID)) {
    inputs <- data$Input[data$ID == id]
    expect_equal(length(unique(inputs)), 1, info = paste("ID =", id))
  }
})

test_that("simu_db_kernel: group effect is diff_group per group level", {
  set.seed(1)
  data <- simu_db_kernel(nb_id = 200, nb_group = 2, nb_sample = 1,
                          diff_group = 10, var_sample = 0.01, kernel = make_kernel())
  mean_g1 <- mean(data$Output[data$Group == 1])
  mean_g2 <- mean(data$Output[data$Group == 2])
  expect_equal(mean_g2 - mean_g1, 10, tolerance = 1)
})

test_that("simu_db_kernel is reproducible with set.seed", {
  set.seed(123)
  d1 <- simu_db_kernel(nb_id = 6, nb_group = 2, nb_sample = 3, kernel = make_kernel())
  set.seed(123)
  d2 <- simu_db_kernel(nb_id = 6, nb_group = 2, nb_sample = 3, kernel = make_kernel())
  expect_equal(d1, d2)
})

test_that("simu_db_kernel: exposes ground truth (mu_true, base_input) via attributes", {
  data <- simu_db_kernel(nb_id = 4, nb_group = 3, nb_sample = 2, kernel = make_kernel(), diff_group = 5)
  mu_true <- attr(data, "mu_true")
  expect_equal(names(mu_true), c("1", "2", "3"))
  expect_equal(length(mu_true[["1"]]), 4)
  # deterministic (mu_random = FALSE): every id shares the same baseline within a
  # group, shifted by diff_group across groups.
  expect_equal(length(unique(mu_true[["1"]])), 1)
  expect_equal(unname(mu_true[["2"]][1] - mu_true[["1"]][1]), 5)

  base_input <- attr(data, "base_input")
  expect_equal(names(base_input), paste0("ID_", 1:4))
})

test_that("simu_db_kernel: mu_true matches the prior draw when mu_random = TRUE", {
  set.seed(13)
  data <- simu_db_kernel(nb_id = 5, nb_group = 1, nb_sample = 1, kernel = make_kernel(),
                          mu_random = TRUE, mu_0 = 10, lambda_0 = 1)
  mu_true <- attr(data, "mu_true")[["1"]]
  expect_equal(length(unique(mu_true)), 5)  # not collapsed to a single constant
})

test_that("simu_db_kernel: replicates are independent draws (variance matches kernel + var_sample)", {
  # Each replicate y_n | mu ~ N(mu, Sigma_theta + var_sample*I) independently, so with
  # nb_id = 1, Sigma_theta is the 1x1 kernel self-covariance (variance_se) and the
  # per-replicate variance across many samples should match variance_se + var_sample.
  set.seed(1)
  ker  <- make_kernel(hp = c(4, 5))  # variance_se = 4
  data <- simu_db_kernel(nb_id = 1, nb_group = 1, nb_sample = 5000,
                          var_sample = 2, kernel = ker)
  expect_equal(stats::var(data$Output), 4 + 2, tolerance = 0.5)
})

test_that("simu_db_kernel: Sigma_theta reflects kernel distance across ids, not within-id repeats", {
  # Two ids close together (small distance) should be more correlated than two ids far
  # apart, for a kernel with positive length scale -- this is the cross-id structure the
  # new generative model relies on (the old model had no such structure at all).
  set.seed(4)
  ker <- methods::new("SEKernel")
  ker <- keRnel::set_hyperparameters(ker, c(variance_se = 5, length_scale_se = 10))
  input_mat <- matrix(c(0, 1, 100), ncol = 1)
  K <- keRnel::pairwise_kernel(ker, input_mat, input_mat)
  expect_gt(K[1, 2], K[1, 3])  # id at distance 1 more correlated than id at distance 100
})

test_that("simu_db_kernel works with nb_id = nb_group = nb_sample = 1", {
  data <- simu_db_kernel(nb_id = 1, nb_group = 1, nb_sample = 1, kernel = make_kernel())
  expect_equal(nrow(data), 1)
  expect_named(data, c("ID", "Group", "Sample", "Input", "Output"))
})

test_that("simu_db_kernel works with large nb_id and nb_sample", {
  expect_no_error(
    simu_db_kernel(nb_id = 50, nb_group = 3, nb_sample = 10, kernel = make_kernel())
  )
})

test_that("chol_inv_jitter_diag produces a valid, symmetric, PD matrix from a constant raw matrix", {
  # Generic sanity check of the helper itself (used by simu_db_kernel internally),
  # independent of any particular Input configuration.
  ker <- make_kernel(hp = c(4, 5))
  input_mat <- matrix(rep(10, 4), ncol = 1)
  raw_cov <- 3 * keRnel::pairwise_kernel(ker, input_mat, input_mat)
  jittered_cov <- BayesOmics:::chol_inv_jitter_diag(raw_cov, 1e-6)
  expect_equal(dim(jittered_cov), c(4, 4))
  expect_true(isSymmetric(jittered_cov))
  expect_no_error(chol(jittered_cov))
  expect_equal(unname(raw_cov[1, 1]), 3 * 4, tolerance = 1e-8)
})

test_that("simu_db_kernel works with nb_group = 1 and several samples", {
  expect_no_error(
    simu_db_kernel(nb_id = 5, nb_group = 1, nb_sample = 6, kernel = make_kernel())
  )
})

# ── simu_db_kernel: Sigma_theta shared across groups ────────────────────────

test_that("simu_db_kernel: every group shares the same kernel matrix (same Input set)", {
  # Required for multi_posterior_mean() + calculate_group_overlaps() to work across groups.
  data <- simu_db_kernel(nb_id = 5, nb_group = 4, nb_sample = 2, kernel = make_kernel())
  posterior <- multi_posterior_mean(data, make_kernel())
  keys <- vapply(posterior$groups, function(g) g$kernel_key, character(1))
  expect_equal(length(unique(keys)), 1)
  expect_no_error(calculate_group_overlaps(posterior))
})

# ── simu_db_kernel: mu_random option ────────────────────────────────────────

test_that("simu_db_kernel: mu_random = FALSE (default) is deterministic given the same base profile", {
  set.seed(9)
  d1 <- simu_db_kernel(nb_id = 4, nb_group = 2, nb_sample = 2, kernel = make_kernel())
  set.seed(9)
  d2 <- simu_db_kernel(nb_id = 4, nb_group = 2, nb_sample = 2, kernel = make_kernel())
  expect_equal(d1, d2)
})

test_that("simu_db_kernel: mu_random = TRUE requires a valid lambda_0", {
  expect_error(simu_db_kernel(kernel = make_kernel(), mu_random = TRUE, lambda_0 = 0), "lambda_0")
  expect_error(simu_db_kernel(kernel = make_kernel(), mu_random = TRUE, lambda_0 = -1), "lambda_0")
  expect_no_error(simu_db_kernel(kernel = make_kernel(), mu_random = TRUE, lambda_0 = 1))
})

test_that("simu_db_kernel: mu_random validates its type", {
  expect_error(simu_db_kernel(kernel = make_kernel(), mu_random = "yes"), "mu_random")
  expect_error(simu_db_kernel(kernel = make_kernel(), mu_random = NA), "mu_random")
})

test_that("simu_db_kernel: mu_random = TRUE makes the group mean vary with mu_0/lambda_0 prior spread", {
  # With a much larger lambda_0 (tighter prior), the drawn group means should land
  # closer to mu_0 on average across many simulated datasets than with a small lambda_0.
  set.seed(11)
  ker <- methods::new("SEKernel")
  ker <- keRnel::set_hyperparameters(ker, c(variance_se = 25, length_scale_se = 5))
  draw_group1_mean <- function(lambda_0) {
    data <- simu_db_kernel(nb_id = 6, nb_group = 1, nb_sample = 1, kernel = ker,
                            var_sample = 0.01, mu_random = TRUE, mu_0 = 20, lambda_0 = lambda_0)
    mean(data$Output)
  }
  tight  <- replicate(50, draw_group1_mean(lambda_0 = 100))
  loose  <- replicate(50, draw_group1_mean(lambda_0 = 0.1))
  expect_lt(stats::sd(tight), stats::sd(loose))
})

# ── simu_db_kernel: input validation ────────────────────────────────────────

test_that("simu_db_kernel requires a valid AbstractKernel object", {
  expect_error(simu_db_kernel(), "kernel")
  expect_error(simu_db_kernel(kernel = list(pairwise_kernel = function(x) diag(length(x)))), "kernel")
  expect_error(simu_db_kernel(kernel = "not a kernel"), "kernel")
})

test_that("simu_db_kernel validates nb_id, nb_group, nb_sample", {
  expect_error(simu_db_kernel(nb_id = 0, kernel = make_kernel()), "nb_id")
  expect_error(simu_db_kernel(nb_id = 2.5, kernel = make_kernel()), "nb_id")
  expect_error(simu_db_kernel(nb_group = 0, kernel = make_kernel()), "nb_group")
  expect_error(simu_db_kernel(nb_sample = -1, kernel = make_kernel()), "nb_sample")
})

test_that("simu_db_kernel validates range_input and range_output", {
  expect_error(simu_db_kernel(range_output = c(10, 1), kernel = make_kernel()), "range_output")
  expect_error(simu_db_kernel(range_input = c(1, 2, 3), kernel = make_kernel()), "range_input")
})

test_that("simu_db_kernel validates diff_group and var_sample", {
  expect_error(simu_db_kernel(diff_group = c(1, 2), kernel = make_kernel()), "diff_group")
  expect_error(simu_db_kernel(var_sample = 0, kernel = make_kernel()), "var_sample")
  expect_error(simu_db_kernel(var_sample = -1, kernel = make_kernel()), "var_sample")
})

test_that("simu_db_kernel validates pen_diag", {
  expect_error(simu_db_kernel(pen_diag = -1, kernel = make_kernel()), "pen_diag")
})

# ── chol_inv_jitter_diag: bounded recursion ─────────────────────────────────

test_that("chol_inv_jitter_diag warns when a large jitter (relative to pen_diag) was needed", {
  mat <- diag(c(-50, 5, 5))  # one large negative eigenvalue forces many jitter doublings
  expect_warning(
    BayesOmics:::chol_inv_jitter_diag(mat, pen_diag = 1e-6),
    "ill-conditioned"
  )
})

test_that("chol_inv_jitter_diag does not warn when no/little jitter was needed", {
  mat <- diag(c(4, 3, 2))  # already positive-definite
  expect_warning(BayesOmics:::chol_inv_jitter_diag(mat, pen_diag = 1e-6), NA)
})

test_that("chol_inv_jitter_diag stops with an informative error instead of recursing forever", {
  mat <- matrix(NaN, 2, 2)
  expect_error(
    BayesOmics:::chol_inv_jitter_diag(mat, pen_diag = 1e-6, max_tries = 5),
    "could not be made positive-definite"
  )
})

test_that("chol_inv_jitter_diag returns the matrix unchanged when pen_diag = 0 and already PD", {
  mat <- diag(2) * 4
  expect_equal(BayesOmics:::chol_inv_jitter_diag(mat, pen_diag = 0), mat)
})
