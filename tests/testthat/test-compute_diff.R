# -- calculate_group_overlaps: input validation --------------------------------

test_that("calculate_group_overlaps errors when top-level structure is wrong", {
  expect_error(calculate_group_overlaps(list(foo = 1)), "kernels.*groups|groups.*kernels")
})

test_that("calculate_group_overlaps errors when a group is missing required elements", {
  bad <- list(
    kernels = list(),
    groups  = list(
      g1 = list(muk = c(0), id_to_input = c(a = 1), kernel_key = "k", scale = 1),
      g2 = list(muk = c(0))
    )
  )
  expect_error(calculate_group_overlaps(bad), "muk.*id_to_input|id_to_input.*muk")
})

# -- calculate_group_overlaps: return structure ---------------------------------

test_that("calculate_group_overlaps returns a square matrix", {
  res <- make_simple_results(n_groups = 3, n_ids = 1)
  mat <- calculate_group_overlaps(res)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), ncol(mat))
})

test_that("calculate_group_overlaps has correct dimension", {
  n_groups <- 4
  res <- make_simple_results(n_groups = n_groups, n_ids = 1)
  mat <- calculate_group_overlaps(res)
  expect_equal(dim(mat), c(n_groups, n_groups))
})

test_that("calculate_group_overlaps has row and column names equal to group names", {
  res <- make_simple_results(n_groups = 3, n_ids = 1)
  mat <- calculate_group_overlaps(res)
  expect_equal(rownames(mat), names(res$groups))
  expect_equal(colnames(mat), names(res$groups))
})

# -- calculate_group_overlaps: numeric properties ------------------------------

test_that("calculate_group_overlaps: diagonal is all 1", {
  res <- make_simple_results(n_groups = 3, n_ids = 1)
  mat <- calculate_group_overlaps(res)
  expect_equal(unname(diag(mat)), rep(1, 3))
})

test_that("calculate_group_overlaps: matrix is symmetric", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  mat <- calculate_group_overlaps(res)
  expect_equal(mat, t(mat))
})

test_that("calculate_group_overlaps: all off-diagonal values in [0, 1]", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  mat <- calculate_group_overlaps(res)
  off_diag <- mat[row(mat) != col(mat)]
  expect_true(all(off_diag >= 0 & off_diag <= 1))
})

test_that("calculate_group_overlaps: identical groups yield overlap = 1", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 5.0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 5.0), sigma = matrix(1, 1, 1))
  ))
  mat <- calculate_group_overlaps(res)
  expect_equal(mat["g1", "g2"], 1)
})

test_that("calculate_group_overlaps: very separated groups yield overlap near 0", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = -100.0), sigma = matrix(0.01, 1, 1)),
    g2 = list(muk = c(ID_1 =  100.0), sigma = matrix(0.01, 1, 1))
  ))
  mat <- calculate_group_overlaps(res)
  expect_lt(mat["g1", "g2"], 0.05)
})

test_that("calculate_group_overlaps: overlap decreases as group means diverge", {
  res_close <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0.0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 1.0), sigma = matrix(1, 1, 1))
  ))
  res_far <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0.0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 5.0), sigma = matrix(1, 1, 1))
  ))
  ov_close <- calculate_group_overlaps(res_close)["g1", "g2"]
  ov_far   <- calculate_group_overlaps(res_far)["g1", "g2"]
  expect_gt(ov_close, ov_far)
})

# -- calculate_group_overlaps: edge cases --------------------------------------

test_that("calculate_group_overlaps with single group returns 1x1 identity matrix", {
  res <- make_simple_results(n_groups = 1, n_ids = 2)
  mat <- calculate_group_overlaps(res)
  expect_equal(dim(mat), c(1, 1))
  expect_equal(mat[1, 1], 1)
})

test_that("calculate_group_overlaps with two groups returns 2x2 matrix", {
  res <- make_simple_results(n_groups = 2, n_ids = 2)
  mat <- calculate_group_overlaps(res)
  expect_equal(dim(mat), c(2, 2))
})

test_that("calculate_group_overlaps works with multivariate posteriors", {
  n_ids <- 5
  res   <- make_simple_results(n_groups = 2, n_ids = n_ids)
  expect_no_error(calculate_group_overlaps(res))
})

# -- calculate_group_overlaps: ID alignment (fix for the former silent
# misalignment / unclear-error gaps) ----------------------------------------

test_that("ID alignment is name-based: swapped ID labels between groups change the overlap", {
  # g1 and g2 have the same muk values but for swapped IDs; a correct
  # ID-aware comparison must realign by name and therefore must NOT report a
  # perfect overlap here, since the IDs don't actually correspond to each
  # other across groups.
  res <- make_custom_results(list(
    A = list(muk = c(ID_1 = 0, ID_2 = 5), sigma = diag(2)),
    B = list(muk = c(ID_2 = 0, ID_1 = 5), sigma = diag(2))
  ))
  mat <- calculate_group_overlaps(res)
  expect_lt(mat["A", "B"], 1)
})

test_that("mismatched number of IDs between groups raises a clear, package-level error", {
  res <- make_custom_results(list(
    A = list(muk = c(ID_1 = 0, ID_2 = 5),          sigma = diag(2)),
    B = list(muk = c(ID_1 = 0, ID_2 = 5, ID_3 = 1), sigma = diag(3))
  ))
  expect_error(calculate_group_overlaps(res), "do not share the same set of IDs")
})

test_that("mismatched ID sets (same count, different labels) raises a clear error", {
  res <- make_custom_results(list(
    A = list(muk = c(ID_1 = 0, ID_2 = 5), sigma = diag(2)),
    B = list(muk = c(ID_3 = 0, ID_4 = 5), sigma = diag(2))
  ))
  expect_error(calculate_group_overlaps(res), "do not share the same set of IDs")
})

test_that("calculate_group_overlaps does not crash with near-singular (very small variance) Sigma", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0.0), sigma = matrix(1e-8, 1, 1)),
    g2 = list(muk = c(ID_1 = 0.0), sigma = matrix(1e-8, 1, 1))
  ))
  expect_no_error(mat <- calculate_group_overlaps(res))
  expect_true(is.finite(mat["g1", "g2"]))
})

test_that("calculate_group_overlaps is deterministic (no randomness)", {
  res  <- make_simple_results(n_groups = 3, n_ids = 2)
  mat1 <- calculate_group_overlaps(res)
  mat2 <- calculate_group_overlaps(res)
  expect_equal(mat1, mat2)
})

# -- calculate_group_overlaps: unequal scale (c != 1) closed-form branch -------
# Sigma1 = raw/scale1, Sigma2 = raw/scale2, c = scale1/scale2 != 1. Reference
# values below were cross-checked against direct numerical integration of
# pmin(dnorm(x, mu1, sd1), dnorm(x, mu2, sd2)) over the real line (1D, so the
# OVL coefficient is an unambiguous ground truth independent of this package).

test_that("calculate_group_overlaps: c != 1 matches numerical integration (c > 1)", {
  # raw = 4 (sd1 = 1), scale1 = 1 -> Sigma1 = 4; scale2 = 4 -> Sigma2 = 1
  # i.e. sd1 = 2, sd2 = 1, mu1 = 0, mu2 = 2 -> c = Sigma2/Sigma1 = 1/4
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(4, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1), scale = 4)
  ))
  mat <- calculate_group_overlaps(res)
  expect_equal(mat["g1", "g2"], 0.453388, tolerance = 1e-5)
})

test_that("calculate_group_overlaps: c != 1 matches numerical integration (c < 1)", {
  # raw = 1, scale1 = 1 -> Sigma1 = 1 (sd1 = 1); scale2 = 0.25 -> Sigma2 = 4 (sd2 = 2)
  # mu1 = 0, mu2 = 2 -> same configuration as above with groups swapped
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(1, 1, 1), scale = 0.25)
  ))
  mat <- calculate_group_overlaps(res)
  expect_equal(mat["g1", "g2"], 0.453388, tolerance = 1e-5)
})

test_that("calculate_group_overlaps: c != 1, identical means, different scale", {
  # mu1 = mu2 = 0, sd1 = 1 (scale = 1), sd2 = 3 (scale = 1/9)
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1), scale = 1 / 9)
  ))
  mat <- calculate_group_overlaps(res)
  expect_equal(mat["g1", "g2"], 0.515672, tolerance = 1e-5)
})

test_that("calculate_group_overlaps: c != 1 result is symmetric and in [0, 1]", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 1, ID_2 = -2), sigma = diag(c(2, 3)), scale = 1),
    g2 = list(muk = c(ID_1 = -1, ID_2 = 1), sigma = diag(c(2, 3)), scale = 2.5)
  ))
  mat <- calculate_group_overlaps(res)
  expect_equal(mat["g1", "g2"], mat["g2", "g1"])
  expect_true(mat["g1", "g2"] >= 0 && mat["g1", "g2"] <= 1)
})

# -- calculate_group_overlaps: near-equal scale (avoids 1/(1-c)^2 blow-up) ----
# When c_ratio = scale1/scale2 is very close to 1 but not exactly equal, the
# c != 1 closed-form divides by (1 - c_ratio)^2, which blows up numerically.
# calculate_group_overlaps() should detect this (relative tolerance 1e-6),
# warn, and fall back to the c = 1 formula using the smaller-scale (larger,
# more conservative) posterior covariance.

test_that("near-equal but unequal scales warn and fall back to the c = 1 formula", {
  # raw sigma is shared (same kernel_key); scale1 < scale2 by less than the
  # 1e-6 relative tolerance, so Sigma_eq must come from g1 (the smaller scale).
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(4, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1), scale = 1 + 1e-9)
  ))
  expect_warning(
    mat <- calculate_group_overlaps(res),
    "nearly identical"
  )
  # Sigma_eq = raw / scale1 = 4 -> sd = 2; delta = 2 -> D2 = delta^2/Sigma_eq = 1
  expected <- 2 * stats::pnorm(-sqrt(1) / 2)
  expect_equal(mat["g1", "g2"], expected, tolerance = 1e-5)
})

test_that("near-equal scale fallback picks the smaller-scale group regardless of order", {
  # Same as above but with g1/g2 roles swapped: now g2 has the smaller scale.
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(4, 1, 1), scale = 1 + 1e-9),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1), scale = 1)
  ))
  expect_warning(
    mat <- calculate_group_overlaps(res),
    "nearly identical"
  )
  expected <- 2 * stats::pnorm(-sqrt(1) / 2)
  expect_equal(mat["g1", "g2"], expected, tolerance = 1e-5)
})

test_that("exactly equal scales still use the c = 1 formula without warning", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(1, 1, 1), scale = 1)
  ))
  expect_warning(mat <- calculate_group_overlaps(res), NA)
  expect_equal(mat["g1", "g2"], 2 * stats::pnorm(-1), tolerance = 1e-5)
})

test_that("scales far enough apart still use the c != 1 formula without warning", {
  # Same configuration as the existing c < 1 closed-form test (ratio = 1/4,
  # well outside the 1e-6 near-equal tolerance).
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(1, 1, 1), scale = 0.25)
  ))
  expect_warning(mat <- calculate_group_overlaps(res), NA)
  expect_equal(mat["g1", "g2"], 0.453388, tolerance = 1e-5)
})

# -- calculate_group_overlaps: O(G^2 * d^3) cost warning ----------------------

test_that("calculate_group_overlaps warns when the number of groups exceeds max_groups_warn", {
  res <- make_simple_results(n_groups = 5, n_ids = 1)
  expect_warning(
    calculate_group_overlaps(res, max_groups_warn = 4),
    "O\\(G\\^2"
  )
})

test_that("calculate_group_overlaps warns when the shared dimension exceeds max_dim_warn", {
  res <- make_simple_results(n_groups = 2, n_ids = 5)
  expect_warning(
    calculate_group_overlaps(res, max_dim_warn = 4),
    "O\\(G\\^2"
  )
})

test_that("calculate_group_overlaps does not warn below the default thresholds", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  expect_warning(calculate_group_overlaps(res), NA)
})

# -- calculate_group_overlaps: different kernel_key (no shared raw Sigma) -----

test_that("calculate_group_overlaps errors when groups don't share the same kernel_key", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 0), sigma = matrix(2, 1, 1))
  ))
  expect_error(calculate_group_overlaps(res), "kernel matrix|kernel_key")
})

# -- compute_multi_diff: input validation --------------------------------------

test_that("compute_multi_diff errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(compute_multi_diff(bad), "following columns are missing")
})

test_that("compute_multi_diff errors with fewer than two groups", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:4, nrow = 4, dimnames = list(NULL, "ID_1"))
  ))
  expect_error(compute_multi_diff(sd), "at least two groups")
})

test_that("compute_multi_diff errors when two groups don't share the same set of IDs", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_3")))
  ))
  expect_error(compute_multi_diff(sd), "do not share the same set of IDs")
})

test_that("compute_multi_diff errors when an ID contributes a different number of draws within a group", {
  sd <- rbind(
    data.frame(ID = "ID_1", Group = "G1", Sample = c(1, 2, 3, 4)),
    data.frame(ID = "ID_2", Group = "G1", Sample = c(1, 2, 3)),
    make_aligned_sample_distrib(list(
      G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
    ))
  )
  expect_error(compute_multi_diff(sd), "different numbers of posterior draws")
})

test_that("compute_multi_diff errors when the two groups being compared have different draw counts", {
  sd <- rbind(
    make_aligned_sample_distrib(list(
      G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
    )),
    make_aligned_sample_distrib(list(
      G2 = matrix(1:6, nrow = 3, dimnames = list(NULL, c("ID_1", "ID_2")))
    ))
  )
  expect_error(compute_multi_diff(sd), "different numbers of posterior draws")
})

test_that("compute_multi_diff errors when 'results' groups don't match 'sample_distrib' groups", {
  sd  <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
  ))
  res <- make_simple_results(n_groups = 2, n_ids = 2)
  names(res$groups) <- c("Other1", "Other2")
  expect_error(compute_multi_diff(sd, results = res), "'results' groups do not match")
})

# -- compute_multi_diff: return structure ---------------------------------------

test_that("compute_multi_diff returns Diff_proba/Diff_mean with the expected columns", {
  multi_diff <- make_multi_diff(nb_id = 3, nb_group = 2, n = 50)
  expect_true(all(c("Group1", "Group2", "Nb_id", "Proba", "Cumul_proba") %in% names(multi_diff$Diff_proba)))
  expect_true(all(c("ID", "Group", "Mean") %in% names(multi_diff$Diff_mean)))
})

test_that("compute_multi_diff omits Overlap_coef when results is NULL", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
  ))
  res <- compute_multi_diff(sd)
  expect_false("Overlap_coef" %in% names(res))
})

test_that("compute_multi_diff attaches Overlap_coef with the expected columns when results is given", {
  multi_diff <- make_multi_diff(nb_id = 3, nb_group = 2, n = 50)
  expect_true("Overlap_coef" %in% names(multi_diff))
  expect_true(all(c("Group1", "Group2", "Overlap_coef") %in% names(multi_diff$Overlap_coef)))
})

test_that("compute_multi_diff: Diff_proba has n_shared_ids + 1 rows per pair", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:12, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2", "ID_3"))),
    G2 = matrix(1:12, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2", "ID_3")))
  ))
  res <- compute_multi_diff(sd)
  expect_equal(nrow(res$Diff_proba), 4)
})

# -- compute_multi_diff: numeric properties -------------------------------------

test_that("compute_multi_diff: Proba sums to 1 per pair and Cumul_proba is its cumulative sum", {
  multi_diff <- make_multi_diff(nb_id = 4, nb_group = 3, n = 200)
  by_pair <- split(multi_diff$Diff_proba, paste(multi_diff$Diff_proba$Group1, multi_diff$Diff_proba$Group2))
  for (pair in by_pair) {
    expect_equal(sum(pair$Proba), 1, tolerance = 1e-8)
    expect_equal(pair$Cumul_proba, cumsum(pair$Proba))
    expect_equal(pair$Cumul_proba[length(pair$Cumul_proba)], 1, tolerance = 1e-8)
  }
})

test_that("compute_multi_diff: exact Nb_id/Proba/Cumul_proba on a hand-verifiable example", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(c(0, 10, 2.5, 100, 6, 6, 6, 6), nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
  ))
  res <- compute_multi_diff(sd)
  expect_equal(res$Diff_proba$Nb_id, c(0, 1, 2))
  expect_equal(res$Diff_proba$Proba, c(0.25, 0.5, 0.25))
  expect_equal(res$Diff_proba$Cumul_proba, c(0.25, 0.75, 1))
})

# -- compute_multi_diff: integration / determinism ------------------------------

test_that("compute_multi_diff works end-to-end via the real pipeline with more than two groups", {
  multi_diff <- make_multi_diff(nb_id = 4, nb_group = 3, n = 100)
  expect_equal(length(unique(paste(multi_diff$Diff_proba$Group1, multi_diff$Diff_proba$Group2))), choose(3, 2))
  expect_true(all(multi_diff$Overlap_coef$Overlap_coef >= 0 & multi_diff$Overlap_coef$Overlap_coef <= 1))
})

test_that("compute_multi_diff is deterministic given identical inputs", {
  sd <- make_aligned_sample_distrib(list(
    G1 = matrix(1:8, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2"))),
    G2 = matrix(8:1, nrow = 4, dimnames = list(NULL, c("ID_1", "ID_2")))
  ))
  res1 <- compute_multi_diff(sd)
  res2 <- compute_multi_diff(sd)
  expect_equal(res1, res2)
})
