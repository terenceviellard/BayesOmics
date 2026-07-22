# -- compute_group_diff: input validation --------------------------------------

test_that("compute_group_diff errors when top-level structure is wrong", {
  expect_error(compute_group_diff(list(foo = 1), mahalanobis_metric()), "kernels.*groups|groups.*kernels")
})

test_that("compute_group_diff errors when a group is missing required elements", {
  bad <- list(
    kernels = list(),
    groups  = list(
      g1 = list(muk = c(0), id_to_input = c(a = 1), kernel_key = "k", scale = 1),
      g2 = list(muk = c(0))
    )
  )
  expect_error(compute_group_diff(bad, mahalanobis_metric()), "muk.*id_to_input|id_to_input.*muk")
})

test_that("compute_group_diff errors when 'metric' is not a distance_metric", {
  res <- make_simple_results(n_groups = 2, n_ids = 1)
  expect_error(compute_group_diff(res, "not a metric"), "distance_metric")
})

# -- compute_group_diff: return structure --------------------------------------

test_that("compute_group_diff returns a square matrix with group dimnames", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  mat <- compute_group_diff(res, mahalanobis_metric())
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(3, 3))
  expect_equal(rownames(mat), names(res$groups))
  expect_equal(colnames(mat), names(res$groups))
})

test_that("compute_group_diff with a single group returns a 1x1 matrix", {
  res <- make_simple_results(n_groups = 1, n_ids = 2)
  mat <- compute_group_diff(res, mahalanobis_metric())
  expect_equal(dim(mat), c(1, 1))
  expect_equal(mat[1, 1], 0)
})

# -- calculate_group_overlaps delegates to compute_group_diff(..., ovl_metric()) --

test_that("calculate_group_overlaps matches compute_group_diff(..., ovl_metric())", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  expect_equal(calculate_group_overlaps(res), compute_group_diff(res, ovl_metric()))
})

# -- diagonal: computed generically via self-comparison, not hardcoded ---------

test_that("ovl_metric diagonal is 1 (self-overlap)", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  mat <- compute_group_diff(res, ovl_metric())
  expect_equal(unname(diag(mat)), rep(1, 3))
})

test_that("distance-like metrics have a diagonal near 0 (self-distance)", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  for (m in list(mahalanobis_metric(), kl_metric(), jeffreys_metric(),
                 bhattacharyya_metric(), hellinger_metric(), wasserstein_metric())) {
    mat <- compute_group_diff(res, m)
    expect_true(all(abs(diag(mat)) < 1e-3), info = paste(class(m)[1], "diagonal not near 0"))
  }
})

# -- mahalanobis_metric: hand-verified values ----------------------------------

test_that("mahalanobis_metric matches hand-computed 1D value", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(1, 1, 1))
  ))
  mat <- compute_group_diff(res, mahalanobis_metric())
  expect_equal(mat["g1", "g2"], 2, tolerance = 1e-5)
})

test_that("mahalanobis_metric is symmetric", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 1, ID_2 = -2), sigma = diag(c(2, 3))),
    g2 = list(muk = c(ID_1 = -1, ID_2 = 1), sigma = diag(c(2, 3)))
  ))
  mat <- compute_group_diff(res, mahalanobis_metric())
  expect_equal(mat["g1", "g2"], mat["g2", "g1"])
})

# -- kl_metric: hand-verified, asymmetric ---------------------------------------

test_that("kl_metric matches hand-computed 1D values in both directions", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1))
  ))
  mat12 <- compute_group_diff(res, kl_metric("1to2"))
  mat21 <- compute_group_diff(res, kl_metric("2to1"))
  expect_equal(mat12["g1", "g2"], 0.8181472, tolerance = 1e-5)
  expect_equal(mat21["g1", "g2"], 2.806853, tolerance = 1e-5)
})

test_that("kl_metric is not symmetric and is_symmetric_metric() reflects that", {
  expect_false(is_symmetric_metric(kl_metric()))
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1))
  ))
  mat <- compute_group_diff(res, kl_metric("1to2"))
  expect_false(isTRUE(all.equal(mat["g1", "g2"], mat["g2", "g1"])))
})

test_that("kl_metric() only accepts its documented directions", {
  expect_error(kl_metric("bogus"))
})

# -- jeffreys_metric: sum of both KL directions --------------------------------

test_that("jeffreys_metric equals KL(1to2) + KL(2to1)", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1))
  ))
  jef  <- compute_group_diff(res, jeffreys_metric())["g1", "g2"]
  kl12 <- compute_group_diff(res, kl_metric("1to2"))["g1", "g2"]
  kl21 <- compute_group_diff(res, kl_metric("2to1"))["g1", "g2"]
  expect_equal(jef, kl12 + kl21, tolerance = 1e-8)
  expect_equal(jef, 3.625, tolerance = 1e-5)
})

test_that("jeffreys_metric is symmetric", {
  expect_true(is_symmetric_metric(jeffreys_metric()))
})

# -- bhattacharyya_metric / hellinger_metric: hand-verified, bounded -----------

test_that("bhattacharyya_metric and hellinger_metric match hand-computed 1D values", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1))
  ))
  db <- compute_group_diff(res, bhattacharyya_metric())["g1", "g2"]
  h  <- compute_group_diff(res, hellinger_metric())["g1", "g2"]
  expect_equal(db, 0.3115718, tolerance = 1e-5)
  expect_equal(h, 0.5174021, tolerance = 1e-5)
  expect_equal(h, sqrt(1 - exp(-db)), tolerance = 1e-8)
})

test_that("hellinger_metric stays in [0, 1] for well-separated groups", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = -100), sigma = matrix(0.01, 1, 1)),
    g2 = list(muk = c(ID_1 =  100), sigma = matrix(0.01, 1, 1))
  ))
  h <- compute_group_diff(res, hellinger_metric())["g1", "g2"]
  expect_true(h >= 0 && h <= 1)
})

# -- wasserstein_metric: shortcut path and general path agree ------------------

test_that("wasserstein_metric general-path (different kernel_key) matches hand-computed value", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1))
  ))
  w2 <- compute_group_diff(res, wasserstein_metric())["g1", "g2"]
  expect_equal(w2, 2.236068, tolerance = 1e-6)
})

test_that("wasserstein_metric shortcut path (shared kernel_key, different scale) matches the general-path value", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(4, 1, 1), scale = 1),
    g2 = list(muk = c(ID_1 = 2), sigma = matrix(4, 1, 1), scale = 4)
  ))
  w2 <- compute_group_diff(res, wasserstein_metric())["g1", "g2"]
  expect_equal(w2, 2.236068, tolerance = 1e-6)
})

test_that("wasserstein_metric is symmetric", {
  expect_true(is_symmetric_metric(wasserstein_metric()))
})

# -- per_feature_metric: exact division identity -------------------------------

test_that("per_feature_metric divides the base metric's value by d^power exactly", {
  res <- make_simple_results(n_groups = 2, n_ids = 4)
  base_mat <- compute_group_diff(res, mahalanobis_metric())
  norm_mat <- compute_group_diff(res, per_feature_metric(mahalanobis_metric(), power = 0.5))
  expect_equal(norm_mat["g1", "g2"], base_mat["g1", "g2"] / 4^0.5, tolerance = 1e-10)
})

test_that("per_feature_metric forwards requires_shared_kernel/is_symmetric_metric to its base", {
  expect_true(requires_shared_kernel(per_feature_metric(ovl_metric(), power = 1)))
  expect_false(requires_shared_kernel(per_feature_metric(mahalanobis_metric(), power = 1)))
  expect_false(is_symmetric_metric(per_feature_metric(kl_metric(), power = 1)))
  expect_true(is_symmetric_metric(per_feature_metric(mahalanobis_metric(), power = 1)))
})

test_that("per_feature_metric() rejects a base that isn't a distance_metric", {
  expect_error(per_feature_metric("not a metric", power = 1), "distance_metric")
})

# -- marginal_metric: exact averaging identity ---------------------------------

test_that("marginal_metric equals the mean of the base metric evaluated per shared ID", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0, ID_2 = 0), sigma = diag(c(1, 4))),
    g2 = list(muk = c(ID_1 = 1, ID_2 = 3), sigma = diag(c(1, 4)))
  ))
  marg <- compute_group_diff(res, marginal_metric(mahalanobis_metric()))["g1", "g2"]

  d1 <- evaluate_metric(mahalanobis_metric(), c(ID_1 = 0), c(ID_1 = 1), matrix(1, 1, 1), matrix(1, 1, 1))
  d2 <- evaluate_metric(mahalanobis_metric(), c(ID_2 = 0), c(ID_2 = 3), matrix(4, 1, 1), matrix(4, 1, 1))
  expect_equal(marg, mean(c(d1, d2)), tolerance = 1e-10)
})

test_that("marginal_metric(ovl_metric()) resists the joint-overlap collapse relative to plain OVL", {
  # A scenario with several shared IDs and a real per-feature effect: joint
  # OVL requires every dimension to overlap simultaneously and collapses
  # toward 0, while the marginal average of per-ID OVLs does not.
  set.seed(42)
  data <- simu_db(nb_id = 10, nb_group = 2, nb_sample = 5, diff_group = 3)
  kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
  posterior <- multi_posterior_mean(data, kern)

  joint_ovl    <- calculate_group_overlaps(posterior)["1", "2"]
  marginal_ovl <- compute_group_diff(posterior, marginal_metric(ovl_metric()))["1", "2"]
  expect_gt(marginal_ovl, joint_ovl)
})

test_that("marginal_metric forwards requires_shared_kernel/is_symmetric_metric to its base", {
  expect_true(requires_shared_kernel(marginal_metric(ovl_metric())))
  expect_false(is_symmetric_metric(marginal_metric(kl_metric())))
})

# -- significant_fraction_metric: exact proportion -----------------------------

test_that("significant_fraction_metric returns the exact proportion of IDs exceeding the threshold", {
  # z-scores: (0-1)/sqrt(1+1) = -0.707 (below threshold), (0-5)/sqrt(1+1) = -3.54 (above)
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0, ID_2 = 0), sigma = diag(c(1, 1))),
    g2 = list(muk = c(ID_1 = 1, ID_2 = 5), sigma = diag(c(1, 1)))
  ))
  frac <- compute_group_diff(res, significant_fraction_metric(threshold = 1.96))["g1", "g2"]
  expect_equal(frac, 0.5, tolerance = 1e-10)
})

test_that("significant_fraction_metric stays in [0, 1]", {
  res <- make_simple_results(n_groups = 2, n_ids = 5)
  frac <- compute_group_diff(res, significant_fraction_metric())["g1", "g2"]
  expect_true(frac >= 0 && frac <= 1)
})

# -- ovl_metric: requires_shared_kernel enforced by the driver -----------------

test_that("compute_group_diff errors for ovl_metric when groups don't share the same kernel_key", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 0), sigma = matrix(2, 1, 1))
  ))
  expect_error(compute_group_diff(res, ovl_metric()), "kernel matrix|kernel_key")
})

test_that("compute_group_diff does not require a shared kernel_key for metrics that don't need it", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1, 1, 1)),
    g2 = list(muk = c(ID_1 = 0), sigma = matrix(2, 1, 1))
  ))
  expect_no_error(compute_group_diff(res, mahalanobis_metric()))
})

test_that("requires_shared_kernel is TRUE only for ovl_metric among the base metrics", {
  expect_true(requires_shared_kernel(ovl_metric()))
  for (m in list(mahalanobis_metric(), kl_metric(), jeffreys_metric(),
                 bhattacharyya_metric(), hellinger_metric(), wasserstein_metric())) {
    expect_false(requires_shared_kernel(m), info = class(m)[1])
  }
})

# -- edge cases -----------------------------------------------------------------

test_that("mismatched ID sets raise a clear, package-level error", {
  res <- make_custom_results(list(
    A = list(muk = c(ID_1 = 0, ID_2 = 5), sigma = diag(2)),
    B = list(muk = c(ID_3 = 0, ID_4 = 5), sigma = diag(2))
  ))
  expect_error(compute_group_diff(res, mahalanobis_metric()), "do not share the same set of IDs")
})

test_that("compute_group_diff does not crash with near-singular Sigma", {
  res <- make_custom_results(list(
    g1 = list(muk = c(ID_1 = 0), sigma = matrix(1e-8, 1, 1)),
    g2 = list(muk = c(ID_1 = 0), sigma = matrix(1e-8, 1, 1))
  ))
  for (m in list(mahalanobis_metric(), kl_metric(), bhattacharyya_metric(), wasserstein_metric())) {
    expect_no_error(mat <- compute_group_diff(res, m))
    expect_true(is.finite(mat["g1", "g2"]), info = class(m)[1])
  }
})

test_that("compute_group_diff is deterministic (no randomness) for analytic metrics", {
  res <- make_simple_results(n_groups = 3, n_ids = 2)
  mat1 <- compute_group_diff(res, wasserstein_metric())
  mat2 <- compute_group_diff(res, wasserstein_metric())
  expect_equal(mat1, mat2)
})

test_that("compute_group_diff warns when the number of groups exceeds max_groups_warn", {
  res <- make_simple_results(n_groups = 5, n_ids = 1)
  expect_warning(compute_group_diff(res, mahalanobis_metric(), max_groups_warn = 4), "O\\(G\\^2")
})
