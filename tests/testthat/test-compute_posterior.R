# -- posterior_mean: input validation ------------------------------------

test_that("posterior_mean errors when Group or Output is missing", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(posterior_mean(data[, c("ID", "Output", "Input")], kern), "Group")
  expect_error(posterior_mean(data[, c("Group", "ID", "Input")], kern), "Output")
})

test_that("posterior_mean warns and switches to univariate mode when both ID and Input are missing", {
  data <- make_data()
  kern <- make_kernel()
  expect_warning(posterior_mean(data[, c("Group", "Output")], kern), "univariate")
})

test_that("posterior_mean errors when only one of ID/Input is missing (ambiguous, not univariate)", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(posterior_mean(data[, c("Group", "Output", "Input")], kern), "one of 'ID'/'Input'")
  expect_error(posterior_mean(data[, c("Group", "ID", "Output")], kern), "one of 'ID'/'Input'")
})

test_that("posterior_mean errors on non-kernel argument", {
  data <- make_data()
  expect_error(posterior_mean(data, list()), "keRnel")
})

test_that("posterior_mean errors when lambda_0 <= 0", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(posterior_mean(data, kern, lambda_0 = 0),  "lambda_0")
  expect_error(posterior_mean(data, kern, lambda_0 = -1), "lambda_0")
})

test_that("posterior_mean errors on NA in Group or ID", {
  data <- make_data()
  kern <- make_kernel()
  data_na <- data
  data_na$Group[1] <- NA
  expect_error(posterior_mean(data_na, kern), "NA")
})

test_that("posterior_mean errors when mu_0 is not numeric", {
  data <- make_data()
  kern <- make_kernel()
  expect_error(posterior_mean(data, kern, mu_0 = "0"), "mu_0")
})

test_that("posterior_mean errors when data$Input contains NaN", {
  data <- make_data()
  kern <- make_kernel()
  data$Input[1] <- NaN
  expect_error(posterior_mean(data, kern), "NaN|Inf|NA")
})

test_that("posterior_mean errors when data$Output contains Inf", {
  data <- make_data()
  kern <- make_kernel()
  data$Output[1] <- Inf
  expect_error(posterior_mean(data, kern), "NaN|Inf|NA")
})

test_that("posterior_mean errors when data$Output contains NA", {
  data <- make_data()
  kern <- make_kernel()
  data$Output[1] <- NA_real_
  expect_error(posterior_mean(data, kern), "NaN|Inf|NA")
})

test_that("posterior_mean errors when an ID maps to several distinct Input values within a group", {
  # Two different IDs sharing the same single Input value in a group: 2
  # distinct IDs but only 1 distinct Input position -> muk/id_to_input length mismatch.
  data <- data.frame(
    ID     = c("ID_1", "ID_2"),
    Group  = "A",
    Output = c(0, 0),
    Input  = c(5, 5)
  )
  kern <- make_kernel()
  expect_error(posterior_mean(data, kern), "distinct Input position")
})

# -- posterior_mean: duplicated Input_ID rows within a replicate
# (regression for the silently-wrong muk bug found by dev/multi_agent_audit)

test_that("posterior_mean errors when Input_ID rows are duplicated within a (Group, ID, Sample) replicate", {
  # Every row of Sample 1 duplicated 2x: n_distinct(Input_ID) and n() both
  # double at the (Group, ID) level, so the OLD "is lengths an integer?"
  # check couldn't tell this apart from clean data (2 replicates, 1 row
  # each, vs. 1 replicate duplicated into 2 rows) -- it silently returned
  # muk = 12 instead of the correct 10 for the case below.
  data <- data.frame(
    ID     = c("ID_1", "ID_1", "ID_1"),
    Group  = "A",
    Sample = c(1, 2, 2),
    Output = c(10, 14, 14),
    Input  = c(5, 5, 5)
  )
  kern <- make_kernel()
  expect_error(
    posterior_mean(data, kern),
    "Duplicated.*Input_ID|Input_ID.*[Dd]uplicat"
  )
})

test_that("posterior_mean does not false-positive on genuinely balanced replicates", {
  data <- data.frame(
    ID     = rep(c("ID_1", "ID_2"), each = 3),
    Group  = "A",
    Sample = rep(1:3, times = 2),
    Output = c(10, 11, 9, 20, 21, 19),
    Input  = rep(c(5, 6), each = 3)
  )
  kern <- make_kernel()
  expect_no_error(posterior_mean(data, kern))
})

test_that("posterior_mean: muk is correct (not inflated) on clean replicated data", {
  # Direct numeric check the duplication bug's fix direction targets:
  # lambda_0 = 1, mu_0 = 0, 2 genuine replicates (Output = 10, 14) -> muk =
  # (0 + 10 + 14) / (2 + 1) = 8, not 12 (what the old code silently returned
  # under a 2x row duplication of this exact same data).
  data <- data.frame(
    ID     = c("ID_1", "ID_1"),
    Group  = "A",
    Sample = c(1, 2),
    Output = c(10, 14),
    Input  = c(5, 5)
  )
  kern <- make_kernel()
  res <- posterior_mean(data, kern, mu_0 = 0, lambda_0 = 1)
  expect_equal(unname(res$groups[["A"]]$muk["ID_1"]), 8, tolerance = 1e-10)
})

# -- posterior_mean: muk/sigmak alignment (regression for the former
# alphabetical-ID vs Input-value-order misalignment bug) ---------------------

test_that("posterior_mean: muk and sigmak (via get_sigmak) are correctly aligned by ID", {
  ker <- keRnel::variance_kernel(variance = 10) * keRnel::se_kernel(length_scale = 200)

  # dist(ID_1, ID_2) = 90, dist(ID_1, ID_3) = 100, dist(ID_2, ID_3) = 10
  data <- data.frame(
    ID     = c("ID_1", "ID_2", "ID_3"),
    Group  = "A",
    Output = c(0, 0, 0),
    Input  = c(100, 10, 0)
  )
  res <- posterior_mean(data, ker)
  g   <- res$groups[["A"]]
  sig <- BayesOmics:::get_sigmak(g, res$kernels)

  pk <- function(xi, xj) as.numeric(keRnel::evaluate(ker, as.matrix(xi), as.matrix(xj))) / 2
  true_cov_12 <- pk(100, 10)  # cov(ID_1, ID_2)
  true_cov_13 <- pk(100, 0)   # cov(ID_1, ID_3)
  true_cov_23 <- pk(10, 0)    # cov(ID_2, ID_3)

  ids <- names(g$muk)
  expect_equal(unname(sig[ids == "ID_1", ids == "ID_2"]), true_cov_12, tolerance = 1e-8)
  expect_equal(unname(sig[ids == "ID_1", ids == "ID_3"]), true_cov_13, tolerance = 1e-8)
  expect_equal(unname(sig[ids == "ID_2", ids == "ID_3"]), true_cov_23, tolerance = 1e-8)
})

# -- posterior_mean: return structure -------------------------------------

test_that("posterior_mean returns a list with kernels and groups", {
  res <- make_posteriors()
  expect_true(is.list(res))
  expect_named(res, c("kernels", "groups"))
  expect_named(res$groups)
})

test_that("posterior_mean result has one group entry per group", {
  data <- make_data(nb_group = 3)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  expect_length(res$groups, 3)
})

test_that("posterior_mean each group entry has muk, id_to_input, kernel_key, scale", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    expect_true(all(c("muk", "id_to_input", "kernel_key", "scale") %in% names(res$groups[[g]])))
  }
})

test_that("posterior_mean: muk is a named numeric vector", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    expect_true(is.numeric(res$groups[[g]]$muk))
    expect_true(!is.null(names(res$groups[[g]]$muk)))
  }
})

test_that("posterior_mean: get_sigmak returns a square matrix", {
  res <- make_posteriors()
  for (g in names(res$groups)) {
    sk <- BayesOmics:::get_sigmak(res$groups[[g]], res$kernels)
    expect_true(is.matrix(sk))
    expect_equal(nrow(sk), ncol(sk))
  }
})

test_that("posterior_mean: muk length equals sigmak dimension", {
  res <- make_posteriors(nb_id = 6)
  for (g in names(res$groups)) {
    sk <- BayesOmics:::get_sigmak(res$groups[[g]], res$kernels)
    expect_equal(length(res$groups[[g]]$muk), nrow(sk))
  }
})

test_that("posterior_mean: muk names match ID values", {
  data <- make_data(nb_id = 4)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  ids  <- sort(unique(data$ID))
  for (g in names(res$groups)) {
    expect_equal(sort(names(res$groups[[g]]$muk)), ids)
  }
})

# -- posterior_mean: numerics ---------------------------------------------

test_that("posterior_mean: larger lambda_0 pulls muk toward mu_0", {
  data <- make_data(nb_id = 10)
  kern <- make_kernel()
  mu_0 <- 0
  res_small <- posterior_mean(data, kern, mu_0 = mu_0, lambda_0 = 0.01)
  res_large <- posterior_mean(data, kern, mu_0 = mu_0, lambda_0 = 100)
  muk_small <- unlist(lapply(res_small$groups, `[[`, "muk"))
  muk_large <- unlist(lapply(res_large$groups, `[[`, "muk"))
  # with large lambda_0, muk is pulled toward mu_0=0
  expect_lt(mean(abs(muk_large)), mean(abs(muk_small)) + 0.1 + abs(mean(muk_small)))
})

test_that("posterior_mean: kernel cache is shared (by reference) across groups with identical inputs", {
  data <- make_data(nb_id = 5, nb_group = 2)
  # Force both groups to share the exact same inputs
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, each = 2)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  expect_equal(res$groups[["1"]]$kernel_key, res$groups[["2"]]$kernel_key)
  expect_length(res$kernels, 1)
  sig1 <- BayesOmics:::get_sigmak(res$groups[["1"]], res$kernels)
  sig2 <- BayesOmics:::get_sigmak(res$groups[["2"]], res$kernels)
  expect_equal(unname(sig1), unname(sig2), tolerance = 1e-10)
})

test_that("posterior_mean: kernel cache is shared (by reference) across three groups with identical inputs", {
  data <- make_data(nb_id = 5, nb_group = 3)
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, length.out = nrow(data))
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  keys <- vapply(res$groups, `[[`, character(1), "kernel_key")
  expect_length(unique(keys), 1)
  expect_length(res$kernels, 1)
})

test_that("posterior_mean: get_sigmak reconstructs a fresh matrix that does not mutate the shared cache", {
  data <- make_data(nb_id = 5, nb_group = 2)
  shared_inputs <- sort(unique(data$Input))[seq_len(5)]
  data$Input <- rep(shared_inputs, each = 2)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)

  sig1 <- BayesOmics:::get_sigmak(res$groups[["1"]], res$kernels)
  sig1[1, 1] <- -999  # mutate the value returned to the caller

  sig2 <- BayesOmics:::get_sigmak(res$groups[["2"]], res$kernels)
  expect_false(any(sig2 == -999))
})

test_that("posterior_mean: Group converted to character if integer", {
  data <- make_data(nb_id = 4, nb_group = 2)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  expect_true(all(vapply(names(res$groups), is.character, logical(1))))
})

# -- posterior_mean: edge cases ------------------------------------------

test_that("posterior_mean works with nb_id = 1", {
  data <- make_data(nb_id = 1, nb_group = 2)
  kern <- make_kernel()
  expect_no_error(posterior_mean(data, kern))
})

test_that("posterior_mean works with a single group", {
  data <- make_data(nb_id = 5, nb_group = 1)
  kern <- make_kernel()
  res  <- posterior_mean(data, kern)
  expect_length(res$groups, 1)
})

# -- sample_posterior: validation ----------------------------------------------

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

# -- sample_posterior: return structure -----------------------------------------

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

# -- posterior_mean: univariate mode / closed-form kernel (kern = NULL) ------

test_that("posterior_mean univariate mode matches the manual conjugate update", {
  set.seed(1)
  cpg <- data.frame(
    Group  = rep(c("A", "B", "C"), each = 6),
    Sample = rep(1:6, 3),
    Output = c(rnorm(6, 0, 1.5), rnorm(6, 3, 1.5), rnorm(6, -2, 1.5))
  )
  mu_0 <- mean(cpg$Output)
  post <- suppressWarnings(posterior_mean(cpg, mu_0 = mu_0, lambda_0 = 1))

  cell_means <- tapply(cpg$Output, cpg$Group, mean)
  ssr <- sum((cpg$Output - cell_means[cpg$Group])^2)
  closed_var <- ssr / (nrow(cpg) - length(unique(cpg$Group)))

  for (g in unique(cpg$Group)) {
    y_g <- cpg$Output[cpg$Group == g]
    manual_mu  <- (1 * mu_0 + sum(y_g)) / (length(y_g) + 1)
    manual_var <- closed_var / (length(y_g) + 1)
    sigmak <- BayesOmics:::get_sigmak(post$groups[[g]], post$kernels)
    expect_equal(dim(sigmak), c(1L, 1L))
    expect_equal(unname(post$groups[[g]]$muk), manual_mu)
    expect_equal(as.numeric(sigmak), manual_var)
  }
})

test_that("posterior_mean warns exactly once and switches to univariate mode without ID/Input", {
  cpg <- data.frame(Group = rep(c("A", "B"), each = 5), Sample = rep(1:5, 2), Output = rnorm(10))
  expect_warning(posterior_mean(cpg), "univariate")
  post <- suppressWarnings(posterior_mean(cpg))
  expect_equal(names(post$groups[["A"]]$muk), "1")
})

test_that("posterior_mean pooled = TRUE shares one kernel_key across groups, pooled = FALSE does not", {
  set.seed(2)
  cpg <- data.frame(
    Group  = rep(c("A", "B", "C"), each = 6),
    Sample = rep(1:6, 3),
    Output = rnorm(18)
  )
  post_pooled <- suppressWarnings(posterior_mean(cpg, pooled = TRUE))
  keys_pooled <- vapply(post_pooled$groups, function(g) g$kernel_key, character(1))
  expect_length(unique(keys_pooled), 1)

  post_unpooled <- suppressWarnings(posterior_mean(cpg, pooled = FALSE))
  keys_unpooled <- vapply(post_unpooled$groups, function(g) g$kernel_key, character(1))
  expect_length(unique(keys_unpooled), 3)
})

test_that("posterior_mean pooled = FALSE: ovl_metric() refuses to compare groups, other metrics still work", {
  set.seed(3)
  cpg <- data.frame(
    Group  = rep(c("A", "B"), each = 6),
    Sample = rep(1:6, 2),
    Output = c(rnorm(6, 0, 1), rnorm(6, 3, 4))
  )
  post <- suppressWarnings(posterior_mean(cpg, pooled = FALSE))
  expect_error(compute_group_diff(post, ovl_metric()), "kernel matrix")
  d <- compute_group_diff(post, mahalanobis_metric())
  expect_equal(dim(d), c(2L, 2L))
})

test_that("posterior_mean closed-form kernel (kern = NULL) is diagonal for real multi-ID data", {
  data <- make_data(nb_id = 6, nb_group = 2, nb_sample = 5)
  post <- posterior_mean(data)
  sigmak <- BayesOmics:::get_sigmak(post$groups[["1"]], post$kernels)
  expect_equal(dim(sigmak), c(6L, 6L))
  expect_true(all(sigmak[upper.tri(sigmak)] == 0))
  expect_true(all(diag(sigmak) > 0))
})

test_that("posterior_mean errors when closed-form degrees of freedom are <= 0", {
  cpg <- data.frame(Group = c("A", "B"), Sample = c(1, 1), Output = c(1, 5))
  expect_error(suppressWarnings(posterior_mean(cpg)), "degrees of freedom")
})

test_that("posterior_mean warns when closed-form degrees of freedom are below df_warn", {
  # ID/Input supplied explicitly (single dummy feature) so only the
  # degrees-of-freedom warning fires, not the separate univariate-mode one.
  cpg <- data.frame(
    Group = rep(c("A", "B"), each = 3), ID = "1", Input = 0,
    Sample = rep(1:3, 2), Output = rnorm(6)
  )
  expect_warning(posterior_mean(cpg), "residual degree")
})

test_that("posterior_mean kern = NULL still requires an explicit non-kernel argument to error", {
  data <- make_data()
  expect_error(posterior_mean(data, list()), "keRnel")
})

# -- fit_block_posterior(): partition strategies -----------------------------

test_that("resolve_partition() resolves each strategy correctly", {
  ids <- paste0("ID_", 1:6)
  manual <- stats::setNames(rep(c("A", "B", "C"), 2), ids)
  expect_equal(unname(resolve_partition(partition_by_id(manual), ids)), unname(manual[ids]))

  expect_equal(unname(resolve_partition(NULL, ids)), rep("ALL", length(ids)))
  expect_equal(unname(resolve_partition(partition_dense(), ids)), rep("ALL", length(ids)))
  expect_equal(length(unique(resolve_partition(partition_diagonal(ids), ids))), length(ids))

  annot <- data.frame(ID = ids, Block = rep(c("X", "Y"), 3))
  res_annot <- resolve_partition(partition_from_annotation(annot), ids)
  expect_equal(unname(res_annot[annot$ID]), annot$Block)
})

test_that("resolve_partition() errors clearly on missing/ambiguous input", {
  ids <- paste0("ID_", 1:4)
  expect_error(resolve_partition(partition_by_id(c(ID_1 = "A")), ids), "missing")

  data_2d <- data.frame(
    ID = rep(ids, 2), Input_ID = rep(1:2, each = 4), Input = rnorm(8)
  )
  expect_error(resolve_partition(partition_by_range(n_bins = 2), ids, data_2d), "input_id")
  expect_error(resolve_partition(partition_dbscan(eps = 1, minPts = 3), ids), "not implemented")
})

# -- fit_block_posterior(): numerical equivalence with direct posterior_mean() --

test_that("fit_block_posterior() with nothing specified stays agnostic to the non-partitioned baseline", {
  data <- make_data(nb_id = 8, nb_group = 2, nb_sample = 5)
  data$Group <- as.character(data$Group)
  fit <- fit_block_posterior(data, partition_dense())
  direct <- posterior_mean(data)
  for (g in unique(data$Group)) {
    expect_equal(fit$block_results[["ALL"]]$groups[[g]]$muk, direct$groups[[g]]$muk)
  }
})

test_that("fit_block_posterior() pooled = FALSE matches posterior_mean() called directly per (Group, block)", {
  data <- make_data(nb_id = 9, nb_group = 2, nb_sample = 5)
  data$Group <- as.character(data$Group)
  ids <- unique(data$ID)
  partition <- partition_by_id(stats::setNames(rep(c("b1", "b2", "b3"), length.out = length(ids)), ids))
  fit <- fit_block_posterior(data, partition, kern = NULL, pooled = FALSE)

  for (g in unique(data$Group)) {
    for (b in fit$blocks) {
      ids_b <- names(fit$block_of_id)[fit$block_of_id == b]
      sub_direct <- data[data$ID %in% ids_b & data$Group == g, , drop = FALSE]
      if (nrow(sub_direct) == 0) next
      direct <- posterior_mean(sub_direct, kern = NULL, pooled = FALSE)
      from_wrapper <- fit$block_results[[paste(b, g, sep = "__")]]
      expect_equal(from_wrapper$groups[[g]]$muk[ids_b], direct$groups[[g]]$muk[ids_b])
    }
  }
})

test_that("fit_block_posterior() pooled = TRUE matches posterior_mean(pooled = TRUE) called directly per block", {
  data <- make_data(nb_id = 9, nb_group = 2, nb_sample = 5)
  data$Group <- as.character(data$Group)
  ids <- unique(data$ID)
  partition <- partition_by_id(stats::setNames(rep(c("b1", "b2"), length.out = length(ids)), ids))
  fit <- fit_block_posterior(data, partition, kern = NULL, pooled = TRUE)

  for (b in fit$blocks) {
    ids_b <- names(fit$block_of_id)[fit$block_of_id == b]
    sub_direct <- data[data$ID %in% ids_b, , drop = FALSE]
    direct <- posterior_mean(sub_direct, kern = NULL, pooled = TRUE)
    for (g in unique(sub_direct$Group)) {
      expect_equal(fit$block_results[[b]]$groups[[g]]$muk[ids_b], direct$groups[[g]]$muk[ids_b])
    }
  }
})

test_that("partition_diagonal() gives the same result as one posterior_mean() call per ID", {
  data <- make_data(nb_id = 5, nb_group = 2, nb_sample = 5)
  data$Group <- as.character(data$Group)
  ids <- unique(data$ID)
  fit <- fit_block_posterior(data, partition_diagonal(ids), kern = NULL, pooled = FALSE)

  id <- ids[1]; g <- unique(data$Group)[1]
  sub_direct <- data[data$ID == id & data$Group == g, , drop = FALSE]
  direct <- posterior_mean(sub_direct, kern = NULL, pooled = FALSE)
  from_wrapper <- fit$block_results[[paste(id, g, sep = "__")]]
  expect_equal(unname(from_wrapper$groups[[g]]$muk[id]), unname(direct$groups[[g]]$muk[id]))
})

test_that("fit_block_posterior() works with an un-coerced integer Group column (simu_db()'s raw output)", {
  data <- make_data(nb_id = 6, nb_group = 2, nb_sample = 5)
  expect_true(is.integer(data$Group) || is.numeric(data$Group))
  ids <- unique(data$ID)
  partition <- partition_by_id(stats::setNames(rep(c("b1", "b2"), length.out = length(ids)), ids))

  fit <- fit_block_posterior(data, partition, kern = NULL, pooled = FALSE)
  for (g in unique(data$Group)) {
    for (b in fit$blocks) {
      from_wrapper <- fit$block_results[[paste(b, g, sep = "__")]]
      expect_true(as.character(g) %in% names(from_wrapper$groups))
    }
  }
})

# -- fit_block_posterior(): downstream consumers -----------------------------

test_that("get_sigmak_block()/get_muk_block()/as_block_diag_matrix() reconstruct a valid block-diagonal covariance", {
  data <- make_data(nb_id = 8, nb_group = 2, nb_sample = 5)
  data$Group <- as.character(data$Group)
  ids <- unique(data$ID)
  partition <- partition_by_id(stats::setNames(rep(c("b1", "b2"), length.out = length(ids)), ids))
  fit <- fit_block_posterior(data, partition, kern = NULL, pooled = FALSE)

  g <- unique(data$Group)[1]
  Sigma_blocks <- get_sigmak_block(fit, g)
  muk <- get_muk_block(fit, g)
  expect_equal(length(muk), length(ids))
  expect_equal(sort(names(muk)), sort(ids))

  Sigma_dense <- as_block_diag_matrix(Sigma_blocks)
  expect_equal(dim(Sigma_dense), c(length(ids), length(ids)))
  n1 <- nrow(Sigma_blocks[[1]])
  n_total <- nrow(Sigma_dense)
  expect_true(all(Sigma_dense[1:n1, (n1 + 1):n_total] == 0))
  expect_true(all(diag(Sigma_dense) > 0))
})
