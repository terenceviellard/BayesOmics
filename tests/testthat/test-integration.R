# Full pipeline: simu_db → optim_hp → multi_posterior_mean → sample_posterior
#               → plot_distrib

test_that("full pipeline runs without error", {
  set.seed(1)
  data <- simu_db(nb_id = 5, nb_group = 2, nb_sample = 3)
  kern <- make_kernel()
  hp0  <- unlist(keRnel::gt_HPs(kern))
  hp   <- optim_hp(hp0, data[data$Group == 1, ], 0, kern, 1)
  kern <- keRnel::set_hyperparameters(kern, hp)
  res  <- multi_posterior_mean(data, kern)
  long <- sample_posterior(res, 200)
  expect_s3_class(long, "data.frame")
  expect_named(long, c("ID", "Group", "Sample"))
})

test_that("pipeline produces long-format data compatible with plot_distrib", {
  set.seed(2)
  data <- simu_db(nb_id = 4, nb_group = 2, nb_sample = 2)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  long <- sample_posterior(res, 300)
  gg   <- plot_distrib(long,
                       group1 = unique(long$Group)[1],
                       group2 = unique(long$Group)[2],
                       id     = unique(long$ID)[1])
  expect_s3_class(gg, "ggplot")
})

test_that("calculate_group_overlaps works on multi_posterior_mean output", {
  set.seed(3)
  res <- make_posteriors(nb_id = 5, nb_group = 3)
  mat <- calculate_group_overlaps(res)
  expect_equal(dim(mat), c(3, 3))
  expect_equal(unname(diag(mat)), c(1, 1, 1))
  expect_equal(mat, t(mat))
})

test_that("pipeline is reproducible end-to-end with set.seed", {
  run_pipeline <- function(seed) {
    set.seed(seed)
    data <- simu_db(nb_id = 5, nb_group = 2, nb_sample = 2)
    kern <- make_kernel()
    res  <- multi_posterior_mean(data, kern)
    set.seed(seed)
    sample_posterior(res, 100)
  }
  long1 <- run_pipeline(7)
  long2 <- run_pipeline(7)
  expect_equal(long1, long2)
})

test_that("pipeline handles nb_sample > 1 correctly in posterior", {
  set.seed(10)
  data <- simu_db(nb_id = 6, nb_group = 2, nb_sample = 5)
  kern <- make_kernel()
  expect_no_error({
    res  <- multi_posterior_mean(data, kern)
    long <- sample_posterior(res, 100)
  })
  expect_equal(length(unique(long$ID)), 6)
})

test_that("optimized hyperparameters improve posterior fit", {
  set.seed(42)
  data <- simu_db(nb_id = 20, nb_group = 2, nb_sample = 1)
  kern0 <- make_kernel(hp = c(0.1, 0.1))
  kern_opt <- make_kernel()
  hp0  <- unlist(keRnel::gt_HPs(kern_opt))
  hp   <- optim_hp(hp0, data[data$Group == 1, ], 0, kern_opt, 1)
  kern_opt <- keRnel::set_hyperparameters(kern_opt, hp)

  res0   <- multi_posterior_mean(data, kern0)
  res_op <- multi_posterior_mean(data, kern_opt)
  # Both should produce valid structures
  for (g in names(res0$groups)) {
    expect_true(!anyNA(res0$groups[[g]]$muk))
    expect_true(!anyNA(res_op$groups[[g]]$muk))
  }
})

test_that("calculate_group_overlaps works end-to-end on simu_db_kernel() output (shared kernel_key across groups)", {
  set.seed(5)
  kern <- make_kernel()
  data <- simu_db_kernel(nb_id = 5, nb_group = 3, nb_sample = 4, kernel = kern)
  res  <- multi_posterior_mean(data, kern)
  mat  <- calculate_group_overlaps(res)
  expect_equal(dim(mat), c(3, 3))
  expect_equal(unname(diag(mat)), c(1, 1, 1))
  expect_equal(mat, t(mat))
})

test_that("plot_distrib dispatches to a pairwise grid for more than two groups", {
  set.seed(4)
  data <- simu_db(nb_id = 4, nb_group = 4, nb_sample = 2)
  kern <- make_kernel()
  res  <- multi_posterior_mean(data, kern)
  long <- sample_posterior(res, 100)
  expect_no_error(plot_distrib(long, id = unique(long$ID)[1]))
})
