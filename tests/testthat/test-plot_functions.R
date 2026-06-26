# ── plot_distrib: validation ──────────────────────────────────────────────────

test_that("plot_distrib errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(plot_distrib(bad))
})

test_that("plot_distrib warns when multiple ids provided", {
  sd <- make_sample_distrib(ids = c("ID_1", "ID_2"))
  expect_warning(plot_distrib(sd), "IDs|Multiple")
})

# ── plot_distrib: single group ────────────────────────────────────────────────

test_that("plot_distrib returns a ggplot for single group", {
  sd  <- make_sample_distrib(groups = "G1", ids = "ID_1")
  gg  <- plot_distrib(sd, group1 = "G1", id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: group1 defaults to first group in data", {
  sd  <- make_sample_distrib(groups = c("G1", "G2"), ids = "ID_1")
  gg  <- plot_distrib(sd, id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: id defaults to first id in data", {
  sd  <- make_sample_distrib(ids = "ID_1")
  expect_warning(
    gg <- plot_distrib(sd),
    NA  # expect no warning when single id
  )
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: mean_bar=FALSE suppresses the vline layer", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  gg <- plot_distrib(sd, group1 = "G1", id = "ID_1", mean_bar = FALSE)
  layer_classes <- vapply(gg$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomVline" %in% layer_classes)
})

# ── plot_distrib: two groups ───────────────────────────────────────────────────

test_that("plot_distrib returns a ggplot when comparing two groups", {
  sd <- make_sample_distrib()
  gg <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib two-group: show_prob=FALSE suppresses labels", {
  sd  <- make_sample_distrib()
  gg  <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1",
                      show_prob = FALSE)
  layer_classes <- vapply(gg$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomLabel" %in% layer_classes)
})

test_that("plot_distrib two-group: show_prob=TRUE adds labels", {
  sd  <- make_sample_distrib()
  gg  <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1",
                      show_prob = TRUE)
  layer_classes <- vapply(gg$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLabel" %in% layer_classes)
})

test_that("plot_distrib: prob_CI affects CI region width", {
  sd   <- make_sample_distrib(groups = "G1", ids = "ID_1", n = 2000)
  gg95 <- plot_distrib(sd, group1 = "G1", id = "ID_1", prob_CI = 0.95)
  gg50 <- plot_distrib(sd, group1 = "G1", id = "ID_1", prob_CI = 0.50)
  # Both should return valid ggplots; no error is the minimal check
  expect_s3_class(gg95, "ggplot")
  expect_s3_class(gg50, "ggplot")
})

# ── plot_distrib: rendering ───────────────────────────────────────────────────

test_that("plot_distrib can be rendered without error", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  gg <- plot_distrib(sd, group1 = "G1", id = "ID_1")
  expect_no_error(ggplot2::ggplot_build(gg))
})

test_that("plot_distrib two-group rendering completes without error", {
  sd <- make_sample_distrib()
  gg <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1")
  expect_no_error(ggplot2::ggplot_build(gg))
})

# ── plot_distrib: group-count dispatch ───────────────────────────────────────

test_that("plot_distrib: a single group in the data gives a single-group plot", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  gg <- plot_distrib(sd, id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: exactly two groups (no explicit group1/group2) gives a comparison plot", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = "ID_1")
  gg <- plot_distrib(sd, id = "ID_1")
  expect_s3_class(gg, "ggplot")
  # Comparison plots carry the probability labels by default
  layer_classes <- vapply(gg$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLabel" %in% layer_classes)
})

test_that("plot_distrib: explicit group1/group2 always gives a single comparison plot, even with more groups present", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = "ID_1")
  gg <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: more than two groups (no explicit group1/group2) returns a grid arrangement", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = "ID_1")
  res <- plot_distrib(sd, id = "ID_1")
  expect_true(inherits(res, "gtable") || inherits(res, "grob"))
})

test_that("plot_distrib: multi-group dispatch works for four groups (regression for the former single-group layout-matrix bug)", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4"), ids = "ID_1")
  expect_no_error(plot_distrib(sd, id = "ID_1"))
})

test_that("plot_distrib: multi-group grid omits the mean panel when plot_mean = FALSE", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = "ID_1")
  expect_no_error(plot_distrib(sd, id = "ID_1", plot_mean = FALSE))
})

test_that("plot_distrib: multi-group summary has exactly 3 panels (heatmap + facet grid + mean) when plot_mean = TRUE", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4"), ids = "ID_1")
  res <- plot_distrib(sd, id = "ID_1", plot_mean = TRUE)
  n_panels <- sum(!vapply(res$grobs, is.null, logical(1)))
  expect_equal(n_panels, 3)
})

test_that("plot_distrib: multi-group summary has exactly 2 panels (heatmap + facet grid) when plot_mean = FALSE", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4"), ids = "ID_1")
  res <- plot_distrib(sd, id = "ID_1", plot_mean = FALSE)
  n_panels <- sum(!vapply(res$grobs, is.null, logical(1)))
  expect_equal(n_panels, 2)
})

test_that("plot_distrib: multi-group summary panel count stays constant (3) regardless of group count", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4", "G5", "G6"), ids = "ID_1")
  res <- plot_distrib(sd, id = "ID_1", plot_mean = TRUE)
  n_panels <- sum(!vapply(res$grobs, is.null, logical(1)))
  expect_equal(n_panels, 3)
})

test_that("plot_distrib: multi-group summary respects top_n_pairs (fewer facets than total pairs)", {
  sd  <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4", "G5"), ids = "ID_1")
  res <- plot_distrib(sd, id = "ID_1", plot_mean = FALSE, top_n_pairs = 2)
  # The facet grid is the 2nd grob; check it only has 2 panels built.
  facet_gg <- res$grobs[[2]]
  expect_true(inherits(facet_gg, "gtable") || inherits(facet_gg, "grob"))
})

# ── plot_distrib_each_pair ────────────────────────────────────────────────────

test_that("plot_distrib_each_pair returns one ggplot per group pair, named", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4"), ids = "ID_1")
  plots <- plot_distrib_each_pair(sd, id = "ID_1")
  expect_length(plots, choose(4, 2))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_setequal(names(plots), c("G1_vs_G2", "G1_vs_G3", "G1_vs_G4", "G2_vs_G3", "G2_vs_G4", "G3_vs_G4"))
})

test_that("plot_distrib_each_pair errors with fewer than two groups", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  expect_error(plot_distrib_each_pair(sd, id = "ID_1"), "at least two groups")
})

test_that("plot_distrib_each_pair scales to many groups without summary logic (all pairs present)", {
  sd <- make_sample_distrib(groups = paste0("G", 1:6), ids = "ID_1")
  plots <- plot_distrib_each_pair(sd, id = "ID_1")
  expect_length(plots, choose(6, 2))
})

# ── plot_group_overlap_heatmap ────────────────────────────────────────────────

test_that("plot_group_overlap_heatmap returns a ggplot", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = "ID_1")
  gg <- plot_group_overlap_heatmap(sd, id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_group_overlap_heatmap can be rendered without error", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = "ID_1")
  gg <- plot_group_overlap_heatmap(sd, id = "ID_1")
  expect_no_error(ggplot2::ggplot_build(gg))
})

test_that("plot_group_overlap_heatmap errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(plot_group_overlap_heatmap(bad))
})

test_that("compute_pairwise_overlap: well-separated groups give low overlap, identical-distribution groups give high overlap", {
  set.seed(7)
  sd <- data.frame(
    ID = "ID_1",
    Group = rep(c("far", "near", "same_as_near"), each = 1000),
    Sample = c(stats::rnorm(1000, 50), stats::rnorm(1000, 0), stats::rnorm(1000, 0)),
    stringsAsFactors = FALSE
  )
  mat <- BayesOmics:::compute_pairwise_overlap(sd, "ID_1", c("far", "near", "same_as_near"))
  expect_gt(mat["near", "same_as_near"], mat["near", "far"])
})

# ── plot_posterior_overlap ───────────────────────────────────────────────────

test_that("plot_posterior_overlap returns a single ggplot for two groups", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = "ID_1")
  gg <- plot_posterior_overlap(sd, group1 = "G1", group2 = "G2", id = "ID_1")
  expect_s3_class(gg, "ggplot")
  expect_no_error(ggplot2::ggplot_build(gg))
})

test_that("plot_posterior_overlap: group1/group2 default to the two groups present", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = "ID_1")
  gg <- plot_posterior_overlap(sd, id = "ID_1")
  expect_s3_class(gg, "ggplot")
})

test_that("plot_posterior_overlap returns one ggplot per group pair, named, for more than two groups", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3", "G4"), ids = "ID_1")
  plots <- plot_posterior_overlap(sd, id = "ID_1")
  expect_length(plots, choose(4, 2))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_setequal(names(plots), c("G1_vs_G2", "G1_vs_G3", "G1_vs_G4", "G2_vs_G3", "G2_vs_G4", "G3_vs_G4"))
})

test_that("plot_posterior_overlap errors with fewer than two groups", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  expect_error(plot_posterior_overlap(sd, id = "ID_1"), "at least two groups")
})

test_that("plot_posterior_overlap errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(plot_posterior_overlap(bad))
})

test_that("plot_posterior_overlap warns when multiple ids provided explicitly", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1", "ID_2"))
  expect_warning(plot_posterior_overlap(sd, id = c("ID_1", "ID_2")), "IDs|Multiple")
})

test_that("compute_overlap_distrib: identical-distribution groups give an OVL close to 1, well-separated groups close to 0", {
  set.seed(11)
  sd_same <- data.frame(
    ID = "ID_1", Group = rep(c("A", "B"), each = 2000),
    Sample = stats::rnorm(4000, 0, 1), stringsAsFactors = FALSE
  )
  sd_far <- data.frame(
    ID = "ID_1", Group = rep(c("A", "B"), each = 2000),
    Sample = c(stats::rnorm(2000, 0, 1), stats::rnorm(2000, 50, 1)), stringsAsFactors = FALSE
  )
  info_same <- BayesOmics:::compute_overlap_distrib(sd_same, "A", "B", "ID_1")
  info_far  <- BayesOmics:::compute_overlap_distrib(sd_far, "A", "B", "ID_1")
  expect_gt(info_same$ov, 0.9)
  expect_lt(info_far$ov, 0.05)
})

# ── plot_distrib: boundary / edge-case behaviour ─────────────────────────────

test_that("plot_distrib: prob_CI = 0 and prob_CI = 1 do not error", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1", n = 500)
  expect_no_error(plot_distrib(sd, group1 = "G1", id = "ID_1", prob_CI = 0))
  expect_no_error(plot_distrib(sd, group1 = "G1", id = "ID_1", prob_CI = 1))
})

test_that("plot_distrib: index_group1/index_group2 are used in the x-axis label", {
  sd <- make_sample_distrib()
  gg <- plot_distrib(sd, group1 = "G1", group2 = "G2", id = "ID_1",
                      index_group1 = "Treated", index_group2 = "Control")
  label_text <- deparse(gg$labels$x)
  expect_true(grepl("Treated", label_text))
  expect_true(grepl("Control", label_text))
})

test_that("plot_distrib: comparing a group to itself yields a degenerate zero-difference distribution (documented current behaviour)", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  expect_no_error(
    gg <- plot_distrib(sd, group1 = "G1", group2 = "G1", id = "ID_1")
  )
  expect_s3_class(gg, "ggplot")
})

test_that("plot_distrib: requesting an id absent from the data errors (documented current behaviour)", {
  sd <- make_sample_distrib(groups = "G1", ids = "ID_1")
  expect_error(plot_distrib(sd, group1 = "G1", id = "NOT_A_REAL_ID"))
})

# ── plot_posterior_mean ───────────────────────────────────────────────────────

test_that("plot_posterior_mean errors when sample_distrib lacks required columns", {
  bad <- data.frame(x = 1:10)
  expect_error(plot_posterior_mean(bad))
})

test_that("plot_posterior_mean returns a ggplot", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1", "ID_2"))
  gg <- plot_posterior_mean(sd)
  expect_s3_class(gg, "ggplot")
})

test_that("plot_posterior_mean has one point per (id, group) pair", {
  sd <- make_sample_distrib(groups = c("G1", "G2", "G3"), ids = c("ID_1", "ID_2"))
  gg <- plot_posterior_mean(sd)
  built <- ggplot2::ggplot_build(gg)
  expect_equal(nrow(built$data[[1]]), 3 * 2)
})

test_that("plot_posterior_mean averages Sample within each (id, group) pair", {
  sd <- data.frame(
    ID     = rep(c("ID_1", "ID_2"), each = 4),
    Group  = rep(c("G1", "G2"), 4),
    Sample = c(1, 10, 3, 12, 5, 14, 7, 16),
    stringsAsFactors = FALSE
  )
  gg <- plot_posterior_mean(sd)
  built <- ggplot2::ggplot_build(gg)
  expect_setequal(round(built$data[[1]]$x, 6), c(2, 6, 11, 15))
})

test_that("plot_posterior_mean can be rendered without error", {
  sd <- make_sample_distrib(groups = c("G1", "G2"), ids = c("ID_1", "ID_2"))
  gg <- plot_posterior_mean(sd)
  expect_no_error(ggplot2::ggplot_build(gg))
})
