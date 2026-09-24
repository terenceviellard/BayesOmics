# Troubleshooting: common errors and how to fix them

``` r

knitr::opts_chunk$set(warning = FALSE, message = FALSE)
library(keRnel)
library(BayesOmics)
```

This page lists the errors and warnings you are most likely to run into,
what they mean, and how to fix them. Each one is reproduced below with
[`try()`](https://rdrr.io/r/base/try.html) so the message matches
exactly what you would see in your own console.

## Optimization did not converge

``` r

data <- simu_db(nb_id = 20)
kern <- variance_kernel(variance = 0.01) * se_kernel(length_scale = 0.01)
opt  <- fit_kernel(kern, data, prior_mean = 0, prior_cov = 1, max_iter = 1)
#> Warning in fit_kernel(kern, data, prior_mean = 0, prior_cov = 1, max_iter = 1):
#> L-BFGS-B did not converge (code 1). Results may be unreliable. Consider
#> increasing max_iter or adjusting pen_diag.
```

**Cause:** L-BFGS-B stopped before satisfying its convergence criterion
– usually because `max_iter` is too low, or the data/kernel combination
is poorly conditioned (`pen_diag` too small, or the kernel’s initial
construction values far from a good optimum).

**Fix:** increase `max_iter`, try different initial values when
constructing `kern`, or increase `pen_diag`. You don’t need
`verbose = TRUE` to check: the returned vector always carries the
optimizer’s status as attributes:

``` r

attr(opt, "convergence")  # 0 means success
#> [1] 1
attr(opt, "value")        # final (negative log-likelihood) objective value
#> [1] 1025.449
```

## Groups do not share the same set of IDs

``` r

bad_results <- list(
  kernels = list(k = diag(2)),
  groups  = list(
    A = list(muk = c(ID_1 = 0, ID_2 = 1), id_to_input = c(ID_1 = 1, ID_2 = 2), kernel_key = "k", scale = 1),
    B = list(muk = c(ID_1 = 0, ID_3 = 1), id_to_input = c(ID_1 = 1, ID_3 = 2), kernel_key = "k", scale = 1)
  )
)
group_diff(bad_results)
#> Error in `kern_mat[input_keys, input_keys, drop = FALSE]`:
#> ! no 'dimnames' attribute for array
```

**Cause:**
[`group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/group_diff.md)
(like every group-comparison function) compares two groups ID by ID; it
requires every group to be measured on the exact same set of IDs.

**Fix:** restrict the comparison to groups that share all their IDs, or
make sure your input data frame has one row per (ID, Group, Sample) for
every ID in every group before calling
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md).

## Groups do not share the same kernel matrix (no longer an error for OVL)

Two groups can have different sets of `Input` values (e.g. different
genomic positions), or be fitted with genuinely different kernels (e.g.
`posterior_mean(pooled = FALSE)`), so their posterior covariances don’t
derive from the same raw kernel matrix.
[`ovl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
(and therefore
[`group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/group_diff.md)/[`calculate_group_overlaps()`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md))
handles this automatically: it uses the exact closed-form formula when
the two groups do share one kernel matrix, and falls back to a Monte
Carlo/KDE estimate otherwise – no error, no metric to switch by hand.
See
[`?ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
for how that fallback works and its accuracy tradeoff (noisier, and
increasingly unreliable as the number of shared IDs grows).

A custom metric can still opt into the strict behavior by defining its
own
[`requires_shared_kernel()`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
method returning `TRUE`;
[`compute_group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
then errors for that metric specifically when the two groups don’t share
a kernel matrix, with no fallback.

See
[`vignette("03_pooled_vs_nonpooled")`](https://terenceviellard.github.io/BayesOmics/articles/03_pooled_vs_nonpooled.md)
for when `pooled = FALSE` produces this on purpose, and
[`vignette("05_metric_choice")`](https://terenceviellard.github.io/BayesOmics/articles/05_metric_choice.md)
for why the OVL fallback exists and its accuracy/cost tradeoff in more
detail.

## Each ID must map to a single Input value within its group

``` r

data_bad <- data.frame(
  ID     = c("ID_1", "ID_1"),
  Group  = "A",
  Output = c(0, 0),
  Input  = c(5, 6)  # same ID, two different Input values
)
posterior_mean(data_bad, se_kernel(length_scale = 1))
#> Error in `purrr::map()`:
#> ℹ In index: 1.
#> Caused by error in `input_matrix_by_id()`:
#> ! Duplicated 'Input_ID' value '1' for observation 'ID_1' (key columns: ID).
```

**Cause:** within a single group, `ID` must uniquely determine `Input`
(e.g. one CpG site cannot sit at two different genomic positions in the
same condition). This usually comes from a data-cleaning issue –
duplicated rows, or an `ID` reused across two physically different
positions.

**Fix:** check for duplicate/inconsistent `(Group, ID)` rows in your
input data frame with `dplyr::distinct(data, Group, ID, Input)` and
resolve the inconsistency before calling
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md).

## Group has IDs with different numbers of observations

``` r

data_uneven <- data.frame(
  ID     = c("ID_1", "ID_1", "ID_2"),
  Group  = "A",
  Sample = c(1, 2, 1),
  Output = c(0, 1, 0),
  Input  = c(5, 5, 6)
)
posterior_mean(data_uneven, se_kernel(length_scale = 1))
#> Error in `purrr::map_int()`:
#> ℹ In index: 1.
#> Caused by error in `.f()`:
#> ! Group 'A' has IDs with different numbers of observations (1, 2); posterior_mean() requires every ID within a group to have the same number of observations.
```

**Cause:**
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
uses a single observation count per group (the posterior `scale` is
`n_obs + lambda_0`); it requires every ID within a group to have the
same number of replicate samples.

**Fix:** make sure every ID in a group has the same number of `Sample`
rows – e.g. by filling in missing replicates or subsetting to a balanced
design. See
[`vignette("07_unequal_sample_size")`](https://terenceviellard.github.io/BayesOmics/articles/07_unequal_sample_size.md)
for designs with different replicate counts *across* groups (which is
supported), as opposed to across IDs *within* the same group (which is
not).

## Matrix could not be made positive-definite

**Cause:** internally, `chol_inv_jitter()` adds a small jitter to a
covariance matrix’s diagonal to make it invertible, increasing the
jitter geometrically up to 20 times. If it still fails, the matrix is
likely severely ill-conditioned – e.g. near-duplicate `Input` values
combined with a very large kernel length scale.

**Fix:** increase `pen_diag`, reduce the kernel’s length scale, or check
for near-duplicate `Input` values in your data.

## Empty groups

**Cause:** a `Group` value in your data has no `ID` or no `Input` value
associated with it – usually an artefact of filtering the data frame
before calling
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md).

**Fix:** drop unused factor levels / filter out the offending rows
before fitting,
e.g. `data <- droplevels(data[data$Group %in% kept_groups, ])`.

## Go further

| Vignette | What you will find |
|----|----|
| [`vignette("01_basic_pipeline")`](https://terenceviellard.github.io/BayesOmics/articles/01_basic_pipeline.md) | The complete two-group walkthrough – kernel choice, hyperparameter fitting, posterior computation, [`group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/group_diff.md). Start here for the full model pipeline. |
| [`vignette("02_univariate")`](https://terenceviellard.github.io/BayesOmics/articles/02_univariate.md) | Comparing two groups with no feature axis at all (a single scalar per replicate). |
| [`vignette("03_pooled_vs_nonpooled")`](https://terenceviellard.github.io/BayesOmics/articles/03_pooled_vs_nonpooled.md) | Fitting one shared kernel vs. one independent kernel per group – and exactly when that makes two groups stop sharing a kernel matrix (see above). |
| [`vignette("04_multi_group")`](https://terenceviellard.github.io/BayesOmics/articles/04_multi_group.md) | Reading a [`group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/group_diff.md) matrix and [`plot_multi_diff()`](https://terenceviellard.github.io/BayesOmics/reference/plot_multi_diff.md) output for more than two groups. |
| [`vignette("05_metric_choice")`](https://terenceviellard.github.io/BayesOmics/articles/05_metric_choice.md) | Why the joint OVL saturates as the number of features grows, when to prefer per-feature Wasserstein instead, and the logic [`group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/group_diff.md) automates. |
| [`vignette("06_kernel_choice")`](https://terenceviellard.github.io/BayesOmics/articles/06_kernel_choice.md) | Choosing a kernel other than Squared Exponential. |
| [`vignette("07_unequal_sample_size")`](https://terenceviellard.github.io/BayesOmics/articles/07_unequal_sample_size.md) | How unbalanced designs (different replicate counts *per group*) affect posterior width and the final comparison, and why it remains valid even with strongly unequal counts. |
| [`vignette("08_multi_dim_input")`](https://terenceviellard.github.io/BayesOmics/articles/08_multi_dim_input.md) | `Input` with more than one dimension (`Input_ID`), and why the fit code needs no change. |
| [`vignette("09_block_diagonal")`](https://terenceviellard.github.io/BayesOmics/articles/09_block_diagonal.md) | Scaling up to many features via a block-diagonal partition ([`fit_block_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)). |
| [`vignette("10_real_data")`](https://terenceviellard.github.io/BayesOmics/articles/10_real_data.md) | Reshaping a real dataset (not a simulation) into the BayesOmics long format, including a ready-to-use AI-assistant prompt. |
