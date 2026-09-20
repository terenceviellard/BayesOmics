# Fit Independent Posteriors per Block

Splits `data` according to `partition`, then calls the unchanged
[`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)/[`optim_hp`](https://terenceviellard.github.io/BayesOmics/reference/optim_hp.md)
pipeline independently per sub-problem (per block, or per
`(Group, block)` when `pooled = FALSE`), collecting the results into a
list consultable via
[`get_sigmak_block`](https://terenceviellard.github.io/BayesOmics/reference/get_sigmak_block.md)/[`get_muk_block`](https://terenceviellard.github.io/BayesOmics/reference/get_muk_block.md).
No result fusion happens here – each sub-problem stays a genuine
`bayesomics_posterior` object, so no numeric primitive is reimplemented.

Passing
`partition = `[`partition_dense()`](https://terenceviellard.github.io/BayesOmics/reference/partition_dense.md)
(or, equivalently, `NULL`) collapses this to a single sub-problem
covering the whole dataset – with the same `pooled` default as
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
(`TRUE`), this reproduces `posterior_mean(data)` exactly, so calling
`fit_block_posterior()` with nothing else specified stays agnostic to
the non-partitioned baseline.

## Usage

``` r
fit_block_posterior(
  data,
  partition,
  kern = NULL,
  pooled = TRUE,
  mu_0 = 1,
  lambda_0 = 1,
  obs_noise = 0,
  df_warn = 8,
  prior_mean = 0,
  prior_cov = 1,
  pen_diag = 1e-06
)
```

## Arguments

- data:

  A data frame, same requirements as
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md).

- partition:

  A `partition_strategy` object (e.g.
  [`partition_by_id`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_id.md),
  [`partition_by_range`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_range.md),
  [`partition_from_annotation`](https://terenceviellard.github.io/BayesOmics/reference/partition_from_annotation.md),
  [`partition_diagonal`](https://terenceviellard.github.io/BayesOmics/reference/partition_diagonal.md)),
  or `NULL` /
  [`partition_dense()`](https://terenceviellard.github.io/BayesOmics/reference/partition_dense.md)
  for a single block.

- kern:

  `NULL` (closed form, as in
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)),
  or an UNFITTED kernel object (template) – in that case
  [`optim_hp()`](https://terenceviellard.github.io/BayesOmics/reference/optim_hp.md)
  is called once per sub-problem before building its posterior
  (mirroring how
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
  never fits `kern` itself when it is supplied, see `R/optim_kernel.R`).

- pooled:

  A single parameter threaded to both branches: when `kern = NULL`,
  passed straight through to
  `posterior_mean(kern = NULL, pooled = pooled)`; when `kern` is an
  unfitted template, `TRUE` fits one `optim_hp(group_col = "Group")`
  shared across every Group within a block, `FALSE` fits an independent
  [`optim_hp()`](https://terenceviellard.github.io/BayesOmics/reference/optim_hp.md)
  per `(Group, block)`. Defaults to `TRUE`, matching
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)'s
  own default.

- mu_0, lambda_0, obs_noise, df_warn:

  Forwarded to
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
  on every sub-problem; see its documentation.

- prior_mean, prior_cov, pen_diag:

  Forwarded to
  [`optim_hp()`](https://terenceviellard.github.io/BayesOmics/reference/optim_hp.md)
  when `kern` is an unfitted template; see its documentation.

## Value

A `block_posterior_fit` object: a list with `block_results` (named list
of per-sub-problem `bayesomics_posterior` objects), `block_index` (named
list per Group, mapping block label -\> `block_results` key),
`block_of_id` (the resolved ID -\> block assignment), `blocks` (block
labels, in natural-sort order), and `pooled`.

## Examples

``` r
data <- simu_db(nb_id = 12, nb_group = 2, nb_sample = 5)
ids <- unique(data$ID)
partition <- partition_by_id(stats::setNames(rep(c("A", "B"), length.out = length(ids)), ids))
fit <- fit_block_posterior(data, partition, kern = NULL, pooled = TRUE)
fit$blocks
#> [1] "A" "B"
```
