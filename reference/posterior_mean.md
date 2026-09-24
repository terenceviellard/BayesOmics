# Compute Posterior Means for Different Groups

This function computes the posterior means for different groups within a
dataset. It expects the input data frame to contain a column named
'Group' that identifies the groups, as well as columns 'ID', 'Output',
and 'Input' – and, for a multi-dimensional Input, an additional
'Input_ID' column (see `normalize_input_cols`): one row per (Group, ID,
Sample, Input_ID), with 'Output' repeated identically across the rows of
one observation. A legacy data frame with just a single 'Input' column
(no 'Input_ID') keeps working unchanged.

## Usage

``` r
posterior_mean(
  data,
  kern = NULL,
  mu_0 = 1,
  lambda_0 = 1,
  obs_noise = 0,
  pooled = TRUE,
  df_warn = 8
)
```

## Arguments

- data:

  A data frame containing the data to be analyzed. Must include columns
  'Group' and 'Output'. 'ID' and 'Input' are optional, but must be
  supplied together or omitted together: if BOTH are missing,
  `posterior_mean()` runs in **univariate mode** (a warning is issued) –
  each group is treated as a single feature (a dummy constant
  'ID'/'Input' is added internally), with no cross-feature correlation
  structure, exactly the `nb_id = 1` case validated in
  `dev/univariate/NOTES_univariate.md`. Supplying only one of the two is
  ambiguous and errors instead. Supply both 'ID' and 'Input' (plus
  'Input_ID' for a multi-dimensional Input) for the general
  multi-feature case.

- kern:

  A kernel object (from the keRnel package) used to compute pairwise
  covariances, or a named list of kernel objects (one entry per group in
  `data`, e.g. from
  [`fit_kernel`](https://terenceviellard.github.io/BayesOmics/reference/fit_kernel.md)`(..., group_col = "Group", pooled = FALSE)`)
  to give each group its own, independently fitted kernel – every group
  then gets its own `kernel_key`, so metrics with
  `requires_shared_kernel() == TRUE` (e.g.
  [`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md))
  cannot compare them (same caveat as `pooled = FALSE` below). Defaults
  to `NULL`, in which case a diagonal
  ([`keRnel::white_noise_kernel()`](https://rdrr.io/pkg/keRnel/man/white_noise_kernel.html))
  kernel is built automatically from a closed-form residual-variance
  estimate – the exact REML minimizer
  `fit_kernel(..., group_col = "Group")` would converge to numerically,
  computed directly instead (see `resolve_closed_form_kernel()` in
  `R/optim_kernel.R` and `dev/univariate/NOTES_univariate.md` for the
  derivation). This is the natural default in univariate mode, but also
  works with a real 'ID'/'Input' design (every feature is then treated
  as independent – no spatial structure – unlike a real kernel fit via
  [`fit_kernel()`](https://terenceviellard.github.io/BayesOmics/reference/fit_kernel.md),
  which must be supplied explicitly via `kern` for that).

- mu_0:

  Prior mean parameter.

- lambda_0:

  Prior precision parameter.

- obs_noise:

  Fixed, known observation-noise variance to add back into the posterior
  covariance (see `get_sigmak()`'s `@details`). Needed whenever `kern`'s
  hyperparameters were fit with `fit_kernel(..., prior_cov = )` treating
  that same value as a separate additive nugget rather than composing it
  into `kern` itself (e.g. via a `NoiseKernel()` term) – otherwise the
  reported credible interval only reflects uncertainty about whatever
  correlated structure the kernel captures BEYOND that fixed noise
  floor, which collapses toward zero whenever there is little such
  structure to find (the exact mechanism behind BayesOmics's severe
  under-coverage on ProteoBayes's own univariate/multivariate paper
  scenarios,
  `dev/benchmark_server/E11_proteobayes_paper_replay_server.R`).
  Defaults to 0 (previous behavior: `kern`'s own HPs assumed to already
  capture the full marginal variance).

- pooled:

  Only used when `kern = NULL`. If `TRUE` (default), a single noise
  variance is estimated and shared by every group (more residual degrees
  of freedom, assumes homogeneous noise across groups). If `FALSE`, each
  group gets its own independently-estimated variance (heteroscedastic);
  groups then never share a `kernel_key`, so metrics with
  `requires_shared_kernel() == TRUE` (e.g.
  [`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md))
  cannot compare them. Ignored when `kern` is supplied directly.

- df_warn:

  Only used when `kern = NULL`. A
  [`warning()`](https://rdrr.io/r/base/warning.html) is issued whenever
  the residual degrees of freedom backing a closed-form variance
  estimate fall below this (the estimate runs but is noisy); an
  `error()` is always raised at 0 or fewer degrees of freedom (the
  estimate is undefined). Defaults to `8` (see
  `dev/optim_exploration/14_group_and_id_demean/README.md`, point 3:
  empirically, the fitted variance's coefficient of variation drops
  below 0.5 around 8 residual degrees of freedom).

## Value

A list with two elements:

- `kernels`:

  A named list of kernel/correlation matrices, one per distinct set of
  Input positions found across groups. Each matrix has `dimnames` set to
  a string key per position (see `row_input_key()`), and is shared by
  reference across every group with that same set of positions (no
  duplication).

- `groups`:

  A named list (one entry per group) with: `muk`, a named vector of
  posterior means keyed by ID; `id_to_input`, an n_ids x D matrix
  (`rownames = ID`) giving that group's ID -\> Input position mapping;
  `kernel_key`, which entry of `kernels` to use; `scale`, the divisor
  (`n_obs + lambda_0`) applied to that kernel matrix to get the
  posterior covariance; and `obs_noise` (this call's value, reused by
  `get_sigmak()`). The (internal) `get_sigmak()` helper reconstructs the
  actual (ID-aligned) posterior covariance matrix for a group.

## Details

The posterior distribution is given by:

\\p(\mathbf{\mu} \mid y_1, \dots, y_N, \Sigma\_{\hat{\theta}}) =
\mathcal{N}\left(\mathbf{\mu}; \\ \dfrac{\lambda_0 \mu_0 +
\sum\_{n=1}^{N} y_n}{N + \lambda_0}, \dfrac{1}{N + \lambda_0}
\Sigma\_{\hat{\theta}}\right)\\

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
posterior$groups[["1"]]$muk
#>      ID_1      ID_2      ID_3      ID_4      ID_5      ID_6      ID_7      ID_8 
#> 21.443325 25.540041  1.661871 39.638712 21.827400 35.635310 36.030173 38.510365 

# Univariate mode (single CpG/feature): no 'ID'/'Input' columns, no kernel
# object needed -- a warning is issued and the noise variance is estimated
# in closed form (pooled across every group by default):
cpg <- data.frame(
  Group  = rep(c("A", "B"), each = 5),
  Sample = rep(1:5, 2),
  Output = c(rnorm(5, 0, 1), rnorm(5, 3, 1))
)
uni_posterior <- posterior_mean(cpg)
#> Warning: posterior_mean(): no 'ID'/'Input' columns found in 'data' -- running in univariate mode (each group treated as a single feature, no cross-feature correlation structure).
```
