# Optimize Hyperparameters for a kernel (with additive kernel support)

Also the single entry point that unifies univariate and multivariate
fitting: `kern = NULL` delegates to the exact same closed-form fit
`posterior_mean(kern = NULL)` uses internally (works with or without
`ID`/`Input` – univariate mode is triggered exactly as it is there, see
`inject_univariate_dummy_cols()`), while a real, unfitted kernel object
fits it via L-BFGS-B. Both branches respect `pooled`/`group_col`
identically in spirit (one shared fit vs. one independent fit per group)
and return the same shape (a single kernel object, or a named list keyed
by group), so
`posterior_mean(data, kern = fit_kernel(kern, data, pooled = pooled, ...))`
is a uniform two-step pipeline regardless of which branch ran.

## Usage

``` r
fit_kernel(
  kern = NULL,
  db,
  prior_mean = NULL,
  prior_cov,
  pen_diag = 1e-06,
  verbose = FALSE,
  max_iter = 1000,
  factr = 1e+07,
  pgtol = 0,
  track_trace = FALSE,
  group_col = NULL,
  pooled = TRUE,
  df_warn = 8
)
```

## Arguments

- kern:

  `NULL` (default) for the closed-form fit (see `@description`), or a
  kernel object inheriting from keRnel's `kernel` class for a real
  kernel fit. In the latter case, its construction values (e.g.
  `se_kernel(length_scale = 2)`) are the optimization's starting point –
  there is no separate `hp` argument. Optimizing in
  `keRnel::get_free_params(kern)`'s free/unconstrained space (rather
  than by natural-space hyperparameter name) is what lets two sibling
  sub-kernels sharing a bare hyperparameter name (e.g.
  `se_kernel(length_scale=1) + se_kernel(length_scale=3)`) be given
  genuinely different starting values and be fit independently – a case
  [`keRnel::kupdate()`](https://rdrr.io/pkg/keRnel/man/kupdate.html)
  cannot express, since it updates every occurrence of a bare name at
  once.

- db:

  The dataset used for optimization. Must contain columns `Group` and
  `Output`; `ID`/`Input` (plus `Input_ID` for a multi-dimensional Input)
  may be omitted together to trigger univariate mode (only meaningful
  when `kern = NULL` – a real kernel has nothing to fit a spatial
  structure to with one dummy position per group). When `kern` is a real
  kernel object, `db` must contain columns `Input` and `Output` (plus
  `Input_ID` for a multi-dimensional Input; see `normalize_input_cols`).
  If it also contains a `Sample` column and the same `Input` positions
  are repeated across samples (i.e. replicated measurements), the
  function automatically uses the correct replicated likelihood: it
  builds a pxp kernel matrix on the *unique* input positions and sums
  independent log N(y_s; mu, K_p) terms over each sample s. Without this
  correction, same-position observations in different samples would be
  treated as perfectly correlated, biasing all hyperparameter estimates.

- prior_mean:

  Prior mean: either a scalar or a vector of length equal to the number
  of unique `Input` positions in `db`. Required unless `group_col` is
  supplied (in which case it is ignored – see `group_col`).

- prior_cov:

  Prior covariance: either a scalar (diagonal value) or a square matrix
  of size equal to the number of unique `Input` positions in `db`.

- pen_diag:

  Jitter added to the diagonal for numerical stability. Defaults to
  `1e-6`.

- verbose:

  If `FALSE` (default), returns the optimized parameter vector; if
  `TRUE`, returns the full `optim` result list (with the fitted kernel
  object attached as `result$kern`).

- max_iter:

  Maximum number of L-BFGS-B iterations. Defaults to `1000`.

- factr:

  L-BFGS-B relative convergence tolerance, passed straight through to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html)'s
  `control$factr`. Smaller values demand tighter convergence (more
  iterations); defaults to `1e7` (the value previously hardcoded).

- pgtol:

  L-BFGS-B projected-gradient convergence tolerance, passed straight
  through to [`stats::optim()`](https://rdrr.io/r/stats/optim.html)'s
  `control$pgtol`. Defaults to `0` (R's own
  [`optim()`](https://rdrr.io/r/stats/optim.html) default, previously
  left unset).

- track_trace:

  If `TRUE`, records every objective/gradient evaluation during the
  optimization (hyperparameter values, NLL or gradient norm, elapsed
  time) and attaches it as a `trace` data.frame. Defaults to `FALSE`, in
  which case nothing is recorded and the return value is identical to
  before this parameter existed.

- group_col:

  Name of a column in `db` (e.g. `"Group"`) identifying which
  experimental group each observation belongs to. When supplied, pools
  replicates from every group into a single fit – centering each
  `(group, ID)` cell on its own empirical mean first (`prior_mean` is
  then ignored, forced to `0`; requires an `ID` column) – and adds the
  exact closed-form REML correction for the degrees of freedom spent
  estimating those group means: \\\text{NLL}\_\text{REML}(\theta) =
  \text{NLL}\_\text{demeaned}(\theta) - (G/2)\log\|\Sigma'\_\theta\|\\,
  \\G\\ = number of distinct groups (NOT groups x ids – the correction
  tracks the rank of the replicate space the mean was projected out of,
  unchanged whether a per-group mean is a scalar or a full per-id
  vector). Without the demeaning step, pooling multiple groups under one
  shared `prior_mean` lets an unmodeled between-group mean shift bias
  the fitted kernel variance; demeaning by the full `(group, ID)` cell
  (rather than by group alone) is required for that correction to hold
  for a realistic, non-uniform differential pattern (only some ids
  shifted) – group-alone demeaning only exactly cancels a shift that is
  the SAME on every id, and otherwise leaves a residual that inflates
  the fitted variance and roughly doubles the package's own OVL metric
  on a realistic partial-shift pattern (confirmed empirically, not just
  a theoretical concern). See `dev/optim_exploration/NOTES_math.md`,
  `dev/optim_exploration/14_group_and_id_demean/README.md` and
  `demeaning_likelihood.tex` for the full derivation and empirical
  validation (HP recovery, end-to-end OVL, joint coverage). Requires a
  balanced Group x Sample design (the same features observed for every
  group/sample combination); errors otherwise. Known pathological
  regimes (not yet guarded against, see the README above): very small
  `N` per `(group, ID)` cell, and extreme `mu_random` prior spread when
  `mu_random = TRUE` was used to simulate the data. Defaults to `NULL`
  (no group handling, behavior identical to before this parameter
  existed).

- pooled:

  When `kern = NULL`: passed straight to `resolve_closed_form_kernel()`
  – `TRUE` (default) shares one estimate across every group, `FALSE`
  gives each group its own (see `resolve_closed_form_kernel()`'s own
  documentation; this is exactly
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)'s
  own `pooled` argument). When `kern` is a real kernel object: only used
  when `group_col` is also supplied. `TRUE` (default) shares a single
  kernel fit across every group (see `group_col`'s own documentation for
  how that pooled fit is built). If `FALSE`, `db` is split by
  `group_col` and `fit_kernel()` is called independently on each group's
  own subset (with that group's `prior_mean`, `group_col = NULL` – no
  demeaning, a genuine per-group MLE/REML fit of the real kernel),
  returning a named list keyed by group instead of a single result.
  Ignored (with no effect) when `kern` is a real kernel object and
  `group_col` is `NULL`.

- df_warn:

  Only used when `kern = NULL`. Forwarded to
  `resolve_closed_form_kernel()`/`check_closed_form_df()`: a
  [`warning()`](https://rdrr.io/r/base/warning.html) is issued below
  this many residual degrees of freedom (the estimate runs but is
  noisy), and an `error()` at 0 or fewer (undefined). Defaults to `8` –
  see
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)'s
  own `df_warn` argument.

## Value

When `kern = NULL`: a single kernel object (`pooled = TRUE`) or a named
list of kernel objects keyed by group (`pooled = FALSE`) – see
`resolve_closed_form_kernel()`; `verbose`/`max_iter`/`factr`/`pgtol`/
`track_trace` are not used on this path. Otherwise (a real kernel
object), if `verbose` is `FALSE`, a named vector of optimized
hyperparameters in natural (constrained) space (via
[`keRnel::get_trainable_params()`](https://rdrr.io/pkg/keRnel/man/get_trainable_params.html)),
with the optimizer's `convergence` code and final objective `value`
attached as attributes (`attr(result, "convergence")`,
`attr(result, "value")`) so convergence can be checked without
re-running with `verbose = TRUE`; a `convergence` of `0` means success.
If `track_trace = TRUE`, a `trace` data.frame (one row per
objective/gradient evaluation) is also attached as
`attr(result, "trace")`. Otherwise (`verbose = TRUE`) the full list
returned by [`stats::optim()`](https://rdrr.io/r/stats/optim.html), with
`result$kern` (the fitted kernel object) and `result$trace` (when
`track_trace = TRUE`) added.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 1, nb_sample = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
fit_kernel(kern, data, prior_mean = 0, prior_cov = 1)
#>     variance length_scale 
#>    850.23062      0.89126 
#> attr(,"convergence")
#> [1] 0
#> attr(,"value")
#> [1] 186.8494

# Pooling replicates from multiple groups into a single fit, without an
# unmodeled group mean shift biasing the fitted variance:
data2 <- simu_db(nb_id = 8, nb_group = 3, nb_sample = 5, diff_group = 4)
fit_kernel(kern, data2, prior_cov = 1, group_col = "Group")
#>     variance length_scale 
#>   2.51962703   0.05659641 
#> attr(,"convergence")
#> [1] 0
#> attr(,"value")
#> [1] 218.6767

# Two independent fits (one per group), each its own real, non-demeaned
# kernel HP optimization -- pooled = FALSE:
fits <- fit_kernel(kern, data2, prior_mean = 0, prior_cov = 1,
                    group_col = "Group", pooled = FALSE, verbose = TRUE)
names(fits)
#> [1] "1" "2" "3"

# kern = NULL: the same pooled/non-pooled choice, but for the closed-form
# fit -- works identically in univariate mode (no ID/Input at all) and in
# multivariate mode (a real ID/Input design), and returns a kernel (or
# named list of kernels) ready for posterior_mean(kern = ...), exactly
# like the real-kernel branch above:
cpg <- data.frame(Group = rep(c("A", "B"), each = 5),
                   Output = c(rnorm(5, 0, 1), rnorm(5, 3, 1)))
fit_kernel(NULL, cpg, pooled = TRUE)   # one shared closed-form kernel
#> Warning: fit_kernel(): no 'ID'/'Input' columns found in 'data' -- running in univariate mode (each group treated as a single feature, no cross-feature correlation structure).
#> white_noise_kernel(noise = 1.011)
fit_kernel(NULL, cpg, pooled = FALSE)  # one closed-form kernel per group
#> Warning: fit_kernel(): no 'ID'/'Input' columns found in 'data' -- running in univariate mode (each group treated as a single feature, no cross-feature correlation structure).
#> Warning: posterior_mean(): only 4 residual degree(s) of freedom for group 'A' (below df_warn = 8); the closed-form variance estimate will be noisy -- consider more replicates, or pooled = TRUE to share degrees of freedom across groups.
#> Warning: posterior_mean(): only 4 residual degree(s) of freedom for group 'B' (below df_warn = 8); the closed-form variance estimate will be noisy -- consider more replicates, or pooled = TRUE to share degrees of freedom across groups.
#> $A
#> white_noise_kernel(noise = 1.075)
#> 
#> $B
#> white_noise_kernel(noise = 0.9469)
#> 
```
