# Compare Two Groups' Posteriors, with a Sensible Default Metric

A convenience wrapper around
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md).
When `metric = NULL` (the default), the metric is chosen automatically
from the number of shared IDs `d`:
[`ovl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
for `d <= id_threshold` (still interpretable, not yet
dimension-collapsed), or
[`per_feature_metric`](https://terenceviellard.github.io/BayesOmics/reference/per_feature_metric.md)`(`[`wasserstein_metric`](https://terenceviellard.github.io/BayesOmics/reference/wasserstein_metric.md)`(), power = 0.5)`
otherwise – a per-feature-normalized distance that keeps growing
informatively with `d` instead of collapsing toward a degenerate value
(see `dev/30_package_dev/distance_metrics/distance-metrics.Rmd` for the
empirical study behind this choice, and the "Choice of metric" article).
Supplying `metric` explicitly bypasses this choice entirely and behaves
exactly like
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md).

The returned matrix carries the metric actually used (as a label) and
prints it above the matrix, so an auto-selected default is never
silently invisible even if the result is stored and printed again later.

## Usage

``` r
group_diff(results, metric = NULL, id_threshold = 6, ...)
```

## Arguments

- results:

  A list, typically from
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md);
  see
  [`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md).

- metric:

  A `distance_metric` object, or `NULL` (default) to auto-select based
  on the number of shared IDs (see `id_threshold`).

- id_threshold:

  The number of shared IDs at or below which
  [`ovl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
  is auto-selected (above it,
  `per_feature_metric(wasserstein_metric(), power = 0.5)` is used
  instead). Only used when `metric = NULL`. Defaults to `6`.

- ...:

  Forwarded to
  [`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
  (e.g. `max_groups_warn`, `max_dim_warn`).

## Value

A `group_diff_result` object: the same matrix
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
returns (still usable as a plain matrix, e.g. `result["g1", "g2"]`),
with an extra `"metric_label"` attribute and its own `print` method.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
group_diff(posterior)                # d = 8 > 6 -> per_feature_metric(wasserstein_metric(), power = 0.5)
#> Metric: per_feature_metric(wasserstein_metric(), power = 0.5) 
#> 
#>          1        2
#> 1 0.000000 4.436319
#> 2 4.436319 0.000000
group_diff(posterior, ovl_metric())  # explicit metric, bypasses auto-selection
#> Metric: ovl_metric(n_mc = 2000) 
#> 
#>              1            2
#> 1 1.000000e+00 3.920899e-49
#> 2 3.920899e-49 1.000000e+00
```
