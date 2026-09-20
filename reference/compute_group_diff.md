# Compute a Pairwise Group-Difference Metric

Generalizes
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md)
to any `distance_metric`: for every pair of groups in `results`, aligns
their posterior means/covariances by shared ID and evaluates `metric`
via
[`evaluate_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md).
Every metric requires the two groups being compared to share the same
set of IDs (checked via `setequal`); metrics for which
[`requires_shared_kernel`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
is `TRUE` additionally require the two groups to share the same
`kernel_key`.

## Usage

``` r
compute_group_diff(results, metric, max_groups_warn = 50, max_dim_warn = 500)
```

## Arguments

- results:

  A list, typically from
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md),
  with elements `kernels` and `groups` (one entry per group, each with
  `muk`, `id_to_input`, `kernel_key`, `scale`).

- metric:

  A `distance_metric` object, e.g.
  [`ovl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md),
  [`kl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/kl_metric.md),
  [`wasserstein_metric()`](https://terenceviellard.github.io/BayesOmics/reference/wasserstein_metric.md).

- max_groups_warn:

  Emit a warning about the \\O(G^2 d^3)\\ cost of the underlying matrix
  inversions if the number of groups exceeds this. Defaults to `50`.

- max_dim_warn:

  Emit the same warning if the number of shared IDs \\d\\ exceeds this.
  Defaults to `500`.

## Value

A square matrix of metric values, one row/column per group. The diagonal
is each group's self-comparison value (computed via `evaluate_metric`,
so it is `1` for
[`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
and `0` for distance-like metrics, without special-casing). For metrics
where
[`is_symmetric_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
is `FALSE` (e.g.
[`kl_metric`](https://terenceviellard.github.io/BayesOmics/reference/kl_metric.md)),
the matrix itself is not symmetric.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
compute_group_diff(posterior, mahalanobis_metric())
#>          1        2
#> 1  0.00000 30.46497
#> 2 30.46497  0.00000
```
