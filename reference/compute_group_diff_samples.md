# Compute a Pairwise Draw-Based Group-Difference Metric

The draw-based counterpart to
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md):
for every pair of groups present in `sample_distrib`, aligns their
posterior draws by shared ID (via the same melt-order invariant
documented for
[`compute_multi_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_multi_diff.md))
and evaluates `metric` via
[`evaluate_sample_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_sample_metric.md).
All `sample_metric`s are symmetric, and the diagonal is fixed at `0` by
convention (a group's distance to itself), rather than computed from
finite-sample draws, which would be a slightly-biased nonzero estimate
rather than exactly `0`.

## Usage

``` r
compute_group_diff_samples(sample_distrib, metric)
```

## Arguments

- sample_distrib:

  A data frame, typically from
  [`sample_posterior`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md),
  with columns `ID`, `Group`, `Sample`.

- metric:

  A `sample_metric` object, e.g.
  [`energy_metric()`](https://terenceviellard.github.io/BayesOmics/reference/energy_metric.md),
  [`mmd_metric()`](https://terenceviellard.github.io/BayesOmics/reference/mmd_metric.md),
  [`sliced_wasserstein_metric()`](https://terenceviellard.github.io/BayesOmics/reference/sliced_wasserstein_metric.md).

## Value

A square, symmetric matrix of metric values, one row/column per group,
diagonal 0.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 300)
compute_group_diff_samples(samples, energy_metric())
#>          1        2
#> 1  0.00000 19.81038
#> 2 19.81038  0.00000
```
