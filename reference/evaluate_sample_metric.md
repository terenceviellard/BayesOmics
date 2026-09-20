# Draw-Based Group-Difference Metrics

A `sample_metric` is a lightweight S3 object describing which draw-based
(rather than closed-form Gaussian) distance to compute between two
groups' posterior samples. Unlike `distance_metric` objects, these don't
need `mu`/`Sigma` at all – they work directly from posterior draws (e.g.
from
[`sample_posterior`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)),
which is also why they're driven by
[`compute_group_diff_samples`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff_samples.md)
rather than
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md).

## Usage

``` r
evaluate_sample_metric(metric, draws1, draws2, ...)
```

## Arguments

- metric:

  A `sample_metric` object.

- draws1, draws2:

  Numeric matrices of posterior draws for the two groups (rows = draws,
  columns = shared IDs, same column order). Row counts may differ
  between the two groups.

- ...:

  Unused; included for extensibility.

## Value

`evaluate_sample_metric()` returns a single numeric value.
