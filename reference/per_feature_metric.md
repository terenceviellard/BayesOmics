# Per-Feature-Normalized Metric

Divides a base metric's value by \\d^{power}\\ (\\d\\ = number of shared
IDs), turning a joint quantity into a per-feature rate that stays
roughly constant (rather than growing/collapsing) as \\d\\ changes under
a constant per-feature effect. Forwards
[`requires_shared_kernel`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)/[`is_symmetric_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
to the wrapped metric.

## Usage

``` r
per_feature_metric(base, power)
```

## Arguments

- base:

  A `distance_metric` object to normalize.

- power:

  The exponent applied to \\d\\ (e.g. `0.5` for \\\sqrt d\\, `1` for
  \\d\\).

## Value

A `distance_metric` object.
