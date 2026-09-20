# Marginal-Averaged Metric

Averages a base metric's *univariate* value over each shared ID
individually, rather than evaluating it jointly on the full
\\d\\-dimensional posterior. This sidesteps metrics like
[`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
requiring *every* dimension to overlap simultaneously, at the cost of
losing joint (cross-ID correlation) information. Forwards
[`requires_shared_kernel`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)/
[`is_symmetric_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
to the wrapped metric.

## Usage

``` r
marginal_metric(base)
```

## Arguments

- base:

  A `distance_metric` object to marginalize.

## Value

A `distance_metric` object.
