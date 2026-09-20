# Symmetrized Jeffreys Divergence Metric

\\J = \mathrm{KL}(\mathcal N_1\\\mathcal N_2) + \mathrm{KL}(\mathcal
N_2\\\mathcal N_1)\\, computed by summing two
[`kl_metric`](https://terenceviellard.github.io/BayesOmics/reference/kl_metric.md)
evaluations rather than re-deriving the formula.

## Usage

``` r
jeffreys_metric()
```

## Value

A `distance_metric` object.
