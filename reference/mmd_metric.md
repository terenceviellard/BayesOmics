# Maximum Mean Discrepancy (MMD) Metric

Squared MMD between two groups' posterior draws under a Gaussian (RBF)
kernel. When `bandwidth` is `NULL` (default), it is set via the median
heuristic (the median pairwise distance in the pooled sample of both
groups).

## Usage

``` r
mmd_metric(bandwidth = NULL)
```

## Arguments

- bandwidth:

  The Gaussian kernel bandwidth, or `NULL` for the median-heuristic
  default.

## Value

A `sample_metric` object.
