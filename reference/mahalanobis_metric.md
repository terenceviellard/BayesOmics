# Mahalanobis Distance Metric

The Mahalanobis distance between two groups' posterior means under
\\\Sigma_1\\: \\D =
\sqrt{(\mu_1-\mu_2)^\top\Sigma_1^{-1}(\mu_1-\mu_2)}\\. This is the same
quantity
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md)
computes internally before converting it to an overlap coefficient.
Asymmetric whenever \\\Sigma_1 \neq \Sigma_2\\ – swapping the two groups
uses \\\Sigma_2^{-1}\\ instead, generally a different value – so
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
computes both directions separately rather than mirroring one across the
diagonal (see `is_symmetric_metric.mahalanobis_metric` below).

## Usage

``` r
mahalanobis_metric()
```

## Value

A `distance_metric` object.
