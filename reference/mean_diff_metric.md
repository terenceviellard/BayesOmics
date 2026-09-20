# Mean Absolute Difference of Means Metric

The mean absolute difference between the two groups' posterior means
over shared IDs, evaluated marginally (per ID) and ignoring the
covariance structure entirely: \\\frac{1}{d}\sum_k
\|\mu\_{1k}-\mu\_{2k}\|\\. `Sigma1`/`Sigma2` are accepted (for interface
consistency with
[`evaluate_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md))
but unused.

## Usage

``` r
mean_diff_metric()
```

## Value

A `distance_metric` object.
