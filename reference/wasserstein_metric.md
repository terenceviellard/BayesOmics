# Wasserstein-2 (Frechet) Distance Metric

The 2-Wasserstein distance between two multivariate Gaussian posteriors.
When the two groups share the same `kernel_key` (so \\\Sigma_2 =
c\Sigma_1\\), a closed-form shortcut is used that avoids a matrix square
root entirely; the general eigendecomposition-based formula is used
otherwise.

## Usage

``` r
wasserstein_metric()
```

## Value

A `distance_metric` object.
