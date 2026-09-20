# Gaussian Overlapping Coefficient (OVL) Metric

The exact closed-form overlapping coefficient between two groups'
posteriors, as used by
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md).
Requires the two groups to share the same `kernel_key` (so \\\Sigma_2 =
c\Sigma_1\\ for a scalar \\c = \mathrm{scale}\_1/\mathrm{scale}\_2\\);
there is no general fallback for this metric.

## Usage

``` r
ovl_metric()
```

## Value

A `distance_metric` object.
