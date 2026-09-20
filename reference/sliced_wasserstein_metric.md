# Sliced Wasserstein Distance Metric

Average of the 1D Wasserstein-1 distance (via sorted order
statistics/interpolated quantiles) between two groups' posterior draws,
projected onto `n_projections` random unit directions in ID-space.
Explicitly designed to avoid the curse of dimensionality that afflicts
the general (non-sliced) multivariate optimal-transport distance.

## Usage

``` r
sliced_wasserstein_metric(n_projections = 50)
```

## Arguments

- n_projections:

  Number of random projection directions to average over. Defaults to
  `50`.

## Value

A `sample_metric` object.
