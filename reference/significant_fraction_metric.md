# Fraction of Individually Significant Features

For each shared ID, computes the z-score
\\(\mu\_{1k}-\mu\_{2k})/\sqrt{\Sigma\_{1,kk}+\Sigma\_{2,kk}}\\ and
returns the proportion of IDs whose absolute z-score exceeds `threshold`
– a nonparametric, per-feature summary that stays a proportion in \\\[0,
1\]\\ regardless of \\d\\.

## Usage

``` r
significant_fraction_metric(threshold = 1.96)
```

## Arguments

- threshold:

  The absolute z-score threshold above which an ID counts as
  individually significant. Defaults to `1.96`.

## Value

A `distance_metric` object.
