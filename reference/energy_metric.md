# Energy Distance Metric

The (U-statistic) energy distance between two groups' posterior draws:
\\2\\\overline{\\X-Y\\} - \overline{\\X-X'\\} - \overline{\\Y-Y'\\}\\,
with within-group terms excluding self-pairs. Unlike
[`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md),
it does not require the joint density to overlap in every dimension
simultaneously, so it does not collapse toward a degenerate value as the
number of shared IDs grows.

## Usage

``` r
energy_metric()
```

## Value

A `sample_metric` object.
