# Hellinger Distance Metric

\\H = \sqrt{1 - e^{-D_B}}\\, derived from the
[`bhattacharyya_metric`](https://terenceviellard.github.io/BayesOmics/reference/bhattacharyya_metric.md)
distance \\D_B\\. Bounded in \\\[0, 1\]\\.

## Usage

``` r
hellinger_metric()
```

## Value

A `distance_metric` object.
