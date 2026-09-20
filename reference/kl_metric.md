# Kullback-Leibler Divergence Metric

The (asymmetric) KL divergence between two groups' Gaussian posteriors,
\\\mathrm{KL}(\mathcal N_1\\\mathcal N_2)\\ (`direction = "1to2"`) or
\\\mathrm{KL}(\mathcal N_2\\\mathcal N_1)\\ (`"2to1"`).

## Usage

``` r
kl_metric(direction = c("1to2", "2to1"))
```

## Arguments

- direction:

  Which direction of the (asymmetric) divergence to compute.

## Value

A `distance_metric` object.
