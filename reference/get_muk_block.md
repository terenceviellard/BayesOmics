# Concatenate Per-Block Posterior Means for a Group

Concatenate Per-Block Posterior Means for a Group

## Usage

``` r
get_muk_block(fit, group)
```

## Arguments

- fit:

  A `block_posterior_fit` object from
  [`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md).

- group:

  Which group's blocks to concatenate.

## Value

A single named numeric vector (posterior mean per ID), concatenated
across every block that group participates in.
