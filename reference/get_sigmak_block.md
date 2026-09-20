# Reconstruct Per-Block Posterior Covariances for a Group

Calls the real, unchanged (internal) `get_sigmak()` on each of a group's
sub-problem results, collecting one (ID-aligned) covariance matrix per
block.

## Usage

``` r
get_sigmak_block(fit, group)
```

## Arguments

- fit:

  A `block_posterior_fit` object from
  [`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md).

- group:

  Which group's blocks to reconstruct.

## Value

A named list of covariance matrices, one per block that group
participates in.
