# Assemble a Dense Block-Diagonal Covariance Matrix

Assembles a list of per-block covariance matrices (as returned by
[`get_sigmak_block`](https://terenceviellard.github.io/BayesOmics/reference/get_sigmak_block.md))
into one dense block-diagonal matrix – e.g. to feed into the dense
[`evaluate_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
family (such as
[`wasserstein_metric`](https://terenceviellard.github.io/BayesOmics/reference/wasserstein_metric.md))
when a metric has no block-decoupled shortcut of its own.
Off-block-diagonal entries are exactly `0`, not merely small.

## Usage

``` r
as_block_diag_matrix(sigma_blocks)
```

## Arguments

- sigma_blocks:

  A named list of covariance matrices, e.g. from
  [`get_sigmak_block`](https://terenceviellard.github.io/BayesOmics/reference/get_sigmak_block.md).

## Value

A dense `p x p` matrix (`p` = total IDs across blocks), with `dimnames`
set to the IDs in block order.
