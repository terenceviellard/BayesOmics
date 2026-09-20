# Explicit Manual ID-to-Block Assignment

Explicit Manual ID-to-Block Assignment

## Usage

``` r
partition_by_id(assignment)
```

## Arguments

- assignment:

  A named character (or coercible) vector, `assignment[id]` giving that
  ID's block label. Must cover every ID being partitioned.

## Value

A `partition_strategy` object usable as
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)'s
`partition` argument.
