# Diagonal Partition (One Block per ID)

The opposite extreme from
[`partition_dense`](https://terenceviellard.github.io/BayesOmics/reference/partition_dense.md):
every ID is its own singleton block, so every block is a `nb_id = 1`
problem.

## Usage

``` r
partition_diagonal(ids)
```

## Arguments

- ids:

  Character vector of IDs to assign one-per-block.

## Value

A `partition_strategy` object usable as
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)'s
`partition` argument.
