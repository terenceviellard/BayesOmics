# Partition IDs from an External Annotation

Partition IDs from an External Annotation

## Usage

``` r
partition_from_annotation(mapping, by = "ID")
```

## Arguments

- mapping:

  Either a data frame with a `by` column and a `Block` column, or a
  named vector (ID -\> block label).

- by:

  Name of the ID column in `mapping` when it is a data frame. Defaults
  to `"ID"`.

## Value

A `partition_strategy` object usable as
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)'s
`partition` argument.
