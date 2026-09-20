# Resolve a Partition Strategy into a Block Assignment

Generic dispatching on the S3 class of `strategy` (a
`partition_strategy` object, e.g. from
[`partition_by_id`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_id.md)).
Not normally called directly –
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)
calls it internally – but exported so custom partition strategies can be
tested in isolation, or implemented by defining a new
`resolve_partition.*` method.

## Usage

``` r
resolve_partition(strategy, ids, db = NULL)
```

## Arguments

- strategy:

  A `partition_strategy` object, or `NULL` (one single dense block, see
  [`partition_dense`](https://terenceviellard.github.io/BayesOmics/reference/partition_dense.md)).

- ids:

  A character vector of all IDs to assign.

- db:

  The full data frame being partitioned (only needed by strategies that
  inspect `Input`/annotation columns, e.g.
  [`partition_by_range`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_range.md)).

## Value

A named character vector, `block_of_id[id]`, giving each ID's block
label.
