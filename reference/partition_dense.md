# Dense Partition (a Single Block)

The "no partitioning at all" extreme: every ID goes into one block, so
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)
degrades exactly to a single call to
[`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
on the whole dataset. Equivalent to passing `partition = NULL`.

## Usage

``` r
partition_dense()
```

## Value

`NULL` (resolved by
[`resolve_partition()`](https://terenceviellard.github.io/BayesOmics/reference/resolve_partition.md)'s
default method).
