# Partition IDs by Slicing Input into Ranges

Cuts `Input` values into intervals (via
[`cut`](https://rdrr.io/r/base/cut.html)); IDs whose Input falls in the
same interval go into the same block.

## Usage

``` r
partition_by_range(
  breaks = NULL,
  n_bins = NULL,
  input_col = "Input",
  input_id = NULL
)
```

## Arguments

- breaks:

  Explicit interval breakpoints, passed to
  [`cut()`](https://rdrr.io/r/base/cut.html). If `NULL` (default),
  `n_bins` equal-width bins spanning the observed range of `Input` are
  used instead.

- n_bins:

  Number of equal-width bins when `breaks` is `NULL`. Defaults to `3`.

- input_col:

  Name of the Input column to slice. Defaults to `"Input"`.

- input_id:

  When `db` has a multi-dimensional Input (an `Input_ID` column with
  more than one distinct value), which axis to slice on. Required in
  that case – without it, guessing an axis could silently cut on the
  wrong one, so
  [`resolve_partition`](https://terenceviellard.github.io/BayesOmics/reference/resolve_partition.md)
  errors instead.

## Value

A `partition_strategy` object usable as
[`fit_block_posterior`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)'s
`partition` argument.
