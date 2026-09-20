# Print a BayesOmics Posterior Object

Pretty-prints the result of
[`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md):
for each group, the posterior mean vector (`muk`) and the reconstructed
(ID-aligned) posterior covariance matrix, obtained via the internal
`get_sigmak()` helper.

## Usage

``` r
# S3 method for class 'bayesomics_posterior'
print(x, digits = 3, ...)
```

## Arguments

- x:

  A list returned by
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md).

- digits:

  Number of significant digits used when rounding the displayed mean
  vector and covariance matrix. Defaults to `3`.

- ...:

  Unused, included for S3 consistency.

## Value

`x`, invisibly.
