# Sample from a Normal multivariate distribution

Sample n elements from the posterior distribution of each group, and
reshape the result directly into the long-format data frame expected by
[`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md).

## Usage

``` r
sample_posterior(results, n)
```

## Arguments

- results:

  A list returned by
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md),
  with elements `kernels` and `groups`.

- n:

  A number indicating the number of samples

## Value

A data frame with columns `ID`, `Group`, and `Sample` (one row per draw,
per ID, per group).

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 100)
head(samples)
#>     ID Group   Sample
#> 1 ID_1     1 21.09615
#> 2 ID_1     1 21.40886
#> 3 ID_1     1 21.31273
#> 4 ID_1     1 21.92299
#> 5 ID_1     1 21.72313
#> 6 ID_1     1 20.97156
```
