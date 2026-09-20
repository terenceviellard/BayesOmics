# Plot the posterior mean as a function of id

Display, for every id present in `sample_distrib`, the posterior mean of
its distribution (averaged over the drawn samples), coloured by group.
This gives a region-wide view of the posterior profile (e.g. one point
per CpG site in a methylation analysis), complementing the per-id
distribution plots produced by
[`plot_distrib()`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md).

## Usage

``` r
plot_posterior_mean(sample_distrib)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  function, containing the following columns: `ID`, `Group` and
  `Sample`. This argument should contain the empirical posterior
  distributions to be summarized.

## Value

A `ggplot` object, with one point per (id, group) pair, the id on the
y-axis, the posterior mean on the x-axis, and colour indicating the
group.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
plot_posterior_mean(samples)
```
