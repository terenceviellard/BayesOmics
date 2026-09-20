# Plot overlaid posterior distributions and their overlap

For a given id, plot both groups' posterior density curves on the *same*
panel (unlike
[`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md),
which plots the distribution of their *difference*), with the region
where the two densities overlap shaded and labelled with the OVL
coefficient – i.e. a direct visualization of the quantity
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md)
computes.

Follows the same pairwise ("two by two") logic as
[`plot_distrib_each_pair`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib_each_pair.md):
for exactly two groups (given explicitly via `group1`/`group2`, or the
only two present in `sample_distrib`), a single `ggplot` is returned;
for more than two groups (with neither `group1` nor `group2` given),
every pairwise comparison is returned as a named list of `ggplot`
objects, one per pair.

## Usage

``` r
plot_posterior_overlap(sample_distrib, group1 = NULL, group2 = NULL, id = NULL)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  function, containing the following columns: `ID`, `Group` and
  `Sample`.

- group1:

  A character string, the first group to compare. If NULL (default) and
  `group2` is also NULL, the groups are inferred from `sample_distrib`
  (see Description).

- group2:

  A character string, the second group to compare. If NULL (default),
  see Description.

- id:

  A character string, the id to plot. If NULL (default), only the first
  id appearing in `sample_distrib` is used.

## Value

Either a single `ggplot` (two groups), or a named list of `ggplot`
objects, one per group pair, named `"<group1>_vs_<group2>"` (more than
two groups).

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
plot_posterior_overlap(samples, group1 = "1", group2 = "2", id = unique(samples$ID)[1])
```
