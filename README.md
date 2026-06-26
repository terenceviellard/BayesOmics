
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesOmics

<!-- badges: start -->

<!-- badges: end -->

BayesOmics is a Bayesian differential analysis method for omics data
that involve capturing the correlation structure through a kernel
function. In this approach, users are required to select an appropriate
kernel, optimize its hyperparameters, and subsequently utilize the
correlation matrix derived from this kernel to compute the posterior.

## See it in action

Given a long-format dataset (one row per id/group/sample), BayesOmics
fits a single kernel-structured posterior per group and turns it into a
constant-size overlap heatmap and a per-id posterior-mean profile,
however many groups are present:

``` r
library(keRnel)
library(BayesOmics)

set.seed(42)
data <- simu_db(nb_id = 25, nb_group = 4, nb_sample = 3, diff_group = 8)

kern <- new("SEKernel")
kern <- set_hyperparameters(kern, c(1.0, 1.0))
opt  <- optim_hp(c(1.0, 1.0), data, rep(0, nrow(data)), kern, diag(0.1, nrow(data)))
kern <- set_hyperparameters(kern, opt)

posterior <- multi_posterior_mean(data, kern)
samples   <- sample_posterior(posterior, n = 2000)
```

``` r
plot_group_overlap_heatmap(samples, id = unique(samples$ID)[1])
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

``` r
plot_posterior_mean(samples)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="" width="100%" />

## Installation

You can install the development version of BayesOmics from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("terenceviellard/BayesOmics")
```

## Bayes Omics in a nutshell

This is a basic example which shows you how to solve a common problem:

Generate a synthetic dataset with the correct format:

``` r
set.seed(123)
data <- simu_db(nb_id = 400,
                nb_group = 2,
                nb_sample = 1)
head(data)
#>     ID Group Sample     Input   Output
#> 1 ID_1     1      1 49.302715 14.23176
#> 2 ID_1     2      1 49.302715 15.04157
#> 3 ID_2     1      1  6.853374 38.14576
#> 4 ID_2     2      1  6.853374 42.35757
#> 5 ID_3     1      1 45.265479 21.79024
#> 6 ID_3     2      1 45.265479 20.14775
```

We choose a Squared Exponential (SE) kernel, given by:

$$
K_{\text{SE}}(x, x') = \sigma^2 \exp\left(-\frac{\|x - x'\|^2}{2\ell^2}\right)
$$

where:

- $K_{\text{SE}}(x, x')$ is the covariance between inputs $x$ and $x'$

- $\sigma^2$ is the signal variance.

- $\ell$ is the length scale.

- $\|x - x'\|$ is the Euclidean distance between $x$ and $x'$.

``` r
SEKernel <- new("SEKernel")
```

We initialize the hyperparameters with $\sigma^2$ = 1 and $\ell$ = 1.

``` r
hp <- c(1.0, 1.0)
SEKernel <- set_hyperparameters(SEKernel, hp)
```

We set identical priors for the different groups.

``` r
n <- dim(data)[1]
prior_mean_nul <- rep(0, n)
prior_cov_nul <- diag(0.1, n)
```

Then we need to optimize $\sigma^2$ and $\ell$ :

``` r
opt= optim_hp(
  hp = hp,
  db = data,
  prior_mean = prior_mean_nul,
  kern = SEKernel,
  prior_cov = prior_cov_nul
)
opt
#> [1] 900.575825   0.000001
#> attr(,"convergence")
#> [1] 0
#> attr(,"value")
#> [1] 18142.41
```

We create a new kernel with optimal hyperparameters.

``` r
SEKernelopt <- set_hyperparameters(SEKernel, opt)
```

The parameters of all posterior distributions can be computed. Each
group’s posterior mean follows a Normal-Normal conjugate update:

$$
p(\mathbf{\mu} \mid y_1, \dots, y_N, \Sigma_{\hat{\theta}}) = \mathcal{N}\left(\mathbf{\mu}; \ \dfrac{\lambda_0 \mu_0 + \sum_{n=1}^{N} y_n}{N + \lambda_0}, \dfrac{1}{N + \lambda_0} \Sigma_{\hat{\theta}}\right)
$$

where $\Sigma_{\hat{\theta}}$ is the correlation matrix built from the
optimized kernel, $\mu_0$ and $\lambda_0$ are the prior mean and prior
precision, and $N$ is the number of observations.

``` r
posterior <- multi_posterior_mean(data, SEKernelopt)
```

Once parameters of the posterior distributions are available, we can
calculate the overlapping coefficient (OVL) between the two groups, used
here as our differential analysis statistic:

$$
OVL = 2 \Phi\!\left(-\frac{\delta}{2}\right)
$$

where $\delta$ is the Mahalanobis distance between the two groups’
posterior means under their pooled covariance, and $\Phi$ is the
standard normal CDF. $OVL = 1$ means the two posterior distributions are
identical (no differential signal); $OVL
\approx 0$ means they are fully separated (strong differential signal).

``` r
calculate_group_overlaps(posterior)
#>           1         2
#> 1 1.0000000 0.3434029
#> 2 0.3434029 1.0000000
```

## Visualizing results

Once a posterior is computed, draw samples from it with
`sample_posterior()` and feed them to any of the `plot_*()` functions
below. They all share the same long-format input (`ID`, `Group`,
`Sample`), so they compose freely.

``` r
samples <- sample_posterior(posterior, n = 2000)
```

`plot_posterior_overlap()` overlays the two groups’ posterior densities
for one id and shades the probability mass used to compute the OVL
coefficient above:

``` r
plot_posterior_overlap(samples, id = unique(samples$ID)[1])
```

<img src="man/figures/README-unnamed-chunk-14-1.png" alt="" width="100%" />

`plot_posterior_mean()` gives a region-wide view: one point per id,
coloured by group, so you can spot which ids drive the differential
signal:

``` r
plot_posterior_mean(samples)
```

<img src="man/figures/README-unnamed-chunk-15-1.png" alt="" width="100%" />

`plot_group_overlap_heatmap()` summarizes every pairwise group
comparison for one id as a single tile grid – the same plot already
shown above, useful here as it stays readable regardless of how many
groups are compared:

``` r
plot_group_overlap_heatmap(samples, id = unique(samples$ID)[1])
```

<img src="man/figures/README-unnamed-chunk-16-1.png" alt="" width="100%" />

## Bayes Omics for multigroup analysis

We can use Bayes Omics to compare more than 2 groups. Generate a
synthetic dataset with 3 groups with the correct format:

``` r
data <- simu_db(nb_id = 200,
                nb_group = 3,
                nb_sample = 1)
```

``` r
SEKernel <- new("SEKernel")

hp <- c(1.0, 1.0)
SEKernel <- set_hyperparameters(SEKernel, hp)

n <- dim(data)[1]

prior_mean_nul <- rep(0, n)
prior_cov_nul <- diag(0.1, n)

opt <- optim_hp(
  hp = hp,
  db = data,
  prior_mean = prior_mean_nul,
  kern = SEKernel,
  prior_cov = prior_cov_nul
)
opt
#> [1] 1059.832891    0.000001
#> attr(,"convergence")
#> [1] 0
#> attr(,"value")
#> [1] 27284.05
```

We create a new kernel with optimal hyperparameters.

``` r
SEKernelopt <- set_hyperparameters(SEKernel, opt)
```

The parameters of the posterior distributions for the 3 groups can be
computed using the same Normal-Normal conjugate update as above:

``` r
posterior <- multi_posterior_mean(data, SEKernelopt)
```

Once parameters of the posterior distributions are available,
`calculate_group_overlaps()` returns a full pairwise OVL matrix, with
one row/column per group:

``` r
calculate_group_overlaps(posterior)
#>           1         2         3
#> 1 1.0000000 0.5288475 0.3025926
#> 2 0.5288475 1.0000000 0.5266329
#> 3 0.3025926 0.5266329 1.0000000
```

As with the 2-group case, `plot_group_overlap_heatmap()` gives a
constant-size summary of every pair at once – unlike a full grid of one
density panel per pair, which grows as the number of groups squared and
quickly becomes illegible:

``` r
samples <- sample_posterior(posterior, n = 2000)
plot_group_overlap_heatmap(samples, id = unique(samples$ID)[1])
```

<img src="man/figures/README-unnamed-chunk-22-1.png" alt="" width="100%" />

## Go further

| Vignette | What you will find |
|----|----|
| `vignette("BayesOmics")` | Step-by-step model walkthrough: kernel selection, hyperparameter optimization, posterior computation, and every plot function, on both a two-group and a multi-group dose-response example. |
| `vignette("sample_size_scenarios")` | How unbalanced designs (different replicate counts per group) behave: each group’s posterior precision depends only on its own sample size, and the comparison remains valid even with strongly unequal counts. |
| `vignette("troubleshooting")` | Exact error messages you may encounter, what causes each of them, and the targeted fix. |
