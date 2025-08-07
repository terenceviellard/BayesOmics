
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesOmics

<!-- badges: start -->

<!-- badges: end -->

Bayesian differential analysis for omics data.The correlation structure
is captured by a kernel. Users need to choose a kernel, optimize the
hyperparameters, and then use the correlation matrix computed by this
kernel in the analysis. \## Installation

You can install the development version of BayesOmics from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("terenceviellard/BayesOmics")
```

## Bayes Omics in a nutshell

This is a basic example which shows you how to solve a common problem:

``` r
library(BayesOmics)
```

Generate a synthetic dataset with the correct format:

``` r
set.seed(123)
data <- simu_db(nb_id = 400,
                nb_group = 2,
                nb_sample = 1)
head(data)
#>     ID Group Sample     Input    Output
#> 1 ID_1     1      1 14.378876 28.682116
#> 2 ID_1     2      1 39.415257  9.008909
#> 3 ID_2     1      1 20.448846 16.073155
#> 4 ID_2     2      1 44.150870 24.392639
#> 5 ID_3     1      1 47.023364 42.020415
#> 6 ID_3     2      1  2.277825  5.049311
```

``` r
colnames(data)
#> [1] "ID"     "Group"  "Sample" "Input"  "Output"
```

You need to have 5 columns in long format like this

We choose a Squared Exponential (SE) kernel, given by:

$$
K_{\text{SE}}(x, x') = \sigma^2 \exp\left(-\frac{\|x - x'\|^2}{2\ell^2}\right)
$$

where: - $K_{\text{SE}}(x, x')$ is the covariance between inputs $x$ and
$x'$. - $\sigma^2$ is the signal variance. - $\ell$ is the length
scale. - $\|x - x'\|$ is the Euclidean distance between $x$ and $x'$.
Then we need to optimize $\sigma^2$ and $\ell$ :

``` r
SEKernel <- new("SEKernel")

hp <- c(1.0, 1.0)
SEKernel <- set_hyperparameters(SEKernel, hp)

n <- dim(data)[1]

mean_nul <- rep(0, n)
post_cov_nul <- diag(0.1, n)

opt <- optim_hp(
  hp = hp,
  db = data,
  mean = mean_nul,
  kern = SEKernel,
  post_cov = post_cov_nul
)
opt
#> [1] 1.052419e+03 1.134949e-03
```

We create a new kernel with optimal HpS

``` r
SEKernelopt <- set_hyperparameters(SEKernel, opt)
```

The parameters of all posterior distributions can be computed thanks to:

``` r
posterior <- multi_posterior_mean(data, SEKernelopt)
```

Once parameters of the posterior distributions are available, we can
calculate the overlapping coefficient for our differential analysis.

``` r
calculate_group_overlaps(posterior)
#>           Group1    Group2
#> Group1 1.0000000 0.9066489
#> Group2 0.9066489 1.0000000
```

## Bayes Omics for multigroup analysis

We can use Bayes Omics to compare more than 2 groups. Generate a
synthetic dataset with 5 groups with the correct format:

``` r
data <- simu_db(nb_id = 200,
                nb_group = 5,
                nb_sample = 1)
```

We also choose a SEKernel so we need to optimize Hps $\sigma^2$ and
$\ell$ :

``` r
SEKernel <- new("SEKernel")

hp <- c(1.0, 1.0)
SEKernel <- set_hyperparameters(SEKernel, hp)

n <- dim(data)[1]

mean_nul <- rep(0, n)
post_cov_nul <- diag(0.1, n)

opt <- optim_hp(
  hp = hp,
  db = data,
  mean = mean_nul,
  kern = SEKernel,
  post_cov = post_cov_nul
)
opt
#> [1]  1.457344e+03 -3.114536e-05
```

We create a new kernel with optimal HpS

``` r
SEKernelopt <- set_hyperparameters(SEKernel, opt)
```

The parameters of all posterior distributions can be computed thanks to:

``` r
posterior <- multi_posterior_mean(data, SEKernelopt)
```

Once parameters of the posterior distributions are available.

``` r
library(bayestestR)
#> Warning: le package 'bayestestR' a été compilé avec la version R 4.5.1
calculate_group_overlaps(posterior)
#>           Group1    Group2    Group3    Group4    Group5
#> Group1 1.0000000 0.8343609 0.8964618 0.7706932 0.8590161
#> Group2 0.8343609 1.0000000 0.8137074 0.9100066 0.7991936
#> Group3 0.8964618 0.8137074 1.0000000 0.7767861 0.9316448
#> Group4 0.7706932 0.9100066 0.7767861 1.0000000 0.7954788
#> Group5 0.8590161 0.7991936 0.9316448 0.7954788 1.0000000
```
