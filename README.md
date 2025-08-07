
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

## Example

This is a basic example which shows you how to solve a common problem:

### DATA

``` r
library(BayesOmics)
```

``` r
data <- simu_db()
head(data)
#>     ID Group Sample     Input    Output
#> 1 ID_1     1      1 21.789283 45.173005
#> 2 ID_1     1      2  6.568757 36.764282
#> 3 ID_1     1      3 49.113050 13.280864
#> 4 ID_1     1      4 40.326233 45.576368
#> 5 ID_1     1      5 18.248666  6.405751
#> 6 ID_1     2      1 40.349434 46.062372
```

``` r
colnames(data)
#> [1] "ID"     "Group"  "Sample" "Input"  "Output"
```

You need to have 6 columns in long format like this

### Kernel

First, choose a kernel, in our example, kernel SE defined as:

$$
k_{\text{SE}}(x, x') = \sigma^2 \exp\left(-\frac{\ x - x'\|^2}{2\ell^2}\right)
$$

where: - $\sigma^2$ is the signal variance, - $\ell$ is the length
scale, - $x$ and $x'$ are input vectors.

Then we need to optimize Hps $\sigma^2$ and $\ell$ :

``` r
SEKernel <- new("SEKernel")

hp <- c(1.0, 1.0)
SEKernel <- set_hyperparameters(SEKernel, hp)

n <- dim(data)[1]

mean_nul <- rep(0, n)
post_cov_nul <- diag(0.1, n)
pen_diag <- 1e-6



opt <- optim_hp(
  hp = hp,
  db = df,
  mean = mean_nul,
  kern = SEKernel,
  post_cov = post_cov_nul,
  pen_diag = pen_diag
)
opt$par
```

We create a new kernel with optimal HpS

``` r
# En attendant de reparer optimHP :
SEKernel <- new("SEKernel")
opt <- list()
opt$par <- c(1, 1)
SEKernelopt <- set_hyperparameters(SEKernel, opt$par)
```

### Posterior

``` r
posterior <- multi_posterior_mean(data, SEKernelopt)
```

### DIfferential analtsis

``` r
calculate_group_overlaps(posterior)
#>           Group1    Group2
#> Group1 1.0000000 0.3279821
#> Group2 0.3279821 1.0000000
```
