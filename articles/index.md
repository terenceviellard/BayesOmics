# Articles

### Start here

The model in one page, before any code.

- [Concepts: the model in one
  page](https://terenceviellard.github.io/BayesOmics/articles/concepts.md):

  What BayesOmics models, how the kernel enters, and what comparing two
  posteriors means.

### Beginner

Your first comparison, from simulated data to your own dataset.

- [01 · Basic pipeline: comparing two groups with a spatial
  kernel](https://terenceviellard.github.io/BayesOmics/articles/01_basic_pipeline.md):

  The full two-group pipeline on simulated data: kernel, hyperparameter
  fit, posteriors, group comparison.

- [02 · Univariate mode: comparing two groups with no feature
  axis](https://terenceviellard.github.io/BayesOmics/articles/02_univariate.md):

  Compare two groups when there is no feature axis: one scalar per
  replicate, fitted in closed form.

- [03 · Working from real data instead of a
  simulation](https://terenceviellard.github.io/BayesOmics/articles/03_real_data.md):

  Turn your own dataset into the BayesOmics long format, with a
  ready-to-use AI prompt.

### Intermediate

Modelling choices and reading results with several groups.

- [04 · Pooled vs. non-pooled: one kernel for every group, or one per
  group](https://terenceviellard.github.io/BayesOmics/articles/04_pooled_vs_nonpooled.md):

  One shared kernel for every group, or one independent kernel per
  group: what changes.

- [05 · Comparing more than two
  groups](https://terenceviellard.github.io/BayesOmics/articles/05_multi_group.md):

  Read a pairwise matrix and the multi-group overview plot for three or
  more groups.

- [06 · Choosing a metric: OVL vs.
  Wasserstein](https://terenceviellard.github.io/BayesOmics/articles/06_metric_choice.md):

  Why the joint OVL saturates as features grow, and when to use
  per-feature Wasserstein.

- [07 · Comparing groups with different sample
  sizes](https://terenceviellard.github.io/BayesOmics/articles/07_unequal_sample_size.md):

  Compare groups with different replicate counts and see how it affects
  posterior width.

### Advanced

Other kernels, multi-dimensional inputs and scaling up.

- [08 · Using a kernel other than Squared
  Exponential](https://terenceviellard.github.io/BayesOmics/articles/08_kernel_choice.md):

  Beyond the Squared Exponential kernel: choosing and composing other
  kernels.

- [09 · Input with more than one
  dimension](https://terenceviellard.github.io/BayesOmics/articles/09_multi_dim_input.md):

  Features located by several covariates at once, and why the kernel
  factorizes on a grid.

- [10 · Scaling up with a block-diagonal
  partition](https://terenceviellard.github.io/BayesOmics/articles/10_block_diagonal.md):

  Split many features into independent blocks to fit and compare at
  scale.

- [11 · Performance at a
  glance](https://terenceviellard.github.io/BayesOmics/articles/11_performance.md):

  Timings versus feature count, and when to switch to block-diagonal
  fitting.

### Gallery

The key figures at a glance, each linking to its article.

- [Gallery](https://terenceviellard.github.io/BayesOmics/articles/gallery.md):

  The key figures of the package at a glance, each linking to the
  article that explains it.

### Troubleshooting

- [Troubleshooting: common errors and how to fix
  them](https://terenceviellard.github.io/BayesOmics/articles/troubleshooting.md):
