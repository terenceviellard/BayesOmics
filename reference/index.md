# Package index

## All functions

- [`as_block_diag_matrix()`](https://terenceviellard.github.io/BayesOmics/reference/as_block_diag_matrix.md)
  : Assemble a Dense Block-Diagonal Covariance Matrix
- [`bhattacharyya_metric()`](https://terenceviellard.github.io/BayesOmics/reference/bhattacharyya_metric.md)
  : Bhattacharyya Distance Metric
- [`calculate_group_overlaps()`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md)
  : Compute Overlapping Coefficient between Groups
- [`compute_group_diff()`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
  : Compute a Pairwise Group-Difference Metric
- [`compute_group_diff_samples()`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff_samples.md)
  : Compute a Pairwise Draw-Based Group-Difference Metric
- [`compute_multi_diff()`](https://terenceviellard.github.io/BayesOmics/reference/compute_multi_diff.md)
  : Compute the Multivariate Distribution of Group Differences
- [`energy_metric()`](https://terenceviellard.github.io/BayesOmics/reference/energy_metric.md)
  : Energy Distance Metric
- [`evaluate_metric()`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
  [`requires_shared_kernel()`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
  [`is_symmetric_metric()`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
  : Group-Difference Metrics for BayesOmics Posteriors
- [`evaluate_sample_metric()`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_sample_metric.md)
  : Draw-Based Group-Difference Metrics
- [`fit_block_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/fit_block_posterior.md)
  : Fit Independent Posteriors per Block
- [`get_muk_block()`](https://terenceviellard.github.io/BayesOmics/reference/get_muk_block.md)
  : Concatenate Per-Block Posterior Means for a Group
- [`get_sigmak_block()`](https://terenceviellard.github.io/BayesOmics/reference/get_sigmak_block.md)
  : Reconstruct Per-Block Posterior Covariances for a Group
- [`hellinger_metric()`](https://terenceviellard.github.io/BayesOmics/reference/hellinger_metric.md)
  : Hellinger Distance Metric
- [`jeffreys_metric()`](https://terenceviellard.github.io/BayesOmics/reference/jeffreys_metric.md)
  : Symmetrized Jeffreys Divergence Metric
- [`kl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/kl_metric.md)
  : Kullback-Leibler Divergence Metric
- [`mahalanobis_metric()`](https://terenceviellard.github.io/BayesOmics/reference/mahalanobis_metric.md)
  : Mahalanobis Distance Metric
- [`marginal_metric()`](https://terenceviellard.github.io/BayesOmics/reference/marginal_metric.md)
  : Marginal-Averaged Metric
- [`mean_diff_metric()`](https://terenceviellard.github.io/BayesOmics/reference/mean_diff_metric.md)
  : Mean Absolute Difference of Means Metric
- [`mmd_metric()`](https://terenceviellard.github.io/BayesOmics/reference/mmd_metric.md)
  : Maximum Mean Discrepancy (MMD) Metric
- [`optim_hp()`](https://terenceviellard.github.io/BayesOmics/reference/optim_hp.md)
  : Optimize Hyperparameters for a kernel (with additive kernel support)
- [`ovl_metric()`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
  : Gaussian Overlapping Coefficient (OVL) Metric
- [`partition_by_id()`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_id.md)
  : Explicit Manual ID-to-Block Assignment
- [`partition_by_range()`](https://terenceviellard.github.io/BayesOmics/reference/partition_by_range.md)
  : Partition IDs by Slicing Input into Ranges
- [`partition_dbscan()`](https://terenceviellard.github.io/BayesOmics/reference/partition_dbscan.md)
  : Partition IDs by Density-Based Clustering (Reserved, Not
  Implemented)
- [`partition_dense()`](https://terenceviellard.github.io/BayesOmics/reference/partition_dense.md)
  : Dense Partition (a Single Block)
- [`partition_diagonal()`](https://terenceviellard.github.io/BayesOmics/reference/partition_diagonal.md)
  : Diagonal Partition (One Block per ID)
- [`partition_from_annotation()`](https://terenceviellard.github.io/BayesOmics/reference/partition_from_annotation.md)
  : Partition IDs from an External Annotation
- [`per_feature_metric()`](https://terenceviellard.github.io/BayesOmics/reference/per_feature_metric.md)
  : Per-Feature-Normalized Metric
- [`plot_distrib()`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md)
  : Plot the posterior distribution(s) of the difference of means
- [`plot_distrib_each_pair()`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib_each_pair.md)
  : Plot each pairwise group comparison separately
- [`plot_group_overlap_heatmap()`](https://terenceviellard.github.io/BayesOmics/reference/plot_group_overlap_heatmap.md)
  : Plot a heatmap of pairwise group overlap coefficients
- [`plot_multi_diff()`](https://terenceviellard.github.io/BayesOmics/reference/plot_multi_diff.md)
  : Plot a multivariate summary of group differences
- [`plot_posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/plot_posterior_mean.md)
  : Plot the posterior mean as a function of id
- [`plot_posterior_overlap()`](https://terenceviellard.github.io/BayesOmics/reference/plot_posterior_overlap.md)
  : Plot overlaid posterior distributions and their overlap
- [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
  : Compute Posterior Means for Different Groups
- [`print(`*`<bayesomics_posterior>`*`)`](https://terenceviellard.github.io/BayesOmics/reference/print.bayesomics_posterior.md)
  : Print a BayesOmics Posterior Object
- [`resolve_partition()`](https://terenceviellard.github.io/BayesOmics/reference/resolve_partition.md)
  : Resolve a Partition Strategy into a Block Assignment
- [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  : Sample from a Normal multivariate distribution
- [`significant_fraction_metric()`](https://terenceviellard.github.io/BayesOmics/reference/significant_fraction_metric.md)
  : Fraction of Individually Significant Features
- [`simu_db()`](https://terenceviellard.github.io/BayesOmics/reference/simu_db.md)
  : Generate a Synthetic Dataset for BayesOmics
- [`simu_db_kernel()`](https://terenceviellard.github.io/BayesOmics/reference/simu_db_kernel.md)
  : Generate a Synthetic Dataset with Kernel-Structured Covariance
- [`sliced_wasserstein_metric()`](https://terenceviellard.github.io/BayesOmics/reference/sliced_wasserstein_metric.md)
  : Sliced Wasserstein Distance Metric
- [`tvd_metric()`](https://terenceviellard.github.io/BayesOmics/reference/tvd_metric.md)
  : Total Variation Distance
- [`wasserstein_metric()`](https://terenceviellard.github.io/BayesOmics/reference/wasserstein_metric.md)
  : Wasserstein-2 (Frechet) Distance Metric
