#' @importFrom ggplot2 ggplot aes geom_ribbon geom_area geom_vline geom_label geom_point
#' @importFrom ggplot2 geom_line geom_tile geom_text facet_wrap
#' @importFrom ggplot2 scale_fill_manual scale_colour_manual scale_fill_gradient
#' @importFrom ggplot2 theme_classic theme element_line element_rect element_text ylab xlab labs
#' @importFrom tibble tibble
#' @importFrom gridExtra grid.arrange
#' @importFrom stats density quantile
#' @importFrom rlang .data

## Shared colour palette (ColorBrewer "Set2" -- a colourblind-safe
## qualitative set; https://colorbrewer2.org) so every plot_*() function
## draws from the same family instead of each picking its own (the previous
## mix of a custom pastel pink/blue, discrete viridis, and continuous
## viridis read as several unrelated plots). `.bo_centre`/`.bo_tail` are the
## two-colour pair used for every "two things being compared" plot (the
## within-/outside-credible-interval shading, and the two groups' density
## curves in plot_posterior_overlap()); `.bo_palette` is the full 8-colour
## set used when more than two categories need distinct colours (e.g.
## plot_posterior_mean() with more than 2 groups).
.bo_palette <- c(
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
)
.bo_centre <- .bo_palette[3]   # blue-purple: within-CI / first group of a pair
.bo_tail   <- .bo_palette[2]   # salmon: outside-CI / second group of a pair
.bo_line   <- "grey35"         # neutral reference line / axis colour

## Discrete colour scale for an arbitrary number of categories, recycling
## .bo_palette and extending it (via interpolation) past its 8 native
## colours rather than erroring -- plot_posterior_mean() has no fixed
## number of groups.
#' @noRd
bo_discrete_colours <- function(n) {
  if (n <= length(.bo_palette)) {
    .bo_palette[seq_len(n)]
  } else {
    grDevices::colorRampPalette(.bo_palette)(n)
  }
}

#' @noRd
theme_bayesomics <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "plain", size = 11, colour = "grey25"),
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = .bo_line),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(face = "plain", colour = "grey25")
    )
}

## Shared computation behind every density-of-the-difference panel: extracts
## the (group1) or (group1 - group2) sample vector for `id`, its density, and
## the credible-interval shading data frame. Factored out so the single-pair
## plot (build_single_distrib) and the multi-pair facets
## (build_faceted_pairs_distrib) compute it identically.
#' @noRd
compute_pair_distrib <- function(sample_distrib, group1, group2, id, prob_CI) {
  db <- sample_distrib %>%
    dplyr::filter(.data$ID %in% id) %>%
    dplyr::filter(.data$Group == group1) %>%
    dplyr::pull(.data$Sample)

  bar <- mean(db)

  if (!is.null(group2)) {
    db2 <- sample_distrib %>%
      dplyr::filter(.data$ID %in% id) %>%
      dplyr::filter(.data$Group == group2) %>%
      dplyr::pull(.data$Sample)

    db <- db - db2
    bar <- 0
  }

  dens <- stats::density(db, n = 5000)
  CI <- stats::quantile(db, prob = c((1 - prob_CI) / 2, (1 + prob_CI) / 2))

  db_plot <- tibble::tibble(x = dens$x, y = dens$y) %>%
    dplyr::mutate(quant = factor(findInterval(.data$x, CI)))

  list(db = db, bar = bar, db_plot = db_plot)
}

## Empirical overlapping coefficient between two 1-D samples: the area of
## intersection of their density estimates, evaluated on a shared grid. Used
## to rank group pairs by how differentiated they are for a given id.
#' @noRd
estimate_overlap_1d <- function(x, y, n = 512) {
  rng <- range(c(x, y))
  dx <- stats::density(x, from = rng[1], to = rng[2], n = n)
  dy <- stats::density(y, from = rng[1], to = rng[2], n = n)
  overlap_y <- pmin(dx$y, dy$y)
  area <- sum((overlap_y[-1] + overlap_y[-n]) / 2 * diff(dx$x))
  min(max(area, 0), 1)
}

#' @noRd
compute_pairwise_overlap <- function(sample_distrib, id, groups) {
  k <- length(groups)
  mat <- diag(k)
  rownames(mat) <- colnames(mat) <- groups
  if (k < 2) return(mat)

  samples_by_group <- stats::setNames(
    lapply(groups, function(g) {
      sample_distrib %>%
        dplyr::filter(.data$ID %in% id) %>%
        dplyr::filter(.data$Group == g) %>%
        dplyr::pull(.data$Sample)
    }),
    groups
  )

  for (i in seq_len(k - 1)) {
    for (j in seq(i + 1, k)) {
      ov <- estimate_overlap_1d(samples_by_group[[i]], samples_by_group[[j]])
      mat[i, j] <- ov
      mat[j, i] <- ov
    }
  }
  mat
}

#' @noRd
build_single_distrib <- function(
    sample_distrib,
    group1,
    group2,
    id,
    prob_CI,
    show_prob,
    mean_bar,
    index_group1,
    index_group2
){

  ## Retrieve the name of the first group in 'sample_distrib' if needed
  if(is.null(group1)){
    group1 = sample_distrib$Group[1]
  }

  info <- compute_pair_distrib(sample_distrib, group1, group2, id, prob_CI)
  db <- info$db
  bar <- info$bar
  db_plot <- info$db_plot

  ## Define the name of index for the label of group1
  if(is.null(index_group1)){
    index_group1 = group1
  }
  ## Define the name of index for the label of group2
  if(is.null(index_group2)){
    index_group2 = group2
  }

  gg = ggplot2::ggplot(db_plot) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x = .data$x,
                   y = .data$y,
                   ymin=0,
                   ymax=.data$y,
                   fill = .data$quant)
    ) +
    ggplot2::ylab('Density') +
    ggplot2::scale_fill_manual(values=c(.bo_tail, .bo_centre, .bo_tail)) +
    theme_bayesomics() +
    ggplot2::theme(legend.position="none")

  if(mean_bar){
    gg = gg + ggplot2::geom_vline(xintercept = bar, color = .bo_line)
  }

  if(is.null(group2)){
    ## Add the adequate label for the x-axis
    gg = gg + ggplot2::xlab( bquote(mu[.(index_group1)]) )
  } else {
    ## Add the adequate label for the x-axis
    gg = gg + ggplot2::xlab(bquote(mu[.(index_group1)] - mu[.(index_group2)]))
    ## Add probabilities of the group comparison if required
    if( show_prob == TRUE ){
      p_inf = (sum(db<0)/length(db)) %>% round(2) %>% as.character()
      p_sup = (sum(db>0)/length(db)) %>% round(2) %>% as.character()
      exp_l = bquote(P(mu[.(index_group1)] <= mu[.(index_group2)]) == .(p_inf))
      exp_r = bquote(P(mu[.(index_group1)] >= mu[.(index_group2)]) == .(p_sup))

      ## Anchored at the panel corners (-Inf/Inf) rather than at `bar` --
      ## anchoring at `bar` clips the right-hand label whenever the two
      ## groups are well separated and `bar` (0) sits near the panel edge.
      gg = gg +
        ggplot2::geom_label(
          data = tibble::tibble(x = -Inf),
          ggplot2::aes(
            x = .data$x,
            y = Inf,
            label = deparse(exp_l)
          ),
          parse = TRUE,
          size = 4,
          hjust=0,
          vjust=1) +
        ggplot2::geom_label(
          data = tibble::tibble(x = Inf),
          ggplot2::aes(
            x = .data$x,
            y = Inf,
            label = deparse(exp_r)
          ),
          parse = TRUE,
          size = 4,
          hjust=1,
          vjust=1)
    }
  }
  return(gg)
}

## One ggplot, faceted by group pair, for a curated *subset* of pairs (e.g.
## the top_n_pairs most differentiated ones). Unlike pasting together N
## independent build_single_distrib() panels with gridExtra::grid.arrange(),
## this shares a single legend/scale across all facets, which is both cleaner
## and avoids the panel count growing as O(k^2) being mistaken for "more
## information" -- it stays a single, constant-size ggplot object.
#' @noRd
build_faceted_pairs_distrib <- function(sample_distrib, pairs, id, prob_CI, show_prob, mean_bar) {
  panels <- lapply(pairs, function(p) {
    info <- compute_pair_distrib(sample_distrib, p[1], p[2], id, prob_CI)
    db_plot <- info$db_plot
    db_plot$pair <- paste0(p[1], " vs ", p[2])
    db <- info$db
    list(
      db_plot = db_plot,
      pair    = db_plot$pair[1],
      p_inf   = (sum(db < 0) / length(db)) %>% round(2) %>% as.character(),
      p_sup   = (sum(db > 0) / length(db)) %>% round(2) %>% as.character(),
      group1  = p[1],
      group2  = p[2]
    )
  })

  combined <- do.call(rbind, lapply(panels, `[[`, "db_plot"))

  gg <- ggplot2::ggplot(combined, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = .data$y, fill = .data$quant)) +
    ggplot2::facet_wrap(~pair, scales = "free") +
    ggplot2::scale_fill_manual(values = c(.bo_tail, .bo_centre, .bo_tail)) +
    ggplot2::ylab("Density") +
    ggplot2::xlab("Difference of posterior means (group1 - group2)") +
    theme_bayesomics() +
    ggplot2::theme(legend.position = "none")

  if (mean_bar) {
    gg <- gg + ggplot2::geom_vline(xintercept = 0, color = .bo_line)
  }

  if (show_prob) {
    label_df <- do.call(rbind, lapply(panels, function(pp) {
      exp_l <- bquote(P(mu[.(pp$group1)] <= mu[.(pp$group2)]) == .(pp$p_inf))
      exp_r <- bquote(P(mu[.(pp$group1)] >= mu[.(pp$group2)]) == .(pp$p_sup))
      tibble::tibble(pair = pp$pair, label_l = deparse(exp_l), label_r = deparse(exp_r))
    }))

    gg <- gg +
      ggplot2::geom_label(
        data = label_df, ggplot2::aes(x = -Inf, y = Inf, label = .data$label_l),
        parse = TRUE, size = 3, hjust = 0, vjust = 1, inherit.aes = FALSE
      ) +
      ggplot2::geom_label(
        data = label_df, ggplot2::aes(x = Inf, y = Inf, label = .data$label_r),
        parse = TRUE, size = 3, hjust = 1, vjust = 1, inherit.aes = FALSE
      )
  }

  gg
}

#' @title Plot the posterior mean as a function of id
#'
#' @description
#' Display, for every id present in \code{sample_distrib}, the posterior mean
#' of its distribution (averaged over the drawn samples), coloured by group.
#' This gives a region-wide view of the posterior profile (e.g. one point per
#' CpG site in a methylation analysis), complementing the per-id distribution
#' plots produced by \code{plot_distrib()}.
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}. This argument should
#'    contain the empirical posterior distributions to be summarized.
#'
#' @return A \code{ggplot} object, with one point per (id, group) pair, the
#'    id on the y-axis, the posterior mean on the x-axis, and colour
#'    indicating the group.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' plot_posterior_mean(samples)
plot_posterior_mean <- function(sample_distrib){

  required_cols <- c("ID", "Group", "Sample")
  if(!all(required_cols %in% names(sample_distrib))){
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  mean_db = sample_distrib %>%
    dplyr::group_by(.data$ID, .data$Group) %>%
    dplyr::summarise(Mean = mean(.data$Sample), .groups = "drop")

  n_groups <- length(unique(mean_db$Group))

  ggplot2::ggplot(mean_db) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$Mean, y = .data$ID, col = factor(.data$Group))
    ) +
    ggplot2::xlab('Posterior Mean') +
    theme_bayesomics() +
    ggplot2::scale_colour_manual(values = bo_discrete_colours(n_groups), name = "Group")
}

#' @title Plot a heatmap of pairwise group overlap coefficients
#'
#' @description
#' For a given id, compute and display an empirical overlapping coefficient
#' (the area of intersection of the two groups' posterior density estimates)
#' between every pair of groups present in \code{sample_distrib}, as a
#' group x group heatmap. This gives a single, constant-size overview of how
#' differentiated every pair of groups is for that id, regardless of how many
#' groups are present -- unlike a full grid of one density panel per pair,
#' which grows as the number of groups squared and becomes illegible.
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}.
#' @param id A character string, the id for which pairwise overlaps are
#'    computed. If NULL (default), only the first id appearing in
#'    \code{sample_distrib} is used.
#' @param digits Number of decimal digits used when displaying the overlap
#'    coefficient on each tile. Defaults to \code{2}.
#'
#' @return A \code{ggplot} object: a tiled heatmap of groups x groups, filled
#'    by the estimated overlapping coefficient (1 = identical distributions,
#'    0 = fully separated).
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 4, nb_sample = 5, diff_group = 4)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' plot_group_overlap_heatmap(samples, id = unique(samples$ID)[1])
plot_group_overlap_heatmap <- function(sample_distrib, id = NULL, digits = 2){

  required_cols <- c("ID", "Group", "Sample")
  if(!all(required_cols %in% names(sample_distrib))){
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  if(is.null(id)){
    id = sample_distrib$ID %>% unique() %>% utils::head(1)
  }
  if(length(id) > 1){
    warning("Multiple IDs provided to plot_group_overlap_heatmap; only '", id[1], "' will be displayed.")
    id <- id[1]
  }

  groups <- sample_distrib$Group %>% unique() %>% sort()
  mat <- compute_pairwise_overlap(sample_distrib, id, groups)

  long <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(long) <- c("Group1", "Group2", "Overlap")

  ggplot2::ggplot(long, ggplot2::aes(x = .data$Group1, y = .data$Group2, fill = .data$Overlap)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf(paste0("%.", digits, "f"), .data$Overlap)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F2F3F8", high = "#5E72A4", limits = c(0, 1), name = "Overlap") +
    ggplot2::labs(title = paste0("Pairwise group overlap (id: ", id, ")"), x = NULL, y = NULL) +
    theme_bayesomics() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' @title Plot each pairwise group comparison separately
#'
#' @description
#' Return one difference-of-means posterior plot per pair of groups present
#' in \code{sample_distrib}, as a named list of \code{ggplot} objects -- with
#' no summary/selection logic, regardless of how many groups (and therefore
#' pairs) there are. Use this when you want to inspect or save every pairwise
#' comparison individually (e.g. in a report with one figure per pair),
#' as opposed to \code{plot_distrib()}'s automatic top-N summary view for
#' more than two groups.
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}.
#' @param id A character string, the id to plot. If NULL (default), only the
#'    first id appearing in \code{sample_distrib} is used.
#' @param prob_CI A number, between 0 and 1, the level of the Credible
#'    Interval. See \code{\link{plot_distrib}}.
#' @param show_prob A boolean, whether to display the probability labels.
#'    See \code{\link{plot_distrib}}.
#' @param mean_bar A boolean, whether to display the vertical bar at 0.
#'    See \code{\link{plot_distrib}}.
#'
#' @return A named list of \code{ggplot} objects, one per group pair, named
#'    \code{"<group1>_vs_<group2>"}.
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 4, nb_sample = 5, diff_group = 4)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' plots <- plot_distrib_each_pair(samples, id = unique(samples$ID)[1])
#' names(plots)
plot_distrib_each_pair <- function(sample_distrib, id = NULL, prob_CI = 0.95,
                                    show_prob = TRUE, mean_bar = TRUE){

  required_cols <- c("ID", "Group", "Sample")
  if(!all(required_cols %in% names(sample_distrib))){
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  if(is.null(id)){
    id = sample_distrib$ID %>% unique() %>% utils::head(1)
  }
  if(length(id) > 1){
    warning("Multiple IDs provided to plot_distrib_each_pair; only '", id[1], "' will be displayed.")
    id <- id[1]
  }

  groups <- sample_distrib$Group %>% unique() %>% sort()
  if(length(groups) < 2){
    stop("plot_distrib_each_pair() requires at least two groups in 'sample_distrib'.")
  }

  # TODO: materializes one full ggplot object per pair (O(choose(G, 2))) with
  # no cap and no warning, unlike calculate_group_overlaps()'s
  # max_groups_warn/max_dim_warn. Fine for the documented use case (2-5
  # groups) but a direct call with many groups (e.g. 30+ dose levels) builds
  # a large list of plots simultaneously in memory. Consider documenting in
  # @details that plot_group_overlap_heatmap() is the safer alternative for
  # many groups, and/or adding a similar guard-rail warning here.
  pairs <- utils::combn(groups, 2, simplify = FALSE)
  plots <- lapply(pairs, function(p){
    build_single_distrib(sample_distrib, p[1], p[2], id, prob_CI, show_prob, mean_bar, NULL, NULL)
  })
  names(plots) <- vapply(pairs, function(p) paste0(p[1], "_vs_", p[2]), character(1))
  plots
}

## Shared computation behind build_overlap_distrib(): both groups' posterior
## densities on a common grid, plus the pointwise minimum (the overlap region
## that calculate_group_overlaps()'s OVL coefficient is the area of).
#' @noRd
compute_overlap_distrib <- function(sample_distrib, group1, group2, id, n = 512) {
  x1 <- sample_distrib %>%
    dplyr::filter(.data$ID %in% id, .data$Group == group1) %>%
    dplyr::pull(.data$Sample)
  x2 <- sample_distrib %>%
    dplyr::filter(.data$ID %in% id, .data$Group == group2) %>%
    dplyr::pull(.data$Sample)

  rng <- range(c(x1, x2))
  d1 <- stats::density(x1, from = rng[1], to = rng[2], n = n)
  d2 <- stats::density(x2, from = rng[1], to = rng[2], n = n)
  overlap_y <- pmin(d1$y, d2$y)
  area <- sum((overlap_y[-1] + overlap_y[-n]) / 2 * diff(d1$x))
  ov <- min(max(area, 0), 1)

  dens_df <- tibble::tibble(
    x = c(d1$x, d2$x),
    y = c(d1$y, d2$y),
    Group = rep(c(group1, group2), each = n)
  )
  overlap_df <- tibble::tibble(x = d1$x, y = overlap_y)

  list(dens_df = dens_df, overlap_df = overlap_df, ov = ov)
}

## One ggplot: both groups' posterior densities overlaid for `id`, with the
## overlap region (pointwise minimum of the two densities) shaded -- the same
## quantity calculate_group_overlaps()'s OVL coefficient measures the area of.
#' @noRd
build_overlap_distrib <- function(sample_distrib, group1, group2, id) {
  info <- compute_overlap_distrib(sample_distrib, group1, group2, id)
  pair_colours <- c(.bo_centre, .bo_tail)

  ggplot2::ggplot() +
    ## Each group's *full* density, lightly filled -- not just the shaded
    ## overlap intersection -- so both whole distributions are readable at
    ## a glance, not only where they intersect.
    ggplot2::geom_area(
      data = info$dens_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$Group),
      alpha = 0.18, position = "identity"
    ) +
    ggplot2::geom_area(
      data = info$overlap_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "grey45", alpha = 0.55
    ) +
    ggplot2::geom_line(
      data = info$dens_df,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$Group),
      linewidth = 0.9
    ) +
    ggplot2::geom_label(
      data = tibble::tibble(x = Inf, y = Inf, label = sprintf("OVL = %.2f", info$ov)),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 1, vjust = 1, size = 4, linewidth = 0
    ) +
    ggplot2::labs(title = paste0(group1, " vs ", group2), x = "Posterior value", y = "Density") +
    ggplot2::scale_fill_manual(values = pair_colours, guide = "none") +
    ggplot2::scale_colour_manual(values = pair_colours, name = NULL) +
    theme_bayesomics()
}

#' @title Plot overlaid posterior distributions and their overlap
#'
#' @description
#' For a given id, plot both groups' posterior density curves on the *same*
#' panel (unlike \code{\link{plot_distrib}}, which plots the distribution of
#' their *difference*), with the region where the two densities overlap shaded
#' and labelled with the OVL coefficient -- i.e. a direct visualization of the
#' quantity \code{\link{calculate_group_overlaps}} computes.
#'
#' Follows the same pairwise ("two by two") logic as
#' \code{\link{plot_distrib_each_pair}}: for exactly two groups (given
#' explicitly via \code{group1}/\code{group2}, or the only two present in
#' \code{sample_distrib}), a single \code{ggplot} is returned; for more than
#' two groups (with neither \code{group1} nor \code{group2} given), every
#' pairwise comparison is returned as a named list of \code{ggplot} objects,
#' one per pair.
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}.
#' @param group1 A character string, the first group to compare. If NULL
#'    (default) and \code{group2} is also NULL, the groups are inferred from
#'    \code{sample_distrib} (see Description).
#' @param group2 A character string, the second group to compare. If NULL
#'    (default), see Description.
#' @param id A character string, the id to plot. If NULL (default), only the
#'    first id appearing in \code{sample_distrib} is used.
#'
#' @return Either a single \code{ggplot} (two groups), or a named list of
#'    \code{ggplot} objects, one per group pair, named
#'    \code{"<group1>_vs_<group2>"} (more than two groups).
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' plot_posterior_overlap(samples, group1 = "1", group2 = "2", id = unique(samples$ID)[1])
plot_posterior_overlap <- function(sample_distrib, group1 = NULL, group2 = NULL, id = NULL) {

  required_cols <- c("ID", "Group", "Sample")
  if (!all(required_cols %in% names(sample_distrib))) {
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  if (is.null(id)) {
    id <- sample_distrib$ID %>% unique() %>% utils::head(1)
  }
  if (length(id) > 1) {
    warning("Multiple IDs provided to plot_posterior_overlap; only '", id[1], "' will be displayed.")
    id <- id[1]
  }

  if (is.null(group1) && is.null(group2)) {
    groups <- sample_distrib$Group %>% unique() %>% sort()
    if (length(groups) < 2) {
      stop("plot_posterior_overlap() requires at least two groups in 'sample_distrib'.")
    }
    if (length(groups) > 2) {
      pairs <- utils::combn(groups, 2, simplify = FALSE)
      plots <- lapply(pairs, function(p) build_overlap_distrib(sample_distrib, p[1], p[2], id))
      names(plots) <- vapply(pairs, function(p) paste0(p[1], "_vs_", p[2]), character(1))
      return(plots)
    }
    group1 <- groups[1]
    group2 <- groups[2]
  }

  build_overlap_distrib(sample_distrib, group1, group2, id)
}

## Constant-size multi-group overview: a pairwise-overlap heatmap (covering
## every pair, however many groups there are) plus a faceted detail grid
## limited to the top_n_pairs most differentiated pairs, plus (optionally)
## the id-mean panel. Unlike the former approach of pasting one
## build_single_distrib() panel per pair (O(k^2) panels), this always
## produces 2 or 3 panels regardless of k.
#' @noRd
plot_distrib_multi_group <- function(
    sample_distrib,
    groups,
    id,
    prob_CI,
    show_prob,
    mean_bar,
    plot_mean,
    top_n_pairs
){
  heatmap_gg <- plot_group_overlap_heatmap(sample_distrib, id)

  overlap_mat <- compute_pairwise_overlap(sample_distrib, id, groups)
  all_pairs <- utils::combn(groups, 2, simplify = FALSE)
  overlap_vals <- vapply(all_pairs, function(p) overlap_mat[p[1], p[2]], numeric(1))
  ## Ascending overlap = most differentiated pairs first.
  ord <- order(overlap_vals)
  n_show <- min(top_n_pairs, length(all_pairs))
  selected_pairs <- all_pairs[ord[seq_len(n_show)]]

  grid_gg <- build_faceted_pairs_distrib(
    sample_distrib, selected_pairs, id, prob_CI, show_prob, mean_bar
  )

  gg_list <- list(heatmap_gg, grid_gg)
  if(plot_mean){
    gg_list[[length(gg_list) + 1]] <- plot_posterior_mean(sample_distrib)
  }

  gridExtra::grid.arrange(grobs = gg_list, ncol = 1)
}

#' @title Plot the posterior distribution(s) of the difference of means
#'
#' @description
#' Display the posterior distribution of the difference of means between
#' groups for a specific id. Behaviour depends on how many groups are
#' involved:
#' \itemize{
#'   \item a single group (\code{group2} is \code{NULL} and only one group is
#'     present): the posterior distribution of the mean for that group;
#'   \item exactly two groups (\code{group1}/\code{group2} given explicitly,
#'     or exactly two groups present in \code{sample_distrib}): the posterior
#'     distribution of the difference of means, with reference at 0 on the
#'     x-axis and probabilities of \code{group1} > \code{group2} (and
#'     conversely);
#'   \item more than two groups (and neither \code{group1} nor \code{group2}
#'     given): a constant-size summary -- a heatmap of the empirical pairwise
#'     overlap between every group (see \code{\link{plot_group_overlap_heatmap}}),
#'     plus a faceted detail grid limited to the \code{top_n_pairs} most
#'     differentiated pairs, plus (if \code{plot_mean = TRUE}) the id-mean
#'     panel. To inspect *every* pairwise comparison individually instead of
#'     this automatic summary, use \code{\link{plot_distrib_each_pair}}.
#' }
#'
#' @param sample_distrib A data frame, typically coming from the
#'    \code{sample_posterior()} function, containing the following columns:
#'    \code{ID}, \code{Group} and \code{Sample}. This argument should
#'    contain the empirical posterior distributions to be displayed.
#' @param group1 A character string, corresponding to the name of the group
#'    for which we plot the posterior distribution of the mean. If NULL
#'    (default) and \code{group2} is also NULL, the groups are inferred from
#'    \code{sample_distrib} (see Description).
#' @param group2 A character string, corresponding to the name of the group
#'    we want to compare to \code{group1}. If NULL (default), see
#'    Description.
#' @param id A character string, corresponding to the name of the id
#'    for which we plot the posterior distribution of the mean. If NULL
#'    (default), only the first appearing in \code{sample_distrib} is displayed.
#' @param prob_CI A number, between 0 and 1, corresponding the level of the
#'    Credible Interval (CI), represented as side regions (in red) of the
#'    posterior distribution. The default value (0.95) display the 95% CI,
#'    meaning that the central region (in blue) contains 95% of the probability
#'    distribution of the mean.
#' @param show_prob A boolean, indicating whether we display the label of
#'    probability comparisons between two groups (ignored for single-group
#'    plots).
#' @param mean_bar A boolean, indicating whether we display the vertical bar
#'    corresponding to 0 on the x-axis (when comparing two groups), of the mean
#'    value of the distribution (when displaying a unique group).
#' @param index_group1 A character string, used as the index of \code{group1} in
#'    the legends. If NULL (default), \code{group1} is used. Only used for the
#'    single-group and two-group plots.
#' @param index_group2 A character string, used as the index of \code{group2} in
#'    the legends. If NULL (default), \code{group2} is used. Only used for the
#'    two-group plot.
#' @param plot_mean A boolean, indicating whether the panel showing the
#'    posterior mean of every id (coloured by group) should be added when
#'    more than two groups are involved.
#' @param top_n_pairs An integer, the number of most-differentiated group
#'    pairs to detail in the facet grid when more than two groups are
#'    involved (ignored otherwise). Defaults to 3.
#'
#' @return Either a single \code{ggplot} (one or two groups), or the result of
#'    \code{gridExtra::grid.arrange} (more than two groups: the overlap
#'    heatmap, the top-\code{top_n_pairs} facet grid, and optionally the
#'    id-mean panel).
#' @export
#'
#' @examples
#' data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5)
#' kern <- methods::new("SEKernel")
#' kern <- keRnel::set_hyperparameters(kern, c(variance_se = 1, length_scale_se = 1))
#' posterior <- multi_posterior_mean(data, kern)
#' samples <- sample_posterior(posterior, n = 500)
#' plot_distrib(samples, group1 = "1", group2 = "2", id = unique(samples$ID)[1])
plot_distrib = function(
    sample_distrib,
    group1 = NULL,
    group2 = NULL,
    id = NULL,
    prob_CI = 0.95,
    show_prob = TRUE,
    mean_bar = TRUE,
    index_group1 = NULL,
    index_group2 = NULL,
    plot_mean = TRUE,
    top_n_pairs = 3
){

  required_cols <- c("ID", "Group", "Sample")
  if(!all(required_cols %in% names(sample_distrib))){
    stop(paste0("The following columns are missing: ",
                paste(setdiff(required_cols, names(sample_distrib)), collapse = ", ")))
  }

  ## Retrieve the name of ids in the 'sample_distrib' argument if needed
  if(is.null(id)){
    id = sample_distrib$ID %>% unique()
  }

  ## If we have multiple values in 'id', warn the user and keep only first
  if(length(id) > 1){
    warning("Multiple IDs provided to plot_distrib; only '", id[1], "' will be displayed.")
    id <- id[1]
  }

  ## If neither group1 nor group2 are given explicitly, dispatch based on the
  ## number of distinct groups actually present in 'sample_distrib'.
  if(is.null(group1) && is.null(group2)){
    all_groups = sample_distrib$Group %>% unique() %>% sort()

    if(length(all_groups) > 2){
      return(plot_distrib_multi_group(
        sample_distrib, all_groups, id,
        prob_CI, show_prob, mean_bar, plot_mean, top_n_pairs
      ))
    }

    group1 = all_groups[1]
    if(length(all_groups) == 2){
      group2 = all_groups[2]
    }
  }

  build_single_distrib(
    sample_distrib, group1, group2, id,
    prob_CI, show_prob, mean_bar, index_group1, index_group2
  )
}
