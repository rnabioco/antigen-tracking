
#' Create two-dimensional scatter plot
plot_features <- function(sobj_in, x = "UMAP_1", y = "UMAP_2", feature, data_slot = "data", 
                          split_id = NULL, pt_size = 0.25, pt_outline = NULL, plot_cols = NULL,
                          feat_levels = NULL, split_levels = NULL, min_pct = NULL, max_pct = NULL, 
                          calc_cor = F, lm_line = F, lab_size = 3.7, lab_pos = c(0.8, 0.9), ...) {
  
  # Format imput data
  counts <- sobj_in
  
  if ("Seurat" %in% class(sobj_in)) {
    vars <- c(x, y, feature)
    
    if (!is.null(split_id)) {
      vars <- c(vars, split_id)
    }
    
    counts <- sobj_in %>%
      FetchData(vars = unique(vars), slot = data_slot) %>%
      as_tibble(rownames = "cell_ids")
  }
  
  # Rename features
  if (!is.null(names(feature))) {
    counts <- counts %>%
      rename(!!!syms(feature))
    
    feature <- names(feature)
  }
  
  if (!is.null(names(x))) {
    counts <- counts %>%
      rename(!!!syms(x))
    
    x <- names(x)
  }
  
  if (!is.null(names(y))) {
    counts <- counts %>%
      rename(!!!syms(y))
    
    y <- names(y)
  }
  
  # Set min and max values for feature
  if (!is.null(min_pct) || !is.null(max_pct)) {
    counts <- counts %>%
      mutate(
        pct_rank = percent_rank(!!sym(feature)),
        max_val  = ifelse(pct_rank > max_pct, !!sym(feature), NA),
        max_val  = min(max_val, na.rm = T),
        min_val  = ifelse(pct_rank < min_pct, !!sym(feature), NA),
        min_val  = max(min_val, na.rm = T),
        !!sym(feature) := if_else(!!sym(feature) > max_val, max_val, !!sym(feature)),
        !!sym(feature) := if_else(!!sym(feature) < min_val, min_val, !!sym(feature))
      )
  }
  
  # Set feature order
  if (!is.null(feat_levels)) {
    counts <- counts %>%
      mutate(!!sym(feature) := fct_relevel(!!sym(feature), feat_levels))
  }
  
  # Set facet order
  if (!is.null(split_id) && length(split_id) == 1) {
    counts <- counts %>%
      mutate(split_id = !!sym(split_id))
    
    if (!is.null(split_levels)) {
      counts <- counts %>%
        mutate(split_id = fct_relevel(split_id, split_levels))
    }
  }
  
  # Calculate correlation
  if (calc_cor) {
    if (!is.null(split_id)) {
      counts <- counts %>%
        group_by(split_id)
    }
    
    counts <- counts %>%
      mutate(
        r       = tidy(cor.test(!!sym(x), !!sym(y)))$estimate,
        r       = round(r, digits = 2),
        pval    = tidy(cor.test(!!sym(x), !!sym(y)))$p.value,
        cor_lab = str_c("r = ", r, ", p = ", format(pval, digits = 2))
        
        # cor_lab = cor(!!sym(x), !!sym(y)),
        # cor_lab = round(cor_lab, digits = 2),
        # cor_lab = str_c("r = ", cor_lab),
      )
    
    if (lab_pos != "strip") {
      counts <- counts %>%
        mutate(
          min_x = min(!!sym(x)),
          max_x = max(!!sym(x)),
          min_y = min(!!sym(y)),
          max_y = max(!!sym(y)),
          lab_x = (max_x - min_x) * lab_pos[1] + min_x,
          lab_y = (max_y - min_y) * lab_pos[2] + min_y
        )
    }
  }
  
  # Create scatter plot
  # To add outline for each cluster create separate layers
  res <- counts %>%
    arrange(!!sym(feature))
  
  if (!is.null(pt_outline)) {
    
    if (!is.numeric(counts[[feature]])) {
      res <- res %>%
        ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature), fill = !!sym(feature)))
      
      feats <- counts[[feature]] %>%
        unique()
      
      if (!is.null(feat_levels)) {
        feats <- feat_levels[feat_levels %in% feats]
      }
      
      for (feat in feats) {
        f_counts <- counts %>%
          filter(!!sym(feature) == feat)
        
        res <- res +
          geom_point(data = f_counts, aes(fill = !!sym(feature)), size = pt_outline, color = "black", show.legend = F) +
          geom_point(data = f_counts, size = pt_size)
      }
      
    } else {
      res <- res %>%
        ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature))) +
        geom_point(aes(fill = !!sym(feature)), size = pt_outline, color = "black", show.legend = F) +
        geom_point(size = pt_size)
    }
    
  } else {
    res <- res %>%
      ggplot(aes(!!sym(x), !!sym(y), color = !!sym(feature))) +
      geom_point(size = pt_size)
  }
  
  # Add regression line
  if (lm_line) {
    res <- res +
      geom_smooth(method = "lm", se = F, color = "black", size = 0.5, linetype = 2)
  }
  
  # Add correlation coefficient label
  if (calc_cor && lab_pos != "strip") {
    res <- res +
      geom_text(
        aes(lab_x, lab_y, label = cor_lab),
        color = "black",
        size  = lab_size,
        check_overlap = T, 
        show.legend = F
      )
  }
  
  # Set feature colors
  if (!is.null(plot_cols)) {
    if (is.numeric(counts[[feature]])) {
      res <- res +
        scale_color_gradientn(colors = plot_cols)
      
    } else {
      res <- res +
        scale_color_manual(values = plot_cols) +
        scale_fill_manual(values = plot_cols)
    }
  }
  
  # Split plot into facets
  cor_labeller <- function(labels) {
    labels %>%
      map(~ {
        cor_labs <- counts %>%
          ungroup() %>%
          select(!!sym(feature), cor_lab) %>%
          unique()
        
        cor_labs <- set_names(cor_labs$cor_lab, cor_labs[[feature]])
        
        str_c(.x, "\n", cor_labs[.x])
      })
  }
  
  if (!is.null(split_id)) {
    if (length(split_id) == 1) {
      
      if (calc_cor && lab_pos == "strip") {
        my_labs <- cor_labeller
        
      } else {
        my_labs <- "label_value"
      }
      
      res <- res +
        facet_wrap(~ split_id, labeller = my_labs, ...)
      
    } else if (length(split_id) == 2) {
      eq <- str_c(split_id[1], " ~ ", split_id[2])
      
      res <- res +
        facet_grid(as.formula(eq), ...)
    }
  }
  
  res
}

#' Plot percentage of cells in given groups
plot_cell_count <- function(sobj_in, group_id, split_id = NULL, group_order = NULL, fill_id, 
                            plot_cols = NULL, x_lab = "Cell type", y_lab = "Fraction of cells",
                            bar_pos = "fill", order_count = T, bar_line = 0, ...) {
  
  res <- sobj_in
  
  if ("Seurat" %in% class(sobj_in)) {
    res <- sobj_in@meta.data %>%
      rownames_to_column("cell_id")
  }
  
  res <- res %>%
    mutate(
      group_id := !!sym(group_id),
      fill_id  := !!sym(fill_id)
    )
  
  if (!is.null(group_order)) {
    res <- res %>%
      mutate(group_id = fct_relevel(group_id, group_order))
  }
  
  if (!is.null(split_id)) {
    res <- res %>%
      mutate(split_id := !!sym(split_id))
  }
  
  if (order_count) {
    res <- res %>%
      mutate(fill_id = fct_reorder(fill_id, cell_id, n_distinct))
  }
  
  res <- res %>%
    ggplot(aes(group_id, fill = fill_id)) +
    geom_bar(position = bar_pos, size = bar_line, color = "black") +
    labs(x = x_lab, y = y_lab)
  
  if (!is.null(plot_cols)) {
    res <- res +
      scale_fill_manual(values = plot_cols)
  }
  
  if (!is.null(split_id)) {
    res <- res +
      facet_wrap(~ split_id, ...)
  }
  
  res
}

#' Create reference UMAP for comparisons
#' 
#' @param sobj_in Seurat object.
#' @param pt_size Point size.
#' @param pt_outline Point outline size.
#' @param pt_mtplyr Point size multiplier.
#' @param color_guide Legend guide.
#' @param ... Additional arguments for plotting.
#' @return ggplot object
#' @export
create_ref_umap <- function(sobj_in, pt_size = 0.1, pt_outline = NULL, pt_mtplyr = 1, color_guide, ...) {
  
  if (is.null(pt_outline)) {
    pt_outline <- pt_size * pt_mtplyr + 0.3
  }
  
  res <- sobj_in %>%
    plot_features(
      pt_size     = pt_size * pt_mtplyr,
      pt_outline  = pt_outline,
      ...
    ) +
    guides(color = color_guide) +
    umap_theme +
    theme(
      legend.position = "top",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10)
    )
  
  res
}

#' Create UMAPs showing marker gene signal
#' 
#' @param sobj_in Seurat object.
#' @param input_markers Marker genes to plot.
#' @param low_col Gradient color for low signal.
#' @param high_col Gradient color for high signal.
#' @param pt_size Point size.
#' @param pt_outline Point outline size.
#' @param pt_mtplyr Point size multiplier.
#' @param low_col_mtplyr Multiplier to adjust gradient scale for low color.
#' @param high_col_mtplyr Multiplier to adjust gradient scale for high color.
#' @return ggplot object
#' @export
create_marker_umaps <- function(sobj_in, input_markers, low_col = "#fafafa", high_col = NULL, pt_size = 0.1, 
                                pt_outline = NULL, pt_mtplyr = 1, low_col_mtplyr = 1, high_col_mtplyr = 1) {
  
  # Set UMAP point size
  pt_size <- pt_size * pt_mtplyr
  
  # Set UMAP colors
  if (!is.null(high_col)) {
    input_markers <- set_names(
      x  = rep(high_col, length(input_markers)),
      nm = input_markers
    )
  }
  
  # Set color multipliers
  if (length(high_col_mtplyr) == 1) {
    high_col_mtplyr <- rep(high_col_mtplyr, length(input_markers))
  }
  
  if (length(low_col_mtplyr) == 1) {
    low_col_mtplyr <- rep(low_col_mtplyr, length(input_markers))
  }
  
  # Plot parameters
  umap_args <- list(
    unname(input_markers),
    names(input_markers),
    low_col_mtplyr,
    high_col_mtplyr
  )
  
  # Create UMAPs
  res <- umap_args %>%
    pmap(~ {
      sobj_in %>%
        plot_features(
          feature    = ..2,
          plot_cols  = c(rep(low_col, ..3), rep(..1, ..4)),
          pt_outline = pt_outline,
          pt_size    = pt_size
        ) +
        ggtitle(.y) +
        umap_theme +
        theme(
          plot.title        = element_text(size = 13),
          legend.position   = "bottom",
          legend.title      = element_blank(),
          legend.text       = element_text(size = 8),
          legend.key.height = unit(0.1, "cm"),
          legend.key.width  = unit(0.3, "cm"),
          axis.title.y      = element_text(size = 13, color = "white"),
          axis.text.y       = element_text(size = 8, color = "white")
        )
    })
  
  res
}


# Figure helpers ----

#' Set equal x-axis scales
#' 
#' @param gg_list_in List of ggplot objects.
#' @param log_tran Log-transform x-axis.
#' @param ... Additional arguments to pass to scale_x_log10
#' @return List of ggplot objects
#' @export
set_equal_x <- function(gg_list_in, log_tran = TRUE, ...) {
  
  set_lims <- function(gg_in, min_x, max_x, log_tran, ...) {
    res <- gg_in +
      coord_cartesian(xlim = c(min_x, max_x))
    
    if (log_tran) {
      res <- res +
        scale_x_log10(labels = trans_format("log10", math_format(10^.x)), ...)
    }
    
    res
  }
  
  gg_ranges <- gg_list_in %>%
    map(~ ggplot_build(.x)$layout$panel_scales_x[[1]]$range$range)
  
  min_val <- gg_ranges %>%
    map_dbl(~ .x[1]) %>%
    min()
  
  max_val <- gg_ranges %>%
    map_dbl(~ .x[2]) %>%
    max()
  
  res <- gg_list_in %>%
    map(
      set_lims,
      min_x    = min_val, 
      max_x    = max_val, 
      log_tran = log_tran,
      breaks   = 10^(-4:4),
      ...
    )
  
  res
}

#' Create GMM figure panels
#' 
#' @param sobj_in Seurat object.
#' @param gmm_column meta.data column containing GMM groups.
#' @param ova_cols Colors for GMM groups.
#' @param bar_cols Colors for bar graphs.
#' @param pt_size UMAP point size.
#' @param pt_outline Point outline size.
#' @param show_bars Include bar graphs.
#' @param plot_labs Labels for figure panels.
#' @param legd_pos Legend position.
#' @param umap_margin Margin adjustment for UMAP.
#' @param ... Additional arguments to pass to plot_grid.
#' @return ggplot object
#' @export
create_gmm_panels <- function(sobj_in, gmm_column = "GMM_grp", ova_cols, bar_cols, pt_size = 0.00001,
                              pt_outline = 0.4, show_bars = TRUE, plot_labs = c(letters[1:3], ""), legd_pos = "top",
                              umap_margin = unit(c(0.2, 1, 0.2, 1.5), "cm"), ...) {
  
  # Theme elements
  guide <- legd_guide(size = cir_size, reverse = TRUE)
  
  # Data for ova group UMAP
  ova_order <- names(ova_cols)
  
  data_df <- sobj_in@meta.data %>%
    as_tibble(rownames = "cell_id") %>%
    group_by(!!sym(gmm_column)) %>%
    mutate(cell_count = n_distinct(cell_id)) %>%
    ungroup() %>%
    mutate(!!sym(gmm_column) := fct_relevel(!!sym(gmm_column), ova_order)) %>%
    arrange(!!sym(gmm_column)) %>%
    mutate(
      cell_count = str_c(!!sym(gmm_column), "\n(n = ", cell_count, ")"),
      cell_count = fct_inorder(cell_count)
    )
  
  # Set ova group colors
  names(ova_order) <- levels(data_df$cell_count)
  
  cols_df <- tibble(
    grp        = ova_order,
    cell_count = names(ova_order)
  )
  
  cols_df <- cols_df %>%
    mutate(color = ova_cols[grp])
  
  umap_cols <- set_names(
    x  = cols_df$color,
    nm = cols_df$cell_count
  )
  
  # Create ova group UMAP
  ova_labeller <- function(labels) {
    labels %>%
      map(str_remove, "^[a-zA-Z0-9 \\-\\+]+-")
  }
  
  ova_grp_umap <- data_df %>%
    plot_features(
      feature     = "cell_count",
      pt_size     = pt_size,
      pt_outline  = pt_outline,
      plot_cols   = umap_cols,
      feat_levels = names(ova_order)
    ) +
    scale_color_manual(values = umap_cols, labels = ova_labeller) +
    scale_fill_manual(values = umap_cols, labels = ova_labeller) +
    guides(color = guide, fill = guide) +
    umap_theme +
    theme(
      plot.margin       = umap_margin,
      legend.position   = legd_pos,
      legend.key.height = unit(1, "cm"),
      legend.title      = element_blank()
    )
  
  # ova hist
  # Add +1 pseudo count to adt_ovalbumin since using log scale
  ova_hist_order <- names(ova_cols) %>%
    grep("ova (low|high)", ., value = T)
  
  ova_hist <- sobj_in@meta.data %>%
    filter(!!sym(gmm_column) != "Other") %>%
    group_by(!!sym(gmm_column)) %>%
    mutate(ave_ova= mean(adt_ovalbumin)) %>%
    ungroup() %>%
    mutate(!!sym(gmm_column) := fct_relevel(!!sym(gmm_column), ova_hist_order)) %>%
    
    ggplot(aes(adt_ovalbumin + 1, after_stat(density), fill = !!sym(gmm_column))) +
    stat_density(geom = "point", position = "identity", size = 0, color = "#ffffff") +  # for square legend symbols
    geom_density(size = 0.3, alpha = 0.8, show.legend = F) +
    
    guides(fill = legd_guide(shape = 22)) +                                             # for square legend symbols
    
    geom_vline(aes(xintercept = ave_ova), size = 0.5, linetype = 2, color = "grey35") +
    coord_cartesian(ylim = c(0, 1.7)) +
    scale_fill_manual(values = ova_cols, labels = ova_labeller) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x)), breaks = 10^(-4:4)) +
    labs(x = "ova signal (UMI counts + 1)", y = "Density") +
    
    theme_minimal_hgrid() +
    base_theme +
    theme(
      plot.margin        = unit(c(0.2, 1, 0.2, 1), "cm"),
      legend.position    = c(0.05, 0.92),
      legend.key.height  = unit(0.5, "cm"),
      legend.title       = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, color = ln_col)
    )
  
  # ova subtype bar graphs
  type_bar <- sobj_in@meta.data %>%
    as_tibble(rownames = "cell_id") %>%
    filter(!!sym(gmm_column) != "Other") %>%
    mutate(
      subtype = fct_reorder(subtype, cell_id, n_distinct),
      !!sym(gmm_column) := fct_relevel(!!sym(gmm_column), c("ova low", "ova high"))
    ) %>%
    
    ggplot(aes(!!sym(gmm_column), fill = subtype)) +
    stat_count(position = "fill", geom = "point", size = 0, color = "#ffffff") +    # for square legend symbols
    geom_bar(position = "fill", size = 0.25, color = "black", show.legend = F) +
    
    guides(fill = legd_guide(ncol = 1, shape = 22)) +                               # for square legend symbols
    
    scale_fill_manual(values = bar_cols) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(y = "Fraction of cells") +
    
    theme_minimal_hgrid() +
    base_theme +
    theme(
      plot.margin        = unit(c(0.2, 6.5, 0.2, 0.2), "cm"),
      legend.title       = element_blank(),
      legend.key.height  = unit(0.45, "cm"),
      axis.title.x       = element_blank(),
      axis.text.x        = element_text(hjust = c(0.6, 0.4)),
      panel.grid.major.y = element_blank()
    )
  
  # Create final figure
  if (!show_bars) {
    type_bar <- ggplot() +
      theme_void()
  }
  
  res <- plot_grid(
    ova_grp_umap, ova_hist, type_bar,
    labels         = plot_labs,
    label_fontface = "plain",
    label_size     = 27,
    ...
  )
  
  res
}

#' Create UMAPs showing feature signal
#' 
#' @param sobj_in Seurat object.
#' @param feats Features to plot.
#' @param ref_cols Colors for reference UMAP showing cell subtypes.
#' @param ref_size Point size for reference UMAP.
#' @param pt_size Point size for feature UMAPs.
#' @param pt_outline Point outline size.
#' @param panels_n_row Number of rows for final figure.
#' @param low_col Gradient color for low signal.
#' @param low_col_mtplyr Multiplier to adjust gradient scale for low color.
#' @param high_col_mtplyr Multiplier to adjust gradient scale for high color.
#' @param on_top Cell subtypes for reference UMAP that should be plotted on
#' top.
#' @param return_list Return a list of plots. If set to FALSE, plots will be
#' be combined with plot_grid().
#' @param ... Additional arguments to adjust the reference UMAP theme.
#' @return ggplot object, or list of ggplot objects
#' @export
create_feat_umap_panels <- function(sobj_in, feats, ref_cols, ref_size = 0.00001, pt_size = ref_size, pt_outline = 0.4,
                                    panels_n_row = 2, low_col = "#ffffff", low_col_mtplyr = 2, high_col_mtplyr = 1,
                                    on_top = NULL, return_list = FALSE, ...) {
  
  # Reference UMAP
  if (!is.null(on_top)) {
    on_top <- ref_cols[on_top]
    
    ref_cols <- ref_cols[!names(ref_cols) %in% names(on_top)]
    ref_cols <- c(ref_cols, on_top)
  }
  
  ref_umap <- sobj_in %>%
    create_ref_umap(
      feature     = "subtype",
      color_guide = legd_guide(size = cir_size),
      plot_cols   = ref_cols,
      pt_size     = ref_size,
      pt_outline  = pt_outline,
      feat_levels = names(ref_cols)
    ) +
    umap_theme +
    theme(
      legend.position   = "left",
      legend.key.height = unit(0.45, "cm"),
      legend.title      = element_blank(),
      legend.text       = element_text(size = txt_pt1)
    ) +
    theme(...)
  
  # Create list of feature UMAPs
  feat_umaps <- sobj_in %>%
    create_marker_umaps(
      pt_size         = pt_size,
      pt_outline      = pt_outline,
      input_markers   = feats,
      low_col         = low_col,
      low_col_mtplyr  = low_col_mtplyr,
      high_col_mtplyr = high_col_mtplyr
    ) %>%
    map(~ {
      .x +
        guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.2)) +
        theme(plot.title = element_text(size = ttl_pt1))
    })
  
  # Return plot list
  if (return_list) {
    return(append(list(ref_umap), feat_umaps))
  }
  
  # Top panel of feature UMAPs
  n_umaps <- (length(feat_umaps) + 2) / panels_n_row
  n_top <- n_umaps - 2
  
  umaps <- append(list(ref_umap), feat_umaps[1:n_top]) %>%
    plot_grid(
      plotlist   = .,
      rel_widths = c(1, rep(0.5, n_top)),
      nrow       = 1
    ) %>%
    list()
  
  panel_heights <- 1
  
  # Bottom panel of feature UMAPs
  if (panels_n_row > 1) {
    bot_umaps <- feat_umaps[(n_top + 1):length(feat_umaps)] %>%
      plot_grid(
        plotlist = .,
        nrow     = panels_n_row - 1,
        align    = "h",
        axis     = "tb"
      )
    
    umaps <- append(umaps, list(bot_umaps))
    panel_heights <- c(1, panels_n_row - 1)
  }
  
  # Final feature UMAP figure
  res <- plot_grid(
    plotlist    = umaps,
    rel_heights = panel_heights,
    ncol        = 1,
    align       = "vh",
    axis        = "trbl"
  ) + 
    theme(plot.margin = unit(c(1, 0.2, 1.5, 0.2), "cm"))
  
  res
}

#' Create boxplots showing feature signal
#' 
#' @param sobj_in Seurat object.
#' @param feats Features to plot.
#' @param grp_column meta.data column containing groups to plot.
#' @param box_cols Colors to use for boxplots.
#' @param median_pt Size of point used to mark median.
#' @param panels_n_col Number of columns for final figure.
#' @param ... Additional arguments to adjust the boxplot theme.
#' @return ggplot object
#' @export
create_feat_boxes <- function(sobj_in, feats, grp_column = "subtype", box_cols,
                              median_pt = 1, panels_n_col = 1, ...) {
  
  # Boxplot data
  # What to arrange boxplots by median and 3rd quartile
  box_data <- sobj_in %>%
    FetchData(c(feats, grp_column)) %>%
    as_tibble(rownames = "cell_id") %>%
    pivot_longer(cols = c(-cell_id, -!!sym(grp_column))) %>%
    mutate(type_name = str_c(!!sym(grp_column), "_", name)) %>%
    group_by(type_name) %>%
    mutate(
      up_qt = boxplot.stats(value)$stats[4],
      med   = median(value)
    ) %>%
    ungroup() %>%
    arrange(desc(med), desc(up_qt)) %>%
    mutate(
      name = fct_relevel(name, feats),
      type_name = fct_inorder(type_name)
    )
  
  # Create boxplots
  feat_boxes <- box_data %>%
    ggplot(aes(type_name, value, fill = !!sym(grp_column))) +
    
    facet_wrap(~ name, ncol = panels_n_col, scales = "free_x") +
    
    scale_fill_manual(values = box_cols) +
    scale_color_manual(values = box_cols) +
    labs(y = "Counts") +
    
    theme_minimal_hgrid() +
    base_theme +
    theme(
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.key.height  = unit(0.45, "cm"),
      axis.title.x       = element_blank(),
      axis.text.x        = element_blank(),
      axis.line.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    theme(...)
  
  res <- feat_boxes +
    stat_summary(geom = "point", shape = 22, fun = median, size = 0) +
    stat_summary(geom = "point", shape = 22, fun = median, size = median_pt + 1, color = "black") +
    geom_boxplot(
      size           = 0.3,
      width          = 0.6,
      fatten         = 0,
      outlier.colour = "grey85",
      outlier.alpha  = 1,
      outlier.size   = 0.1,
      show.legend    = F,
      coef           = 0,
    ) +
    stat_summary(
      aes(color = subtype),
      geom        = "point",
      shape       = 22,
      fun         = median,
      size        = median_pt,
      stroke      = 1,
      fill        = "#ffffff",
      show.legend = F
    ) +
    guides(fill = legd_guide(ncol = 2, shape = 22))
  
  res
}

#' Create violin plots showing DNA-tag counts
#' 
#' @param sobj_in Seurat object.
#' @param plot_cols Colors to use for violin plots.
#' @param box_counts Data to use for violin plots.
#' @param type_column meta.data column containing cell subtypes.
#' @param plot_title Title to add to plot.
#' @param control_types Control subtypes to plot last.
#' @return ggplot object
#' @export
create_count_violins <- function(sobj_in, plot_cols, box_counts, type_column = "subtype",
                                 plot_title = NULL, control_types = c("B cell", "T cell")) {
  
  # Violin plot helper
  create_violins <- function(df_in, filt = NULL, color_id = type_column, plot_cols,
                             gg_title = NULL, plot_levels, as_swarm = F) {
    
    if (!is.null(filt)) {
      df_in <- df_in %>%
        filter(name == filt)
    }
    
    res <- df_in %>%
      mutate(!!sym(color_id) := fct_relevel(!!sym(color_id), plot_levels)) %>%
      ggplot(aes(counts + 1, !!sym(color_id))) +
      facet_wrap(~ name) +
      theme_minimal_vgrid() +
      theme(
        legend.position    = "none",
        plot.margin        = unit(c(0.2, 0.2, 1, 0.2), "cm"),
        plot.title         = element_text(size = 12, face = "plain"),
        strip.text         = element_text(size = 10),
        axis.title.y       = element_blank(),
        axis.title         = element_text(size = 10),
        axis.text          = element_text(size = 8),
        axis.ticks.x       = element_line(size = 0.1),
        panel.grid.major.x = element_line(size = 0.1)
      )
    
    if (as_swarm) {
      res <- res +
        geom_quasirandom(aes(color = !!sym(color_id)), size = 0.5, groupOnX = F) +
        scale_color_manual(values = plot_cols)
      
    } else {
      res <- res +      
        geom_violin(aes(fill = !!sym(color_id)), size = 0.3, draw_quantiles = c(0.25, 0.75), alpha = 0.75) +
        stat_summary(geom = "point", color = "black", fun = median) +
        scale_fill_manual(values = plot_cols)
    }
    
    if (!is.null(gg_title)) {
      res <- res +
        ggtitle(gg_title)
      
    }
    
    res
  }
  
  # Data for violin and beeswarm plots
  box_data <- sobj_in %>%
    FetchData(slot = "counts", vars = c(type_column, box_counts)) %>%
    as_tibble(rownames = "cell_id") %>%
    rename(all_of(box_counts)) %>%
    pivot_longer(cols = c(-cell_id, -!!sym(type_column)), values_to = "counts")
  
  # Plot levels
  box_levels <- box_data %>%
    filter(name == names(box_counts[1])) %>%
    mutate(!!sym(type_column) := fct_reorder(!!sym(type_column), counts, median)) %>%
    pull(type_column) %>%
    levels()
  
  box_levels <- box_levels[!box_levels %in% control_types] %>%
    c(control_types, .)
  
  # Create violin plots
  plot_args <- list(
    filt       = names(box_counts),
    as_swarm   = c(F, T, T),
    gg_title   = list(plot_title, NULL, NULL)
  )
  
  vlns <- plot_args %>%
    pmap(
      create_violins,
      df_in       = box_data,
      plot_cols   = plot_cols,
      plot_levels = box_levels
    )
  
  # Set equal x-axis scales
  vlns <- vlns %>%
    set_equal_x(log_tran = T)
  
  # Create final figure
  res <- plot_grid(
    plotlist = vlns,
    nrow       = 1,
    align      = "vh",
    axis       = "trbl"
  )
  
  res
}

#' Create correlation heatmaps for cell types
#' 
#' @param sobj_in Seurat object.
#' @param ref_in Cell type refence to use for clustifyr.
#' @param type_column meta.data column containing cell subtypes to compare.
#' @param plot_cols Gradient color to use for high values.
#' @param rename_ref Named vector to use for changing reference labels.
#' @param exclude_ref Regular expression to use for filtering cell type references.
#' @param plot_title Title to add to plot.
#' @param return_mat Return a correlation matrix instead of a heatmap.
#' @return ggplot object or correlation matrix
#' @export
create_ref_heatmap <- function(sobj_in, ref_in, type_column = "subtype", plot_cols = "red", rename_ref = NULL,
                               exclude_ref = "^Endothelial cells ", plot_title = NULL, return_mat = FALSE) {
  
  # Get clustifyr matrix
  cor_mat <- sobj_in %>%
    clustify(
      cluster_col = type_column,
      ref_mat     = ref_in,
      seurat_out  = FALSE
    )
  
  if (return_mat) {
    return(cor_mat)
  }
  
  # Heatmaps
  cor_df <- cor_mat %>%
    as_tibble(rownames = "Cell type") %>%
    rename(all_of(rename_ref)) %>%
    pivot_longer(cols = -`Cell type`, names_to = "Reference type", values_to = "r") %>%
    mutate(
      `Reference type` = fct_reorder(`Reference type`, r, mean, .desc = T),
      `Cell type` = fct_reorder(`Cell type`, r, mean)
    )
  
  # Exclude select cell types from reference
  if (!is.null(exclude_ref)) {
    cor_df <- cor_df %>%
      filter(!grepl(exclude_ref, `Reference type`))
  }
  
  # Create heatmap
  res <- cor_df %>%
    ggplot(aes(`Reference type`, `Cell type`, fill = r)) +
    geom_tile(color = "#ffffff", size = 1) +
    scale_fill_gradient(low = "#ffffff", high = plot_cols) +
    guides(fill = guide_colorbar(barwidth = unit(0.2, "cm"), frame.colour = "black", frame.linewidth = 0.1)) +
    base_theme +
    theme(
      plot.title   = element_text(size = 12, face = "plain"),
      strip.text   = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8),
      axis.title   = element_text(size = 10),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text    = element_text(size = 8),
      axis.line.x  = element_blank(),
      axis.line.y  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  
  if (!is.null(plot_title)) {
    res <- res +
      ggtitle(plot_title)
  }
  
  res
}

#' Calculate pairwise p-values for cell subtypes
#' 
#' @param sobj_in Seurat object.
#' @param sample_name Sample name to add to results table.
#' @param data_column meta.data column containing data to use for test.
#' @param type_column meta.data column containing cell subtypes to test.
#' @return tibble
#' @export
calc_p_vals <- function(sobj_in, sample_name, data_column, type_column) {
  
  # Helper to calc p-values for give cell types
  calc_p <- function(x, y, df_in, data_column, type_column) {
    
    x_dat <- df_in %>%
      filter(!!sym(type_column) == x) %>%
      pull(data_column)
    
    y_dat <- df_in %>%
      filter(!!sym(type_column) == y) %>%
      pull(data_column)
    
    res <- wilcox.test(
      x        = x_dat, 
      y        = y_dat, 
      conf.int = T
    ) %>%
      tidy() %>%
      mutate(
        group1 = x,
        group2 = y
      )
    
    res
  }
  
  # Pull data from gg object
  p_data <- sobj_in@meta.data %>%
    as_tibble(rownames = "cell_id") %>%
    mutate("cell_type" = !!sym(type_column))

  # Calculate additional stats for table
  p_stats <- p_data %>%
    group_by(cell_type) %>%
    summarize(
      n_cells = n_distinct(cell_id),
      med     = median(!!sym(data_column)),
      .groups = "drop"
    ) %>%
    mutate(frac_cells = n_cells / sum(n_cells))
  
  # Create data.frame with all combinations of cell types for comparison
  # Flip group columns for B and T cells so group1 includes all cell types
  c_types <- p_data %>%
    pull(type_column) %>%
    unique()
  
  c_comps <- combinations(
    n = length(c_types),
    r = 2,
    v = c_types
  ) %>%
    as_tibble(.name_repair = "unique") %>%
    rename(x = ...1, y = ...2) %>%
    mutate(
      x = if_else(x == "B cell" & y == "T cell", "T cell", x),
      y = if_else(x == "T cell" & y == "T cell", "B cell", y)
    )
  
  # Calculate pairwise wilcox test for all combinations of cell types
  res <- c_comps %>%
    pmap_dfr(
      calc_p,
      df_in       = p_data,
      data_column = data_column,
      type_column = type_column
    )
  
  # Apply multiple testing correction
  res <- res %>%
    mutate(p_adj = p.adjust(p.value, method = "bonf"))
  
  # Add medians to data.frame
  res <- p_stats %>%
    rename(med_1 = med) %>%
    right_join(res, by = c(cell_type = "group1"))
  
  res <- p_stats %>%
    select(cell_type, med_2 = med) %>%
    right_join(res, by = c("cell_type" = "group2"))
  
  # Format final table
  res <- res %>%
    mutate(Sample = sample_name) %>%
    select(
      Sample,
      `Cell type 1`      = cell_type.y,
      `n cells 1`        = n_cells,
      `frac cells 1`     = frac_cells,
      `Median rel ova 1` = med_1,
      `Cell type 2`      = cell_type,
      `Median rel ova 2` = med_2,
      method,
      alternative,
      statistic,
      estimate,
      conf.low,
      conf.high,
      p.value,
      p_adj
    ) %>%
    arrange(`Cell type 1`, p_adj)
  
  res
}

#' Helper to write excel files
#' 
#' @param df_in data.frame to add to excel document.
#' @param file_out Path to output file.
#' @param sheet_id Name of sheet.
#' @return Excel file
#' @export
write_xlsx <- function(df_in, file_out, sheet_id, ...) {
  
  if (file.exists(file_out)) {
    stop("File ", file_out, " exists, nothing written.")
  }
  
  samples <- df_in %>%
    pull(sheet_id) %>%
    unique()
  
  samples <- samples %>%
    str_replace_all("[ \\/]", "_") %>%
    str_remove_all("[\\(\\)]") %>%
    str_remove_all("Cells_") %>%
    set_names(
      x  = samples,
      nm = .
    )
  
  samples %>%
    iwalk(~ {
      df_in %>%
        filter(!!sym(sheet_id) == .x) %>%
        as.data.frame() %>%
        write.xlsx(
          file      = file_out,
          sheetName = .y,
          row.names = F,
          append    = T
        )
    })
}

#' Create files for Cell Browser
#' 
#' @param sobj_in Seurat object.
#' @param markers Marker genes to include in browser.
#' @param browser_dir Path to Cell Browser directory.
#' @param browser_name Name for dataset.
#' @param cols Colors to use for Cell Browser.
#' @param browser_title Title to use for Cell Browser description.
#' @param abstract Abstract to use for Cell Browser description.
#' @param methods_desc Methods to use for Cell Browser description.
#' @param ... Additional arguments to pass to ExportToCellbrowser().
#' @return Cell Browser files
#' @export
create_browser <- function(sobj_in, markers, browser_dir, browser_name, cols = NULL,
                           browser_title = NULL, abstract = NULL, methods_desc = NULL, ...) {
  
  # Set file paths
  brow_path   <- file.path(browser_dir, browser_name)  
  marker_path <- file.path(brow_path, "markers.tsv")
  
  dir.create(brow_path, recursive = T)
  
  # Create tsv file for marker genes
  markers %>%
    select(
      cluster = group,
      gene,
      score = logFC
    ) %>%
    write_tsv(marker_path)
  
  # Change adt_ovalbumin column name so raw counts are pulled
  sobj_in@meta.data <- sobj_in@meta.data %>%
    rename(ova_counts = adt_ovalbumin)
  
  # Create Cell Browser files for gene expression data
  sobj_in %>%
    ExportToCellbrowser(
      markers.file   = marker_path,
      dir            = brow_path,
      reductions     = "umap",
      dataset.name   = browser_name,
      subtype        = "Cell type",
      GMM_grp        = "Cluster",
      type_GMM_grp_2 = "Cell type ova groups",
      nCount_RNA     = "RNA UMI count",
      nFeature_RNA   = "Gene count",
      Percent_mito   = "Percent mito",
      ova_counts     = "ova counts",
      ova_fc         = "Relative ova signal",
      ...
    )
  
  # Helper to write desc file
  write_desc <- function(str_in, key) {
    if (!is.null(str_in)) {
      str_in %>%
        str_trim() %>%
        str_replace_all("\\]", " ") %>%
        str_remove_all("[\\[\\]]") %>%
        str_c("\"", ., "\"\n") %>%
        str_c(key, "=", .) %>%
        cat(file = desc_path, append = T)
    }
  }
  
  # Write desc.conf
  desc_path <- file.path(brow_path, "desc.conf")
  
  write_desc(browser_title, "title")
  write_desc(abstract, "abstract")
  write_desc(methods_desc, "methods")
  
  # Save color.tsv
  if (!is.null(cols)) {
    cols %>%
      write_tsv(file.path(brow_path, "colors.tsv"), col_names = F)
  }
}


# Create figure panels ----

#' Create panels for figure 3
#' 
#' @param sobj_in Seurat object.
#' @param cols_in Colors to use for cell subtypes.
#' @param subtype_column meta.data column containing cell subtypes.
#' @param data_slot Slot to pull data.
#' @param ova_cols Colors to use for plotting signal on UMAP.
#' @param box_counts Data to use for violin plots.
#' @param umap_counts Data to use for plotting signal on UMAP.
#' @param plot_title Title to use for plot.
#' @param pt_size Point size for cell subtype UMAP.
#' @param pt_outline Point outline size for cell subtype UMAP.
#' @param pt_size_2 Point size for plotting signal on UMAP.
#' @param pt_outline_2 Point outline size for plotting signal on UMAP.
#' @param box_cell_count For violin plots include cell count label.
#' @param control_types Control cell subtypes that should be plotted last.
#' @param on_top Cell subtypes that should be plotted first.
#' @param ... Additional arguments to adjust the cell type UMAP theme.
#' @return List of ggplot objects
#' @export
create_fig3 <- function(sobj_in, cols_in, subtype_column = "subtype", data_slot = "counts", ova_cols = c("#fafafa", "#d7301f"),
                        box_counts = c("Relative ova signal" = "ova_fc"), umap_counts = c("ova counts" = "adt_ovalbumin"), 
                        plot_title = NULL, pt_size = 0.1, pt_outline = 0.4, pt_size_2 = 0.3, pt_outline_2 = 0.5, box_cell_count = TRUE,
                        control_types = c("B cell", "T cell"), on_top = NULL, ...) {
  
  box_column <- umap_column <- "cell_type"
  box_cols <- umap_cols <- cols_in
  
  # Fetch plotting data
  data_df <- sobj_in %>%
    FetchData(c(subtype_column, box_counts, umap_counts, "UMAP_1", "UMAP_2"), slot = data_slot) %>%
    as_tibble(rownames = "cell_id")
  
  if (!is.null(names(box_counts))) {
    data_df <- data_df %>%
      rename(!!box_counts)
    
    box_counts <- names(box_counts)
  }
  
  if (!is.null(names(umap_counts))) {
    data_df <- data_df %>%
      rename(!!umap_counts)
    
    umap_counts <- names(umap_counts)
  }
  
  # Set subtype order
  # Move select cell types to front of order
  # Rename subtype_column to 'cell_type'
  data_df <- data_df %>%
    mutate(
      cell_type = !!sym(subtype_column),
      cell_type = fct_reorder(cell_type, !!sym(box_counts), median)
    )
  
  type_order <- levels(data_df$cell_type)
  
  if (!is.null(control_types)) {
    control_types <- control_types[control_types %in% type_order]
    type_order    <- type_order[!type_order %in% control_types]
    type_order    <- c(control_types, type_order)
  }
  
  # Count cells for each subtype
  data_df <- data_df %>%
    group_by(cell_type) %>%
    mutate(cell_count = n_distinct(cell_id)) %>%
    ungroup() %>%
    mutate(cell_type = fct_relevel(cell_type, type_order)) %>%
    arrange(cell_type) %>%
    mutate(
      cell_count = str_c(cell_type, "\n(n = ", cell_count, ")"),
      cell_count = fct_inorder(cell_count)
    )
  
  # Set cell type colors
  names(type_order) <- levels(data_df$cell_count)
  
  cols_df <- tibble(
    cell_type = type_order,
    cell_count = names(type_order)
  )
  
  cols_df <- cols_df %>%
    mutate(color = cols_in[cell_type])
  
  # Subtype UMAP
  # Reorder layers of UMAP
  if (!is.null(on_top)) {
    on_top <- umap_cols[on_top]
    
    umap_cols <- umap_cols[!names(umap_cols) %in% names(on_top)]
    umap_cols <- c(umap_cols, on_top)
  }
  
  umap <- data_df %>%
    plot_features(
      feature     = umap_column,
      pt_size     = pt_size,
      pt_outline  = pt_outline,
      plot_cols   = umap_cols,
      feat_levels = names(umap_cols)
    ) +
    guides(color = guide_legend(override.aes = list(size = cir_size))) +
    ggtitle(plot_title) +
    umap_theme +
    theme(
      plot.title = element_text(size = 12),
      legend.position = "none"
    ) +
    theme(...)
  
  # ova UMAP
  if (!is.null(umap_counts)) {
    ova_umap <- data_df %>%
      plot_features(
        feature    = umap_counts,
        plot_cols  = ova_cols,
        pt_size    = pt_size_2,
        pt_outline = pt_outline_2,
        min_pct    = 0.01,
        max_pct    = 0.99
      ) +
      guides(color = guide_colorbar(frame.colour = "black", frame.linewidth = 0.2)) +
      umap_theme +
      theme(
        plot.title        = element_text(size = 10, hjust = 0.5),
        legend.position   = "right",
        legend.key.width  = unit(0.15, "cm"),
        legend.key.height = unit(0.30, "cm"),
        legend.title      = element_text(size = 10),
        legend.text       = element_text(size = 8)
      )
  }
  
  # ova boxes
  if (box_cell_count) {
    box_column <- "cell_count"
    box_cols <- set_names(
      x  = cols_df$color,
      nm = cols_df$cell_count
    )
  }
  
  boxes <- data_df %>%
    ggplot(aes(!!sym(box_counts), !!sym(box_column), fill = !!sym(box_column))) +
    geom_violin(size = 0.3, draw_quantiles = c(0.25, 0.75), alpha = 0.75) +
    stat_summary(geom = "point", color = "black", fun = median) +
    
    scale_color_manual(values = box_cols) +
    scale_fill_manual(values = box_cols) +
    theme_minimal_vgrid() +
    theme(
      legend.position    = "none",
      axis.title.y       = element_blank(),
      axis.title         = element_text(size = 10),
      axis.text          = element_text(size = 8),
      axis.ticks.x       = element_line(size = 0.1),
      panel.grid.major.x = element_line(size = 0.1)
    )
  
  res <- list(umap, boxes)
  
  if (!is.null(umap_counts)) {
    res <- append(res, list(ova_umap))
  }
  
  names(res) <- rep(plot_title, length(res))
  
  res
}

#' Create panels for figure 4
#' 
#' @param sobj_in Seurat object.
#' @param feat_cols Named vector of colors with features to plot as names.
#' @param ref_cols Colors to use for cell subtypes.
#' @param ova_cols Colors to use for ova-low and -high groups.
#' @param pt_size Point size to use for UMAPs.
#' @param pt_outline Point outline size to use for UMAPs.
#' @param median_pt Size of point used to mark median.
#' @param panels_n_col Number of columns for feature UMAPs.
#' @param low_col_mtplyr Multiplier for feature UMAPs to adjust gradient scale
#' for low color.
#' @param high_col_mtplyr Multiplier for feature UMAPs to adjust gradient scale
#' for high color.
#' @param on_top Cell subtypes that should be plotted first.
#' @param gmm_column meta.data column containing GMM groups.
#' @param gmm_grps GMM groups to plot, groups not included will be labeled as
#' 'Other'.
#' @param panel_heights Relative heights for figure panels.
#' @param plot_labs Labels to use for each panel (e.g. 'A', 'B', ...).
#' @param box_pad Additional padding to add to the right of the feature boxplots.
#' @param umap_legd_pad Additional padding to add to the top of the legend for
#' the cell type reference UMAP.
#' @param ... Additional arguments to use for create feature UMAPs.
#' @return ggplot object
#' @export
create_fig4 <- function(sobj_in, feat_cols, ref_cols, ova_cols, pt_size = 0.00001, pt_outline = 0.4, median_pt = 1.5,
                        panels_n_col = 4, low_col_mtplyr = 1, high_col_mtplyr = 1, on_top = NULL, gmm_column = "GMM_grp",
                        gmm_grps = c("ova low", "ova high"), panel_heights = c(0.5, 1, 0.4), plot_labs, box_pad = 0.2, 
                        umap_legd_pad = 2, ...) {
  
  # Set GMM groups in Seurat object
  sobj_in@meta.data <- sobj_in@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(!!sym(gmm_column) := ifelse(
      !(!!sym(gmm_column)) %in% gmm_grps,
      "Other",
      !!sym(gmm_column)
    )) %>%
    column_to_rownames("cell_id")
  
  # GMM summary plots
  gmm_panels <- create_gmm_panels(
    sobj_in         = sobj_in,
    gmm_column      = gmm_column,
    ova_cols        = ova_cols,
    bar_cols        = ref_cols,
    pt_size         = pt_size,
    pt_outline      = pt_outline,
    show_bars       = F,
    legd_pos        = "left",
    umap_margin     = unit(c(0.2, 0.2, 0.2, 0.4), "cm"),
    plot_labs       = c(plot_labs[1:2], "", ""), 
    rel_widths      = c(0.86, 1, 0.7, 0.2),
    nrow            = 1,
    align           = "h",
    axis            = "tb"
  )
  
  # Feature UMAPs
  feat_umaps <- create_feat_umap_panels(
    sobj_in         = sobj_in,
    ref_cols        = ref_cols,
    feats           = feat_cols,
    pt_size         = pt_size,
    pt_outline      = pt_outline,
    on_top          = on_top,
    low_col_mtplyr  = low_col_mtplyr,
    high_col_mtplyr = high_col_mtplyr,
    return_list     = TRUE,
    legend.margin   = margin(0.2, 0, 0.2, 2, "cm"),
    ...
  )
  
  umap_legd <- get_legend(feat_umaps[[1]]) %>%
    plot_grid(
      ., NULL,
      rel_heights = c(1, umap_legd_pad),
      ncol = 1
    )
  
  feat_umaps[[1]] <- feat_umaps[[1]] +
    theme(legend.position = "none")
  
  feat_umaps <- feat_umaps %>%
    plot_grid(
      plotlist = .,
      ncol     = panels_n_col,
      align    = "vh",
      axis     = "trbl"
    )
  
  feat_umaps <- feat_umaps %>%
    plot_grid(
      umap_legd, .,
      rel_widths = c(0.1, 1),
      nrow = 1
    )
  
  # Feature boxplots
  feat_boxes <- sobj_in %>%
    create_feat_boxes(
      feats           = names(feat_cols),
      grp_column      = "subtype",
      box_cols        = ref_cols,
      median_pt       = median_pt,
      panels_n_col    = length(feat_cols),
      legend.position = "left",
      legend.margin   = margin(0.8, 0.2, 0.2, 0.3, "cm")
    ) +
    guides(fill = legd_guide(ncol = 1, shape = 22))
  
  feat_boxes <- plot_grid(
    feat_boxes, NULL,
    rel_widths = c(1, box_pad)
  )
  
  # Create final figure
  # Add spacing between panels
  panel_heights <- c(
    panel_heights[1], 0.05,
    panel_heights[2], 0.05,
    panel_heights[3], 0.15
  )
  
  plot_list <- list(
    gmm_panels, NULL,
    feat_umaps, NULL,
    feat_boxes, NULL
  )
  
  res <- plot_grid(
    plotlist       = plot_list,
    rel_heights    = panel_heights,
    ncol           = 1,
    labels         = c("", "", plot_labs[3], "", plot_labs[4]),
    label_fontface = "plain",
    label_size     = 27,
    label_y        = 1.05
  )
  
  res
}

#' Create panels for figure 5
#' 
#' @param sobj_in Seurat object.
#' @param feat_cols Named vector of colors with features to plot as names.
#' @param ref_cols Colors to use for cell subtypes.
#' @param ova_cols Colors to use for ova-low and -high groups.
#' @param ref_size Point size for cell type reference UMAP.
#' @param pt_size Point size to use for UMAPs.
#' @param pt_outline Point outline size to use for UMAPs.
#' @param median_pt Size of point used to mark median.
#' @param low_col Gradient color for low signal.
#' @param low_col_mtplyr Multiplier for feature UMAPs to adjust gradient scale
#' for low color.
#' @param high_col_mtplyr Multiplier for feature UMAPs to adjust gradient scale
#' for high color.
#' @param on_top Cell subtypes that should be plotted first.
#' @param gmm_column meta.data column containing GMM groups.
#' @param gmm_grps GMM groups to plot, groups not included will be labeled as
#' 'Other'.
#' @return ggplot object
#' @export
create_fig5 <- function(sobj_in, feat_cols, ref_cols, ova_cols, ref_size = 0.00001, pt_size = ref_size, pt_outline = 0.4, 
                        median_pt = 1.5, low_col = "#ffffff", low_col_mtplyr = 2, high_col_mtplyr = 1, on_top = NULL, 
                        gmm_column = "GMM_grp", gmm_grps = c("ova low", "ova high")) {
  
  # Set GMM groups in Seurat object
  sobj_in@meta.data <- sobj_in@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(!!sym(gmm_column) := ifelse(
      !(!!sym(gmm_column)) %in% gmm_grps,
      "Other",
      !!sym(gmm_column)
    )) %>%
    column_to_rownames("cell_id")
  
  # GMM summary plots
  gmm_panels <- create_gmm_panels(
    sobj_in     = sobj_in,
    gmm_column  = gmm_column,
    ova_cols    = ova_cols,
    bar_cols    = ref_cols,
    pt_size     = pt_size,
    pt_outline  = pt_outline,
    show_bars   = T,
    umap_margin = unit(c(0.2, 2, 0.2, 1.5), "cm"),
    plot_labs   = letters[1:3],
    rel_heights = c(1, 0.65, 0.7),
    nrow        = 3,
    align       = "v",
    axis        = "l"
  )
  
  # Feature UMAPs
  feat_umaps <- create_feat_umap_panels(
    sobj_in         = sobj_in,
    ref_cols        = ref_cols,
    feats           = feat_cols,
    ref_size        = ref_size,
    pt_size         = pt_size,
    pt_outline      = pt_outline,
    low_col_mtplyr  = low_col_mtplyr,
    high_col_mtplyr = high_col_mtplyr,
    on_top          = on_top,
    return_list     = T,
    legend.margin   = margin(0.1, 0, 0.1, 1.5, "cm")
  )
  
  umap_legd <- get_legend(feat_umaps[[1]]) %>%
    plot_grid(
      ., NULL,
      rel_heights = c(0.4, 1),
      ncol = 1
    )
  
  feat_umaps[[1]] <- feat_umaps[[1]] +
    theme(legend.position = "none")
  
  feat_umaps <- feat_umaps %>%
    plot_grid(
      plotlist = .,
      ncol     = 2,
      align    = "vh",
      axis     = "trbl"
    )
  
  feat_umaps <- feat_umaps %>%
    plot_grid(
      umap_legd, .,
      rel_widths = c(0.1, 1),
      nrow = 1
    )
  
  feat_boxes <- sobj_in %>%
    create_feat_boxes(
      feats       = names(feat_cols),
      grp_column   = "subtype",
      box_cols     = ref_cols,
      median_pt    = median_pt,
      panels_n_col = 1,
      plot.margin  = unit(c(0.1, 4, 0.1, 2.25), "cm")
    ) +
    labs(y = "log normalized counts")
  
  # Create final figure
  res <- plot_grid(
    gmm_panels, feat_umaps, feat_boxes,
    rel_widths     = c(0.8, 1, 0.6),
    nrow           = 1,
    labels         = c("", "d", "e"),
    label_fontface = "plain",
    label_size     = 27
  )
  
  res
}

#' Create panels for figure S6
#' 
#' Generates scatter plots comparing ova counts and mRNA counts.
#' 
#' @param sobj_in Seurat object.
#' @param x Feature to plot on the x-axis.
#' @param y Feature to plot on the y-axis.
#' @param feat Feature to use for splitting plots.
#' @param data_slot Slot to pull data.
#' @param cols_in Colors to use for scatter plots.
#' @param plot_title Title to use for plot.
#' @param add_y_pseudo Pseudo count to add to the y-axis data.
#' @param ... Additional arguments to use for generating plots.
#' @return ggplot object
#' @export
create_figS6 <- function(sobj_in, x, y, feat, data_slot, cols_in, plot_title,
                         add_y_pseudo = 1, ...) {
  
  res <- sobj_in %>%
    FetchData(c(x, y, feat), slot = data_slot) %>%
    mutate(
      !!sym(x) := log10(!!sym(x)),
      !!sym(y) := log10(!!sym(y) + add_y_pseudo)
    ) %>%
    plot_features(
      x         = x,
      y         = y,
      feature   = feat,
      data_slot = "counts",
      plot_cols = cols_in,
      lab_pos   = "strip",
      lab_size  = 5,
      calc_cor  = T,
      lm_line   = T,
      ...
    ) +
    ggtitle(plot_title) +
    guides(color = guide_legend(override.aes = list(size = cir_size))) +
    base_theme +
    theme(legend.title = element_blank())
  
  res
}
