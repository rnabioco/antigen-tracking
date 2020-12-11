
#' Create GO plot
#' 
#' @param GO_df data.frame of GO results.
#' @param plot_colors Plot colors.
#' @param n_terms Number of terms to label.
#' @return ggplot object
#' @export
create_bubbles <- function(GO_df, plot_colors = c("#E69F00", "#56B4E9", "#009E73", "#d7301f"),
                           n_terms = 15) {
  
  # Check for empty inputs
  if (is_empty(GO_df) || nrow(GO_df) == 0) {
    res <- ggplot() +
      geom_blank()
    
    return(res)
  }
  
  # Shorten GO terms and database names
  GO_data <- GO_df %>%
    mutate(
      term_id = str_remove(term_id, "(GO|KEGG):"),
      term_id = str_c(term_id, " ", term_name),
      term_id = str_to_lower(term_id),
      term_id = str_trunc(term_id, 40, "right"),
      source  = fct_recode(
        source,
        "Biological\nProcess" = "GO:BP",
        "Cellular\nComponent" = "GO:CC",
        "Molecular\nFunction" = "GO:MF",
        "KEGG"                = "KEGG"
      )
    )
  
  # Reorder database names
  plot_levels <- c(
    "Biological\nProcess",
    "Cellular\nComponent",
    "Molecular\nFunction",
    "KEGG"
  )
  
  GO_data <- GO_data %>%
    mutate(source = fct_relevel(source, plot_levels))
  
  # Extract top terms for each database
  top_GO <- GO_data %>%
    group_by(source) %>%
    arrange(p_value) %>%
    dplyr::slice(1:n_terms) %>%
    ungroup()
  
  # Create bubble plots
  res <- GO_data %>%
    ggplot(aes(1.25, -log10(p_value), size = intersection_size)) +
    geom_point(color = plot_colors, alpha = 0.5, show.legend = T) +
    geom_text_repel(
      aes(2, -log10(p_value), label = term_id),
      data         = top_GO,
      size         = 2.3,
      direction    = "y",
      hjust        = 0,
      segment.size = NA
    ) +
    xlim(1, 8) +
    labs(y = "-log10(p-value)") +
    base_theme +
    theme(
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    ) +
    facet_wrap(~ source, scales = "free", nrow = 1)
  
  res
}

#' Create boxplots showing marker gene signal
create_marker_boxes <- function(sobj_in, input_markers, clust_column, box_cols, group = NULL, include_legend = F,
                                all_boxes = F, all_violins = F, order_boxes = T, clust_regex = "\\-[a-zA-Z0-9_ ]+$",
                                n_boxes = 10, median_pt = 0.75, n_rows = 2, pt_mtplyr = 1, exclude_clust = NULL, ...) {
  
  # Retrieve and format data for boxplots
  input_markers <- input_markers %>%
    head(n_boxes)
  
  box_data <- sobj_in %>%
    FetchData(c(clust_column, input_markers)) %>%
    as_tibble(rownames = "cell_id") %>%
    filter(!(!!sym(clust_column) %in% exclude_clust)) %>%
    mutate(grp = str_remove(!!sym(clust_column), clust_regex))
  
  input_markers <- input_markers %>%
    str_trunc(9)
  
  # Filter based on input group
  if (!is.null(group)) {
    box_data <- box_data %>%
      filter(grp == group)
  }
  
  # Format data for plots
  box_data <- box_data %>%
    pivot_longer(
      cols = c(-cell_id, -grp, -!!sym(clust_column)),
      names_to = "key",
      values_to = "Counts"
    ) %>%
    mutate(
      !!sym(clust_column) := fct_relevel(!!sym(clust_column), names(box_cols)),
      key = str_trunc(key, width = 9, side = "right"),
      key = fct_relevel(key, input_markers)
    )
  
  # Order boxes by mean signal
  if (order_boxes) {
    box_data <- box_data %>%
      mutate(!!sym(clust_column) := fct_reorder(!!sym(clust_column), Counts, mean, .desc = T))
  }
  
  n_clust <- box_data %>%
    pull(clust_column) %>%
    n_distinct()
  
  # Create plots
  n_cols <- ceiling(n_boxes / n_rows)
  
  res <- box_data %>%
    ggplot(aes(!!sym(clust_column), Counts, fill = !!sym(clust_column))) + 
    facet_wrap(~ key, ncol = n_cols) +
    scale_color_manual(values = box_cols) +
    base_theme +
    theme(
      panel.spacing.x  = unit(0.7, "cm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = 13),
      legend.position  = "none",
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 13),
      axis.text.x      = element_blank(),
      axis.text.y      = element_text(size = 8),
      axis.ticks.x     = element_blank(),
      axis.line.x      = element_blank()
    )
  
  # Adjust output plot type
  if (n_clust > 6 || all_boxes) {
    res <- res +
      stat_summary(geom = "point", shape = 22, fun = median, size = 0) +
      stat_summary(geom = "point", shape = 22, fun = median, size = median_pt * 2, color = "black") +
      geom_boxplot(
        color          = "white",
        fill           = "white",
        alpha          = 1,
        size           = 0.3,
        outlier.colour = "white",
        outlier.alpha  = 1,
        outlier.size   = 0.1,
        coef           = 0       # To exclude whiskers
      ) +
      geom_boxplot(
        size           = 0.3,
        outlier.colour = "grey85",
        outlier.alpha  = 1,
        outlier.size   = 0.1,
        show.legend    = F,
        coef           = 0,
        fatten = 0
      ) +
      stat_summary(
        aes(color = !!sym(clust_column)),
        geom        = "point",
        shape       = 22,
        fun         = median,
        size        = median_pt,
        stroke      = 0.75,
        fill        = "white",
        show.legend = F
      ) +
      guides(fill = guide_legend(override.aes = list(size = 3.5, stroke = 0.25))) +
      scale_fill_manual(values = box_cols) +
      theme(
        panel.background = element_rect(color = "#fafafa", fill = "#fafafa"),
        panel.spacing.x  = unit(0.2, "cm")
      ) +
      theme(...)
    
  } else if (all_violins) {
    res <- res +
      geom_violin(aes(fill = !!sym(clust_column)), size = 0.2) +
      stat_summary(
        aes(color = !!sym(clust_column)),
        geom   = "point",
        shape  = 22,
        fun    = median,
        size   = median_pt,
        stroke = 0.75,
        fill   = "white"
      ) +
      scale_fill_manual(values = box_cols) +
      scale_color_manual(values = box_cols) +
      theme(
        panel.background = element_rect(color = "#fafafa", fill = "#fafafa"),
        panel.spacing.x  = unit(0.2, "cm")
      ) +
      theme(...)
    
  } else {
    pt_size <- 0.3 * pt_mtplyr
    
    res <- res +
      geom_quasirandom(size = pt_size) +
      theme(...)
  }
  
  # Add legend
  if (include_legend) {
    res <- res +
      guides(color = col_guide) +
      theme(legend.position = "top")
  }
  
  # Add blank space for missing facets
  n_keys <- n_distinct(box_data$key)
  
  if (n_keys <= n_cols && n_rows > 1) {
    n_keys <- if_else(n_keys == 1, 2, as.double(n_keys))
    n_cols <- floor(n_cols / n_keys)
    
    res <- res %>%
      plot_grid(
        ncol = n_cols,
        nrow = 2
      )
  }
  
  res
}

#' Create figure summarizing marker genes
create_marker_fig <- function(sobj_in, input_markers, input_GO, clust_column, input_umap, umap_color, fig_heights = c(0.46, 0.3, 0.3), 
                              GO_genome = "mmusculus", box_colors, n_boxes = 10, umap_outline = NULL, umap_mtplyr = 1, xlsx_name = NULL, 
                              sheet_name = NULL, ...) {
  
  # Set blank plots
  marks_umap <- marks_boxes <- GO_bubbles <- ggplot() +
    geom_blank() +
    theme_void()
  
  if (nrow(input_markers) == 0) {
    res <- plot_grid(
      marks_umap, marks_boxes, GO_bubbles,
      rel_heights = fig_heights,
      ncol = 1
    )
    
    return(res)
  }
  
  # Create UMAPs showing marker gene signal
  top_marks <- input_markers$feature %>%
    head(n_boxes)
  
  clust_legend <- get_legend(input_umap)
  
  input_umap <- input_umap +
    theme(legend.position = "none")
  
  marks_umap <- sobj_in %>%
    create_marker_umaps(
      input_markers = head(top_marks, 7),
      high_col      = umap_color,
      pt_outline    = umap_outline,
      pt_mtplyr     = umap_mtplyr
    ) %>%
    append(list(input_umap), .)
  
  marks_umap <- plot_grid(
    plotlist = marks_umap,
    ncol     = 4,
    nrow     = 2,
    align    = "vh",
    axis     = "trbl"
  )
  
  marks_umap <- plot_grid(
    clust_legend, marks_umap,
    rel_heights = c(0.2, 0.9),
    nrow = 2
  )
  
  # Create boxplots showing marker gene signal
  marks_boxes <- sobj_in %>%
    create_marker_boxes(
      input_markers = top_marks,
      clust_column  = clust_column,
      box_cols      = box_colors,
      n_boxes       = n_boxes,
      plot.margin   = unit(c(0.8, 0.2, 0.2, 0.2), "cm"),
      ...
    )
  
  # Create GO term plots
  if (nrow(input_GO) == 0) {
    res <- plot_grid(
      marks_umap, marks_boxes, GO_bubbles,
      rel_heights = fig_heights,
      ncol = 1
    )
    
    return(res)  
  }
  
  GO_bubbles <- input_GO %>%
    create_bubbles(plot_colors = umap_color) +
    theme(
      plot.margin      = unit(c(0.8, 0.2, 0.2, 0.2), "cm"),
      strip.background = element_blank(),
      strip.text       = element_text(size = 13),
      axis.title.y     = element_text(size = 13),
      axis.text.y      = element_text(size = 8),
      axis.line.x      = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 8)
    )
  
  # Write GO terms to excel file 
  if (!is.null(xlsx_name)) {
    input_GO %>%
      dplyr::select(
        term_name,  term_id,
        source,     effective_domain_size,
        query_size, intersection_size,
        p_value,    significant 
      ) %>%
      arrange(source, p_value) %>%
      write.xlsx(
        file      = str_c(xlsx_name, "_GO.xlsx"),
        sheetName = sheet_name,
        append    = T
      )
  }
  
  # Write markers to excel file
  if (!is.null(xlsx_name)) {
    input_markers %>%
      write.xlsx(
        file      = str_c(xlsx_name, "_markers.xlsx"),
        sheetName = sheet_name,
        append    = T
      )
  }
  
  # Create final figure
  res <- plot_grid(
    marks_umap, marks_boxes, GO_bubbles,
    rel_heights = fig_heights,
    ncol = 1
  )
  
  res
}

# Filter clusters and set cluster order
set_cluster_order <- function(input_cols, input_marks, n_cutoff = 5) {
  input_marks <- input_marks %>%
    group_by(group) %>%
    filter(n() >= n_cutoff) %>%
    ungroup()
  
  marks <- unique(input_marks$group)
  res   <- names(input_cols)
  res   <- res[res %in% marks]
  
  res
}

# Create v1 panel for marker genes
create_marker_panel_v1 <- function(sobj_in, input_cols, input_umap = NULL, clust_column, order_boxes = T,
                                   color_guide = guide_legend(override.aes = list(size = 3.5, shape = 16)),
                                   uniq = params$uniq_GO, umap_mtplyr = 6, xlsx_name = NULL, exclude_clust = NULL,
                                   groups_use = NULL, ...) {
  
  # Set point size
  # ref_mtplyr <- if_else(umap_mtplyr == 1, umap_mtplyr, umap_mtplyr * 2.5)
  umap_mtplyr <- if_else(ncol(sobj_in) < 500, umap_mtplyr, 1)
  ref_mtplyr <- umap_mtplyr
  
  # Find marker genes
  markers <- sobj_in %>%
    find_markers(
      group_column  = clust_column,
      groups_use    = groups_use,
      exclude_clust = exclude_clust
    )
  
  # Find GO terms
  GO_df <- markers
  
  if (nrow(markers) > 0) {
    GO_df <- markers %>%
      group_by(group) %>%
      do({
        arrange(., desc(logFC)) %>%
          pull(feature) %>%
          run_gprofiler(genome = params$genome)
      }) %>%
      ungroup()
  }
  
  if (uniq && nrow(GO_df) > 0) {
    GO_df <- GO_df %>%
      group_by(term_id) %>%
      filter(n() == 1) %>%
      ungroup()
  }
  
  # Set cluster order based on order of input_cols
  fig_clusters <- input_cols %>%
    set_cluster_order(markers)
  
  fig_clusters <- fig_clusters[!fig_clusters %in% exclude_clust]
  
  # Create figures
  for (i in seq_along(fig_clusters)) {
    cat("\n#### ", fig_clusters[i], "\n", sep = "")
    
    # Filter markers and GO terms
    clust <- fig_clusters[i]
    
    fig_marks <- markers %>%
      filter(group == clust)
    
    fig_GO <- GO_df %>%
      filter(group == clust)
    
    # Create reference umap
    ref_umap <- input_umap
    umap_col <- input_cols[clust]
    
    if (is.null(input_umap)) {
      umap_levels <- input_cols[names(input_cols) != clust]
      umap_levels <- names(c(umap_levels, umap_col))
      
      ref_umap <- sobj_in %>%
        create_ref_umap(
          feature     = clust_column,
          plot_cols   = input_cols,
          feat_levels = umap_levels,
          pt_mtplyr   = ref_mtplyr,
          color_guide = color_guide
        )
    }
    
    # Create panel
    marker_fig <- sobj_in %>%
      create_marker_fig(
        input_markers = fig_marks,
        input_GO      = fig_GO,
        clust_column  = clust_column,
        input_umap    = ref_umap,
        umap_color    = umap_col,
        box_colors    = input_cols,
        order_boxes   = order_boxes,
        umap_mtplyr   = umap_mtplyr,
        xlsx_name     = xlsx_name,
        sheet_name    = clust,
        exclude_clust = exclude_clust,
        ...
      )
    
    # Save excel file meta data
    cat(nrow(fig_marks), "marker genes and", nrow(fig_GO), "GO terms were identified.")
    print(marker_fig)
    cat("\n\n---\n\n<br>\n\n<br>\n\n")
  }
  
  # md5sums for output files
  file_out <- str_c(xlsx_name, "_GO.xlsx")
  
  if (file.exists(file_out)) {
    FILES_OUT <<- file_out %>%
      md5sum() %>%
      c(FILES_OUT)
  }
  
  file_out <- str_c(xlsx_name, "_markers.xlsx")
  
  if (file.exists(file_out)) {
    FILES_OUT <<- file_out %>%
      md5sum() %>%
      c(FILES_OUT)
  }
}

# Create v2 panel that splits plots into groups
create_marker_panel_v2 <- function(sobj_in, input_markers, input_cols, grp_column, clust_column,
                                   color_guide = guide_legend(override.aes = list(size = 3.5, shape = 16)), 
                                   uniq_GO = params$uniq_GO, umap_mtplyr = 6, xlsx_name = NULL, 
                                   clust_regex = "\\-[a-zA-Z0-9_ ]+$", ...) {
  
  # Set point size
  # ref_mtplyr <- if_else(ncol(sobj_in) < 500, umap_mtplyr * 2.5, 1)
  umap_mtplyr <- if_else(ncol(sobj_in) < 500, umap_mtplyr, 1)
  ref_mtplyr <- umap_mtplyr
  
  # Figure colors and order
  fig_clusters <- input_cols %>%
    set_cluster_order(input_markers)
  
  # Find GO terms
  GO_df <- input_markers %>%
    group_by(group) %>%
    do({
      arrange(., desc(logFC)) %>%
        pull(feature) %>%
        run_gprofiler(genome = params$genome)
    }) %>%
    ungroup()
  
  if (uniq_GO && nrow(GO_df) > 0) {
    GO_df <- GO_df %>%
      group_by(term_id) %>%
      filter(n() == 1) %>%
      ungroup()
  }
  
  # Create figures
  for (i in seq_along(fig_clusters)) {
    cat("\n#### ", fig_clusters[i], "\n", sep = "")
    
    # Filter markers and GO terms
    clust <- fig_clusters[i]
    
    fig_marks <- input_markers %>%
      filter(group == clust)
    
    fig_GO <- GO_df %>%
      filter(group == clust)
    
    # Set colors
    umap_col <- input_cols[clust]
    
    group <- clust %>%
      str_remove(clust_regex)
    
    grp_regex <- str_c("^", group, "-") %>%
      str_replace("\\+", "\\\\+")            # include this to escape "+" in names
    
    fig_cols <- input_cols[grepl(grp_regex, names(input_cols))]
    fig_cols <- c( "Other" = "#fafafa", fig_cols)
    ref_cols <- fig_cols[names(fig_cols) != clust]
    ref_cols <- c(ref_cols, umap_col)
    
    # Create reference UMAP
    ref_umap <- sobj_in %>%
      FetchData(c("UMAP_1", "UMAP_2", grp_column, clust_column)) %>%
      as_tibble(rownames = "cell_id") %>%
      mutate(!!sym(clust_column) := if_else(
        !!sym(grp_column) != group, 
        "Other", 
        !!sym(clust_column)
      )) %>%
      create_ref_umap(
        feature     = clust_column,
        plot_cols   = ref_cols,
        feat_levels = names(ref_cols),
        pt_mtplyr   = ref_mtplyr,
        color_guide = color_guide
      )
    
    # Create panel
    marker_fig <- sobj_in %>%
      create_marker_fig(
        input_markers = fig_marks,
        input_GO      = fig_GO,
        clust_column  = clust_column,
        input_umap    = ref_umap,
        umap_color    = umap_col,
        box_colors    = fig_cols,
        group         = group,
        umap_mtplyr   = umap_mtplyr,
        xlsx_name     = xlsx_name,
        sheet_name    = clust,
        ...
      )
    
    cat(nrow(fig_marks), "marker genes were identified.", nrow(fig_GO), "GO terms were identified.")
    print(marker_fig)
    cat("\n\n---\n\n<br>\n\n<br>\n\n")
  }
  
  # md5sums for output files
  file_out <- str_c(xlsx_name, "_GO.xlsx")
  
  if (file.exists(file_out)) {
    FILES_OUT <<- file_out %>%
      md5sum() %>%
      c(FILES_OUT)
  }
  
  file_out <- str_c(xlsx_name, "_markers.xlsx")
  
  if (file.exists(file_out)) {
    FILES_OUT <<- file_out %>%
      md5sum() %>%
      c(FILES_OUT)
  }
}
