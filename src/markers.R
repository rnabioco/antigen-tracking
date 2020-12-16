
#' Find markers with presto
#' 
#' @param sobj_in Seurat object.
#' @param group_column meta.data column containing labels to use for grouping cells.
#' @param exclude_clust Name of cell group to exclude from comparison.
#' @param groups_use Vector of cell groups to use for comparison.
#' @param FC_min Minimum log fold change for markers.
#' @param auc_min Minimum AUC for markers.
#' @param p_max Maximum adjusted p-value for filtering markers.
#' @param pct_in_min Minimum percentage of cells within the group that express each marker.
#' @param pct_out_max Maximum percentage of cells outside the group that express each marker.
#' @param ... Additional arguments to pass to wilcoxauc.
#' @return tibble
#' @export
find_markers <- function(sobj_in, group_column = NULL, exclude_clust = NULL, groups_use = NULL, 
                         FC_min = 0.25, auc_min = 0.5, p_max = 0.05, pct_in_min = 50, 
                         pct_out_max = 100, ...) {
  
  if (!is.null(exclude_clust) && is.null(groups_use)) {
    sobj_in <- sobj_in %>%
      subset(subset = !!sym(group_column) != exclude_clust)
  }
  
  res <- sobj_in %>%
    wilcoxauc(
      group_by = group_column,
      groups_use = groups_use,
      ...
    ) %>%
    as_tibble() %>%
    filter(
      padj    < p_max,
      logFC   > FC_min,
      auc     > auc_min,
      pct_in  > pct_in_min,
      pct_out < pct_out_max
    ) %>%
    arrange(desc(logFC))
  
  res
}

#' Helper to find subtype ova markers
#' 
#' @param sobj_in Seurat object.
#' @param prefix Perfix to add to GMM groups.
#' @param gmm_column meta.data column containing GMM groups.
#' @param so_title Sample title to include in output table.
#' @param p_max Maximum adjusted p-value for filtering markers.
#' @param FC_min Minimum log fold change for markers.
#' @param auc_min Minimum AUC for markers.
#' @param pct_in_min Minimum percentage of cells within the group that express each marker.
#' @param pct_out_max Maximum percentage of cells outside the group that express each marker.
#' @param meta_data Additional meta.data columns to include in output table.
#' @return tibble
#' @export
find_ova_markers <- function(sobj_in, prefix = "", gmm_column, so_title, p_max = 0.05, FC_min = 0.25,
                             auc_min = 0.5, pct_in_min = 50, pct_out_max = 100, meta_data = NULL) {
  
  grps <- c("ova low", "ova high") %>%
    set_names(., .)
  
  if (prefix != "") {
    grps <- grps %>%
      map_chr(~ str_c(prefix, "-", .x))
  }
  
  res <- find_markers(
    sobj_in      = sobj_in,
    group_column = gmm_column,
    groups_use   = grps,
    p_max        = p_max,
    FC_min       = FC_min,
    auc_min      = auc_min,
    pct_in_min   = pct_in_min,
    pct_out_max  = pct_out_max,
    uniq         = FALSE
  )
  
  if (!is.null(meta_data)) {
    meta_data <- sobj_in %>%
      FetchData(c(gmm_column, meta_data)) %>%
      unique() %>%
      as_tibble()
    
    res <- res %>%
      left_join(meta_data, by = c("group" = gmm_column))
  }
  
  res <- res %>%
    mutate(Sample = so_title) %>%
    relocate(Sample, before = feature) %>%
    rename(gene = feature) %>%
    group_by(Sample, group) %>%
    arrange(desc(logFC), .by_group = T) %>%
    ungroup()
  
  res
}

