
# Utility ----

#' Plot color palette
#' 
#' @param cols_in Character vector containing colors to plot
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_color_palette <- function(cols_in, ...) {
  
  col_df <- tibble(color = factor(cols_in))
  
  if (!is.null(names(cols_in))) {
    col_df <- col_df %>%
      mutate(color = factor(names(cols_in)))
  }
  
  res <- col_df %>%
    ggplot(aes(x = "color", fill = color)) +
    geom_bar(...) +
    scale_fill_manual(values = cols_in) +
    theme_void()
  
  res
}

#' Create function to generate color gradient
#' 
#' @param cols_in Character vector containing colors
#' @return Color gradient function
#' @export
create_col_fun <- function(cols_in) {
  
  create_gradient <- function(cols_in, n = NULL) {
    
    if (is.null(n)) {
      n <- length(cols_in)
    }
    
    colorRampPalette(cols_in)(n)
  }
  
  function(n = NULL) {
    create_gradient(cols_in, n)
  }
}

#' Capitalize first character without modifying other characters
#' 
#' @param string Charactor string to modify
#' @param exclude Word to exclude from output
#' @return Character string
#' @export
str_to_title_v2 <- function(string, exclude = "cell") {
  
  cap_first_chr <- function(string, exclude) {
    chrs <- string %>%
      str_split(pattern = "") %>%
      unlist()
    
    if (any(chrs %in% LETTERS) || string == exclude) {
      return(string)
    }
    
    chrs[1] <- str_to_upper(chrs[1])
    
    res <- chrs %>%
      reduce(str_c)
    
    res
  }
  
  res <- string %>%
    map_chr(~ {
      .x %>%
        str_split(pattern = " ") %>%
        unlist() %>%
        map_chr(cap_first_chr, exclude = exclude) %>%
        reduce(str_c, sep = " ")
    })
  
  res
}

#' Assign colors for cell subtypes
#' 
#' @param type_in Cell type to assign colors for. Colors will be assigned for
#' each subtype belonging to the cell type.
#' @param sobjs_in List of Seurat objects. Include all objects to ensure that
#' subtypes shared between multiple objects are given the same colors.
#' @param type_key Named vector for cell types with samples as the names. This
#' is used to determine which objects are for the same cell type.
#' @param type_column meta.data column containing cell subtypes to use for
#' assigning colors.
#' @param cols_in Vector of colors.
#' @param other_cols Vector containing set colors for given cells types.
#' @return Named vector containing colors
#' @export
set_type_cols <- function(type_in, sobjs_in, type_key, type_column = "subtype",
                          cols_in, other_cols) {
  
  set_cols <- function(types_in, cols_in, other_cols) {
    
    types_in <- types_in[!types_in %in% names(other_cols)]
    cols_in <- cols_in[!cols_in %in% other_cols]
    
    names(cols_in) <- types_in
    cols_in <- cols_in[!is.na(names(cols_in))]
    
    res <- c(cols_in, other_cols)
    
    res
  }
  
  sobjs_in <- sobjs_in[type_key[names(sobjs_in)] == type_in]
  
  res <- sobjs_in %>%
    map(~ {
      .x@meta.data %>%
        pull(type_column) %>%
        unique()
    }) %>%
    reduce(c) %>%
    unique() %>%
    set_cols(
      cols_in = cols_in,
      other_cols = other_cols
    )
  
  res
}

#' Modify meta.data for Seurat object
#' 
#' @param sobj_in Seurat object
#' @param ... Arguments to pass to specified function
#' @param .fun Function to use for modifying Seurat object
#' @return Seurat object
#' @export
mutate_metadata <- function(sobj_in, ..., .fun = dplyr::mutate) {
  meta_df <- sobj_in@meta.data
  
  meta_df <- tibble::rownames_to_column(meta_df, "cell_id")
  meta_df <- .fun(meta_df, ...)
  meta_df <- tibble::column_to_rownames(meta_df, "cell_id")
  
  sobj_in@meta.data <- meta_df
  
  sobj_in
}

#' Create labeller function to add cell n labels
#' 
#' @param sobj_in Seurat object or data.frame containing plot data.
#' @param lab_col meta.data column containing cell groups.
#' @param nm Should cell group be included in label.
#' @param sep Separator to use for labels.
#' @return Labeller function
#' @export
get_nlab_fun <- function(sobj_in, lab_col, nm = TRUE, sep = "\n") {
  
  dat <- sobj_in
  
  if ("Seurat" %in% class(sobj_in)) {
    dat <- sobj_in@meta.data
  }
  
  labs <- dat %>%
    group_by(!!sym(lab_col)) %>%
    summarize(
      n       = format(n(), big.mark = ",", scientific = FALSE),
      n       = str_c("n = ", n),
      .groups = "drop"
    ) %>%
    distinct()
  
  if (nm) {
    labs <- labs %>%
      mutate(n = str_c(!!sym(lab_col), sep, "(", n, ")"))
  }
  
  labs <- set_names(
    x  = pull(labs, "n"),
    nm = pull(labs, lab_col)
  )
  
  res <- function(x) labs[x]
  
  res
}

# Processing helpers ----

#' Generate shortened project names from matrix path
#' 
#' @param str_in Input string
#' @param extra_path Portion of path to remove from str_in
#' @return Shortened name
#' @export
shorten_names <- function(str_in, extra_path = "/outs/(filtered|raw)_feature_bc_matrix$") {
  res <- str_in %>%
    str_remove(extra_path) %>%
    basename() %>%
    str_extract("^[a-zA-Z0-9_]+") %>%
    str_remove("^GEX_")
  
  res
}

#' Import matrices and create Seurat object
#' 
#' @param matrix_dir Directory containing matrix generated by Cell Ranger.
#' @param proj_name Project name to include in meta.data table.
#' @param hash_ids Name of cell hashing antibodies included in matrix.
#' @param adt_count_min If CITE-seq was performed, this option will remove
#' antibodies where the sum total counts is less than adt_count_min.
#' @param gene_min Minimum number of detected genes for cell.
#' @param gene_max Maximum number of detected genes for cell.
#' @param mito_max Maximum percentage of mitochondrial reads for cell.
#' @param mt_str String to use for identifying mitochondrial genes.
#' @param rna_assay Name of RNA assay if multiple assays are being added to the
#' object (e.g. if CITE-seq data is included).
#' @param adt_assay Name of ADT assay for Seurat object.
#' @return Seurat object
#' @export
create_sobj <- function(matrix_dir, proj_name = "SeuratProject", hash_ids = NULL, adt_count_min = 0,
                        gene_min = 250, gene_max = 5000, mito_max = 15, mt_str = "^mt-",
                        rna_assay = "Gene Expression", adt_assay = "Antibody Capture") {
  
  # Load matrices
  mat_list <- Read10X(matrix_dir)
  rna_mat  <- mat_list
  
  # Create Seurat object using gene expression data
  if (is_list(mat_list)) {
    rna_mat <- mat_list[[rna_assay]]
  }
  
  res <- rna_mat %>%
    CreateSeuratObject(
      project   = proj_name,
      min.cells = 5
    )
  
  # Add antibody capture data to Seurat object
  if (is_list(mat_list)) {
    adt_mat <- mat_list[[adt_assay]]
    
    # Double check that cells match for both assays
    if (!identical(colnames(res), colnames(adt_mat))) {
      adt_mat <- adt_mat[, colnames(res)]
      
      warning("Not all cells are shared between RNA and ADT assays.")
    }
    
    # Remove ADT features that have low total counts and likely failed or
    # were omitted
    n_feats    <- nrow(adt_mat)
    count_sums <- rowSums(as.matrix(adt_mat))
    
    adt_mat <- adt_mat[count_sums >= adt_count_min, ]
    
    if (n_feats != nrow(adt_mat)) {
      warning("Some ADT features were removed due to low counts (<", adt_count_min, ").")
    }
    
    res[["ADT"]] <- CreateAssayObject(adt_mat)
  }
  
  # Calculate percentage of mitochondrial reads
  res <- res %>%
    PercentageFeatureSet(
      pattern  = mt_str,
      col.name = "Percent_mito"
    )
  
  # Add QC classifications to meta.data
  res <- res %>%
    mutate_metadata(
      qc_class = case_when(
        Percent_mito > mito_max ~ "High mito reads",
        nFeature_RNA > gene_max ~ "High gene count",
        nFeature_RNA < gene_min ~ "Low gene count",
        TRUE ~ "Pass filters"
      )
    )
  
  res
}

#' Filter and normalize Seurat object
#' 
#' @param sobj_in Seurat object.
#' @param rna_assay Name of RNA assay in object.
#' @param adt_assay Name of ADT assay in object.
#' @param cc_scoring Score cell cycle genes.
#' @param regress_vars Variables to regress out when scaling data
#' @param rna_method Method to use with NormalizeData for RNA assay.
#' @param adt_method Method to use with NormalizeData for ADT assay.
#' @return Seurat object
#' @export
norm_sobj <- function(sobj_in, rna_assay = "RNA", adt_assay = "ADT", cc_scoring = FALSE,
                      regress_vars = NULL, rna_method = "LogNormalize", adt_method = "CLR") {
  
  # Normalize counts
  res <- sobj_in %>%
    subset(qc_class == "Pass filters") %>%
    NormalizeData(
      assay                = rna_assay,
      normalization.method = rna_method
    )
  
  # Score cell cycle genes
  # cc.genes comes loaded with Seurat
  if (cc_scoring) {
    s.genes <- cc.genes$s.genes %>%
      str_to_title()
    
    g2m.genes <- cc.genes$g2m.genes %>%
      str_to_title()
    
    res <- res %>%
      CellCycleScoring(
        s.features = s.genes,
        g2m.features = g2m.genes
      )
  }
  
  # Scale data
  # By default variable features will be used
  res <- res %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures        = 2000
    ) %>%
    ScaleData(vars.to.regress = regress_vars)
  
  # Normalize ADT data
  if (adt_assay %in% names(res)) {
    res <- res %>%
      NormalizeData(
        assay                = adt_assay,
        normalization.method = adt_method
      ) %>%
      ScaleData(assay = adt_assay)
  }
  
  res
}

#' Run PCA and UMAP for gene expression data
#' 
#' @param sobj_in Seurat object.
#' @param assay Name of assay in object.
#' @param dims Dimensions to use for UMAP.
#' @param prefix Prefix to add to reduction keys and meta.data columns.
#' @param pca_meta Should PC-1 and PC-2 coordinates be added to meta.data table.
#' @param umap_meta Should UMAP coordinates be added to meta.data table.
#' @param ... Additional arguments to pass to RunPCA.
#' @return Seurat object
#' @export
run_UMAP_RNA <- function(sobj_in, assay = "RNA", dims = 1:40, prefix = "",
                         pca_meta = TRUE, umap_meta = TRUE, ...) {
  
  # Reduction keys
  pca_name  = str_c(prefix, "pca")
  pca_key   = str_c(prefix, "PC_")
  umap_name = str_c(prefix, "umap")
  umap_key  = str_c(prefix, "UMAP_")
  
  # Run PCA and UMAP
  # By default only variable features are used for PCA
  res <- sobj_in %>%  
    RunPCA(
      assay          = assay,
      reduction.name = pca_name,
      reduction.key  = pca_key,
      ...
    ) %>%
    RunUMAP(
      assay          = assay,
      dims           = dims,
      reduction.name = umap_name,
      reduction.key  = umap_key
    )
  
  # Add PCA to meta.data
  if (pca_meta) {
    pca_columns = str_c(pca_key, c(1, 2))
    
    res <- res %>%
      AddMetaData(
        metadata = FetchData(., pca_columns),
        col.name = pca_columns
      )
  }
  
  # Add UMAP to meta.data
  if (umap_meta) {
    umap_columns = str_c(umap_key, c(1, 2))
    
    res <- res %>%
      AddMetaData(
        metadata = Embeddings(., reduction = umap_name),
        col.name = umap_columns
      )
  }
  
  res
}

#' Run PCA, cluster, run UMAP, cluster cells
#' 
#' @param sobj_in Seurat object.
#' @param assay Name of assay in object.
#' @param resolution Resolution to use for clustering.
#' @param dims Dimensions to use for clustering.
#' @param prefix Prefix to add to reduction keys and meta.data columns.
#' @param pca_meta Should PC-1 and PC-2 coordinates be added to meta.data table.
#' @param umap_meta Should UMAP coordinates be added to meta.data table.
#' @param rerun_pca Re-run PCA and UMAP even if there is already a reduction
#' with the same name.
#' @param ... Additional arguments to pass to RunPCA.
#' @return Seurat object
#' @export
cluster_RNA <- function(sobj_in, assay = "RNA", resolution = 0.6, dims = 1:40, prefix = "",
                        pca_meta = TRUE, umap_meta = TRUE, rerun_pca = TRUE, ...) {
  # Use FindNeighbors to construct a K-nearest neighbors graph based on the euclidean distance in 
  # PCA space, and refine the edge weights between any two cells based on the
  # shared overlap in their local neighborhoods (Jaccard similarity).
  # Use FindClusters to apply modularity optimization techniques such as the Louvain algorithm 
  # (default) or SLM, to iteratively group cells together
  
  res <- sobj_in
  
  # Run PCA and UMAP
  # Data must be scaled
  umap_name <- str_c(prefix, "umap")
  
  if (!umap_name %in% names(sobj_in@reductions) || rerun_pca) {
    res <- res %>%
      run_UMAP_RNA(
        assay     = assay,
        prefix    = prefix,
        dims      = dims,
        pca_meta  = pca_meta,
        umap_meta = umap_meta,
        ...
      )
  }
  
  # Create nearest neighbors graph and find clusters
  pca_name <- str_c(prefix, "pca")
  
  res <- res %>%
    FindNeighbors(
      assay     = assay,
      reduction = pca_name,
      dims      = dims
    ) %>%
    FindClusters(
      resolution = resolution,
      verbose    = FALSE
    ) %>%
    AddMetaData(
      metadata = Idents(.),
      col.name = str_c(assay, "_clusters")
    )
  
  res
}

#' Calculate feature fold change using median signal from control group
#' 
#' @param sobj_in Seurat object.
#' @param feat Name of feature to use for calculating fold changes.
#' @param data_slot Slot to pull data from.
#' @param add_pseudo Should a pseudo count be added to results? This is useful
#' if results will be log transformed.
#' @param fc_column Name of new column to add to meta.data table.
#' @param grp_column Name of meta.data column containing cell groups.
#' @param control_grps Name of control groups to use for calculating fold changes.
#' @return Seurat object
#' @export
calc_feat_fc <- function(sobj_in, feat = "adt_ovalbumin", data_slot = "counts", add_pseudo = FALSE,
                         fc_column = "ova_fc", grp_column = "cell_type", control_grps = c("B cell", "T cell")) {
  
  # Add data for calculations to meta.data
  res <- sobj_in %>%
    AddMetaData(FetchData(
      object = ., 
      vars   = feat,
      slot   = data_slot
    ))
  
  # Calculate fold changes
  res@meta.data <- res@meta.data %>%
    rownames_to_column("cell_id") %>%
    
    mutate(con_grp = if_else(                                               # Identify control groups
      !!sym(grp_column) %in% control_grps,
      TRUE,
      FALSE
    )) %>%
    
    group_by(con_grp) %>%                                                   # Calculate median counts for control
    mutate(con_med = ifelse(
      con_grp,
      median(!!sym(feat)),
      NA
    )) %>%
    ungroup() %>%
    
    mutate(
      con_med           = replace_na(con_med, max(con_med, na.rm = TRUE)),  # Fill in NAs with median control counts
      !!sym(fc_column) := !!sym(feat) / con_med
    ) %>%
    
    dplyr::select(-con_grp, -con_med) %>%
    column_to_rownames("cell_id")
  
  # Add pseudo count so fold changes can be log transformed
  # pseudo count is smallest non-zero value divided by 2
  if (0 %in% pull(res@meta.data, fc_column) && add_pseudo) {
    res <- res %>%
      mutate_metadata(
        pseudo            = ifelse(!!sym(fc_column) > 0, !!sym(fc_column), NA), 
        pseudo            = min(pseudo, na.rm = TRUE) * 0.5,
        !!sym(fc_column) := !!sym(fc_column) + pseudo
      )
  }
  
  res
}

#' Subset Seurat objects based on cell type
#' 
#' @param sobj_in Seurat object.
#' @param cell_types Cell types to use for subsetting object.
#' @param type_column meta.data column containing cell types.
#' @param dims Dimensions to use for UMAP.
#' @param cc_scoring Score cell cycle genes.
#' @param regress_vars Variables to regress out when scaling data.
#' @return Seurat object
#' @export
subset_sobj <- function(sobj_in, cell_types, type_column = "cell_type", dims = 1:40,
                        cc_scoring = FALSE, regress_vars = NULL) {
  
  # Filter cells based on input cell type
  res <- sobj_in %>%
    subset(subset = !!sym(type_column) %in% cell_types)
  
  # Score cell cycle genes
  # Use str_to_title for mouse genes
  if (cc_scoring) {
    s.genes <- cc.genes$s.genes %>%
      str_to_title()
    
    g2m.genes <- cc.genes$g2m.genes %>%
      str_to_title()
    
    res <- res %>%
      Seurat::CellCycleScoring(
        s.features   = s.genes,
        g2m.features = g2m.genes
      )
  }
  
  # Re-process object
  res <- res %>%
    Seurat::FindVariableFeatures(
      selection.method = "vst",
      nfeatures        = 2000
    ) %>%
    Seurat::ScaleData(vars.to.regress = regress_vars) %>%
    Seurat::RunPCA() %>%
    Seurat::RunUMAP(dims = dims)
  
  res
}

#' Fit gaussian mixture model
#' 
#' @param sobj_in A Seurat object or data.frame with cell IDs as row names.
#' @param gmm_data Column containing data to use for fitting model.
#' @param data_slot Slot to pull data from.
#' @param gmm_grps Names to use for GMM groups.
#' @param prob Probability cutoff to use for classifying cells. If the
#' probability is >= prob, cell will be classified in the high signal group.
#' @param prefix Prefix to use for results table containing containing GMM groups.
#' @param quiet Suppress output messages.
#' @return List containing results
#' @export
fit_gmm <- function(sobj_in, gmm_data, data_slot = "counts", gmm_grps = c("Low", "High"),
                    prob = 0.5, prefix = "", quiet = FALSE) {
  
  set.seed(42)
  
  data_df <- sobj_in
  
  if ("Seurat" %in% class(sobj_in)) {
    data_df <- sobj_in %>%
      Seurat::FetchData(gmm_data, slot = data_slot)
  }
  
  # Fit GMM for ova signal
  quiet_EM <- quietly(~ normalmixEM(.))
  
  if (!quiet) {
    quiet_EM <- normalmixEM
  }
  
  mixmdl <- data_df %>%
    pull(gmm_data) %>%
    quiet_EM()
  
  if (quiet) {
    mixmdl <- mixmdl$result
  }
  
  # New column names
  comp_names <- c("comp.1", "comp.2")
  
  if (mixmdl$mu[1] > mixmdl$mu[2]) {
    gmm_grps <- rev(gmm_grps)
  }
  
  names(comp_names)    <- gmm_grps
  names(mixmdl$mu)     <- gmm_grps
  names(mixmdl$sigma)  <- gmm_grps
  names(mixmdl$lambda) <- gmm_grps
  
  # Divide into gmm groups
  gmm_col <- str_c(prefix, "GMM_grp")
  
  res <- data.frame(
    cell_id = rownames(data_df),
    data    = data_df[, gmm_data],
    mixmdl$posterior
  ) %>%
    dplyr::rename(
      !!sym(gmm_data) := data,
      all_of(comp_names)
    ) %>%
    # rename(all_of(comp_names)) %>%
    mutate(!!sym(gmm_col) := if_else(
      !!sym(gmm_grps[2]) >= prob,
      gmm_grps[2],
      gmm_grps[1]
    )) %>%
    column_to_rownames("cell_id")
  
  # Check that GMM results match input data
  data_chk <- identical(res[, gmm_data], data_df[, gmm_data])
  cell_chk <- identical(rownames(res), rownames(data_df))
  
  if (!data_chk && cell_chk) {
    stop("Input cells do not match cells in GMM output.")
  }
  
  res <- list(
    res     = res,
    mu      = mixmdl$mu,
    sigma   = mixmdl$sigma,
    lambda  = mixmdl$lambda,
    mix_obj = mixmdl
  )
  
  res
}

#' Classify cells using GMM
#' 
#' @param sobj_in Seurat object.
#' @param gmm_data Variable present in Seurat object to use for classifying cells.
#' @param grp_column meta.data column containing cell labels to use for
#' filtering object prior to fitting GMM.
#' @param filt Cell label to use for filtering object.
#' @param data_slot Slot to pull data from.
#' @param prefix Prefix to add to meta.data columns.
#' @param return_sobj Return a Seurat object. If set to FALSE, a tibble will be
#' returned.
#' @param ... Additional arguments to pass to fit_GMM.
#' @return Seurat object or tibble containing results
#' @export
run_gmm <- function(sobj_in, gmm_data, grp_column = "cell_type", filt = NULL,
                    data_slot = "counts", prefix = "", return_sobj = TRUE) {
  
  # Filter Seurat object
  sobj_filt <- sobj_in
  
  if (!is.null(filt)) {
    sobj_filt <- sobj_filt %>%
      subset(!!sym(grp_column) == filt)
  }
  
  # Fit GMM
  # Split meta.data by grp_column
  # Split meta.data by grp_column
  gmm_df <- sobj_filt %>%
    Seurat::SplitObject(grp_column)
  
  # Fit GMM for each group and bind data.frames
  gmm_grps <- c("Low", "High")
  gmm_col  <- str_c(prefix, "GMM_grp")
  mu_col   <- str_c(prefix, "GMM_mu")
  sig_col  <- str_c(prefix, "GMM_sigma")
  lam_col  <- str_c(prefix, "GMM_lambda")
  
  gmm_df <- gmm_df %>%
    imap_dfr(~ {
      res <- .x %>%
        fit_gmm(
          gmm_data  = gmm_data,
          data_slot = data_slot,
          gmm_grps  = gmm_grps,
          prefix    = prefix,
          quiet     = TRUE
        )
      
      res$res %>%
        select(-all_of(gmm_grps)) %>%
        rownames_to_column() %>%
        mutate(
          !!sym(grp_column) := .y,
          !!sym(mu_col)     := res$mu[!!sym(gmm_col)],
          !!sym(sig_col)    := res$sigma[!!sym(gmm_col)],
          !!sym(lam_col)    := res$lambda[!!sym(gmm_col)]
        ) %>%
        column_to_rownames()
    })
  
  # Return data.frame
  if (!return_sobj) {
    return(gmm_df)
  }
  
  # Add ova groups to meta.data for input object
  gmm_df <- gmm_df %>%
    select(-!!sym(gmm_data), -!!sym(grp_column))
  
  res <- sobj_in %>%
    AddMetaData(gmm_df)
  
  # Check that rownames are the same for input and output
  if (!identical(rownames(sobj_in@meta.data), rownames(res@meta.data))) {
    stop("Input and output cell IDs do not match.")
  }
  
  # Add "other" label for cells not included in the comparison
  res <- res %>%
    mutate_metadata(!!sym(gmm_col) := replace_na(!!sym(gmm_col), "Other"))
  
  res
}

#' Export counts and meta.data tables
#' 
#' @param sobj_in Seurat object
#' @param out_dir Output directory
#' @param file_prefix Prefix to add to output files
#' @return Counts and meta.data tables
#' @export
export_matrices <- function(sobj_in, out_dir, file_prefix) {
  
  # Format count matrices
  rna <- sobj_in %>%
    GetAssayData(assay = "RNA", "counts") %>%
    as_tibble(rownames = "gene_symbol")
  
  adt <- sobj_in %>%
    GetAssayData(assay = "ADT", "counts") %>%
    as_tibble(rownames = "gene_symbol") %>%
    mutate(gene_symbol = str_c("adt_", gene_symbol))
  
  # Write count matrix
  counts_out <- file.path(out_dir, str_c(file_prefix, "_count_matrix.tsv.gz"))
  
  bind_rows(rna, adt) %>%
    write_tsv(counts_out)
  
  # meta.data columns to include  
  meta_cols <- c(
    "cell_id",    "orig.ident",
    "nCount_RNA", "nFeature_RNA",
    "nCount_ADT", "Percent_mito",
    "S.Score",    "G2M.Score", 
    "Phase",      "subtype",
    "GMM_grp",    "type_GMM_grp_2",
    "UMAP_1",     "UMAP_2"
  )
  
  # Write meta.data table
  meta_out <- file.path(out_dir, str_c(file_prefix, "_metadata.tsv.gz"))
  
  sobj_in@meta.data %>%
    as_tibble(rownames = "cell_id") %>%
    select(all_of(meta_cols)) %>%
    write_tsv(meta_out)
}


# Processing ----

#' Create Seurat object, normalize data, run UMAP, cluster
#' 
#' @param path_in Path to matrix generated by Cell Ranger.
#' @param resolution Clustering resolution.
#' @param cc_scoring Score cell cycle genes.
#' @param regress_vars Variables to regress out when scaling data.
#' @param proj_name Project name to include in meta.data table.
#' @return Seurat object
#' @export
create_sobjs_01 <- function(path_in, resolution, cc_scoring = FALSE, regress_vars = NULL,
                            proj_name = "SeuratProject") {
  
  # Create Seurat object
  res <- path_in %>%
    create_sobj(
      proj_name = proj_name,
      adt_count_min = 0
    )
  
  # Normalize and cluster
  res <- res %>%
    norm_sobj(
      rna_assay    = "RNA",
      adt_assay    = "ADT",
      cc_scoring   = cc_scoring, 
      regress_vars = regress_vars
    ) %>%
    cluster_RNA(
      assay      = "RNA",
      resolution = resolution,
      pca_meta   = FALSE,
      umap_meta  = FALSE 
    )
  
  res
}

#' Annotate cell types and calculate ova fold change for each object
#' 
#' @param sobj_in Seurat object.
#' @param ref_mat Reference matrix to use for annotating cell types.
#' @param threshold Threshold to use for clustifyr.
#' @param umap_meta Add UMAP coordinates to object.
#' @return Seurat object
#' @export
clustify_cell_types_02 <- function(sobj_in, ref_mat, threshold, umap_meta = FALSE) {
  
  # Assign cell types
  # By default variable features are used for annotation
  res <- sobj_in %>%
    clustify(
      cluster_col   = "RNA_clusters",
      ref_mat       = ref_mat,
      rename_prefix = "t1",
      seurat_out    = TRUE,
      threshold     = threshold
    )
  
  res <- res %>%
    mutate_metadata(cell_type = t1_type)
  
  if (!umap_meta) {
    res <- res %>%
      mutate_metadata(.fun = select, -UMAP_1, -UMAP_2)
  }
  
  # Calculate ova fold change
  res <- res %>%
    calc_feat_fc(
      feat         = "adt_ovalbumin",
      data_slot    = "counts",
      grp_column   = "cell_type",
      control_grps = c("B cell", "T cell"),
      add_pseudo   = TRUE
    )
  
  res
}

#' Subset objects based on cell type, re-cluster and run clustify to annotate
#' subtypes
#' 
#' @param sobj_in Seurat object.
#' @param type_in Cell type to use for subsetting object.
#' @param ref_mat Reference matrix to use for annotating cell types.
#' @param assay Name of assay in object.
#' @param resolution Resolution to use for clustering cells.
#' @param dims Dimensions to use for clustering cells.
#' @param type_column meta.data column containing cell types.
#' @param subtype_column meta.data column containing cell subtypes to generate
#' type/cluster labels.
#' @param prefix Prefix to add to meta.data columns.
#' @param cc_scoring Score cell cycle genes.
#' @param regress_vars Variables to regress out when scaling data.
#' @param ... Additional arguments to pass to clustifyr.
#' @return Seurat object
#' @export
clustify_subsets_03 <- function(sobj_in, type_in = NULL, ref_mat, assay = "RNA", resolution = 1.8, dims = 1:40,
                                type_column = "cell_type", subtype_column = "subtype", prefix = NULL, 
                                cc_scoring = FALSE, regress_vars = NULL, ...) {
  
  # Subset object, scale data, run PCA, and run UMAP
  if (!is.null(type_in)) {
    sobj_in <- sobj_in %>%
      subset_sobj(
        cell_types   = type_in,
        type_column  = type_column,
        cc_scoring   = cc_scoring,
        regress_vars = regress_vars
      )
  }
  
  # Re-cluster data
  res <- sobj_in %>%
    FindNeighbors(
      assay     = assay,
      reduction = "pca",
      dims      = dims
    ) %>%
    FindClusters(
      resolution = resolution,
      verbose    = FALSE
    ) %>%
    AddMetaData(
      metadata = Idents(.),
      col.name = "type_clusters"
    ) %>%
    clustify(
      ref_mat       = ref_mat,
      cluster_col   = "type_clusters",
      rename_prefix = prefix,
      ...
    )
  
  # Set subtype and subtype_cluster columns
  if (!is.null(subtype_column)) {
    new_col <- "type"
    
    if (!is.null(prefix)) {
      new_col <- str_c(prefix, "_", new_col)
    }
    
    sub_col   <- sym(subtype_column)
    clust_col <- sym(str_c(subtype_column, "_clusters"))
    new_col   <- sym(new_col)
    
    res <- res %>%
      mutate_metadata(
        !!sub_col   := str_to_title_v2(!!new_col),
        !!clust_col := str_c(type_clusters, "-", !!sub_col)
      )
  }
  
  res
}

#' Add subtype assignments back to main objects and split again to include
#' additional cell types (B cells, T cells, etc.) for plotting
#' 
#' @param sobj_in Seurat object.
#' @param subtype_so Seurat object containing cell subtype classifications.
#' @param type_in Cell type to subset object.
#' @param inc_types Additional cell types to include in object.
#' @param cc_scoring Score cell cycle genes.
#' @param regress_vars Variables to regress out when scaling data.
#' @return Seurat object
#' @export
resplit_objects_04 <- function(sobj_in, subtype_so, type_in, inc_types, cc_scoring = FALSE,
                               regress_vars = NULL) {
  
  # Retrieve subtype labels from subtype_so
  type_res <- subtype_so %>%
    FetchData(c(
      "t2_type", "t2_r",
      "subtype", "subtype_clusters"
    ))
  
  # Add subtype labels and subset to add back additional
  # cell types (B cells, T cells, etc)
  res <- sobj_in %>%
    AddMetaData(type_res) %>%
    subset_sobj(
      cell_types   = c(type_in, inc_types),
      type_column  = "cell_type",
      cc_scoring   = cc_scoring,
      regress_vars = regress_vars
    )
  
  # Format meta.data
  res <- res %>%
    AddMetaData(FetchData(., c("UMAP_1", "UMAP_2")))
  
  res <- res %>%
    mutate_metadata(
      subtype          = ifelse(is.na(subtype), cell_type, subtype),
      subtype          = str_to_title_v2(subtype),
      subtype_clusters = ifelse(is.na(subtype_clusters), cell_type, subtype_clusters)
    )
  
  res
  
}

#' Classify cells based on ova counts, do this for each cell type and each
#' subtype
#' 
#' @param sobj_in Seurat object.
#' @param type_in Cell type to use for classifying cells.
#' @param type_column meta.data column containing cell types.
#' @param subtype_column meta.data column containing cell subtypes.
#' @param gmm_data Data to use for classifying cells.
#' @param data_slot Slot to pull data from.
#' @return Seurat object
#' @export
classify_ova_05 <- function(sobj_in, type_in, type_column = "cell_type", subtype_column = "subtype",
                            gmm_data = "adt_ovalbumin", data_slot = "counts") {
  
  # Classify cells for cell type based on ova signal
  res <- sobj_in %>%
    run_gmm(
      grp_column  = type_column,
      filt        = type_in,
      gmm_data    = gmm_data,
      data_slot   = data_slot,
      return_sobj = TRUE
    ) %>%
    mutate_metadata(
      GMM_grp = if_else(GMM_grp != "Other", str_to_lower(GMM_grp), GMM_grp),
      GMM_grp = if_else(GMM_grp != "Other", str_c("ova ", GMM_grp), GMM_grp ),
    )
  
  # Classify cells separately for each subtype
  GMM_res <- res %>%
    run_gmm(
      grp_column  = subtype_column,
      gmm_data    = gmm_data,
      data_slot   = data_slot,
      return_sobj = FALSE
    ) %>%
    rownames_to_column("cell_id")
  
  # Format subtype GMM results
  GMM_res <- GMM_res %>%
    mutate(
      GMM_grp        = str_c("ova ", str_to_lower(GMM_grp)),
      type_GMM_grp_2 = str_c(!!sym(subtype_column), "-", GMM_grp)
    ) %>%
    select(
      cell_id,
      type_GMM_grp    = GMM_grp,
      type_GMM_grp_2,
      type_GMM_mu     = GMM_mu,
      type_GMM_sigma  = GMM_sigma,
      type_GMM_lambda = GMM_lambda
    ) %>%
    column_to_rownames("cell_id")
  
  # Add subtype GMM results back to object
  res <- res %>%
    AddMetaData(GMM_res)
  
  res
}

