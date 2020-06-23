migDC_data <- read_csv("https://sharehost.hms.harvard.edu/immgen/GSE15907/GSE15907_Normalized_Data.csv")
migDC <- migDC_data %>%
  dplyr::select(ProbeSetID, GeneSymbol, "DC.IIhilang-103-11b+.SLN#1", "DC.IIhilang-103-11b+.SLN#2",
                "DC.IIhilang-103-11b+.SLN#3","DC.IIhilang+103+11blo.SLN#1", "DC.IIhilang+103+11blo.SLN#2",
                "DC.IIhilang+103+11blo.SLN#3", "DC.IIhilang+103-11b+.SLN#1","DC.IIhilang+103-11b+.SLN#2",
                "DC.IIhilang+103-11b+.SLN#3")

#sub-classify DC cells
d2_DC_so <- subset(d2_so, subset = cell_type1 == "DC")
d2_DC_all.genes <- rownames(d2_DC_so)
d2_DC_so <- ScaleData(d2_DC_so, features = d2_DC_all.genes)
d2_DC_so <- FindVariableFeatures(object = d2_DC_so)
d2_DC_so <- RunPCA(d2_DC_so, features = VariableFeatures(object = d2_DC_so))
ElbowPlot(d2_DC_so)

d2_DC_so <- FindNeighbors(d2_DC_so, dims = 1:30)
d2_DC_so <- FindClusters(d2_DC_so, resolution = 1.5)
d2_DC_so <- RunUMAP(d2_DC_so, dims = 1:30)
plot <- DimPlot(d2_DC_so, reduction = "umap", group.by = "seurat_clusters", label = T)
saveRDS(d2_DC_so, file = paste(so_folder, "d2_DC_so.rds", sep = ""))

d2DC_resDC <- clustify(
  input = d2_DC_so,
  cluster_col = "seurat_clusters",
  ref_mat = ref_lcDC,
  seurat_out = T, threshold = 0.5, verbose = T)
d2DC_resDC_type <- FetchData(d2DC_resDC, vars = c("seurat_clusters", "type", "r")) %>%
  rownames_to_column(var = "cell_id") %>%
  dplyr::select(-cell_id) %>%
  unique()
d2DC_resDC_type

d2DC_resDC_mig <- clustify(
  input = d2_DC_so,
  cluster_col = "seurat_clusters",
  ref_mat = migDC,
  seurat_out = F, threshold = 0.5, verbose = T)
d2DC_resDC_type <- FetchData(d2DC_resDC, vars = c("seurat_clusters", "type", "r")) %>%
  rownames_to_column(var = "cell_id") %>%
  dplyr::select(-cell_id) %>%
  unique()
d2DC_resDC_type

d2DC.new.cluster.ids <- c("CCR7hi cDC2", "cDC2 Tbet-", "cDC2 Tbet-", "CCR7hi cDC2", "cDC2 Tbet-", "CCR7hi cDC2", "CCR7hi cDC2", "CCR7hi cDC2", "cDC2 Tbet-", "cDC1", "cDC2 Tbet-", "CCR7hi cDC2", "CCR7hi cDC2", "cDC2 Tbet-", "cDC1", "cDC2 Tbet-", "cDC1", "Monocyte")
d2_DC_so@meta.data$cell_typeDC <- d2DC.new.cluster.ids[d2_DC_so$seurat_clusters]
plot <- DimPlot(d2_DC_so, reduction = "umap", group.by = "cell_typeDC")

d2DC.new.cluster.ids2 <- c("0-CCR7hi cDC2", "1-cDC2 Tbet-", "2-cDC2 Tbet-", "3-CCR7hi cDC2", "4-cDC2 Tbet-", "5-CCR7hi cDC2", "6-CCR7hi cDC2", "7-CCR7hi cDC2", "8-cDC2 Tbet-", "9-cDC1", "10-cDC2 Tbet-", "11-CCR7hi cDC2", "12-CCR7hi cDC2", "13-cDC2 Tbet-", "14-cDC1", "15-cDC2 Tbet-", "16-cDC1", "17-Monocyte")
d2_DC_so@meta.data$cell_typeDC2 <- d2DC.new.cluster.ids2[d2_DC_so$seurat_clusters]

FeaturePlot(d2_DC_so, c("Itgax", "Itgam", "Xcr1", "Ccr7", "Tbx21", "Id2"))
plot <- VlnPlot(d2_DC_so, features = c("Itgax", "Ccr7", "Xcr1", "Irf8", "Csf1r", "Fcgr3", "Cd209a", "Itgam", "Tbx21", "Siglech", "Tcf4", "Bst2"), group.by = "cell_typeDC2", ncol = 1)

FetchData(d2_DC_so, c("cell_typeDC", "ovalbumin"), slot = "data") %>%
  rownames_to_column(var = "cell_id") %>%
  plot_activity(adt_ovalbumin, group = cell_typeDC, vertical = F) +
  scale_color_OkabeIto() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "CLR Normalized", y = NULL) +
  scale_x_log10()

FetchData(d2_DC_so, c("seurat_clusters", "ovalbumin"), slot = "data") %>%
  rownames_to_column(var = "cell_id") %>%
  plot_activity(adt_ovalbumin, group = seurat_clusters, vertical = F) +
  #scale_color_OkabeIto() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "counts", y = NULL) +
  scale_x_log10()
#ggsave(paste(plots_folder, "d2/DC_clusters_d2.pdf", sep = ""), plot, width = 6, height = 4, units = c("in"), useDingbats = F)