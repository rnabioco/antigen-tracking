library(clustifyr)
library(tidyverse)
library(usethis)

#migratory DC subtypes ImmGen
migDC_data <- read_csv("https://sharehost.hms.harvard.edu/immgen/GSE15907/GSE15907_Normalized_Data.csv")
migDC <- migDC %>%
  rename(migcDC2_A1 = "DC.IIhilang-103-11b+.SLN#1", migcDC2_A2 = "DC.IIhilang-103-11b+.SLN#2", migcDC2_A3 =
           "DC.IIhilang-103-11b+.SLN#3", migcDC1_1 = "DC.IIhilang+103+11blo.SLN#1", migcDC1_2 =
           "DC.IIhilang+103+11blo.SLN#2", migcDC1_3 = "DC.IIhilang+103+11blo.SLN#3", lang_1 =
           "DC.IIhilang+103-11b+.SLN#1", lang_2 ="DC.IIhilang+103-11b+.SLN#2", lang_3=
           "DC.IIhilang+103-11b+.SLN#3", migcDC2_B1 = "DC.IIhilang-103-11blo.SLN#1", migcDC2_B2 =  
           "DC.IIhilang-103-11blo.SLN#2", migcDC2_B3 = "DC.IIhilang-103-11blo.SLN#3") %>%
  gather(key = "A", value = "migcDC2A", "migcDC2_A1", "migcDC2_A2", "migcDC2_A3") %>%
  gather(key = "B", value = "migcDC1", "migcDC1_1", "migcDC1_2", "migcDC1_3") %>%
  gather(key = "C", value = "lang", "lang_1", "lang_2", "lang_3") %>%
  gather(key = "D", value = "migcDC2B", "migcDC2_B1", "migcDC2_B2", "migcDC2_B3") %>%
  dplyr::select(-A, -B, -C ,-D) %>%
  group_by(GeneSymbol) %>%
  summarise_all(mean)
ref_migDC <- migDC %>%
  dplyr::select(-GeneSymbol) %>%
  as.matrix()
rownames(ref_migDC) <- migDC$GeneSymbol

#sub-classify DC cells
d14_DC_so <- subset(d14_so, subset = cell_type1 == "DC")
d14_DC_all.genes <- rownames(d14_DC_so)
d14_DC_so <- ScaleData(d14_DC_so, features = d14_DC_all.genes)
d14_DC_so <- FindVariableFeatures(object = d14_DC_so)
d14_DC_so <- RunPCA(d14_DC_so, features = VariableFeatures(object = d14_DC_so))
ElbowPlot(d14_DC_so)

d14_DC_so <- FindNeighbors(d14_DC_so, dims = 1:30)
d14_DC_so <- FindClusters(d14_DC_so, resolution = 1.5)
d14_DC_so <- RunUMAP(d14_DC_so, dims = 1:30)
plot <- DimPlot(d14_DC_so, reduction = "umap", group.by = "seurat_clusters", label = T)

d14DC_resDC <- clustify(
  input = d14_DC_so,
  cluster_col = "seurat_clusters",
  ref_mat = ref_lcDC,
  seurat_out = T, threshold = 0.5, verbose = T)
d14DC_resDC_type <- FetchData(d14DC_resDC, vars = c("seurat_clusters", "type", "r")) %>%
  rownames_to_column(var = "cell_id") %>%
  dplyr::select(-cell_id) %>%
  unique()
d14DC_resDC_type

d14DC_resDC_mig <- clustify(
  input = d14_DC_so,
  cluster_col = "seurat_clusters",
  ref_mat = ref_migDC,
  seurat_out = T, threshold = 0.5, verbose = T)
d14DC_resDC_mig_type <- FetchData(d14DC_resDC_mig, vars = c("seurat_clusters", "type", "r")) %>%
  rownames_to_column(var = "cell_id") %>%
  dplyr::select(-cell_id) %>%
  unique()
d14DC_resDC_mig_type
d14DC_resDC_mig

d14DC.new.cluster.ids <- c("cDC2 Tbet-", "cDC2 Tbet-", "cDC1", "CCR7hi", "cDC2 Tbet-", "Monocyte", "CCR7hi", "cDC2 Tbet-", "cDC1", "cDC2 Tbet+", "CCR7hi", "CCR7hi", "cDC1", "Siglec-H", "cDC2 Tbet-", "Monocyte", "Monocyte", "cDC2 Tbet-", "Monocyte", "cDC2 Tbet-")
d14_DC_so@meta.data$cell_typeDC <- d14DC.new.cluster.ids[d14_DC_so$seurat_clusters]
plot <- DimPlot(d14_DC_so, reduction = "umap", group.by = "cell_typeDC")

d14DC.new.cluster.ids2 <- c("0-cDC2 Tbet-", "1-cDC2 Tbet-", "2-cDC1", "3-CCR7hi", "4-cDC2 Tbet-", "5-Monocyte", "6-CCR7hi", "7-cDC2 Tbet-", "8-cDC1", "9-cDC2 Tbet+", "10-CCR7hi", "11-CCR7hi", "12-cDC1", "13-Siglec-H", "14-cDC2 Tbet-", "15-Monocyte", "16-Monocyte", "17-cDC2 Tbet-", "18-Monocyte", "19-cDC2 Tbet-")
d14_DC_so@meta.data$cell_typeDC2 <- d14DC.new.cluster.ids2[d14_DC_so$seurat_clusters]

FeaturePlot(d14_DC_so, c("Itgax", "Itgam", "Xcr1", "Ccr7", "Tbx21", "Id2"))
plot <- VlnPlot(d14_DC_so, features = c("Itgax", "Ccr7", "Xcr1", "Irf8", "Csf1r", "Fcgr3", "Cd209a", "Itgam", "Tbx21", "Siglech", "Tcf4", "Bst2"), group.by = "cell_typeDC2", ncol = 1)

FeaturePlot(d14_DC_so, c("Cd207", "Itgae", "Itgam", "Epcam"))
plot <- VlnPlot(d14_DC_so, features = c("Cd207", "Itgae", "Itgam", "Epcam"), group.by = "cell_typeDC2", ncol = 1)

FetchData(d14_DC_so, c("seurat_clusters", "ovalbumin"), slot = "data") %>%
  rownames_to_column(var = "cell_id") %>%
  plot_activity(adt_ovalbumin, group = seurat_clusters, vertical = F) +
  #scale_color_OkabeIto() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "counts", y = NULL) +
  scale_x_log10()
#ggsave(paste(plots_folder, "d14/migDC_gex_d14.pdf", sep = ""), plot, width = 6, height = 16, units = c("in"), useDingbats = F)
#saveRDS(d14_DC_so, file = paste(so_folder, "d14_DC_so.rds", sep = ""))
#d14_so <- readRDS(file = paste(so_folder, "d14_so.rds", sep = ""))

