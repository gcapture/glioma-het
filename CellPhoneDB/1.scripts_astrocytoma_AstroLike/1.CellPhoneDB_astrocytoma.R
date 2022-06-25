rm(list = ls())

library(Seurat)
library(stringr)
Clusters  <- readRDS("Clusters_microglia_primary_annotated.rds")
#Clusters <- subset(Clusters, idents = "Mg-IFNg TAMs")
Idents(Clusters) <- "sample"
Clusters <- subset(Clusters, idents = "primary_astro")
Command(Clusters)
Clusters@active.assay <- "RNA"
Clusters <- NormalizeData(
  Clusters,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)


seurat <- readRDS("OE0145-IDH_integrated_astrocytoma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
Command(seurat)
label <- readRDS("labels_astro_primary.rds")
seurat$Subclusters <- label
Idents(seurat) <- "Subclusters"

table(seurat$Subclusters)
seurat <- subset(seurat, idents = "Astro-like")
seurat@active.assay <- "RNA"
seurat <- NormalizeData(
  seurat,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

seurat@assays$SCT <- NULL
seurat@assays$integrated <- NULL
seurat@reductions$harmony = NULL
seurat@reductions$pca = NULL
seurat@reductions$tsne = NULL
seurat@reductions$umap = NULL

#cpdb <- merge(seurat, y = Clusters, add.cell.ids = c("Astro", "IFN"), merge.data = TRUE)
cpdb <- merge(seurat, y = Clusters, add.cell.ids = c("Astro", "Microglia"), merge.data = TRUE)

#data corresponds to log normalize counts
#counts corresponds to raw data
counts = GetAssayData(object = cpdb, slot = "data")
colnames(counts) =  str_replace_all(colnames(counts),pattern = "-", replacement = "_")
write.table(counts,
            "counts.txt",
            sep="\t",
            quote=FALSE,
            row.names = TRUE)

cpdb[["Cells"]] = Cells(cpdb)
metadata <- as.data.frame(cpdb@meta.data[,c("Cells","Subclusters")]) 
colnames(metadata) = c("Cell","cell_type")

metadata$cell_type <- str_replace_all(metadata$cell_type,
                                       pattern = " ",
                                       replacement = "_")
metadata$cell_type <- str_replace_all(metadata$cell_type,
                                      pattern = "-",
                                      replacement = "_")
metadata$Cell <- str_replace_all(metadata$Cell,
                                      pattern = "-",
                                      replacement = "_")


# write metadata table
write.table(metadata,
            "metadata.txt",
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

saveRDS(cpdb, file = "astrocytoma.cpdb.rds")

