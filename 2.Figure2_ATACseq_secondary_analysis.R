## Inma Hernandez Lopez and Paula Soler Vila
## inmaculada.hernandez-lopez@ukr.de
## LIT -  Leibniz-Institute für Immunotherapie


## Tutorial: https://satijalab.org/signac/articles/pbmc_vignette.html#integrating-with-scrna-seq-data-1

################################################################################
## Input Robjects:

  #Filtered ATAC objects: 
  #OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar
  #OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar

  #Filtered and annotated RNA objects for oligodendroglioma:
  #OE0145-IDH_integrated_oligodendroglioma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds
  #labels_oligo_primary.rds"

  #Filtered and annotated RNA objects for astrocytoma:
  #
  #

## Output Robjects:


## Steps:
    #A. Oligodendroglioma
      #A1. Integrating with scRNA-seq data: label transfer
      #A2. Motif Analysis
      #A3. Comparison Accessibility/Expression. Code from Paula Soler
      


    #B. Astrocytoma
      #B1. Re-RUN ATAC integrated object astrocytes with less dimensions to avoid 
#microglia split in astrocytoma
      #B2. Integrating with scRNA-seq data: label transfer
      #B3. Motif Analysis
      #B4. Comparison Accessibility/Expression. Code from Paula Soler
      

################################################################################

rm(list = ls())

library(Seurat)
library(Signac)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggpubr)
library(ggplot2)
library(TFBSTools)
library(patchwork)
library(chromVAR)
library(motifmatchr)
library(data.table)
library(chromVARmotifs)
library(dplyr)
library(purrr)
library(readxl)
library(openxlsx)
library(tidyr)
library(pheatmap)
library(factoextra)
library(corrplot)
library(plyr)
library(SCpubr)
library(chromVAR)
library(tidyverse)


#-------------------------------------------------------------------------------
#A. Oligodendroglioma
      #A1. Integrating with scRNA-seq data: label transfer
#-------------------------------------------------------------------------------
seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_oligo_primary.rds")


seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
seurat@assays$RNA <- NULL

Idents(seurat) <- "transferred_labels"

peaks <- readRDS("OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar")

peaks@active.assay <- "ACTIVITY"
peaks@assays$ATAC <- NULL
peaks@assays$RNA <- NULL
peaks@assays$chromvar <- NULL

transfer.anchors <- FindTransferAnchors(
  reference = seurat,
  query = peaks,
  normalization.method = "SCT",
  reduction = 'cca',
  recompute.residuals = FALSE
)

saveRDS(transfer.anchors, file = "transfer.anchors_oligodendroglioma.rds")


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat$transferred_labels,
  weight.reduction = peaks[['lsi']],
  dims = 2:30
)

peaks <- AddMetaData(object = peaks, metadata = predicted.labels)

peaks[["Cells"]] <- Cells(peaks)
annotation_ATAC_oligo <- as.vector(peaks@meta.data[,c("predicted.id")])
names(annotation_ATAC_oligo) <- peaks$Cells
saveRDS(annotation_ATAC_oligo, file = "annotation_ATAC_oligo_primary.rds")

cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  "Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  "Astrocytes" = "#279185",                  # Celadon Green.
                  "Microglia" = "#7EB356",                   # Bud Green.
                  "Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  "Neurons" = "#544B81",                     # Purple Navy.
                  "Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E",                  # Ruby Red.
                  "Undefined" = "#63412C"                   # Van Dyke Brown.
                  
)   

cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  #"Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  #"Astrocytes" = "#279185",                  # Celadon Green.
                  #"Microglia" = "#7EB356",                   # Bud Green.
                  #"Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  #"Neurons" = "#544B81",                     # Purple Navy.
                  #"Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E"                  # Ruby Red.
                  #"Undefined" = "#63412C"                   # Van Dyke Brown.
                  
)   


peaks <- readRDS("OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar.rds")
label <- readRDS("annotation_ATAC_oligo_primary.rds")
peaks$predicted.id <- label 

Idents(peaks) <- "predicted.id"
pdf("Fig2_supplement_DimPlot_ATAC_oligodendroglioma.pdf", width = 7.48, height = 8.48 )
do_DimPlot(peaks, colors.use = cluster_cols, plot.title = "ATAC Oligodendroglioma")
dev.off()

table(peaks$predicted.id)
peaks <- subset(peaks, idents = c("Astro-like", "OPC-like","Gradient","RA"))
peaks[["Cell"]] <- Cells(peaks)

emb = as.data.frame(peaks[["umap"]]@cell.embeddings)
cells = emb[(emb$UMAP_1 > (-7)) & (emb$UMAP_2 > (-10)),]
cells = rownames(cells)
Idents(peaks) <- "Cell"
peaks <- subset(peaks, idents = cells)

Idents(peaks) <- "predicted.id"
pdf("Fig2_main_DimPlot_ATAC_oligodendroglioma.pdf", width = 3.4, height = 3.39 )
do_DimPlot(peaks, colors.use = cluster_cols, plot.title = "ATAC Oligodendroglioma",
           legend.text.size = 14, legend.ncol = 2)
dev.off()



pdf("Fig2_main_DimPlot_ATAC_oligodendroglioma_split.pdf",width = 6.81, height = 3.39 )
do_DimPlot(peaks, split.by = "predicted.id", cols.split = cluster_cols,
           legend = FALSE, plot.title = "ATAC oligodendroglioma",pt.size = 1, ncol = 4)
dev.off()

#-------------------------------------------------------------------------------
#A. Oligodendroglioma
      #A2. Motif Analysis
#-------------------------------------------------------------------------------

seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar.rds")

label <- readRDS("annotation_ATAC_oligo_primary.rds")

seurat$predicted.id <- label

data("human_pwms_v1")
length(human_pwms_v1)

seurat <- AddMotifs(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v1
)

seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(seurat, "OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar_IHL.rds")

  
  path_to_object <- ""
  tumor <- "oligodendroglioma"
  
  # Paths
  path_to_obj <- paste0(
    here::here(path_to_object),"OE0145-IDH_integrated_",
    tumor,"_peaks_activity_chromvar_IHL.rds",
    sep = ""
  )
  
  color_palette <-  c("#fed439","#709ae1","#4caf50","#d2af81",
                      "#d5e4a2","#9c27b0","#f05c3b","#46732e",
                      "#71d0f5","#370335","#075149","#c80813",
                      "#91331f","#1a9993","#fd8cc1")
  
  path_to_save <- paste0(here::here(path_to_object),"OE0145-IDH_integrated_",
                         tumor,"_motifs_ChromVar_FC0_adjpvalue0.1.xls",
                         sep = "")
  
  
  seurat <- readRDS(path_to_obj)
  DefaultAssay(seurat) <- 'chromvar'
  Idents(seurat) <- seurat$predicted.id
  seurat <- subset(seurat, idents = c("Astro-like", "OPC-like","Gradient","RA"))
  
  da_regions <- FindAllMarkers(
    object = seurat,
    only.pos = TRUE,
    min.pct = 0.1,
    test.use = 'LR',
    latent.vars = 'nCount_ATAC'
  )
  
  
  da_regions_selected <- (da_regions[da_regions$p_val_adj < 0.05 & da_regions$avg_log2FC > 0, ])
  #da_regions_selected <- (da_regions[da_regions$p_val_adj < 0.25 & da_regions$avg_log2FC > 0, ])
  
  
  motif_name_entire <- da_regions_selected %>% separate(gene, 
                                                        c("Part1", "Part2","Part3",
                                                          "Part4","Part5"), sep="-")
  
  da_regions_selected$motif_name <- motif_name_entire$Part3
  
  da_regions_selected_sorted <- da_regions_selected[with(da_regions_selected,
                                                         order(cluster, -avg_log2FC)), ]
  
  da_regions_selected_sorted$rank <- ave(da_regions_selected_sorted$avg_log2FC,
                                         da_regions_selected_sorted$cluster, FUN = seq_along)
  
  da_regions_selected_sorted_prepared <- da_regions_selected_sorted %>% 
    group_by(cluster) %>% 
    top_n(20,-rank)
  
  output <- split(da_regions_selected_sorted, da_regions_selected_sorted$cluster)
  wb <- createWorkbook()
  
  for (i in 1:length(output)) {
    addWorksheet(wb, sheetName=names(output[i]))
    writeData(wb, sheet=names(output[i]), x=output[[i]]) # Note [[]]
  }
  
  saveWorkbook(wb, file=path_to_save,  overwrite = TRUE)
  
  cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                    "OPC-like" = "#043362",                    # Prussian Blue.
                    #"Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                    #"Astrocytes" = "#279185",                  # Celadon Green.
                    #"Microglia" = "#7EB356",                   # Bud Green.
                    #"Pericytes" = "#AC70FF",                   # Medium Purple.
                    "Gradient" = "#D6D6D6",                    # Light grey.
                    #"Neurons" = "#544B81",                     # Purple Navy.
                    #"Endothelial" = "#da627d",                 # Blush.
                    "Astro-like" = "#9A031E"                  # Ruby Red.
                    #"Undefined" = "#63412C"                   # Van Dyke Brown.
                    
  )   
  
  motifs = c("ENSG00000130816-LINE1722-DNMT1-D", "ENSG00000113658-LINE19373-SMAD5-I-N1",
             "ENSG00000120693-LINE19389-SMAD9-I-N1","ENSG00000101412-LINE1761-E2F1-D-N1",
             "ENSG00000113658-LINE19362-SMAD5-I-N1","ENSG00000120693-LINE19376-SMAD9-I-N1",
             "ENSG00000101076-LINE2935-HNF4A-D-N1","ENSG00000120738-LINE864-EGR1-D-N8")

  
  for (i in 1:length(motifs)) {
    pdf(paste0("Fig2J_main_OligoMotifs_VlnPlot_",motifs[i], ".pdf", sep = ""))
    print(do_VlnPlot(seurat,features = motifs[i] ,colors.use = cluster_cols, group.by = "predicted.id"))
    dev.off()
  }
  
#-------------------------------------------------------------------------------
#A. Oligodendroglioma
      #A3. Comparison Accessibility/Expression. Code from Paula Soler
      #Paula Soler paula.soler@cnag.crg.eu
#-------------------------------------------------------------------------------

seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_oligo_primary.rds")


seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
Idents(seurat) <- "transferred_labels"
seurat@assays$RNA <- NULL


Genes_DE_RA_expression <- FindMarkers(
  seurat,
  ident.1 = "RA",
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  verbose = TRUE
)

seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar.rds")
label <- readRDS("annotation_ATAC_oligo_primary.rds")

seurat$predicted.id <- label
Idents(seurat) <- "predicted.id"

da_peaks <- FindMarkers(
  object = seurat,
  ident.1 = "RA",
  min.pct = 0,
  test.use = 'wilcox',
  only.pos = TRUE
)

peaks <- rownames(da_peaks[da_peaks$avg_log2FC > 0 & da_peaks$p_val < 0.05, ])
closest_genes_oligo_RA <- ClosestFeature(seurat, regions = peaks)

rownames(closest_genes_oligo_RA) <- closest_genes_oligo_RA$query_region

gene_FC <- merge(da_peaks, closest_genes_oligo_RA, by=0)
gene_FC <- gene_FC[,c("query_region","avg_log2FC","gene_name")]


RA_expression = Genes_DE_RA_expression
RA_expression$gene_name <- rownames(RA_expression)


common_genes_classif <- merge(RA_expression, gene_FC, by = "gene_name")
mean_accessibility <- tapply(common_genes_classif$avg_log2FC.y,
                             common_genes_classif$gene, mean)
mean_expresion <- tapply(common_genes_classif$avg_log2FC.x,
                         common_genes_classif$gene, mean)

avglogFC <- cbind(melt(mean_accessibility),
                  melt(mean_expresion)[,2])
colnames(avglogFC) <- c("gene","Accessibility","Expression")
write.table(avglogFC, file = "SuppFig2_avglogFC_oligodendroglioma.txt", sep = "\t", quote = FALSE)


avglogFC.melt <- melt(avglogFC)
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("CD24", "PDGFRA", "CSPG4", "SOX11", "ETV5",
                  "PCDH15", "LINC00639", "FGF14", "ADARB2")

avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)

avglogFC.melt_interesting$gene <- factor(avglogFC.melt_interesting$gene, levels = c("CD24", "PDGFRA", "CSPG4", "SOX11", "ETV5",
                                                                                    "PCDH15", "LINC00639", "FGF14", "ADARB2"))

pdf("Fig2_main_RA_oligodendroglioma_ggdotchart.pdf", width=3.2, height=4)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 14, face = "bold")) 
dev.off()

pdf("Fig2_Suppl_RA_oligodendroglioma_ggscatter.pdf", width=3, height=4)

ggscatter(avglogFC, x = "Accessibility", y = "Expression",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Accesibility (avg_log2FC)", 
          ylab = "Expression (avg_log2FC)")
dev.off()

##dotplot to check stem like genes are only from RA

features_use <- c("PDGFRA", "CD24","ADARB2", "FGF14","LINC00639",
                  "PCDH15", "SOX11", "ETV5","CSPG4" )

seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_oligo_primary.rds")
seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
Idents(seurat) <- "transferred_labels"
seurat = subset(seurat, idents = c("Astro-like", "Cycling-like", "Gradient", "OPC-like", "RA"))

pdf("Fig2_main_RA_oligodendroglioma_DotPlot.pdf", width=4.8, height=4.5)
do_DotPlot(sample = seurat,
           features = rev(features_use),
           cols = c("grey75", "#3c5b8b"),
           ylab = "Oligodendroglioma",
           legend = TRUE,
           legend.position = "right",
           flip = T) 
dev.off()


#-------------------------------------------------------------------------------
#B. Astrocytoma
#-------------------------------------------------------------------------------
rm(list = ls())

library(Seurat)
library(Signac)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggpubr)
library(ggplot2)
library(TFBSTools)
library(patchwork)
library(chromVAR)
library(motifmatchr)
library(ggpubr)
library(data.table)
library(chromVARmotifs)
library(dplyr)
library(purrr)
library(readxl)
library(openxlsx)
library(tidyr)
library(pheatmap)
library(factoextra)
library(corrplot)
library(plyr)

#-------------------------------------------------------------------------------
#B1. Re-RUN ATAC integrated object astrocytes with less dimensions to avoid 
#microglia split in astrocytoma
#-------------------------------------------------------------------------------
peaks <- readRDS("OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar.rds")
peaks@active.assay <- "ATAC"
peaks@assays$chromvar <- NULL

peaks <- Seurat::RunUMAP(object = peaks,
                         reduction = "harmony",
                         dims = 2:10)

# Find neighbors.
print("Compute clusters.")
peaks <- Seurat::FindNeighbors(object = peaks,
                               reduction = "harmony",
                               dims = 2:10)



# Find clusters. Using algorithm = 3 (SLM algorithm) according to Signac's vignette.
peaks <- Seurat::FindClusters(object = peaks,
                              verbose = FALSE,
                              algorithm = 3)

Idents(peaks) <- "predicted.id"

p.umap <- Seurat::DimPlot(peaks,
                          reduction = "umap",
                          label = TRUE,
                          repel = TRUE,
                          label.box = TRUE,
                          label.color = "black",
                          shuffle = TRUE,
                          pt.size = 0.5) +
  ggpubr::theme_pubr() +
  ggpubr::rremove("legend.title") +
  ggpubr::rremove("legend") +
  Seurat::NoAxes()
p.umap$layers[[2]]$aes_params$fontface <- "bold"


saveRDS(peaks,"OE0145-IDH_integrated_astrocytoma_peaks_activity_newIHL.rds")

#-------------------------------------------------------------------------------
#B. Astrocytoma
#B2. Integrating with scRNA-seq data: label transfer
#-------------------------------------------------------------------------------
#Label transfer with oligodendroglioma scRNAseq samples
seurat <- readRDS("OE0145-IDH_integrated_oligodendroglioma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_oligo_primary.rds")

seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
seurat@assays$RNA <- NULL
Idents(seurat) <- "transferred_labels"

peaks <- readRDS("OE0145-IDH_integrated_astrocytoma_peaks_activity_newIHL.rds")

peaks@active.assay <- "ACTIVITY"
peaks@assays$ATAC <- NULL
peaks@assays$RNA <- NULL
peaks@assays$chromvar <- NULL

transfer.anchors <- FindTransferAnchors(
  reference = seurat,
  query = peaks,
  normalization.method = "SCT",
  reduction = 'cca',
  recompute.residuals = FALSE
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat$transferred_labels,
  weight.reduction = peaks[['lsi']],
  dims = 2:30
)


peaks <- AddMetaData(object = peaks, metadata = predicted.labels)

peaks[["Cells"]] <- Cells(peaks)
annotation_ATAC_astro <- as.vector(peaks@meta.data[,c("predicted.id")])
names(annotation_ATAC_astro) <- peaks$Cells
saveRDS(annotation_ATAC_astro, file = "annotation_ATAC_astro_primary.rds")


cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  "Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  "Astrocytes" = "#279185",                  # Celadon Green.
                  "Microglia" = "#7EB356",                   # Bud Green.
                  "Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  "Neurons" = "#544B81",                     # Purple Navy.
                  "Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E",                  # Ruby Red.
                  "Undefined" = "#63412C",                   # Van Dyke Brown.
                  "Cycling-like" = "#5F0F40")                # Tyrian Purple.
                  
 

cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  #"Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  #"Astrocytes" = "#279185",                  # Celadon Green.
                  #"Microglia" = "#7EB356",                   # Bud Green.
                  #"Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  #"Neurons" = "#544B81",                     # Purple Navy.
                  #"Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E")                 # Ruby Red.
                  #"Undefined" = "#63412C")                   # Van Dyke Brown)   
                  
 


peaks <- readRDS("OE0145-IDH_integrated_astrocytoma_peaks_activity_newIHL.rds")
label <- readRDS("annotation_ATAC_astro_primary.rds")
peaks$predicted.id <- label 


Idents(peaks) <- "predicted.id"
pdf("Fig2I_supplement_DimPlot_ATAC_astrocytoma.pdf", width = 7.48, height = 8.48 )
do_DimPlot(peaks, colors.use = cluster_cols, plot.title = "ATAC Astrocytoma")
dev.off()

table(peaks$predicted.id)

peaks <- DietSeurat(peaks, assays = "ATAC",dimreducs = c("lsi","umap","harmony"))
peaks <- subset(peaks, idents = c("Astro-like", "OPC-like","Gradient","RA"))


Idents(peaks) <- "predicted.id"
pdf("Fig2I_main_DimPlot_ATAC_astrocytoma.pdf", width = 3.4, height = 3.39 )
do_DimPlot(peaks, colors.use = cluster_cols, plot.title = "ATAC Astrocytoma",
           legend.text.size = 14, legend.ncol = 2)
dev.off()



pdf("Fig2I_main_DimPlot_ATAC_astrocytoma_split.pdf",width = 6.81, height = 3.39 )
do_DimPlot(peaks, split.by = "predicted.id", colors.split = cluster_cols,
           legend = FALSE, plot.title = "ATAC Astrocytoma",pt.size = 1, ncol = 4)
dev.off()

#-------------------------------------------------------------------------------
#B. Astrocytoma
#B3. Comparison Accessibility/Expression. Code from Paula Soler
#------------------------------------------------------------------------------

library(Seurat)
library(Signac)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggpubr)
library(ggplot2)
library(TFBSTools)
library(patchwork)
library(chromVAR)
library(motifmatchr)
library(ggpubr)
library(data.table)
library(chromVARmotifs)
library(dplyr)
library(purrr)
library(readxl)
library(openxlsx)
library(tidyr)
library(pheatmap)
library(factoextra)
library(corrplot)
library(plyr)
library(SCpubr)

seurat <- readRDS("OE0145-IDH_integrated_astrocytoma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_astro_primary.rds")

seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
Idents(seurat) <- "transferred_labels"
seurat@assays$RNA <- NULL


Genes_DE_RA_expression <- FindMarkers(
  seurat,
  ident.1 = "RA",
  test.use = "wilcox",
  only.pos = TRUE,
  logfc.threshold = 0.25,
  verbose = TRUE
)

seurat <- readRDS("OE0145-IDH_integrated_astrocytoma_peaks_activity_newIHL.rds")
label <- readRDS("annotation_ATAC_astro_primary.rds")

seurat$predicted.id <- label
Idents(seurat) <- "predicted.id"

da_peaks <- FindMarkers(
  object = seurat,
  ident.1 = "RA",
  min.pct = 0,
  test.use = 'wilcox',
  only.pos = TRUE
)

peaks <- rownames(da_peaks[da_peaks$avg_log2FC > 0 & da_peaks$p_val < 0.05, ])
closest_genes_astro_RA <- ClosestFeature(seurat, regions = peaks)

rownames(closest_genes_astro_RA) <- closest_genes_astro_RA$query_region

gene_FC <- merge(da_peaks, closest_genes_astro_RA, by=0)
gene_FC <- gene_FC[,c("query_region","avg_log2FC","gene_name")]


RA_expression = Genes_DE_RA_expression
RA_expression$gene_name <- rownames(RA_expression)


common_genes_classif <- merge(RA_expression, gene_FC, by = "gene_name")
mean_accessibility <- tapply(common_genes_classif$avg_log2FC.y,
                             common_genes_classif$gene, mean)
mean_expresion <- tapply(common_genes_classif$avg_log2FC.x,
                         common_genes_classif$gene, mean)

avglogFC <- cbind(melt(mean_accessibility),
                  melt(mean_expresion)[,2])
colnames(avglogFC) <- c("gene","Accessibility","Expression")
write.table(avglogFC, file = "SuppFig2_avglogFC_astrocytoma.txt", sep = "\t", quote = FALSE)


avglogFC.melt <- melt(avglogFC)
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("SOX4","HIP1","MMP16","OLIG1","OLIG2","EGFR","PDGFRA","CHD7","FGF14")  

avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)

avglogFC.melt_interesting$gene <- factor(avglogFC.melt_interesting$gene, levels =  c("SOX4","HIP1","MMP16","OLIG1","OLIG2","EGFR","PDGFRA","CHD7","FGF14"))

pdf("Fig2I_main_RA_astrocytoma_ggdotchart.pdf", width=3.2, height=4)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                                                            axis.text.y = ggplot2::element_text(size = 14, face = "bold")) 
dev.off()

pdf("Fig2I_Suppl_RA_astrocytoma_ggscatter.pdf", width=3, height=4)

ggscatter(avglogFC, x = "Accessibility", y = "Expression",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Accesibility (avg_log2FC)", 
          ylab = "Expression (avg_log2FC)")
dev.off()

##dotplot to check stem like genes are only from RA

features_use <- c("SOX4","OLIG2","OLIG1", "HIP1","PDGFRA","FGF14", "MMP16","CHD7","EGFR")  

seurat <- readRDS("OE0145-IDH_integrated_astrocytoma_NMF_labelled_with_metacell_mapping_and_1p_scores.rds")
label <- readRDS("labels_astro_primary.rds")
seurat$transferred_labels <- label
table(seurat$transferred_labels)
seurat@active.assay <- "SCT"
Idents(seurat) <- "transferred_labels"
seurat = subset(seurat, idents =  c("Gradient","Astro-like", "OPC-like","RA","Cycling-like"))

seurat$transferred_labels <- factor(x = seurat$transferred_labels, levels =  c("Gradient","Astro-like", "OPC-like","RA","Cycling-like"))
Idents(seurat) = "transferred_labels"

pdf("Fig2I_main_RA_astrocytoma_DotPlot.pdf", width=4.8, height=4.5)
do_DotPlot(sample = seurat,
           features = rev(features_use),
           cols = c("grey75", "#3c5b8b"),
           ylab = "Astrocytoma",
           legend = TRUE,
           legend.position = "right",
           flip = T) 
dev.off()

#-------------------------------------------------------------------------------
#B. Astrocytoma
#B4. Motif Analysis
#------------------------------------------------------------------------------

seurat <- readRDS("OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar.rds")

label <- readRDS("annotation_ATAC_astro_primary.rds")

seurat$predicted.id <- label

data("human_pwms_v1")
length(human_pwms_v1)

seurat <- AddMotifs(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = human_pwms_v1
)

seurat <- RunChromVAR(
  object = seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(seurat, "OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar_IHL.rds")


path_to_object <- ""
tumor <- "astrocytoma"

# Paths
path_to_obj <- paste0(
  here::here(path_to_object),"OE0145-IDH_integrated_",
  tumor,"_peaks_activity_chromvar_IHL.rds",
  sep = ""
)

color_palette <-  c("#fed439","#709ae1","#4caf50","#d2af81",
                    "#d5e4a2","#9c27b0","#f05c3b","#46732e",
                    "#71d0f5","#370335","#075149","#c80813",
                    "#91331f","#1a9993","#fd8cc1")

path_to_save <- paste0(here::here(path_to_object),"OE0145-IDH_integrated_",
                       tumor,"_motifs_ChromVar_FC0_adjpvalue0.1.xls",
                       sep = "")


seurat <- readRDS(path_to_obj)
seurat = DietSeurat(seurat, assays = "chromvar")
#DefaultAssay(seurat) <- 'chromvar'
Idents(seurat) <- seurat$predicted.id
seurat <- subset(seurat, idents = c("Astro-like", "OPC-like","Gradient","RA"))
seurat$predicted.id <- factor(x = seurat$predicted.id, levels =  c("Gradient","Astro-like", "OPC-like","RA"))

da_regions <- FindAllMarkers(
  object = seurat,
  only.pos = TRUE,
  min.pct = 0.1,
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


da_regions_selected <- (da_regions[da_regions$p_val_adj < 0.05 & da_regions$avg_log2FC > 0, ])
#da_regions_selected <- (da_regions[da_regions$p_val_adj < 0.25 & da_regions$avg_log2FC > 0, ])


motif_name_entire <- da_regions_selected %>% separate(gene, 
                                                      c("Part1", "Part2","Part3",
                                                        "Part4","Part5"), sep="-")

da_regions_selected$motif_name <- motif_name_entire$Part3

da_regions_selected_sorted <- da_regions_selected[with(da_regions_selected,
                                                       order(cluster, -avg_log2FC)), ]

da_regions_selected_sorted$rank <- ave(da_regions_selected_sorted$avg_log2FC,
                                       da_regions_selected_sorted$cluster, FUN = seq_along)

da_regions_selected_sorted_prepared <- da_regions_selected_sorted %>% 
  group_by(cluster) %>% 
  top_n(20,-rank)

output <- split(da_regions_selected_sorted, da_regions_selected_sorted$cluster)
wb <- createWorkbook()

for (i in 1:length(output)) {
  addWorksheet(wb, sheetName=names(output[i]))
  writeData(wb, sheet=names(output[i]), x=output[[i]]) # Note [[]]
}

saveWorkbook(wb, file=path_to_save,  overwrite = TRUE)

cluster_cols <- c("RA" = "#ECA809",                          # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  #"Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  #"Astrocytes" = "#279185",                  # Celadon Green.
                  #"Microglia" = "#7EB356",                   # Bud Green.
                  #"Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  #"Neurons" = "#544B81",                     # Purple Navy.
                  #"Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E"                  # Ruby Red.
                  #"Undefined" = "#63412C"                   # Van Dyke Brown.
                  
)   

motifs = c("ENSG00000028277-LINE2636-POU2F2-D-N1", "ENSG00000171540-LINE2483-OTP-I",
           "ENSG00000136630-LINE2355-HLX-I","ENSG00000172789-LINE2485-AC0125311-I",
           "ENSG00000185610-LINE2537-DBX2-I","ENSG00000184486-LINE2677-POU3F2-D-N1",
           "ENSG00000140836-LINE2376-ZFHX3-D-N1","ENSG00000106261-LINE809-ZKSCAN1-I")


for (i in 1:length(motifs)) {
  pdf(paste0("Fig2J_main_AstroMotifs_VlnPlot_",motifs[i], ".pdf", sep = ""))
  print(do_VlnPlot(seurat,features = motifs[i] ,colors.use = cluster_cols, group.by = "predicted.id"))
  dev.off()
}
