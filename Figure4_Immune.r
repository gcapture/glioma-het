rm (list = ls())
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(future)
library(purrr)
library(harmony)
library(Seurat)
library(tidyverse)
library(plyr)
library(xlsx)
library("ggpubr")
library(SCpubr)

cluster_cols <- c("Mg homeostatic"= "#ECA809",                          # Marigold.
                  "Mg inflammatory TAMs" = "#043362",                    # Prussian Blue.
                  "Mg Activated"= "#9A031E",                  # Ruby Red      
                  "Mg resident-like TAMs"= "#009FF5",                     # Carolina Blue.
                  "Mo TAMs anti-inflammatory"= "#BC5210",            # Burnt Orange.
                  "Mg phagocytic"= "#279185",                  # Celadon Green.
                  "Mg IFNg TAMs"= "#7EB356",                   # Bud Green.
                  "Mg stressed TAMs"= "#AC70FF",                   # Medium Purple.
                  "Mg Inflammatory ICAM+"= "#63412C",                   # Van Dyke Brown.
                  "Mo TAMs Infiltrating"= "#D6D6D6")




Clusters  <- readRDS("Clusters_microglia_primary_annotated.rds")

Clusters$grade <- Clusters$gem_id
Clusters$grade <- revalue(Clusters$grade, c("OE0145-IDH_ACB_AD_540" = "grade3", "OE0145-IDH_ACB_AD_785" = "grade2", "OE0145-IDH_ACB_AD_832" = "grade2", "OE0145-IDH_ACB_AD_809" = "grade2",   
                                            "OE0145-IDH_ACB_AD_865" = "grade3", "OE0145-IDH_ACB_AD_883" = "grade3",  "OE0145-IDH_NCH2018" = "grade3" ,"OE0145-IDH_NCH2111" = "grade3",
                                            "OE0145-IDH_NCH2157" = "grade3"  ,"OE0145-IDH_NCH2164" = "grade3", "OE0145-IDH_NCH536" = "grade2" ,"OE0145-IDH_NCH6341" = "grade2","OE0145-IDH_NCH6702" = "grade3", "OE0145-IDH_NCH781"= "grade3"))

Clusters$Subclusters <- revalue(Clusters$Subclusters, c("Mg-IFNg TAMs" = "Mg IFNg TAMs", "Mo-TAMs anti-inflammatory" = "Mo TAMs anti-inflammatory", "Mo-TAMs Infiltrating"="Mo TAMs Infiltrating"))
Idents(Clusters) <- "Subclusters"
pdf("Fig4A_main_DimPlot_microglia_primary.pdf", width = 8, height = 7)
do_DimPlot(Clusters, colors.use = cluster_cols, plot.title = "Microglia Primary",legend.ncol = 3, legend.position = "bottom", fontsize = 15)
dev.off()


Clusters$sample <- revalue(Clusters$sample, c("primary_astro" = "Astrocytoma", "primary_oligo" = "Oligodendroglioma"))

pdf("Fig4A_main_Barplot_grade.pdf",width = 10, height = 3 )
do_BarPlot(Clusters, features = "sample",group.by = "Subclusters", 
           position = "fill", horizontal = TRUE, colors.use = cluster_cols, legend.position = "bottom", legend.ncol = 3, fontsize = 15)

dev.off()



Idents(Clusters) <- "sample"
scale.subtype <- c("Astrocytoma" = "#b38b14",
                   "Oligodendroglioma" = "#3c5b8b") 

pdf("Fig4A_main_DimPlot_microglia_primary_split_Oligo_Astro.pdf", width = 12, height = 5 )
do_DimPlot(Clusters, split.by = "sample", colors.split = scale.subtype,legend = FALSE, plot.title = "Microglia Primary",shuffle = FALSE)
dev.off()

Idents(Clusters) <- "Subclusters"
markers<- FindAllMarkers(Clusters,
                         test.use = "wilcox",
                         only.pos = TRUE,
                         logfc.threshold = 0.25,
                         verbose = TRUE
)
markers <- markers %>% filter(p_val_adj < 0.05)

output <- split(markers, markers$cluster)
wb <- createWorkbook()

for (i in 1:length(output)) {
  addWorksheet(wb, sheetName=names(output[i]))
  writeData(wb, sheet=names(output[i]), x=output[[i]], rowNames = TRUE) # Note [[]]
}

saveWorkbook(wb, file="Fig4_Suppl_markers_microglia_primary_samples_adjPvalue0.05.xls",  overwrite = TRUE)

Idents(Clusters) <- "samples"

##4D. ATAC Expression plots##

rm (list = ls())
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(future)
library(purrr)
library(harmony)
library(Seurat)
library(tidyverse)
library(plyr)
library(xlsx)
library("ggpubr")
library(SCpubr)
library(Signac)

peaks <- readRDS("/Users/Inma/Documents/Glioblastoma_Submission/Robjects/rds_figure2/OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar_IHL.rds")
label <- readRDS("/Users/Inma/Documents/Glioblastoma_Submission/Robjects/rds_figure2/annotation_ATAC_astro_primary.rds")
peaks$predicted.id <- label 

Idents(peaks) <- "predicted.id"
peaks <- subset(peaks, ident = "Microglia")

peaks_oligo <- readRDS("/Users/Inma/Documents/Glioblastoma_Submission/Robjects/rds_figure2/OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar.rds")
label <- readRDS("/Users/Inma/Documents/Glioblastoma_Submission/Robjects/rds_figure2/annotation_ATAC_oligo_primary.rds")
peaks_oligo$predicted.id <- label 

Idents(peaks_oligo) <- "predicted.id"

peaks_oligo <- subset(peaks_oligo, ident = "Microglia")

# first add dataset-identifying metadata
peaks$dataset <- "astro"
peaks_oligo$dataset <- "oligo"

# merge
Microglia_ATAC <- merge(peaks, peaks_oligo)
rm(peaks)
rm(peaks_oligo)

# process the combined dataset
Microglia_ATAC <- FindTopFeatures(Microglia_ATAC, min.cutoff = 10)
Microglia_ATAC <- RunTFIDF(Microglia_ATAC)
Microglia_ATAC <- RunSVD(Microglia_ATAC)
Microglia_ATAC <- RunUMAP(Microglia_ATAC, reduction = "lsi", dims = 2:30)

Idents(Microglia_ATAC) <- "dataset"

da_peaks <- FindMarkers(
  object = Microglia_ATAC,
  ident.1 = "astro",
  ident.2 = "oligo",
  min.pct = 0,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

da_peaks_astro = da_peaks[da_peaks$avg_log2FC > 0, ]
da_peaks_oligo = da_peaks[da_peaks$avg_log2FC < 0, ]
rm(da_peaks)

Clusters  <- readRDS("Clusters_microglia_primary_annotated.rds")
Idents(Clusters) <- "sample"
Clusters$sample <- revalue(Clusters$sample, c("primary_astro" = "Astrocytoma", "primary_oligo" = "Oligodendroglioma"))

Idents(Clusters) <- "sample"
Genes_DE_Microglia_primary <- FindMarkers(
  Clusters,
  ident.1 = "Astrocytoma",
  ident.2 = "Oligodendroglioma",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  verbose = TRUE
)

Astro_expression = Genes_DE_Microglia_primary[Genes_DE_Microglia_primary$avg_log2FC > 0 & Genes_DE_Microglia_primary$p_val_adj < 0.05, ]
Oligo_expression = Genes_DE_Microglia_primary[Genes_DE_Microglia_primary$avg_log2FC < 0 & Genes_DE_Microglia_primary$p_val_adj < 0.05, ]
rm(Genes_DE_Microglia_primary)

## Astrocytoma plots
peaks <- rownames(da_peaks_astro[da_peaks_astro$avg_log2FC > 0, ])
closest_genes_astro <- ClosestFeature(Microglia_ATAC, regions = peaks)
rownames(closest_genes_astro) <- closest_genes_astro$query_region
gene_FC <- merge(da_peaks_astro, closest_genes_astro, by=0)
gene_FC <- gene_FC[,c("query_region","avg_log2FC","gene_name")]
Astro_expression$gene_name <- rownames(Astro_expression)


common_genes_classif <- merge(Astro_expression, gene_FC, by = "gene_name")
mean_accessibility <- tapply(common_genes_classif$avg_log2FC.y,
                             common_genes_classif$gene, mean)
mean_expresion <- tapply(common_genes_classif$avg_log2FC.x,
                         common_genes_classif$gene, mean)

avglogFC <- cbind(melt(mean_accessibility),
                  melt(mean_expresion)[,2])
colnames(avglogFC) <- c("gene","Accessibility","Expression")
write.table(avglogFC, file = "SuppFig4_avglogFC_astrocytoma.txt", sep = "\t", quote = FALSE)


avglogFC.melt <- melt(avglogFC)


##1.Mg resident-like TAMs
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("SERPINE1", "RGS1", "INSR", "NLRP3")
avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)


pdf("Fig4D_Astrocytoma_Mg_resident_like_TAMs_ggdotchart.pdf", width=3.5, height=2.5)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "", title = "Mg resident-like TAMs") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                                                            axis.text.y = ggplot2::element_text(size = 14, face = "bold"),
                                                            plot.title = element_text(color="black", size=14, face="bold.italic")) 
dev.off()

##2.Mg IFNg TAMs gene markers
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("SPTLC2", "L3MBTL4", "IFI44L","TGFBR2")
avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)


pdf("Fig4D_Astrocytoma_Mg_IFNg_TAMs_ggdotchart.pdf", width=3.5, height=2.5)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "", title = "Mg IFNg TAMs") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                                                            axis.text.y = ggplot2::element_text(size = 14, face = "bold"),
                                                            plot.title = element_text(color="black", size=14, face="bold.italic")) 
dev.off()


## Oligodendroglioma plots
peaks <- rownames(da_peaks_oligo[da_peaks_oligo$avg_log2FC < 0, ])
closest_genes_oligo <- ClosestFeature(Microglia_ATAC, regions = peaks)
rownames(closest_genes_oligo) <- closest_genes_oligo$query_region
gene_FC <- merge(da_peaks_oligo, closest_genes_oligo, by=0)
gene_FC <- gene_FC[,c("query_region","avg_log2FC","gene_name")]
Oligo_expression$gene_name <- rownames(Oligo_expression)


common_genes_classif <- merge(Oligo_expression, gene_FC, by = "gene_name")
mean_accessibility <- tapply(common_genes_classif$avg_log2FC.y,
                             common_genes_classif$gene, mean)
mean_expresion <- tapply(common_genes_classif$avg_log2FC.x,
                         common_genes_classif$gene, mean)

avglogFC <- cbind(melt(mean_accessibility),
                  melt(mean_expresion)[,2])
colnames(avglogFC) <- c("gene","Accessibility","Expression")
write.table(avglogFC, file = "SuppFig4_avglogFC_Oligodendroglioma.txt", sep = "\t", quote = FALSE)


avglogFC.melt <- melt(avglogFC)


##1.Mo TAMs anti-inflammatory
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("FCGRT", "GADD45B", "MERTK", "HLA-DRB5")
avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)
avglogFC.melt_interesting$value = avglogFC.melt_interesting$value*(-1)

pdf("Fig4D_Oligodendroglioma_Mo_TAMs_anti_inflammatory_ggdotchart.pdf", width=4, height=2.5)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "", title = "Mo TAMs anti-inflammatory") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                                                            axis.text.y = ggplot2::element_text(size = 14, face = "bold"),
                                                            plot.title = element_text(color="black", size=14, face="bold.italic")) 
dev.off()

##2.Mo TAMs infiltrating
avglogFC.melt_interesting <- avglogFC.melt
features_use <- c("LYVE1","THRB","AFF3","RGL1")
avglogFC.melt_interesting = avglogFC.melt_interesting[avglogFC.melt_interesting$gene %in% features_use,] %>% arrange(gene)
avglogFC.melt_interesting$value = avglogFC.melt_interesting$value*(-1)

pdf("Fig4D_Oligodendroglioma_Mo_TAMs_infiltrating_ggdotchart.pdf", width=3.5, height=2.5)
ggdotchart(avglogFC.melt_interesting, x="gene", 
           y="value", color="variable", 
           add = "segments", ylab = "mean avglog2FC", xlab = "", title = "Mo TAMs infiltrating") + 
  coord_flip() + facet_wrap(~variable) + NoLegend() + theme(axis.text.x = ggplot2::element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
                                                            axis.text.y = ggplot2::element_text(size = 14, face = "bold"),
                                                            plot.title = element_text(color="black", size=14, face="bold.italic")) 
dev.off()



##boxplot

table <- read.table("Microglia_onlyPrimary_counts_percent.txt", sep="\t")


for(i in 2:ncol(table)) {    
  
  p <- ggplot(data=table, aes(y = table[,i],x = sample, fill= sample))+
    geom_boxplot() 
  p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, fill = "black")+ 
    scale_fill_manual(values=c("#0d798c", "#cb5a28"))+
    stat_compare_means(paired = FALSE,label.y = (max(table[i]) + 1))+
    theme(legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
  #theme(legend.position="none",text = element_text(family = "Univers-Black", size = 10),axis.text=element_text(family = "Univers-Black", size = 10),axis.title=element_text(family = "Univers-Black", size = 10))
  
  p <- p +  ggtitle(colnames(table)[i]) + xlab("") + ylab("%")
  print(p)    
  
  mypath <- file.path(paste("boxplot_",colnames(table)[i], ".pdf", sep = ""))
  
  pdf(file=mypath, width = 2.5,height = 3.5)
  print(p)
  dev.off()
}