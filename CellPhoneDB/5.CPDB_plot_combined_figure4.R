# CPDB docs - plotting functions
# https://github.com/Teichlab/cellphonedb/blob/master/cellphonedb/src/plotters/R/plot_dot_by_column_name.R
rm(list = ls())

library(ggplot2)

#0. Common genes

genes <- c("SIRPA_CD47","TNF_PTPRS","TNF_NOTCH1","BMR1A_ACR2A_BMP2",
           "BMPR1B_BMPR2_BMP2","BMPR1A_BMPR2_BMP2","BMR1A_ACR2A_BMP7",
           "BMPR1B_BMPR2_BMP7","BMPR1A_BMPR2_BMP7","ACVR_1A2A receptor_BMP7",
           "ACVR1_BMPR2_BMP7","WNT5A_FZD3","WNT5A_PTPRK","NRP2_VEGFA",
           "NRP1_VEGFA","SPP1_a9b1 complex","PDGFB_PDGFRB","TNF_PTPRS",
           "TNF_RIPK1","DLL1_NOTCH1","TNF_RIPK1","DLL1_NOTCH1",
           "TNF_NOTCH1","LGALS9_DAG1","TNF_DAG1", "NRP1_PGF","NRP1_VEGFA",
           "NRP1_PGF","NRP1_VEGFA","BMR1A_ACR2A_BMP7","ACVR_1A2A receptor_BMP7",
           "ACVR1_BMPR2_BMP7", "CSF1R_CSF1","WNT5A_FZD3","LGALS9_MRC2",
           "LGALS9_COLEC12","CD40_TNFSF13B","TFRC_TNFSF13B","HLA-DPB1_TNFSF13B",
           "BMR1B_AVR2A_BMP7","SPP1_a9b1 complex")



##1. Analysis cellphoneDB OPC-Like Astrocytoma:

library(ggplot2)

dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = NULL,
                    pvalues_path = NULL,
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  print(head(plot.data))
  print(dim(plot.data))
  return(plot.data)
  
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
}

my_palette <- colorRampPalette(c("blue4", "dodgerblue4", "deepskyblue4", "deepskyblue"), alpha=TRUE)(n=399)

data <- dot_plot(selected_rows = NULL,
                 selected_columns = NULL,
                 filename = 'plot_astrocytoma.pdf',
                 width = 8,
                 height = 10,
                 means_path = "Astrocytoma_OPCLike/means_fixorder.txt",
                 pvalues_path = "Astrocytoma_OPCLike/pvalues.txt",
                 means_separator = '\t',
                 pvalues_separator = '\t',
                 output_extension = '.pdf'
)

#data <- data[data$pvalue < 0.05,]
data <- data[grep("OPC_like",data$clusters),]
data <- data[!(data$clusters %in% "OPC_like|OPC_like"),]
data <- data[grepl("\\|OPC_like",data$clusters),]

data <- data[data$pair %in% genes,]

data$clusters_plot <- paste("Astrocytoma",data$clusters,sep= ":")

ggplot(data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

data_OPC_Like_Astrocytoma <- data


##2. Analysis cellphoneDB OPC-Like Oligodendroglioma:

data <- dot_plot(selected_rows = NULL,
                 selected_columns = NULL,
                 filename = 'plot_oligodendroglioma.pdf',
                 width = 8,
                 height = 10,
                 means_path = "Oligodendroglioma_OPCLike/means_fixorder.txt",
                 pvalues_path = "Oligodendroglioma_OPCLike/pvalues.txt",
                 means_separator = '\t',
                 pvalues_separator = '\t',
                 output_extension = '.pdf'
)

#data <- data[data$pvalue < 0.05,]
data <- data[grep("OPC_like",data$clusters),]
data <- data[!(data$clusters %in% "OPC_like|OPC_like"),]
data <- data[grepl("\\|OPC_like",data$clusters),]

data <- data[data$pair %in% genes,]

data$clusters_plot <- paste("Oligodendroglioma",data$clusters,sep= ":")

ggplot(data,aes(x=clusters_plot,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

data_OPC_Like_Oligodendroglioma <- data


##3. Analysis cellphoneDB Astro-Like Astrocytoma:

data <- dot_plot(selected_rows = NULL,
                 selected_columns = NULL,
                 filename = 'plot_astrocytoma.pdf',
                 width = 8,
                 height = 10,
                 means_path = "Astrocytoma_AstroLike/means_fixorder.txt",
                 pvalues_path = "Astrocytoma_AstroLike/pvalues.txt",
                 means_separator = '\t',
                 pvalues_separator = '\t',
                 output_extension = '.pdf'
)

#data <- data[data$pvalue < 0.05,]
data <- data[grep("Astro_like",data$clusters),]
data <- data[!(data$clusters %in% "Astro_like|Astro_like"),]
data <- data[grepl("\\|Astro_like",data$clusters),]

data <- data[data$pair %in% genes,]

data$clusters_plot <- paste("Astrocytoma",data$clusters,sep= ":")
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

ggplot(data,aes(x=clusters_plot,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

data_Astro_Like_Astrocytoma <- data

##4. Analysis cellphoneDB Astro-Like Oligodendroglioma:

data <- dot_plot(selected_rows = NULL,
                 selected_columns = NULL,
                 filename = 'plot_astrocytoma.pdf',
                 width = 8,
                 height = 10,
                 means_path = "Oligodendroglioma_AstroLike/means_fixorder.txt",
                 pvalues_path = "Oligodendroglioma_AstroLike/pvalues.txt",
                 means_separator = '\t',
                 pvalues_separator = '\t',
                 output_extension = '.pdf'
)

#data <- data[data$pvalue < 0.05,]
data <- data[grep("Astro_like",data$clusters),]
data <- data[!(data$clusters %in% "Astro_like|Astro_like"),]
data <- data[grepl("\\|Astro_like",data$clusters),]

data <- data[data$pair %in% genes,]

data$clusters_plot <- paste("Oligodendroglioma",data$clusters,sep= ":")

ggplot(data,aes(x=clusters_plot,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

data_Astro_Like_Oligodendroglioma <- data


data_plot <- as.data.frame(rbind(data_OPC_Like_Astrocytoma, data_OPC_Like_Oligodendroglioma, data_Astro_Like_Astrocytoma, data_Astro_Like_Oligodendroglioma))
levels <- unique(data_plot$clusters_plot)


head(data_plot, n = 100)
#data_plot$clusters_plot <- as.character(data_plot$clusters_plot)
data_plot <- data_plot  %>%
  mutate(clusters_plot = factor(clusters_plot, levels=levels))

pdf("DotPlot_CellPhoneDB_Figure4.pdf", width=15, height=12)

ggplot(data_plot,aes(x=clusters_plot,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

dev.off()

