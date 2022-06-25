library(magrittr)
library(ComplexHeatmap)
# Metadata
figures_folder <- "/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science"
figure_type <- "main_figures"
figure_path <- sprintf("%s/second_iteration/%s/", figures_folder, figure_type)
project_folder <- "/omics/odcf/analysis/OE0145_projects/idh_gliomas"

sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach")
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach")

#saveRDS(sample.oligo, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach", compress = F)
#saveRDS(sample.astro, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach", compress = F)


# Colors clusters:
cluster_cols <- c("RA" = "#ECA809",
                  "RA_MS" = "#ECA809",                    # Marigold.
                  "OPC-like" = "#043362",                    # Prussian Blue.
                  "Oligo-like" = "#043362",
                  "Oligo program" = "#043362",                    # Prussian Blue.
                  "OPC_MS" = "#043362",                    # Prussian Blue.
                  "T Cells" = "#009FF5",                     # Carolina Blue.
                  "T Cell" = "#009FF5",                     # Carolina Blue.
                  "Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                  "Astrocytes" = "#279185",                  # Celadon Green.
                  "Microglia" = "#7EB356",                   # Bud Green.
                  "Pericytes" = "#AC70FF",                   # Medium Purple.
                  "Mixed" = "#8d5d3f",                   # Van Dyke Brown.
                  "Gradient" = "#D6D6D6",                    # Light grey.
                  "Neurons" = "#544B81",                     # Purple Navy.
                  "Neuron" = "#544B81",                     # Purple Navy.
                  "Endothelial" = "#da627d",                 # Blush.
                  "Astro-like" = "#9A031E",                  # Ruby Red.
                  "Astro_MS" = "#9A031E",                  # Ruby Red.
                  "Astro program" = "#9A031E",                  # Ruby Red.
                  "Excluded" = "#4D6880",                    # Dark Electric Blue.
                  "Cycling" = "#5F0F40",
                  "Cycling_MS" = "#5F0F40")                # Tyrian Purple.
cluster_cols <- cluster_cols[sort(names(cluster_cols))]


scale.grade <- c("2" = "#94d2bd",
                 "3" = "#005F73",
                 "4" = "#626D87")

scale.subtype <- c("Astrocytoma" = "#b38b14",
                   "Oligodendroglioma" = "#3c5b8b")
scale.subtype.short <- c("AS" = "#b38b14",
                   "OD" = "#3c5b8b")

scale.relapse <- c("Primary" = "#8D99AE",
                   "Relapse" = "#1C2839")

scale.ident <- c("IDH_ACB_AD_540" = "#E07A5F",
                 "IDH_ACB_AD_809" = "#F2CC8F",
                 "IDH_ACB_AD_883" = "#DAA82B",
                 "IDH_NCH2111" = "#81B29A",
                 "IDH_NCH536" = "#3E745E",
                 "IDH_NCH6341" = "#5CADC1",
                 "IDH_NCH6702" = "#8B6BB8",
                 "IDH_NCH781" = "#3D405B",
                 "IDH_ACB_AD_785" = "#ae2012",
                 "IDH_ACB_AD_832" = "#ca6702",
                 "IDH_ACB_AD_865" = "#ee9b00",
                 "IDH_NCH2018" = "#e9d8a6",
                 "IDH_NCH2157" = "#457b9d",
                 "IDH_NCH2164" = "#1d3557")




Figure_1A <- function(){
  metadata <- as.data.frame(readxl::read_excel("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/metadata/IDH_gliomas_sample_metadata.xlsx"))
  rownames(metadata) <- metadata$Samples
  metadata$Samples <- NULL
  metadata <- metadata %>% dplyr::filter(Diagnosis != "GBM-like")

  do_heatmap <- function(data, name, col){
    fontsize <- 22
    h <- ComplexHeatmap::Heatmap(data,
                                 name = name,
                                 col = col,
                                 show_row_dend = FALSE,
                                 show_column_dend = FALSE,
                                 column_names_gp = grid::gpar(fontsize = fontsize,
                                                              fontface = "bold"),
                                 row_names_gp = grid::gpar(fontsize = fontsize,
                                                           fontface = "bold"),
                                 row_names_side = "left",
                                 column_title_gp = grid::gpar(fontsize = fontsize,
                                                              fontface = "bold"),
                                 border = TRUE,
                                 heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize,
                                                                                   fontface = "bold"),
                                                             labels_gp = grid::gpar(fontsize = fontsize)),
                                 rect_gp = grid::gpar(col= "white"))
    return(h)
  }
  labels.use <- c("IDH_ACB_AD_809",
                  "IDH_NCH6341",
                  "IDH_NCH536",
                  "IDH_ACB_AD_540",
                  "IDH_NCH6702",
                  "IDH_NCH2111",
                  "IDH_NCH781",
                  "IDH_ACB_AD_883",
                  "IDH_ACB_AD_785",
                  "IDH_ACB_AD_832",
                  "IDH_NCH2018",
                  "IDH_ACB_AD_865",
                  "IDH_NCH2157",
                  "IDH_NCH2164")

  no_color <- "#98b9cd"
  yes_color <- "#1b4965"
  col_vector_snRNAseq <- c(yes_color, no_color)
  names(col_vector_snRNAseq) <- c("Yes", "No")
  h_snRNAseq <- do_heatmap(data = t(as.matrix(metadata[labels.use, "snRNAseq", drop = FALSE])), name = "snRNAseq", col = col_vector_snRNAseq)
  t_h_snRNAseq <- do_heatmap(data = as.matrix(metadata[labels.use, "snRNAseq", drop = FALSE]), name = "snRNAseq", col = col_vector_snRNAseq)


  col_vector_snATACseq <- c(yes_color, no_color)
  names(col_vector_snATACseq) <- c("Yes", "No")
  h_snATACseq <- do_heatmap(data = t(as.matrix(metadata[labels.use, "snATACseq", drop = FALSE])), name = "snATACseq", col = col_vector_snATACseq)
  t_h_snATACseq <- do_heatmap(data = as.matrix(metadata[labels.use, "snATACseq", drop = FALSE]), name = "snATACseq", col = col_vector_snATACseq)


  col_vector_IDH_status<- c(yes_color, no_color)
  names(col_vector_IDH_status) <- c("Yes", "No")
  h_IDH_status <- do_heatmap(data = t(as.matrix(metadata[labels.use, "IDH mutated", drop = FALSE])), name = "IDH mutated", col = col_vector_IDH_status)
  t_h_IDH_status <- do_heatmap(data = as.matrix(metadata[labels.use, "IDH mutated", drop = FALSE]), name = "IDH mutated", col = col_vector_IDH_status)


  col_vector_crh1_chr19_codeletion<- c(yes_color, no_color)
  names(col_vector_crh1_chr19_codeletion) <- c("Yes", "No")
  h_crh1_chr19_codeletion <- do_heatmap(data = t(as.matrix(metadata[labels.use, "1p/19q codeletion", drop = FALSE])), name = "1p/19q codeletion", col = col_vector_crh1_chr19_codeletion)
  t_h_crh1_chr19_codeletion <- do_heatmap(data = as.matrix(metadata[labels.use, "1p/19q codeletion", drop = FALSE]), name = "1p/19q codeletion", col = col_vector_crh1_chr19_codeletion)


  col_vector_MGMT<- c(yes_color, no_color, "grey75")
  names(col_vector_MGMT) <- c("Yes", "No", "NA")
  h_MGMT <- do_heatmap(data = t(as.matrix(metadata[labels.use, "MGMT methylated", drop = FALSE])), name = "MGMT methylated", col = col_vector_MGMT)
  t_h_MGMT <- do_heatmap(data = as.matrix(metadata[labels.use, "MGMT methylated", drop = FALSE]), name = "MGMT methylated", col = col_vector_MGMT)


  col_vector_TERT<- c("#1a936f", "#ca6702", "#8c2f39", "grey75")
  names(col_vector_TERT) <- c("WT", "C228T", "C250T", "NA")
  h_TERT <- do_heatmap(data = t(as.matrix(metadata[labels.use, "TERT status", drop = FALSE])), name = "TERT status", col = col_vector_TERT)
  t_h_TERT <- do_heatmap(data = as.matrix(metadata[labels.use, "TERT status", drop = FALSE]), name = "TERT status", col = col_vector_TERT)


  col_vector_ATRX<- c("#4f6d7a", "#99582a", "grey75")
  names(col_vector_ATRX) <- c("WT", "Loss", "NA")
  h_ATRX <- do_heatmap(data = t(as.matrix(metadata[labels.use, "ATRX status", drop = FALSE])), name = "ATRX status", col = col_vector_ATRX)
  t_h_ATRX <- do_heatmap(data = as.matrix(metadata[labels.use, "ATRX status", drop = FALSE]), name = "ATRX status", col = col_vector_ATRX)


  col_vector_Diagnosis<- c("#3c5b8b", "#b38b14", "grey75")
  names(col_vector_Diagnosis) <- c("Oligodendroglioma", "Astrocytoma", "NA")
  h_Diagnosis <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Diagnosis", drop = FALSE])), name = "Diagnosis", col = col_vector_Diagnosis)
  h_Diagnosis
  t_h_Diagnosis <- do_heatmap(data = as.matrix(metadata[labels.use, "Diagnosis", drop = FALSE]), name = "Diagnosis", col = col_vector_Diagnosis)


  col_vector_Histological_grade <- c("#94d2bd", "#005f73")
  names(col_vector_Histological_grade) <- c("2", "3")
  h_Histological_grade <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Grade", drop = FALSE])), name = "Grade", col = col_vector_Histological_grade)
  t_h_Histological_grade <- do_heatmap(data = as.matrix(metadata[labels.use, "Grade", drop = FALSE]), name = "Grade", col = col_vector_Histological_grade)

  col_vector_gender <- c("#723d46", "#af9d6a")
  names(col_vector_gender) <- c("Male", "Female")
  h_gender <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Sex", drop = FALSE])), name = "Sex", col = col_vector_gender)
  t_h_gender <- do_heatmap(data = as.matrix(metadata[labels.use, "Sex", drop = FALSE]), name = "Sex", col = col_vector_gender)



  ht_list <- h_snRNAseq %v%
    h_snATACseq %v%
    h_IDH_status %v%
    h_crh1_chr19_codeletion %v%
    h_MGMT %v%
    h_TERT %v%
    h_ATRX %v%
    h_Diagnosis %v%
    h_Histological_grade %v%
    h_gender

  t_ht_list <- t_h_snRNAseq +
    t_h_snATACseq +
    t_h_IDH_status +
    t_h_crh1_chr19_codeletion +
    t_h_MGMT +
    t_h_TERT +
    t_h_ATRX +
    t_h_Diagnosis +
    t_h_Histological_grade +
    t_h_gender

  return_list <- list("landscape" = ht_list,
                      "portrait" = t_ht_list)
  return(return_list)

}

Figure_3A <- function(){
  metadata <- as.data.frame(readxl::read_excel("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/metadata/IDH_gliomas_paired_samples_metadata.xlsx"))
  rownames(metadata) <- metadata$Samples
  metadata$Samples <- NULL

  do_heatmap <- function(data, name, col){
    fontsize <- 22
    h <- ComplexHeatmap::Heatmap(data,
                                 name = name,
                                 col = col,
                                 cluster_rows = F,
                                 column_split = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5", "6", "6"),
                                 column_title = c("", "", "", "", "", ""),
                                 column_gap = grid::unit(c(2, 2, 2, 2, 2, 2), "mm"),
                                 column_names_gp = grid::gpar(fontsize = fontsize,
                                                              fontface = "bold"),
                                 row_names_gp = grid::gpar(fontsize = fontsize,
                                                           fontface = "bold"),
                                 row_names_side = "left",
                                 column_title_gp = grid::gpar(fontsize = fontsize,
                                                              fontface = "bold"),
                                 border = TRUE,
                                 heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize,
                                                                                   fontface = "bold"),
                                                             labels_gp = grid::gpar(fontsize = fontsize)),
                                 rect_gp = grid::gpar(col= "white"))
    return(h)
  }

  labels.use <- c("IDH_NCH557",
                  "IDH_NCH758a",
                  "IDH_NCH511b",
                  "IDH_NCH678k",
                  "IDH_NCH302",
                  "IDH_NCH645",
                  "IDH_NCH988",
                  "IDH_NCH2375",
                  "IDH_NCH740w",
                  "IDH_NCH2367",
                  "IDH_NCH673d",
                  "IDH_NCH2260")

  no_color <- "#98b9cd"
  yes_color <- "#1b4965"
  na_color <- "grey75"

  col_vector_treatment <- c("white", "#8a817c", "#dda15e", "#3a5a40", "#903738")
  names(col_vector_treatment) <- c("Primary sample", "None", "RT", "TMZ", "TMZ and RT")
  h_treatment <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Treatment", drop = FALSE])), name = "Treatment", col = col_vector_treatment)

  col_vector_TMZ_cycles <- c("white", "#a8dadc", "#457b9d", "#465686")
  names(col_vector_TMZ_cycles) <- c("None", "5", "8", "12")
  h_TMZ_cycles <- do_heatmap(data = t(as.matrix(metadata[labels.use, "TMZ cycles", drop = FALSE])), name = "TMZ cycles", col = col_vector_TMZ_cycles)


  col_vector_IDH_mutated <- c(yes_color, no_color)
  names(col_vector_IDH_mutated) <- c("Yes", "No")
  h_IDH_mutated <- do_heatmap(data = t(as.matrix(metadata[labels.use, "IDH mutated", drop = FALSE])), name = "IDH mutated", col = col_vector_IDH_mutated)

  col_vector_1p19q_codeletion <- c(yes_color, no_color, na_color)
  names(col_vector_1p19q_codeletion) <- c("Yes", "No", "NA")
  h_1p19q_codeletion <- do_heatmap(data = t(as.matrix(metadata[labels.use, "1p/19q", drop = FALSE])), name = "1p/19q codeletion", col = col_vector_1p19q_codeletion)

  col_vector_MGMT <- c(yes_color, no_color, na_color)
  names(col_vector_MGMT) <- c("Yes", "No", "NA")
  h_MGMT <- do_heatmap(data = t(as.matrix(metadata[labels.use, "MGMT methylated", drop = FALSE])), name = "MGMT methylated", col = col_vector_MGMT)

  col_vector_TERT <- c(yes_color, no_color, na_color)
  names(col_vector_MGMT) <- c("Yes", "No", "NA")
  h_TERT <- do_heatmap(data = t(as.matrix(metadata[labels.use, "TERT status", drop = FALSE])), name = "TERT status", col = col_vector_MGMT)

  col_vector_ATRX <- c(yes_color, no_color, na_color)
  names(col_vector_ATRX) <- c("Yes", "No", "NA")
  h_ATRX <- do_heatmap(data = t(as.matrix(metadata[labels.use, "ATRX status", drop = FALSE])), name = "ATRX methylated", col = col_vector_ATRX)

  col_vector_gender <- c(yes_color, no_color, na_color)
  names(col_vector_ATRX) <- c("Yes", "No", "NA")
  h_ATRX <- do_heatmap(data = t(as.matrix(metadata[labels.use, "ATRX status", drop = FALSE])), name = "ATRX methylated", col = col_vector_ATRX)

  col_vector_Diagnosis<- c("#3c5b8b", "#b38b14", "#6F8A2B", "grey75")
  names(col_vector_Diagnosis) <- c("Oligodendroglioma", "Astrocytoma", "sGBM", "NA")
  h_Diagnosis <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Diagnosis", drop = FALSE])), name = "Diagnosis", col = col_vector_Diagnosis)

  col_vector_Histological_grade <- c("#94d2bd", "#005f73", "#626D87")
  names(col_vector_Histological_grade) <- c("2", "3", "4")
  h_Histological_grade <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Grade", drop = FALSE])), name = "Grade", col = col_vector_Histological_grade)

  col_vector_gender <- c("#723d46", "#af9d6a")
  names(col_vector_gender) <- c("Male", "Female")
  h_gender <- do_heatmap(data = t(as.matrix(metadata[labels.use, "Sex", drop = FALSE])), name = "Sex", col = col_vector_gender)

  h <-  h_IDH_mutated %v% h_1p19q_codeletion %v% h_MGMT %v% h_treatment %v% h_TMZ_cycles %v% h_Diagnosis %v% h_Histological_grade %v% h_gender

  return(h)
}

cellphoneDB_output <- function(interactions){
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(magrittr)


  interactions <- dplyr::as_tibble(readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/CellPhoneDB.rds"))
  interactions$pair <- as.character(interactions$pair)
  interactions$group <- "None"
  interactions$group[interactions$pair %in% c("CSF1R_CSF1")] <- "G1"
  interactions$group[interactions$pair %in% c("WNT5A_FZD3",
                                              "WNT5A_PTPRK")] <- "G2"
  interactions$group[interactions$pair %in% c("PDGFB_PDGFRB")] <- "G3"
  interactions$group[interactions$pair %in% c("SIRPA_CD47")] <- "G4"
  interactions$group[interactions$pair %in% c("SPP1_a9b1 complex")] <- "G5"
  interactions$group[interactions$pair %in% c("LGALS9_COLEC12",
                                              "LGALS9_DAG1",
                                              "LGALS9_MRC2")] <- "G6"
  interactions$group[interactions$pair %in% c("TNF_DAG1",
                                              "TNF_PTPRS",
                                              "TNF_RIPK1",
                                              "TNF_NOTCH1",
                                              "DLL1_NOTCH1")] <- "G7"
  interactions$group[interactions$pair %in% c("NRP1_PGF",
                                              "NRP1_VEGFA",
                                              "NRP2_VEGFA")] <- "G8"
  interactions$group[interactions$pair %in% c("CD40_TNFSF13B",
                                              "TFRC_TNFSF13B",
                                              "HLA-DPB1_TNFSF13B")] <- "G9"
  interactions$group[interactions$pair %in% c("BMPR1A_BMPR2_BMP2",
                                              "BMPR1A_BMPR2_BMP7",
                                              "BMPR1B_BMPR2_BMP2",
                                              "BMPR1B_BMPR2_BMP7",
                                              "BMR1A_ACR2A_BMP2",
                                              "BMR1A_ACR2A_BMP7",
                                              "BMR1B_AVR2A_BMP7",
                                              "ACVR_1A2A receptor_BMP7",
                                              "ACVR1_BMPR2_BMP7")] <- "G10"
  interactions$group <- factor(interactions$group, levels = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"))






  replace_last_underscore <- function(x){
    # Reverse the string.
    x <- paste(rev(stringr::str_split(x, pattern = "")[[1]]), collapse = "")
    # Replace
    x <- stringr::str_replace(x, "_", " | ")
    # reverse
    x <- paste(rev(stringr::str_split(x, pattern = "")[[1]]), collapse = "")
    return(x)
  }

  interactions <- interactions %>%
    dplyr::select(pair, pvalue, mean, clusters_plot, group) %>% # Remove unnecessary clusters column.
    dplyr::mutate(pvalue = -log10(pvalue)) %>%
    #dplyr::mutate(mean = log2(mean)) %>%
    dplyr::mutate(lgg_subtype = dplyr::case_when(stringr::str_starts(.$clusters_plot, "^Astrocytoma:")  == TRUE ~ "Astrocytoma", # Reassign the names back to Astrocytoma.
                                                 stringr::str_starts(.$clusters_plot, "^Oligodendroglioma:")  == TRUE ~ "Oligodendroglioma")) %>% # Reassign the names back to Oligodendroglioma.
    dplyr::mutate(target_cluster = dplyr::case_when(stringr::str_ends(.$clusters_plot, "Astro_like$")  == TRUE ~ "Astro-like", # Reassign the names back to targetting Astro-like.
                                                    stringr::str_ends(.$clusters_plot, "OPC_like$")  == TRUE ~ "OPC-like")) %>% # Reassign the names back to targetting OPC-like.
    dplyr::mutate(clusters_plot = stringr::str_remove_all(.$clusters_plot, "^Astrocytoma:|^Oligodendroglioma:|\\|OPC_like$|\\|Astro_like")) %>% # Remove unnecessary long cluster names.
    dplyr::rowwise() %>%
    dplyr::mutate(pair = replace_last_underscore(pair)) %>%
    dplyr::mutate(gene1 = stringr::str_split(pair, pattern = " \\| ", n = 2)[[1]][1]) %>%
    dplyr::mutate(gene2 = stringr::str_split(pair, pattern = " \\| ", n = 2)[[1]][2]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(`-log10(p-value)` = pvalue) %>%
    dplyr::mutate(`log2 mean(Molecule 1, Molecule 2)` = mean)


  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_Activated"] <- "Activated"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_IFNg_TAMs"] <- "IFNg"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_Inflammatory_ICAM+"] <- "Inflammatory ICAM+"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_homeostatic"] <- "Homeostatic"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_inflammatory_TAMs"] <- "Inflammatory"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_phagocytic"] <- "Phagocytic"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_resident_like_TAMs"] <- "Resident-like"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mg_stressed_TAMs"] <- "Stressed"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mo_TAMs_Infiltrating"] <- "Infiltrating"
  interactions$clusters_plot[interactions$clusters_plot %in% "Mo_TAMs_anti_inflammatory"] <- "Anti-inflammatory"
  interactions$clusters_plot <- factor(interactions$clusters_plot, levels = sort(unique(interactions$clusters_plot)))


  row_names <- interactions %>%
    dplyr::select(pair, pvalue, clusters_plot) %>%
    tidyr::pivot_wider(names_from = clusters_plot, values_from = pvalue) %>%
    dplyr::pull(pair)


  col_names <- sort(unique(interactions$clusters_plot))
  input.astrocytoma.astro_like <- interactions %>% dplyr::filter(lgg_subtype == "Astrocytoma" & target_cluster == "Astro-like") %>% dplyr::mutate(pair = factor(pair, levels = row_names)) %>% dplyr::mutate(clusters_plot = factor(clusters_plot, levels = col_names))
  limits <- c(min(interactions$mean), max(interactions$mean))
  color_low <-"#001FA9"
  color_high <- "#00A6A9"
  interacting_color_1 <- "#A90040"
  interacting_color_2 <- "#A95800"
  p4 <- ggplot2::ggplot(input.astrocytoma.astro_like, mapping = ggplot2::aes(x = clusters_plot, y = pair)) +
    ggplot2::geom_point(mapping = ggplot2::aes(size = `-log10(p-value)`, color = `log2 mean(Molecule 1, Molecule 2)`)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Astro-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", color = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(face = "bold")) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_color_gradient(limits = limits, low = color_low, high = color_high) +
    ggplot2::facet_grid(rows = ggplot2::vars(group), scales = "free", space = "free_y")

  input.astrocytoma.opc_like <- interactions %>% dplyr::filter(lgg_subtype == "Astrocytoma" & target_cluster == "OPC-like") %>% dplyr::mutate(pair = factor(pair, levels = row_names)) %>% dplyr::mutate(clusters_plot = factor(clusters_plot, levels = col_names))

  p3 <- ggplot2::ggplot(input.astrocytoma.opc_like, mapping = ggplot2::aes(x = clusters_plot, y = pair)) +
    ggplot2::geom_point(mapping = ggplot2::aes(size = `-log10(p-value)`, color = `log2 mean(Molecule 1, Molecule 2)`)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("OPC-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", color = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(face = "bold")) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_color_gradient(limits = limits, low = color_low, high = color_high) +
    ggplot2::facet_grid(rows = ggplot2::vars(group), scales = "free", space = "free_y")

  input.oligodendroglioma.astro_like <- interactions %>% dplyr::filter(lgg_subtype == "Oligodendroglioma" & target_cluster == "Astro-like") %>% dplyr::mutate(pair = factor(pair, levels = row_names)) %>% dplyr::mutate(clusters_plot = factor(clusters_plot, levels = col_names))

  p2 <- ggplot2::ggplot(input.oligodendroglioma.astro_like, mapping = ggplot2::aes(x = clusters_plot, y = pair)) +
    ggplot2::geom_point(mapping = ggplot2::aes(size = `-log10(p-value)`, color = `log2 mean(Molecule 1, Molecule 2)`)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Astro-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", color = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(face = "bold")) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_color_gradient(limits = limits, low = color_low, high = color_high) +
    ggplot2::facet_grid(rows = ggplot2::vars(group), scales = "free", space = "free_y")

  input.oligodendroglioma.opc_like <- interactions %>% dplyr::filter(lgg_subtype == "Oligodendroglioma" & target_cluster == "OPC-like") %>% dplyr::mutate(pair = factor(pair, levels = row_names)) %>% dplyr::mutate(clusters_plot = factor(clusters_plot, levels = col_names))

  p1 <- ggplot2::ggplot(input.oligodendroglioma.opc_like, mapping = ggplot2::aes(x = clusters_plot, y = pair)) +
    ggplot2::geom_point(mapping = ggplot2::aes(size = `-log10(p-value)`, color = `log2 mean(Molecule 1, Molecule 2)`)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("OPC-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(face = "bold", color = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
                   plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, color = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(face = "bold"),
                   axis.title.y = ggplot2::element_text(face = "bold")) +
    ggplot2::ylab("Gene 1 | Gene 2") +
    ggplot2::xlab("") +
    ggplot2::scale_color_gradient(limits = limits, low = color_low, high = color_high) +
    ggplot2::facet_grid(rows = ggplot2::vars(group), scales = "free", space = "free_y")

  p1 <- p1 + Seurat::NoLegend() + ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                                                 strip.text.y = ggplot2::element_blank(),
                                                 strip.background = ggplot2::element_blank())
  p2 <- p2 + Seurat::NoLegend() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks") + ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0),
                                                                                                          strip.text.y = ggplot2::element_blank(),
                                                                                                          strip.background = ggplot2::element_blank())
  p3 <- p3 + Seurat::NoLegend() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks") + ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                                                                                                          strip.text.y = ggplot2::element_blank(),
                                                                                                          strip.background = ggplot2::element_blank())
  p4 <- p4 + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks") + ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
                                                                                     strip.text.y = ggplot2::element_text(face = "bold", angle = 0),
                                                                                     strip.background = ggplot2::element_blank())
  layout = "ABCDE"
  p <- patchwork::wrap_plots(A = p1,
                             B = p2,
                             C = p3,
                             D = p4,
                             E = patchwork::guide_area(),
                             design = layout)

  return(p)

}


Figure_1D_1E <- function(type){
  scoreSignature <- function(X.center, X.mean, signature, n = 100, verbose = TRUE) {
    if(verbose) {
      message("Number of cells: ", ncol(X.center))
      message("Number of genes: ", nrow(X.center))
      message("Number of genes in signature: ", length(signature))
      message("...\n\n")
    }
    # Compute the signature score: A score per cell for the top 30 genes in this signature.

    s.score <- colMeans(# Compute the mean per column of the output so that you have 1 score per gene.
      do.call(rbind, # Concatenate the output row by row (gene by gene).
              BiocGenerics::lapply(signature, # Apply this function to each gene in the list.
                                   function(gene) {
                                     ## X.mean[gene] - X.mean => Substracts the median expression of the genes to the expression of the queried gene.
                                     ## abs(X.mean[gene] - X.mean) ==> Have it as absolute values.
                                     ## sort(abs(X.mean[gene] - X.mean)) ==> sort the values in ascending order.
                                     ## names(sort(abs(X.mean[gene] - X.mean)) ==> keep only the names of the sorted genes in ascending order with the lowest difference in expression values.
                                     ## names(sort(abs(X.mean[gene] - X.mean))[2:(n+1)]) Keep only the 2 to 101 genes (100 in totol) as control set (the 1st is the queried gene).

                                     g.n <- names(sort(abs(X.mean[gene] - X.mean))[2:(n+1)])
                                     ## X.center[gene, ] => The expression values of the queried gene in each cell.
                                     ## X.center[g.n, ] => The mean expression values of the control set in each cell.
                                     ## X.center[gene, ] - colMeans(X.center[g.n, ]) => Difference in expression between queried gene and the control set per cell.
                                     X.center[gene, ] - colMeans(X.center[g.n, ])
                                   }
              ) # Therefore, this outputs a vector of the difference in expression values between the gene and the control set per cell.
      ) # Therefore, the outputs a matrix where columns are cells and rows are each of the queried genes.
      # Each row will contain the difference of expression between each queried gene and the control set, per cells.
    ) # Finally, this outputs the mean score of all genes belonging to the signature per cell.
    # In other words, it is the  difference expression values between each gene belonging to the signature
    # and the control set of genes, which is defined as the top 100 most similar genes in expression values between each queried
    # gene and the rest of the genes in the centered expression matrix for all tumor cells, averaged for each cell.

    # s.score = named vector with names = each cell and value = signature score for the cells.
    return(s.score)
  }

  # Pretty much the same contept as scoreSignature but with a twist in the end. Instead of doing colMeans to get the average score for each cell for all of the genes in the signature.
  # rowMeans to get the averaged score for each gene for all of the cells in the dataset provided.
  scoreGene <- function(X.center, X.mean, signature, n = 100, verbose = TRUE) {
    if(verbose) {
      message("Number of cells: ", ncol(X.center))
      message("Number of genes: ", nrow(X.center))
      message("Number of genes in signature: ", length(signature))
      message("...\n\n")
    }
    s.score <- rowMeans(do.call(rbind, lapply(signature, function(gene) {
      g.n <- names(sort(abs(X.mean[gene] - X.mean))[2:(n+1)])
      X.center[gene, ] - colMeans(X.center[g.n, ])
    })))
    return(s.score)
  }

  if (type == "oligodendroglioma"){
    scale.ms <- c("Metasignature 1" = "#9A031E",
                  "Metasignature 2" = "#043362",
                  "Metasignature 3" = "#ECA809",
                  "Metasignature 4" = "#5F0F40",
                  "Excluded" = "grey50")
    scale.ident <- c("IDH_ACB_AD_540" = "#E07A5F",
                     "IDH_ACB_AD_809" = "#F2CC8F",
                     "IDH_ACB_AD_883" = "#DAA82B",
                     "IDH_NCH2111" = "#81B29A",
                     "IDH_NCH536" = "#3E745E",
                     "IDH_NCH6341" = "#5CADC1",
                     "IDH_NCH6702" = "#8B6BB8",
                     "IDH_NCH781" = "#3D405B")
    scale.grade <- c("2" = "#94d2bd",
                     "3" = "#005f73")
  } else if (type == "astrocytoma") {
    scale.ms <- c("Metasignature 1" = "#9A031E",
                  "Metasignature 2" = "#043362",
                  "Metasignature 3" = "#5F0F40",
                  "Excluded" = "grey50")
    scale.ident <- c("IDH_ACB_AD_785" = "#ae2012",
                     "IDH_ACB_AD_832" = "#ca6702",
                     "IDH_ACB_AD_865" = "#ee9b00",
                     "IDH_NCH2018" = "#e9d8a6",
                     "IDH_NCH2157" = "#457b9d",
                     "IDH_NCH2164" = "#1d3557")
    scale.grade <- c("2" = "#94d2bd",
                     "3" = "#005f73")
  }

  if (type == "oligodendroglioma"){
    sample.tumor <- sample.oligo[, !(sample.oligo$New_NMF_labelling %in% c("Astrocytes", "Endothelial", "Excluded", "Microglia", "T-Cells", "Neurons", "Oligodendrocytes", "Pericytes", "Undefined"))]
    sample.tumor <- Seurat::ScaleData(sample.tumor)
    tumor_scaled <- sample.tumor@assays$SCT@scale.data
    # Retrieve ribosomal genes.
    ribosomal_genes <- rownames(tumor_scaled)[grepl("^MT", rownames(tumor_scaled))]
    # Filter out ribosomal genes.
    tumor_scaled <- tumor_scaled[!(rownames(tumor_scaled) %in% ribosomal_genes), ]
    # nmf_top30 <- "/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_merged_oligodendroglioma/snRNAseq/10x_3_prime_v3/NMF_analysis/merged_sample/top30_genes_per_signature/"
    # ind_k <- 10
    # nmf.sign <- read.table(paste0(nmf_top30, "all_samples_nmf", ind_k, "_signatures.txt"), sep = "\t")
    # X.mean <- rowMeans(tumor_scaled)
    # sign.scores <- apply(nmf.sign, 2, function(signature){scoreSignature(X.center = tumor_scaled,
    #                                                                      X.mean = X.mean,
    #                                                                      s = signature)})
    # saveRDS(sign.scores, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/main_figures/Figure_1/Figure_1A/oligodendroglioma_sign_scores.rds")
    sign.scores <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/main_figures/Figure_1/Figure_1A/oligodendroglioma_sign_scores.rds")
    corr <- stats::cor(sign.scores)
    range <- max(abs(corr))

    ind_k <- 10
    meta_k <- 6


    h <- pheatmap::pheatmap(corr,
                            color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                            breaks = seq(-range, range, length.out = 100),
                            show_colnames = F,
                            treeheight_col = 0,
                            cutree_rows = meta_k,
                            cutree_cols = meta_k)

    ann_colors <- colortools::wheel("steelblue", meta_k)
    names(ann_colors) <- BiocGenerics::sapply(as.character(seq(1, meta_k)), function(x){paste0("Metasignature_", x)})

    samp_colors <- colortools::wheel("navyblue", length(unique(sample.tumor$orig.ident)))
    names(samp_colors) <- unique(sample.tumor$orig.ident)

    subtype_colors <- viridis::turbo(2, alpha = 0.75)
    names(subtype_colors) <- c("2", "3")

    colors <- list("Metasignature" = ann_colors,
                   "Origin" = samp_colors,
                   "Grade" = subtype_colors)

    metaclust <- stats::cutree(h$tree_row, k = meta_k)
    metaclust.s <- as.character(BiocGenerics::sapply(metaclust, function(x){paste0("Metasignature_", x)}))
    names(metaclust.s) <- names(metaclust)
    metaclust <- metaclust.s

    metaclust.small <- as.character(metaclust)
    names(metaclust.small) <- names(metaclust)

    origin <- BiocGenerics::sapply(colnames(corr), function(x){substring(x, 1, nchar(x)-2)})
    origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "[.]", "-")})
    origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "-$", "")})

    subgroup <- origin
    grade_2 <- unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident)
    grade_3 <- unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident)

    origin <- factor(origin)

    subgroup[subgroup %in% grade_2] <- 2
    subgroup[subgroup %in% grade_3] <- 3

    h.anno <- pheatmap::pheatmap(corr,
                                 color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                                 breaks = seq(-range, range, length.out = 100),
                                 show_colnames = F,
                                 treeheight_col = 0,
                                 annotation_col = data.frame("Metasignature" = metaclust,
                                                             "Origin" = origin,
                                                             "Grade" = subgroup),
                                 annotation_colors = colors,
                                 cellwidth = 10,
                                 cellheight = 10)

    metaclust_small <- metaclust.small

    corr.small <- stats::cor(sign.scores[, names(metaclust_small)])
    range.small <- max(abs(corr.small))

    # OPTIONAL IF THE HEATMAP BREAKS HAVE TO BE REVISED.
    names_sign <- names(metaclust_small)
    assign_sign <- as.character(metaclust_small)

    assign_sign[assign_sign == "Metasignature_6"] <- "Excluded"
    assign_sign[assign_sign == "Metasignature_5"] <- "Metasignature_6"
    assign_sign[assign_sign == "Metasignature_2"] <- "Metasignature_3"
    assign_sign[assign_sign == "Metasignature_1"] <- "Metasignature_2"
    assign_sign[assign_sign == "Metasignature_3"] <- "Metasignature_1"
    assign_sign[assign_sign == "Metasignature_6"] <- "Metasignature_3"


    assign_sign[names(assign_sign) == "OE0145.IDH_NCH781.8"] <- "Excluded"
    assign_sign <- stringr::str_replace_all(assign_sign, "_", " ")


    names(assign_sign) <- names_sign
    metaclust_small <- as.factor(assign_sign)

    origin_small <- BiocGenerics::sapply(colnames(corr.small), function(x){substring(x, 1, nchar(x)-2)})
    origin_small  <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "OE0145", "")})
    origin_small  <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "[.]", "")})
    origin_small <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "-$", "")})

    subgroup_small <- origin_small
    grade_2 <- stringr::str_replace_all(unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident), "OE0145-", "")
    grade_3 <-  stringr::str_replace_all(unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident), "OE0145-", "")

    origin_small <- factor(origin_small)

    subgroup_small[subgroup_small %in% grade_2] <- 2
    subgroup_small[subgroup_small %in% grade_3] <- 3

    subgroup_small <- factor(subgroup_small)

    colors_small <- list("Metasignature" = scale.ms,
                         "Sample" = scale.ident,
                         "Grade" = scale.grade)

    assign_sign[names(assign_sign) == "OE0145.IDH_NCH781.8"] <- "Excluded"
    assign_sign <- stringr::str_replace_all(assign_sign, "_", " ")


    names(assign_sign) <- names_sign
    metaclust_small <- as.factor(assign_sign)

    h.small <- pheatmap::pheatmap(corr.small,
                                  color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                                  breaks = seq(-range.small, range.small, length.out = 100),
                                  show_colnames = F,
                                  show_rownames = F,
                                  treeheight_col = 0,
                                  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                                  legend_labels = c("-1\n", "-0.5", "0", "0.5", "\n1"),
                                  treeheight_row = 0,
                                  fontsize_row = 8,
                                  fontsize_col = 6,
                                  annotation_col = data.frame("Metasignature" = metaclust_small,
                                                              "Sample" = origin_small,
                                                              "Grade" = subgroup_small),
                                  annotation_colors = colors_small,
                                  cellwidth = 5,
                                  cellheight = 5,
                                  fontsize = 12)



    return(h.small)

  } else if (type == "astrocytoma"){
    sample.tumor <- sample.astro[, !(sample.astro$New_NMF_labelling %in% c("Astrocytes", "Endothelial", "Excluded", "Microglia", "T-Cells", "Neurons", "Oligodendrocytes", "Pericytes", "Undefined"))]
    sample.tumor <- Seurat::ScaleData(sample.tumor)
    tumor_scaled <- sample.tumor@assays$SCT@scale.data
    # Retrieve ribosomal genes.
    ribosomal_genes <- rownames(tumor_scaled)[grepl("^MT", rownames(tumor_scaled))]
    # Filter out ribosomal genes.
    tumor_scaled <- tumor_scaled[!(rownames(tumor_scaled) %in% ribosomal_genes), ]

    # nmf_top30 <- "/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/NMF_analysis/merged_sample/top30_genes_per_signature/"
    # ind_k <- 10
    # nmf.sign <- read.table(paste0(nmf_top30, "all_samples_nmf", ind_k, "_signatures.txt"), sep = "\t")
    # X.mean <- rowMeans(tumor_scaled)
    # sign.scores <- apply(nmf.sign, 2, function(signature){scoreSignature(X.center = tumor_scaled,
    #                                                                      X.mean = X.mean,
    #                                                                      s = signature)})
    # saveRDS(sign.scores, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/main_figures/Figure_1/Figure_1A/astrocytomas_sign_scores.rds")
    sign.scores <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/main_figures/Figure_1/Figure_1A/astrocytomas_sign_scores.rds")

    ind_k <- 10
    meta_k <- 4

    corr <- stats::cor(sign.scores)
    range <- max(abs(corr))

    h <- pheatmap::pheatmap(corr,
                            color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                            breaks = seq(-range, range, length.out = 100),
                            show_colnames = F,
                            treeheight_col = 0,
                            cutree_rows = meta_k,
                            cutree_cols = meta_k)

    ann_colors <- colortools::wheel("steelblue", meta_k)
    names(ann_colors) <- BiocGenerics::sapply(as.character(seq(1, meta_k)), function(x){paste0("Metasignature_", x)})

    samp_colors <- colortools::wheel("navyblue", length(unique(sample.tumor$orig.ident)))
    names(samp_colors) <- unique(sample.tumor$orig.ident)

    subtype_colors <- viridis::turbo(2, alpha = 0.75)
    names(subtype_colors) <- c("2", "3")

    colors <- list("Metasignature" = ann_colors,
                   "Origin" = samp_colors,
                   "Grade" = subtype_colors)

    metaclust <- stats::cutree(h$tree_row, k = meta_k)
    metaclust.s <- as.character(BiocGenerics::sapply(metaclust, function(x){paste0("Metasignature_", x)}))
    names(metaclust.s) <- names(metaclust)
    metaclust <- metaclust.s

    metaclust.small <- as.character(metaclust)
    names(metaclust.small) <- names(metaclust)

    origin <- BiocGenerics::sapply(colnames(corr), function(x){substring(x, 1, nchar(x)-2)})
    origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "[.]", "-")})
    origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "-$", "")})

    subgroup <- origin
    grade_2 <- unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident)
    grade_3 <- unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident)

    origin <- factor(origin)

    subgroup[subgroup %in% grade_2] <- 2
    subgroup[subgroup %in% grade_3] <- 3

    h.anno <- pheatmap::pheatmap(corr,
                                 color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                                 breaks = seq(-range, range, length.out = 100),
                                 show_colnames = F,
                                 treeheight_col = 0,
                                 annotation_col = data.frame("Metasignature" = metaclust,
                                                             "Origin" = origin,
                                                             "Grade" = subgroup),
                                 annotation_colors = colors,
                                 cellwidth = 10,
                                 cellheight = 10)

    metaclust_small <- metaclust.small

    corr.small <- stats::cor(sign.scores[, names(metaclust_small)])
    range.small <- max(abs(corr.small))

    # OPTIONAL IF THE HEATMAP BREAKS HAVE TO BE REVISED.
    names_sign <- names(metaclust_small)
    assign_sign <- as.character(metaclust_small)

    assign_sign[assign_sign == "Metasignature_1"] <- "Excluded"
    assign_sign[assign_sign == "Metasignature_2"] <- "Metasignature_1"
    assign_sign[assign_sign == "Metasignature_3"] <- "Metasignature_2"
    assign_sign[assign_sign == "Metasignature_4"] <- "Metasignature_3"


    assign_sign <- stringr::str_replace_all(assign_sign, "_", " ")
    names(assign_sign) <- names_sign
    metaclust_small <- as.factor(assign_sign)

    origin_small <- BiocGenerics::sapply(colnames(corr.small), function(x){substring(x, 1, nchar(x)-2)})
    origin_small  <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "OE0145", "")})
    origin_small  <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "[.]", "")})
    origin_small <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "-$", "")})

    subgroup_small <- origin_small
    grade_2 <- stringr::str_replace_all(unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident), "OE0145-", "")
    grade_3 <-  stringr::str_replace_all(unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident), "OE0145-", "")

    origin_small <- factor(origin_small)

    subgroup_small[subgroup_small %in% grade_2] <- 2
    subgroup_small[subgroup_small %in% grade_3] <- 3

    subgroup_small <- factor(subgroup_small)

    colors_small <- list("Metasignature" = scale.ms,
                         "Sample" = scale.ident,
                         "Grade" = scale.grade)


    h.small <- pheatmap::pheatmap(corr.small,
                                  color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100),
                                  breaks = seq(-range.small, range.small, length.out = 100),
                                  show_colnames = F,
                                  show_rownames = T,
                                  treeheight_col = 0,
                                  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                                  legend_labels = c("-1\n", "-0.5", "0", "0.5", "\n1"),
                                  treeheight_row = 0,
                                  fontsize_row = 8,
                                  fontsize_col = 6,
                                  annotation_col = data.frame("Metasignature" = metaclust_small,
                                                              "Sample" = origin_small,
                                                              "Grade" = subgroup_small),
                                  annotation_colors = colors_small,
                                  cellwidth = 7,
                                  cellheight = 7,
                                  fontsize = 12)
    return(h.small)

  }
}

initial_cluster_assignment <- function(){

  for (type in c("merged", "integrated")){
    if (type == "merged"){
      sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_merged_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_merged_oligodendroglioma_labelled")
      sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_merged_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_merged_astrocytoma_labelled")

    } else if (type == "integrated"){
      sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach")
      sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach")

    }

    suva_folder <- paste0(figures_folder, "/second_iteration/supplementary/suva_markers")
    panglaodb_folder <- paste0(figures_folder, "/second_iteration/supplementary/panglaodb_folder")
    heatmap_folder <- paste0(figures_folder, "/second_iteration/supplementary/heatmap_markers")
    known_markers_folder <- paste0(figures_folder, "/second_iteration/supplementary/known_markers")


    #Known markers:

    features <- c("PTPRC", "VWF", "MBP", "AQP4", "MCAM", "RELN")
    for (feature in features){
      p.oligo <- do_FeaturePlot(sample.oligo, features = feature)
      p.astro <- do_FeaturePlot(sample.astro, features = feature)

      dir.create(paste0(known_markers_folder, "/", feature))
      savefig(plot = p.oligo, figure_path = paste0(known_markers_folder, "/", feature), file_name = paste0(type, "_oligo_", feature), width = 8, heigth = 7)
      savefig(plot = p.astro, figure_path = paste0(known_markers_folder, "/", feature), file_name = paste0(type, "_astro_", feature), width = 8, heigth = 7)
    }



    markers.suva <- lapply(# Omit any NAs in the named list.
      as.list(# Turn tibble to a named list.
        readxl::read_excel("/omics/odcf/analysis/OE0145_projects/idh_gliomas/scripts/main/06_scana/marker_genes/suva_markers.xlsx")),
      function(x) x[!is.na(x)])
    names(markers.suva) <- stringr::str_replace_all(names(markers.suva), " ", "_")
    names(markers.suva) <- stringr::str_replace_all(names(markers.suva), "-", "_")
    names(markers.suva) <- stringr::str_replace_all(names(markers.suva), "/", "_")

    markers.panglaodb <- readr::read_tsv("/omics/odcf/analysis/hipo/hipo_049/ATRT/scripts/data_analysis/06_scana/marker_genes/PanglaoDB_markers_27_Mar_2020.tsv")

    markers.panglaodb <- as.list(markers.panglaodb %>%
                                   dplyr::filter(organ == "Brain") %>%
                                   dplyr::filter(species %in% c("Mm Hs", "Hs")) %>%
                                   dplyr::select(`official gene symbol`, `cell type`) %>%
                                   tidyr::pivot_wider(names_from = `cell type`,
                                                      values_from = `official gene symbol`))
    names(markers.panglaodb) <- stringr::str_replace_all(names(markers.panglaodb), " ", "_")
    names(markers.panglaodb) <- stringr::str_replace_all(names(markers.panglaodb), "-", "_")
    names(markers.panglaodb) <- stringr::str_replace_all(names(markers.panglaodb), "/", "_")

    names.suva <- names(markers.suva)
    names.panglaodb <- names(markers.panglaodb)

    markers <- c(markers.suva, markers.panglaodb)

    for (gene_list in sort(names(markers))){
      message(gene_list)

      if (gene_list %in% names.suva){
        #outfolder <- paste0(suva_folder, "/", gene_list)
        genes <- markers[[gene_list]]
      } else if (gene_list %in% names.panglaodb){
        #outfolder <- paste0(panglaodb_folder, "/", gene_list)
        genes <- markers[[gene_list]][[1]]
      }

      sample.oligo <- Seurat::AddModuleScore(sample.oligo, list(genes), name = gene_list)
      sample.astro <- Seurat::AddModuleScore(sample.astro, list(genes), name = gene_list)

      gene.query <- paste0(gene_list, "1")
      p.oligo <- do_FeaturePlot(sample.oligo, gene.query, plot.title = gene_list)
      p.astro <- do_FeaturePlot(sample.astro, gene.query, plot.title = gene_list)

      #savefig(plot = p.oligo, figure_path = outfolder, file_name = paste0(type, "_oligo_", gene_list), width = 8, heigth = 7)
      #savefig(plot = p.astro, figure_path = outfolder, file_name = paste0(type, "_astro_", gene_list), width = 8, heigth = 7)
    }

    seurat_scoring.oligo <- data.frame("rownames" = unique(sample.oligo$first_labelling))
    seurat_scoring.astro <- data.frame("rownames" = unique(sample.astro$first_labelling))



    # Iterate over each list of marker genes.
    for (celltype in names(markers)){
      # Get the names of the metadata columns for AUCell and Seurat scorings.
      col_name_seurat_scoring <- paste0(celltype, "1")

      # Generate empty vectors for the aggregated scores.
      list_score_seurat.oligo <- c()
      list_score_seurat.astro <- c()


      # Iterate over each cluster.
      for (cluster_name in unique(sample.oligo$first_labelling)){
        # Retrieve which cells are assigned to the cluster.
        scores_seurat <- sample.oligo@meta.data[sample.oligo@meta.data$first_labelling == cluster_name, col_name_seurat_scoring]

        # Append to the vector the mean for each cell type.
        list_score_seurat.oligo <- append(list_score_seurat.oligo, mean(scores_seurat))
      }
      for (cluster_name in unique(sample.astro$first_labelling)){
        # Retrieve which cells are assigned to the cluster.
        scores_seurat <- sample.astro@meta.data[sample.astro@meta.data$first_labelling == cluster_name, col_name_seurat_scoring]

        # Append to the vector the mean for each cell type.
        list_score_seurat.astro <- append(list_score_seurat.astro, mean(scores_seurat))
      }
      # Get the name of the column together with the number of genes used for the enrichment scoring.
      seurat_scoring.oligo[celltype] <- list_score_seurat.oligo
      seurat_scoring.astro[celltype] <- list_score_seurat.astro
    }

    # Remove the rownames column in the object, since you set them to be the rownames of the dataframe.
    rownames(seurat_scoring.oligo) <- seurat_scoring.oligo$rownames
    seurat_scoring.oligo$rownames <- NULL
    rownames(seurat_scoring.astro) <- seurat_scoring.astro$rownames
    seurat_scoring.astro$rownames <- NULL

    # Transform the data frames into a Matrix object.
    seurat_scoring.oligo <- as.matrix(seurat_scoring.oligo)
    seurat_scoring.astro <- as.matrix(seurat_scoring.astro)

    dir.create(heatmap_folder, recursive = T)
    saveRDS(seurat_scoring.oligo, paste0(heatmap_folder, "/", type, "_oligo_matrix_marker_lists.rds"))
    saveRDS(seurat_scoring.astro, paste0(heatmap_folder, "/", type, "_astro_matrix_marker_lists.rds"))

    range.oligo <- max(abs(seurat_scoring.oligo))
    range.astro <- max(abs(seurat_scoring.astro))

    # Taken from: https://github.com/raivokolde/pheatmap/issues/48#issue-402653138
    # use this function to make row or column names bold
    # parameters:
    #   mat: the matrix passed to pheatmap
    #   rc_fun: either rownames or colnames
    #   rc_names: vector of names that should appear in boldface
    make_bold_names <- function(mat, rc_fun, rc_names) {
      bold_names <- rc_fun(mat)
      ids <- rc_names %>% match(rc_fun(mat))
      ids %>%
        purrr::walk(
          function(i)
            bold_names[i] <<-
            bquote(bold(.(rc_fun(mat)[i]))) %>%
            as.expression()
        )
      bold_names
    }

    h.oligo <- pheatmap::pheatmap(seurat_scoring.oligo,
                                  color = viridis::viridis(100),
                                  breaks = seq(0, range.oligo, length.out = 100),
                                  cellwidth = 15,
                                  cellheight = 15,
                                  angle_col = "90",
                                  labels_row = make_bold_names(seurat_scoring.oligo, rownames, rownames(seurat_scoring.oligo)),
                                  labels_col = make_bold_names(seurat_scoring.oligo, colnames, colnames(seurat_scoring.oligo)))
    h.astro <- pheatmap::pheatmap(seurat_scoring.astro,
                                  color = viridis::viridis(100),
                                  breaks = seq(0, range.astro, length.out = 100),
                                  cellwidth = 15,
                                  cellheight = 15,
                                  angle_col = "90",
                                  labels_row = make_bold_names(seurat_scoring.astro, rownames, rownames(seurat_scoring.astro)),
                                  labels_col = make_bold_names(seurat_scoring.astro, colnames, colnames(seurat_scoring.astro)))

    savefig(plot = h.oligo, figure_path = heatmap_folder, file_name = paste0(type, "_oligo_heatmap_marker_lists"))
    savefig(plot = h.astro, figure_path = heatmap_folder, file_name = paste0(type, "_astro_heatmap_marker_lists"))
  }


}

diffusion_maps_analysis <- function(){
  library(destiny)
  exp_name <- "difussion_map_object_with_atrt_all_using_variable_genes"
  diffusion_folder <- paste0(figures_folder, "/second_iteration/main_figures/Figure_2/diffusion/", exp_name)
  dir.create(diffusion_folder, recursive = T)
  name_oligo <- "OE0145-IDH_integrated_oligodendroglioma"
  name_astro <- "OE0145-IDH_integrated_astrocytoma"

  sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach")
  sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach")

  # sample.oligo$New_NMF_labelling[sample.oligo$New_NMF_labelling == "OPC-like"] <- "OPC_MS"
  # sample.oligo$New_NMF_labelling[sample.oligo$New_NMF_labelling == "Astro-like"] <- "Astro_MS"
  # sample.oligo$New_NMF_labelling[sample.oligo$New_NMF_labelling == "Cycling-like"] <- "Cycling_MS"
  # sample.oligo$New_NMF_labelling[sample.oligo$New_NMF_labelling == "RA"] <- "RA_MS"
  Idents(sample.oligo) <- sample.oligo$New_NMF_labelling

  # sample.astro$New_NMF_labelling[sample.astro$New_NMF_labelling == "OPC-like"] <- "OPC_MS"
  # sample.astro$New_NMF_labelling[sample.astro$New_NMF_labelling == "Astro-like"] <- "Astro_MS"
  # sample.astro$New_NMF_labelling[sample.astro$New_NMF_labelling == "Cycling-like"] <- "Cycling_MS"
  # sample.astro$New_NMF_labelling[sample.astro$New_NMF_labelling == "RA"] <- "RA_MS"
  Idents(sample.astro) <- sample.astro$New_NMF_labelling

  cluster_cols <- c("RA" = "#ECA809",
                    "RA_MS" = "#ECA809",                    # Marigold.
                    "OPC-like" = "#043362",                    # Prussian Blue.
                    "Oligo_program" = "#043362",                    # Prussian Blue.
                    "OPC_MS" = "#043362",                    # Prussian Blue.
                    "T-Cells" = "#009FF5",                     # Carolina Blue.
                    "Oligodendrocytes" = "#BC5210",            # Burnt Orange.
                    "Astrocytes" = "#279185",                  # Celadon Green.
                    "Microglia" = "#7EB356",                   # Bud Green.
                    "Pericytes" = "#AC70FF",                   # Medium Purple.
                    "Undefined" = "#63412C",                   # Van Dyke Brown.
                    "Gradient" = "#D6D6D6",                    # Light grey.
                    "Neurons" = "#544B81",                     # Purple Navy.
                    "Endothelial" = "#da627d",                 # Blush.
                    "Astro-like" = "#9A031E",                  # Ruby Red.
                    "Astro_MS" = "#9A031E",                  # Ruby Red.
                    "Astro_program" = "#9A031E",                  # Ruby Red.
                    "Excluded" = "#4D6880",                    # Dark Electric Blue.
                    "Cycling-like" = "#5F0F40",
                    "Cycling_MS" = "#5F0F40")                # Tyrian Purple.
  cluster_cols <- cluster_cols[sort(names(cluster_cols))]

  diff.map.oligo <- readRDS(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/diffusion_maps/OE0145-IDH_integrated_oligodendroglioma_difussion_map_object_with_suva_markers_only_stem"))
  diff.map.astro <- readRDS(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/diffusion_maps/OE0145-IDH_integrated_astrocytoma_difussion_map_object_with_suva_markers_only_stem"))

  sample.tumor.oligo <- sample.oligo[, sample.oligo$New_NMF_labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA", "Cycling-like")]
  sample.tumor.astro <- sample.astro[, sample.astro$New_NMF_labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA", "Cycling-like")]

  rm(sample.oligo)
  rm(sample.astro)

  sample.tumor.oligo@reductions$diffusion <- Seurat::CreateDimReducObject(embeddings = diff.map.oligo@eigenvectors, key = "DC", assay = "SCT")
  sample.tumor.astro@reductions$diffusion <- Seurat::CreateDimReducObject(embeddings = diff.map.astro@eigenvectors, key = "DC", assay = "SCT")

  dc_choice <- 1
  dc_use <- paste0("DC", dc_choice)
  selected.oligo <- names(diff.map.oligo@eigenvectors[, dc_use])
  selected.astro <- names(diff.map.astro@eigenvectors[, dc_use])


  sample.tumor.oligo <- sample.tumor.oligo[, selected.oligo]
  sample.tumor.astro <- sample.tumor.astro[, selected.astro]

  p.dif.oligo <- do_DimPlot(sample = sample.tumor.oligo, reduction = "diffusion", colors.use = cluster_cols[names(cluster_cols) %in% unique(sample.tumor.oligo$New_NMF_labelling)], label = F, plot.title = "Oligodendroglioma")
  p.dif.astro <- do_DimPlot(sample = sample.tumor.astro, reduction = "diffusion", colors.use = cluster_cols[names(cluster_cols) %in% unique(sample.tumor.astro$New_NMF_labelling)], label = F, plot.title = "Astrocytoma")

  # color scale
  colors.use <- cluster_cols[names(cluster_cols) %in% unique(sample.tumor.oligo$New_NMF_labelling)]
  colors.use["Gradient"] <- "black"
  p.cluster.oligo <- do_DimPlot(sample = sample.tumor.oligo, reduction = "diffusion", split.by = "New_NMF_labelling", cols.split = colors.use, label = F, legend = F, ncol = 2)
  p.cluster.astro <- do_DimPlot(sample = sample.tumor.astro, reduction = "diffusion", split.by = "New_NMF_labelling", cols.split = colors.use, label = F, legend = F, ncol = 2)

  diffusion_folder <- paste0(diffusion_folder, "/", "DC_", dc_choice)
  dir.create(diffusion_folder, recursive = T)

  savefig(p.dif.oligo, figure_path = paste0(diffusion_folder, "/", name_oligo, "/UMAP_by_cluster"), file_name = paste0(name_oligo, "_UMAP_by_cluster"), width = 8, heigth = 7)
  savefig(p.dif.astro, figure_path = paste0(diffusion_folder, "/", name_astro, "/UMAP_by_cluster" ), file_name = paste0(name_oligo, "_UMAP_by_cluster"), width = 8, heigth = 7)
  savefig(p.cluster.oligo, figure_path = paste0(diffusion_folder, "/", name_oligo, "/UMAP_split_by_cluster"), file_name = paste0(name_oligo, "_UMAP_split_by_cluster"), width = 8, heigth = 12)
  savefig(p.cluster.astro, figure_path = paste0(diffusion_folder, "/", name_astro, "/UMAP_split_by_cluster"), file_name = paste0(name_oligo, "_UMAP__split_by_cluster"), width = 8, heigth = 12)

  # markers <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/final_marker_genes/IDH_gliomas_final_markers.rds")
  # markers.use <- markers[c("TB_OPC-like", "TB_Astro-like", "TB_Cycling-like", "TB_Underlying", "A_TB_Gradient", "O_TB_Gradient")]
  # names(markers.use) <- stringr::str_replace_all(names(markers.use), "-", ".")
  #
  markers.suva <- lapply(# Omit any NAs in the named list.
    as.list(# Turn tibble to a named list.
      readxl::read_excel("/omics/odcf/analysis/OE0145_projects/idh_gliomas/scripts/main/06_scana/marker_genes/suva_markers.xlsx")),
    function(x) x[!is.na(x)])

  names(markers.suva) <- stringr::str_replace_all(names(markers.suva), "-", ".")

  for (gene_list in "Stemness_program"){
    message(gene_list)
    genes <- markers.suva[[gene_list]]
    sample.tumor.oligo <- Seurat::AddModuleScore(sample.tumor.oligo, list(genes), name = gene_list)
    sample.tumor.astro <- Seurat::AddModuleScore(sample.tumor.astro, list(genes), name = gene_list)
    p.oligo <- do_FeaturePlot(sample.tumor.oligo, reduction = "diffusion", features = paste0(gene_list, "1"), plot.title = gene_list)
    p.astro <- do_FeaturePlot(sample.tumor.astro, reduction = "diffusion", features = paste0(gene_list, "1"), plot.title = gene_list)
  }


  savefig(p.oligo, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment"), width = 10, heigth = 10)
  savefig(p.astro, figure_path = paste0(diffusion_folder, "/", name_astro, "/Stemness_enrichment" ), file_name = paste0(name_oligo, "_Stemness_enrichment"), width = 10, heigth = 10)

  sample.tumor.oligo$New_NMF_labelling <- factor(sample.tumor.oligo$New_NMF_labelling, levels = rev(c("Astro-like", "Cycling-like", "Gradient", "OPC-like", "RA")))
  sample.tumor.astro$New_NMF_labelling <- factor(sample.tumor.astro$New_NMF_labelling, levels = rev(c("Astro-like", "Cycling-like", "Gradient", "OPC-like", "RA")))
  Seurat::Idents(sample.tumor.oligo) <- sample.tumor.oligo$New_NMF_labelling
  Seurat::Idents(sample.tumor.astro) <- sample.tumor.astro$New_NMF_labelling

  ridge.oligo <- Seurat::RidgePlot(sample.tumor.oligo, features = "Stemness_program1", cols = colors.use) + ggtitle("Oligodendroglioma") + Seurat::NoLegend() + ggpubr::rremove("y.title")
  ridge.astro <- Seurat::RidgePlot(sample.tumor.astro, features = "Stemness_program1", cols = colors.use) + ggtitle("Astrocytoma") + Seurat::NoLegend() + ggpubr::rremove("y.title")

  savefig(ridge.oligo, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ridge_plot"), width = 10, heigth = 10)
  savefig(ridge.astro, figure_path = paste0(diffusion_folder, "/", name_astro, "/Stemness_enrichment" ), file_name = paste0(name_oligo, "_Stemness_enrichment_ridge_plot"), width = 10, heigth = 10)


  # DotPlots.
  sample.tumor.oligo$pseudotime_diffmap <- rank(destiny::eigenvectors(diff.map.oligo)[,dc_choice])
  sample.tumor.astro$pseudotime_diffmap <- rank(destiny::eigenvectors(diff.map.astro)[,dc_choice])


  sample.tumor.oligo$New_NMF_labelling <- factor(sample.tumor.oligo$New_NMF_labelling, levels = c("Astro-like", "Cycling-like", "Gradient", "OPC-like", "RA"))
  sample.tumor.astro$New_NMF_labelling <- factor(sample.tumor.astro$New_NMF_labelling, levels = c("Astro-like", "Cycling-like", "Gradient", "OPC-like", "RA"))
  Seurat::Idents(sample.tumor.oligo) <- sample.tumor.oligo$New_NMF_labelling
  Seurat::Idents(sample.tumor.astro) <- sample.tumor.astro$New_NMF_labelling


  p.dotplot.oligo <- ggplot2::ggplot(sample.tumor.oligo@meta.data, mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(New_NMF_labelling), color = New_NMF_labelling)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = cluster_cols[names(cluster_cols) %in% unique(sample.tumor.oligo$New_NMF_labelling)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("DC 1") +
    ggplot2::ggtitle("Oligodendroglioma") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 20,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 20,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 20,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.astro <- ggplot2::ggplot(sample.tumor.astro@meta.data, mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(New_NMF_labelling), color = New_NMF_labelling)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = cluster_cols[names(cluster_cols) %in% unique(sample.tumor.astro$New_NMF_labelling)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::xlab("DC_1") +
    ggplot2::ylab("") +
    ggplot2::ggtitle("Astrocytoma") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 20,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 20,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 20,
                                                      face = "bold",
                                                      hjust = 0.5))
  patch <- p.dotplot.oligo | p.dotplot.astro

  savefig(patch, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering"), width = 11, heigth = 10)
  savefig(p.dotplot.oligo, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_oligo"), width = 8, heigth = 7)
  savefig(p.dotplot.astro, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_astro"), width = 8, heigth = 7)


  # Can we do the same for the enrichment scores?

  sample.tumor.oligo$stemness_rank <- rank(sample.tumor.oligo$Stemness_program1)
  sample.tumor.astro$stemness_rank <- rank(sample.tumor.astro$Stemness_program1)

  p.dotplot.oligo.stemness <- ggplot2::ggplot(sample.tumor.oligo@meta.data, mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(New_NMF_labelling), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr(legend = "right") +
    #ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("Stemness signature Suva") +
    ggplot2::ggtitle("Oligodendroglioma") +
    ggpubr::rremove("legend.title") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 15,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 15,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 15,
                                                      face = "bold",
                                                      hjust = 0.5),
                   legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1))

  p.dotplot.astro.stemness <- ggplot2::ggplot(sample.tumor.astro@meta.data, mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(New_NMF_labelling), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr(legend = "right") +
    #ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("Stemness signature Suva") +
    ggplot2::ggtitle("Astrocytoma") +
    ggpubr::rremove("legend.title") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 15,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 15,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 15,
                                                      face = "bold",
                                                      hjust = 0.5),
                   legend.text = ggplot2::element_text(size = 10, face = "bold", hjust = 1))
  patch <- p.dotplot.oligo.stemness | p.dotplot.astro.stemness
  savefig(p.dotplot.oligo.stemness, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_oligo_by_enrichment_scores"), width = 8, heigth = 7)
  savefig(p.dotplot.astro.stemness, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_astro_by_enrichment_scores"), width = 8, heigth = 7)
  savefig(patch, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_by_enrichment_scores"), width = 10, heigth = 7)

  # And now, for each sample of origin for RA.
  scale.ident <- c("IDH_ACB_AD_540" = "#E07A5F",
                   "IDH_ACB_AD_809" = "#F2CC8F",
                   "IDH_ACB_AD_883" = "#DAA82B",
                   "IDH_NCH2111" = "#81B29A",
                   "IDH_NCH536" = "#3E745E",
                   "IDH_NCH6341" = "#5CADC1",
                   "IDH_NCH6702" = "#8B6BB8",
                   "IDH_NCH781" = "#3D405B",
                   "IDH_ACB_AD_785" = "#ae2012",
                   "IDH_ACB_AD_832" = "#ca6702",
                   "IDH_ACB_AD_865" = "#ee9b00",
                   "IDH_NCH2018" = "#e9d8a6",
                   "IDH_NCH2157" = "#457b9d",
                   "IDH_NCH2164" = "#1d3557")
  sample.tumor.oligo$orig.ident <- stringr::str_remove_all(sample.tumor.oligo$orig.ident, "OE0145-")
  sample.tumor.astro$orig.ident <- stringr::str_remove_all(sample.tumor.astro$orig.ident, "OE0145-")


  sample.tumor.oligo$orig.ident <- factor(sample.tumor.oligo$orig.ident, levels = c("IDH_ACB_AD_540", "IDH_ACB_AD_809", "IDH_ACB_AD_883",  "IDH_NCH2111", "IDH_NCH536", "IDH_NCH6341", "IDH_NCH6702", "IDH_NCH781"))
  sample.tumor.astro$orig.ident <- factor(sample.tumor.astro$orig.ident, levels = c( "IDH_ACB_AD_785", "IDH_ACB_AD_832", "IDH_ACB_AD_865", "IDH_NCH2018", "IDH_NCH2157", "IDH_NCH2164"))


  p.dotplot.oligo.opc <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "OPC-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.oligo$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggplot2::ylab("Oligodendroglioma") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("OPC-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 18,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.oligo.astro <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "Astro-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.oligo$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("Astro-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.oligo.ra <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "RA"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.oligo$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("RA") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.oligo.cycling <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "Cycling-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.oligo$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("Cycling-like") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))


  p.dotplot.oligo.stemness.opc <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "OPC-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("OPC-like") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.oligo.stemness.astro <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "Astro-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("Astro-like") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))
  p.dotplot.oligo.stemness.cycling <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "Cycling-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("Cycling-like") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))

  p.dotplot.oligo.stemness.ra <- ggplot2::ggplot(sample.tumor.oligo@meta.data %>% dplyr::filter(New_NMF_labelling == "RA"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("RA") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))
















  p.dotplot.astro.opc <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "OPC-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.astro$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggplot2::ylab("Astrocytoma") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 18,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.astro.astro <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "Astro-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.astro$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.astro.ra <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "RA"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.astro$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.astro.cycling <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "Cycling-like"), mapping = ggplot2::aes(x = pseudotime_diffmap, y = forcats::fct_rev(orig.ident), color = orig.ident)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    ggplot2::scale_color_manual(values = scale.ident[names(scale.ident) %in% unique(sample.tumor.astro$orig.ident)]) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14,
                                                     face = "bold"),
                   axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))


  p.dotplot.astro.stemness.opc <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "OPC-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14,
                                                      face = "bold"),
                   plot.title = ggplot2::element_text(size = 14,
                                                      face = "bold",
                                                      hjust = 0.5))

  p.dotplot.astro.stemness.astro <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "Astro-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))
  p.dotplot.astro.stemness.cycling <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "Cycling-like"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))

  p.dotplot.astro.stemness.ra <- ggplot2::ggplot(sample.tumor.astro@meta.data %>% dplyr::filter(New_NMF_labelling == "RA"), mapping = ggplot2::aes(x = stemness_rank, y = forcats::fct_rev(orig.ident), color = Stemness_program1)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
    viridis::scale_color_viridis() +
    ggpubr::theme_pubr() +
    ggpubr::rremove("legend") +
    ggpubr::rremove("x.ticks") +
    ggpubr::rremove("x.text") +
    ggpubr::rremove("y.ticks") +
    ggpubr::rremove("y.text") +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::ggtitle("") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14,
                                         face = "bold"),
      plot.title = ggplot2::element_text(size = 14,
                                         face = "bold",
                                         hjust = 0.5))



  p.oligo <- (p.dotplot.oligo.opc | p.dotplot.oligo.stemness.opc | p.dotplot.oligo.astro | p.dotplot.oligo.stemness.astro | p.dotplot.oligo.cycling | p.dotplot.oligo.stemness.cycling | p.dotplot.oligo.ra | p.dotplot.oligo.stemness.ra)
  p.astro <- (p.dotplot.astro.opc | p.dotplot.astro.stemness.opc | p.dotplot.astro.astro | p.dotplot.astro.stemness.astro | p.dotplot.astro.cycling | p.dotplot.astro.stemness.cycling | p.dotplot.astro.ra | p.dotplot.astro.stemness.ra)


  layout <- "A
    B"
  final.patch <- patchwork::wrap_plots(A = p.oligo,
                                       B = p.astro,
                                       design = layout)



  savefig(final.patch, figure_path = paste0(diffusion_folder, "/", name_oligo, "/Stemness_enrichment"), file_name = paste0(name_oligo, "_Stemness_enrichment_ordering_by_sample_and_cell_identity"), width = 15, heigth = 8)

  # list.nebulosa.oligo <- list()
  # list.nebulosa.astro <- list()
  # for (identity in sort(unique(sample.tumor.oligo$New_NMF_labelling))){
  #     message(identity)
  #     sample.tumor.oligo$dummy <- ifelse(sample.tumor.oligo$New_NMF_labelling == identity, TRUE, FALSE)
  #     sample.tumor.astro$dummy <- ifelse(sample.tumor.astro$New_NMF_labelling == identity, TRUE, FALSE)
  #     p.oligo <- do_Nebulosa_plot(sample.tumor.oligo, reduction = "diffusion", features = "dummy", plot.title = identity)
  #     p.astro <- do_Nebulosa_plot(sample.tumor.astro, reduction = "diffusion", features = "dummy", plot.title = identity)
  #     list.nebulosa.oligo[[identity]] <- p.oligo
  #     list.nebulosa.astro[[identity]] <- p.astro
  #     sample.tumor.oligo$dummy <- NULL
  #     sample.tumor.astro$dummy <- NULL
  # }
  #
  # p.oligo <- patchwork::wrap_plots(list.nebulosa.oligo, ncol = 2)
  # p.astro <- patchwork::wrap_plots(list.nebulosa.astro, ncol = 2)



  # list.invididual_plots.oligo <- list()
  # list.invididual_plots.astro <- list()
  # for (gene_list in sort(names(markers.use))){
  #     message(gene_list)
  #     genes <- markers.use[[gene_list]]
  #     p.oligo <- do_FeaturePlot(sample.tumor.oligo, reduction = "diffusion", features = genes[1:5], plot.title = paste0("Top5 genes: ", gene_list), ncol = 5)
  #     p.astro <- do_FeaturePlot(sample.tumor.astro, reduction = "diffusion", features = genes[1:5], plot.title = paste0("Top5 genes: ", gene_list), ncol = 5)
  #     list.invididual_plots.oligo[[gene_list]] <- p.oligo
  #     list.invididual_plots.astro[[gene_list]] <- p.astro
  # }
  #
  # p.oligo <- list.invididual_plots.oligo$TB_Astro.like / list.invididual_plots.oligo$TB_OPC.like / list.invididual_plots.oligo$TB_Cycling.like / list.invididual_plots.oligo$TB_Underlying
  # p.astro <- list.invididual_plots.astro$TB_Astro.like / list.invididual_plots.astro$TB_OPC.like / list.invididual_plots.astro$TB_Cycling.like / list.invididual_plots.astro$TB_Underlying
  #
  #
  #
  #
  # # Plot Eigenvectors.
  # plot(destiny::eigenvalues(diff.map.oligo), ylim = 0:1, pch = 20,
  #      xlab = "Diffusion component (DC)", ylab = "Eigenvalue")
  #
  # plot(destiny::eigenvalues(diff.map.astro), ylim = 0:1, pch = 20,
  #      xlab = "Diffusion component (DC)", ylab = "Eigenvalue")
  #
  #
  # # Generate new clusters.
  # clustering_algorithm <- 1
  # seurat_version <- "v3"
  # clustering_method = "matrix"
  # dims <- 12
  #
  # sample.tumor.oligo <- sample.tumor.oligo %>%
  #     Seurat::FindNeighbors(dims = 1:dims,
  #                           nn.method = ifelse(seurat_version == "v3", "rann", "annoy"),
  #                           reduction = "diffusion") %>%
  #     Seurat::FindClusters(resolution = 0.1,
  #                          random.seed = 42,
  #                          algorithm = clustering_algorithm,
  #                          method = clustering_method)
  #
  # sample.tumor.astro <- sample.tumor.astro %>%
  #     Seurat::FindNeighbors(dims = 1:dims,
  #                           nn.method = ifelse(seurat_version == "v3", "rann", "annoy"),
  #                           reduction = "diffusion") %>%
  #     Seurat::FindClusters(resolution = 0.1,
  #                          random.seed = 42,
  #                          algorithm = clustering_algorithm,
  #                          method = clustering_method)
  #
  # markers.oligo <- Seurat::FindAllMarkers(sample.tumor.oligo, assay = "SCT", only.pos = T, verbose = T,
  #                                         base = ifelse(seurat_version == "v3", exp(1), 2))
  #
  # markers.astro <- Seurat::FindAllMarkers(sample.tumor.astro, assay = "SCT", only.pos = T, verbose = T,
  #                                         base = ifelse(seurat_version == "v3", exp(1), 2))
  #
  # top10_oligo <- markers.oligo %>%
  #     dplyr::group_by(.data$cluster) %>%
  #     dplyr::top_n(n = 10, wt = .data$avg_logFC)
  #
  # top10_astro <- markers.astro %>%
  #     dplyr::group_by(.data$cluster) %>%
  #     dplyr::top_n(n = 10, wt = .data$avg_logFC)
  #
  # cells.to.plot.oligo <- c()
  # for (cluster in levels(sample.tumor.oligo)){
  #     print(cluster)
  #     len <- ifelse(sum(sample.tumor.oligo$seurat_clusters == cluster) >= 200, 200, sum(sample.tumor.oligo$seurat_clusters == cluster))
  #     cells <- colnames(sample.tumor.oligo)[(sample.tumor.oligo$seurat_clusters == cluster)][1:len]
  #     cells.to.plot.oligo <- c(cells.to.plot.oligo, cells)
  #
  # }
  #
  # cells.to.plot.astro <- c()
  # for (cluster in levels(sample.tumor.astro)){
  #     print(cluster)
  #     len <- ifelse(sum(sample.tumor.astro$seurat_clusters == cluster) >= 200, 200, sum(sample.tumor.astro$seurat_clusters == cluster))
  #     cells <- colnames(sample.tumor.astro)[(sample.tumor.astro$seurat_clusters == cluster)][1:len]
  #     cells.to.plot.astro <- c(cells.to.plot.astro, cells)
  #
  # }
  # p.oligo <- Seurat::DoHeatmap(sample.tumor.oligo,
  #                              features = top10_oligo$gene,
  #                              cells = cells.to.plot.oligo,
  #                              angle = 90,
  #                              lines.width = 25) +
  #     viridis::scale_fill_viridis(na.value = "white") +
  #     Seurat::NoLegend()
  #
  # p.astro <- Seurat::DoHeatmap(sample.tumor.astro,
  #                              features = top10_astro$gene,
  #                              cells = cells.to.plot.astro,
  #                              angle = 90,
  #                              lines.width = 25) +
  #     viridis::scale_fill_viridis(na.value = "white") +
  #     Seurat::NoLegend()

}

progeny_analysis <- function(){

  sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach")
  sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach")

  labels.oligo <- colnames(sample.oligo)
  labels.astro <- colnames(sample.astro)

  progeny_scores <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_LGGs/snRNAseq/10x_3_prime_v3/dorothea_and_progeny/RNBR/regressed/progeny/heatmaps/OE0145-IDH_integrated_LGGs_progeny_scores")

  progeny_scores_oligo <- progeny_scores %>% dplyr::filter(Cell %in% labels.oligo)
  progeny_scores_astro <- progeny_scores %>% dplyr::filter(Cell %in% labels.astro)

  Seurat::Idents(sample.oligo) <- sample.oligo$New_NMF_labelling
  Seurat::Idents(sample.astro) <- sample.astro$New_NMF_labelling


  cluster_mapping_oligo <- data.frame(Cell = names(Seurat::Idents(sample.oligo)),
                                      CellType = as.character(Seurat::Idents(sample.oligo)),
                                      stringsAsFactors = FALSE)

  cluster_mapping_astro <- data.frame(Cell = names(Seurat::Idents(sample.astro)),
                                      CellType = as.character(Seurat::Idents(sample.astro)),
                                      stringsAsFactors = FALSE)


  progeny_scores_oligo <- dplyr::inner_join(progeny_scores_oligo, cluster_mapping_oligo)
  progeny_scores_astro <- dplyr::inner_join(progeny_scores_astro, cluster_mapping_astro)


  summarized_progeny_scores_oligo <- progeny_scores_oligo %>%
                                     dplyr::group_by(Pathway, CellType) %>%
                                     dplyr::summarise(avg = mean(Activity), std = stats::sd(Activity))

  summarized_progeny_scores_astro <- progeny_scores_astro %>%
                                     dplyr::group_by(Pathway, CellType) %>%
                                     dplyr::summarise(avg = mean(Activity), std = stats::sd(Activity))


  mat_oligo <- summarized_progeny_scores_oligo %>%
               dplyr::select(-std) %>%
               tidyr::spread(Pathway, avg) %>%
               data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  mat_astro <- summarized_progeny_scores_astro %>%
               dplyr::select(-std) %>%
               tidyr::spread(Pathway, avg) %>%
               data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


  mat_tum_oligo <- as.matrix(mat_oligo[c("Astro-like", "OPC-like", "Gradient", "Cycling", "RA"), ])
  mat_tum_astro <- as.matrix(mat_astro[c("Astro-like", "OPC-like", "Gradient", "Cycling", "RA"), ])

  mat_tum_no_OPC_cycling_oligo <- as.matrix(mat_oligo[c("Astro-like", "OPC-like", "Gradient", "RA"),])
  mat_tum_no_OPC_cycling_astro <- as.matrix(mat_astro[c("Astro-like", "OPC-like", "Gradient", "RA"),])

  range_tum_oligo <- max(abs(mat_tum_oligo))
  range_tum_astro <- max(abs(mat_tum_astro))

  range_tum_no_OPC_cycling_oligo <- max(abs(mat_tum_no_OPC_cycling_oligo))
  range_tum_no_OPC_cycling_astro <- max(abs(mat_tum_no_OPC_cycling_astro))



  make_bold_names <- function(mat, rc_fun, rc_names) {
    bold_names <- rc_fun(mat)
    ids <- rc_names %>% match(rc_fun(mat))
    ids %>%
      purrr::walk(
        function(i)
          bold_names[i] <<-
          bquote(bold(.(rc_fun(mat)[i]))) %>%
          as.expression()
      )
    bold_names
  }

  blue_color <- "#02294b"
  red_color <- "#4b010b"



  h.tumor.oligo <- SCpubr:::heatmap_inner(mat_tum_oligo,
                                          legend_name = " ",
                                          column_title = "",
                                          row_title = "OD",
                                          row_title_side = "left",
                                          row_names_side = "right",
                                          row_title_rotation = 0,
                                          range.data = 1,
                                          cluster_columns = F,
                                          cluster_rows = F,
                                          cell_size = 5,
                                          outlier.data = TRUE)

  h.tumor.oligo.t <- SCpubr:::heatmap_inner(t(mat_tum_oligo),
                                          legend_name = " ",
                                          row_title = "",
                                          column_title = "OD",
                                          column_title_side = "top",
                                          column_names_side = "bottom",
                                          column_title_rotation = 0,
                                          range.data = 1,
                                          cluster_columns = F,
                                          cluster_rows = F,
                                          cell_size = 5,
                                          outlier.data = TRUE)


  h.tumor.astro <- SCpubr:::heatmap_inner(mat_tum_astro,
                                          legend_name = " ",
                                          column_title = "",
                                          row_title = "AS",
                                          row_title_side = "left",
                                          row_names_side = "right",
                                          row_title_rotation = 0,
                                          range.data = 1,
                                          cluster_columns = F,
                                          cluster_rows = F,
                                          cell_size = 5,
                                          outlier.data = TRUE)
  h.tumor.astro.t <- SCpubr:::heatmap_inner(t(mat_tum_astro),
                                            legend_name = " ",
                                            row_title = "",
                                            column_title = "AS",
                                            column_title_side = "top",
                                            column_names_side = "bottom",
                                            column_title_rotation = 0,
                                            range.data = 1,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            cell_size = 5,
                                            outlier.data = TRUE)

  outlier.top.oligo <- mat_tum_oligo[mat_tum_oligo > 1]
  outlier.down.oligo <- mat_tum_oligo[mat_tum_oligo < -1]
  outlier.top.astro <- mat_tum_astro[mat_tum_astro > 1]
  outlier.down.astro <- mat_tum_astro[mat_tum_astro < -1]


  breaks <- c(-1, -0.5, 0, 0.5, 1)
  colors <- c("#023F73", "#809FB9", "#FFFFFF", "#BC8089", "#7A0213")
  col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors)

  lgd_all <- ComplexHeatmap::Legend(at = breaks,
                                    labels = c("-1", "-0.5", "0", "0.5", "1"),
                                    col_fun = col_fun, title = "",
                                    break_dist = c(1, 1, 1, 1),
                                    labels_gp = grid::gpar(fontsize = 12),
                                    legend_height = grid::unit(4, "cm"))
  lgd_all2 <- ComplexHeatmap::Legend(labels = c(">  1", "< -1"),
                                     legend_gp = grid::gpar(fill = c("> 1" = "#4b010b", "< -1" = "#02294b")),
                                     labels_gp = grid::gpar(fontsize = 12))

  pd = ComplexHeatmap::packLegend(list = list(lgd_all, lgd_all2))

  h.tumor <- suppressWarnings({h.tumor.oligo$heatmap %v% h.tumor.astro$heatmap})
  h.tumor.t <- suppressWarnings({h.tumor.oligo.t$heatmap + h.tumor.astro.t$heatmap})

  ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(4, "mm"))
  suppressWarnings({
    grDevices::pdf(NULL)
    h <- ComplexHeatmap::draw(h.tumor, heatmap_legend_list = pd, padding = unit(c(20, 20, 2, 20), "mm"))
    h.t <- ComplexHeatmap::draw(h.tumor.t, heatmap_legend_list = pd, padding = unit(c(20, 20, 2, 20), "mm"))

    dev.off()
  })

  return(list("h" = h, "h.t" = h.t))
}

dorothea_analysis <- function(){

  sample.oligo <- sample.oligo[, sample.oligo$New_NMF_labelling != "Excluded"]

  labels.oligo <- colnames(sample.oligo)
  labels.astro <- colnames(sample.astro)

  dorothea_markers <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_LGGs/snRNAseq/10x_3_prime_v3/dorothea_and_progeny/RNBR/regressed/dorothea/individual_plots/OE0145-IDH_integrated_LGGs_dorothea_markers")
  viper_scores <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_LGGs/snRNAseq/10x_3_prime_v3/dorothea_and_progeny/RNBR/regressed/dorothea/individual_plots/OE0145-IDH_integrated_LGGs_viper_scores_df")

  viper_scores_oligo <- viper_scores[labels.oligo, ]
  viper_scores_astro <- viper_scores[labels.astro, ]


  cell_clusters_oligo <- data.frame(cell = names(Seurat::Idents(sample.oligo)),
                                    cell_type = as.character(Seurat::Idents(sample.oligo)),
                                    stringsAsFactors = FALSE)

  cell_clusters_astro <- data.frame(cell = names(Seurat::Idents(sample.astro)),
                                    cell_type = as.character(Seurat::Idents(sample.astro)),
                                    stringsAsFactors = FALSE)

  viper_scores_clusters_oligo <- viper_scores_oligo  %>%
                                 data.frame() %>%
                                 tibble::rownames_to_column("cell") %>%
                                 tidyr::gather(tf, activity, -cell) %>%
                                 dplyr::inner_join(cell_clusters_oligo)

  viper_scores_clusters_astro <- viper_scores_astro  %>%
                                 data.frame() %>%
                                 tibble::rownames_to_column("cell") %>%
                                 tidyr::gather(tf, activity, -cell) %>%
                                 dplyr::inner_join(cell_clusters_astro)

  summarized_viper_scores_oligo <- viper_scores_clusters_oligo %>%
                                   dplyr::group_by(tf, cell_type) %>%
                                   dplyr::summarise(avg = mean(activity),
                                                    std = stats::sd(activity))

  summarized_viper_scores_astro <- viper_scores_clusters_astro %>%
                                   dplyr::group_by(tf, cell_type) %>%
                                   dplyr::summarise(avg = mean(activity),
                                                    std = stats::sd(activity))

  highly_variable_tfs_oligo <- summarized_viper_scores_oligo %>%
                               dplyr::group_by(tf) %>%
                               dplyr::mutate(var = stats::var(avg))  %>%
                               dplyr::ungroup() %>%
                               dplyr::top_n((50 * length(levels(sample.oligo))), var) %>%
                               dplyr::distinct(tf)

  highly_variable_tfs_astro <- summarized_viper_scores_astro %>%
                               dplyr::group_by(tf) %>%
                               dplyr::mutate(var = stats::var(avg))  %>%
                               dplyr::ungroup() %>%
                               dplyr::top_n((50 * length(levels(sample.astro))), var) %>%
                               dplyr::distinct(tf)

  tfs_use <- highly_variable_tfs_astro$tf[highly_variable_tfs_astro$tf %in% highly_variable_tfs_oligo$tf]
  tfs_use <- tidyr::tibble(tf = tfs_use)

  mat_oligo <- summarized_viper_scores_oligo %>%
               dplyr::semi_join(tfs_use, by = "tf") %>%
               dplyr::select(-std) %>%
               tidyr::spread(tf, avg) %>%
               data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  mat_astro <- summarized_viper_scores_astro %>%
               dplyr::semi_join(tfs_use, by = "tf") %>%
               dplyr::select(-std) %>%
               tidyr::spread(tf, avg) %>%
               data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)



  mat_tum_oligo <- as.matrix(mat_oligo[c("Astro-like", "OPC-like", "Gradient", "Cycling", "RA"), ])
  mat_tum_astro <- as.matrix(mat_astro[c("Astro-like", "OPC-like", "Gradient", "Cycling", "RA"), ])

  mat_tum_no_OPC_cycling_oligo <- as.matrix(mat_oligo[c("Astro-like", "OPC-like", "Gradient", "RA"),])
  mat_tum_no_OPC_cycling_astro <- as.matrix(mat_astro[c("Astro-like", "OPC-like", "Gradient", "RA"),])

  range_tum_oligo <- max(abs(mat_tum_oligo))
  range_tum_astro <- max(abs(mat_tum_astro))

  range_tum_no_OPC_cycling_oligo <- max(abs(mat_tum_no_OPC_cycling_oligo))
  range_tum_no_OPC_cycling_astro <- max(abs(mat_tum_no_OPC_cycling_astro))



  make_bold_names <- function(mat, rc_fun, rc_names) {
    bold_names <- rc_fun(mat)
    ids <- rc_names %>% match(rc_fun(mat))
    ids %>%
      purrr::walk(
        function(i)
          bold_names[i] <<-
          bquote(bold(.(rc_fun(mat)[i]))) %>%
          as.expression()
      )
    bold_names
  }

  red_color <- "#4b010b"



  h.tumor.oligo <- SCpubr:::heatmap_inner(mat_tum_oligo,
                                          legend_name = " ",
                                          column_title = "",
                                          row_title = "OD",
                                          row_title_side = "left",
                                          row_names_side = "right",
                                          row_title_rotation = 0,
                                          range.data = 1,
                                          cluster_columns = F,
                                          cluster_rows = F,
                                          cell_size = 5,
                                          outlier.data = TRUE)

  h.tumor.oligo.t <- SCpubr:::heatmap_inner(t(mat_tum_oligo),
                                            legend_name = " ",
                                            row_title = "",
                                            column_title = "OD",
                                            column_title_side = "top",
                                            column_names_side = "bottom",
                                            column_title_rotation = 0,
                                            range.data = 1,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            cell_size = 5,
                                            outlier.data = TRUE)


  h.tumor.astro <- SCpubr:::heatmap_inner(mat_tum_astro,
                                          legend_name = " ",
                                          column_title = "",
                                          row_title = "AS",
                                          row_title_side = "left",
                                          row_names_side = "right",
                                          row_title_rotation = 0,
                                          range.data = 1,
                                          cluster_columns = F,
                                          cluster_rows = F,
                                          cell_size = 5,
                                          outlier.data = TRUE)
  h.tumor.astro.t <- SCpubr:::heatmap_inner(t(mat_tum_astro),
                                            legend_name = " ",
                                            row_title = "",
                                            column_title = "AS",
                                            column_title_side = "top",
                                            column_names_side = "bottom",
                                            column_title_rotation = 0,
                                            range.data = 1,
                                            cluster_columns = F,
                                            cluster_rows = F,
                                            cell_size = 5,
                                            outlier.data = TRUE)

  outlier.top.oligo <- mat_tum_oligo[mat_tum_oligo > 1]
  outlier.down.oligo <- mat_tum_oligo[mat_tum_oligo < -1]
  outlier.top.astro <- mat_tum_astro[mat_tum_astro > 1]
  outlier.down.astro <- mat_tum_astro[mat_tum_astro < -1]


  breaks <- c(-1, -0.5, 0, 0.5, 1)
  colors <- c("#023F73", "#809FB9", "#FFFFFF", "#BC8089", "#7A0213")
  col_fun <- circlize::colorRamp2(breaks = breaks, colors = colors)

  lgd_all <- ComplexHeatmap::Legend(at = breaks,
                                    labels = c("-1", "-0.5", "0", "0.5", "1"),
                                    col_fun = col_fun, title = "",
                                    break_dist = c(1, 1, 1, 1),
                                    labels_gp = grid::gpar(fontsize = 12),
                                    legend_height = grid::unit(4, "cm"))

  lgd_all2 <- ComplexHeatmap::Legend(labels = c(">  1"),
                                     legend_gp = grid::gpar(fill = c("> 1" = red_color),
                                     labels_gp = grid::gpar(fontsize = 12)))

  pd = ComplexHeatmap::packLegend(list = list(lgd_all, lgd_all2))

  h.tumor <- suppressWarnings({h.tumor.oligo$heatmap %v% h.tumor.astro$heatmap})
  h.tumor.t <- suppressWarnings({h.tumor.oligo.t$heatmap + h.tumor.astro.t$heatmap})

  ComplexHeatmap::ht_opt("HEATMAP_LEGEND_PADDING" = ggplot2::unit(4, "mm"))
  suppressWarnings({
    grDevices::pdf(NULL)
    h <- ComplexHeatmap::draw(h.tumor, heatmap_legend_list = pd, padding = unit(c(20, 20, 2, 20), "mm"))
    h.t <- ComplexHeatmap::draw(h.tumor.t, heatmap_legend_list = pd, padding = unit(c(20, 20, 2, 20), "mm"))

    dev.off()
  })

  return(list("h" = h, "h.t" = h.t))
}


monocle_analysis <- function(){
  sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_oligodendroglioma_relabelled_with_permutation_approach")
  sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_astrocytoma_relabelled_with_permutation_approach")

  monocle_folder <- paste0(figures_folder, "/second_iteration/main_figures/Figure_2/monocle3/")
  name_oligo <- "OE0145-IDH_integrated_oligodendroglioma"
  name_astro <- "OE0145-IDH_integrated_astrocytoma"


  sample.tumor.oligo <- sample.oligo[, sample.oligo$New_NMF_labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA")]
  sample.tumor.astro <- sample.astro[, sample.astro$New_NMF_labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA")]

  rm(sample.oligo)
  rm(sample.astro)

  sample.tumor.oligo[["monocle3_clusters"]] <- sample.tumor.oligo$New_NMF_labelling
  sample.tumor.oligo[["monocle3_partitions"]] <- 1
  cds.oligo <- SeuratWrappers::as.cell_data_set(sample.tumor.oligo)

  sample.tumor.astro[["monocle3_clusters"]] <- sample.tumor.astro$New_NMF_labelling
  sample.tumor.astro[["monocle3_partitions"]] <- 1
  cds.astro <- SeuratWrappers::as.cell_data_set(sample.tumor.astro)

  cds.oligo <- monocle3::learn_graph(cds.oligo)
  cds.astro <- monocle3::learn_graph(cds.astro)

  p.trajectory.oligo <- monocle3::plot_cells(cds.oligo,
                                             label_groups_by_cluster = FALSE,
                                             label_leaves = FALSE,
                                             label_branch_points = FALSE,
                                             label_cell_groups = FALSE,
                                             cell_size = 1) +
    ggplot2::ggtitle("Inferred trajectory monocle 3") +
    ggplot2::scale_color_manual(values = cluster_cols[c("Astro-like", "OPC-like", "Gradient", "RA")]) +
    ggpubr::theme_pubr(legend = "bottom") +
    Seurat::NoAxes() +
    ggplot2::ggtitle("") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.text = ggplot2::element_text(size = 14, face = "bold"),
                   legend.title = ggplot2::element_text(size = 14, face = "bold")) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 5,
                                                  byrow = T,
                                                  override.aes = list(size = 4))) +
    ggpubr::rremove("legend.title")
  p.trajectory.oligo$layers[[2]]$aes_params$size <- 2.5
  p.trajectory.oligo$layers[[2]]$aes_params$colour <- "black"


  p.trajectory.astro <- monocle3::plot_cells(cds.astro,
                                             label_groups_by_cluster = FALSE,
                                             label_leaves = FALSE,
                                             label_branch_points = FALSE,
                                             label_cell_groups = FALSE,
                                             cell_size = 1) +
    ggplot2::ggtitle("Inferred trajectory monocle 3") +
    ggplot2::scale_color_manual(values = cluster_cols[c("Astro-like", "OPC-like", "Gradient", "RA")]) +
    ggpubr::theme_pubr(legend = "bottom") +
    Seurat::NoAxes() +
    ggplot2::ggtitle("") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
                   legend.text = ggplot2::element_text(size = 14, face = "bold"),
                   legend.title = ggplot2::element_text(size = 14, face = "bold")) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = 5,
                                                  byrow = T,
                                                  override.aes = list(size = 4))) +
    ggpubr::rremove("legend.title")
  p.trajectory.astro$layers[[2]]$aes_params$size <- 2.5
  p.trajectory.astro$layers[[2]]$aes_params$colour <- "black"

  savefig(p.trajectory.oligo, figure_path = paste0(monocle_folder, "/", name_oligo), file_name = paste0(name_oligo, "_tumor_cells_ordered"), width = 10, heigth = 10)
  savefig(p.trajectory.astro, figure_path = paste0(monocle_folder, "/", name_astro), file_name = paste0(name_astro, "_tumor_cells_ordered"), width = 10, heigth = 10)

}

# Figure 3A:
fig_3a <- Figure_3A()
SCpubr::save_Plot(plot = fig_3a,
                 figure_path = paste0(figure_path, "/Figure_3/Figure_3A/"),
                 file_name = "Figure_3A_test",
                 height = 12,
                 width = 18,
                 dpi = 300)

# Figure 1A:
return_list <- Figure_1A()
fig_1a <- return_list$portrait


SCpubr::save_Plot(plot = fig_1a,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1A/"),
                 file_name = "Figure_1A",
                 height = 14,
                 width = 14,
                 dpi = 300)


# Figure 1B:
fig_1b <- SCpubr::do_DimPlot(sample = sample.oligo,
                             group.by = "NMF_labelling",
                             colors.use = cluster_cols,
                             plot.title = "Oligodendroglioma",
                             legend.position = "bottom",
                             legend.nrow = 3,
                             fontsize = 22)
fig_1b.raster <- SCpubr::do_DimPlot(sample = sample.oligo,
                                    group.by = "NMF_labelling",
                                    colors.use = cluster_cols,
                                    plot.title = "Oligodendroglioma",
                                    pt.size = 3,
                                    raster = TRUE,
                                    legend.position = "bottom",
                                    legend.nrow = 3,
                                    fontsize = 22)
SCpubr::save_Plot(plot = fig_1b,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1B/"),
                 file_name = "Figure_1B",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = fig_1b.raster,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1B/"),
                 file_name = "Figure_1B_raster",
                 height = 8,
                 width = 8,
                 dpi = 300)


# Figure 1C:
fig_1c <- SCpubr::do_DimPlot(sample = sample.astro,
                             group.by = "NMF_labelling",
                             colors.use = cluster_cols,
                             plot.title = "Astrocyoma",
                             legend.position = "bottom",
                             legend.nrow = 3,
                             fontsize = 22)
fig_1c.raster <- SCpubr::do_DimPlot(sample = sample.astro,
                                    group.by = "NMF_labelling",
                                    colors.use = cluster_cols,
                                    plot.title = "Astrocytoma",
                                    pt.size = 3,
                                    raster = TRUE,
                                    legend.position = "bottom",
                                    legend.nrow = 3,
                                    fontsize = 22)
SCpubr::save_Plot(plot = fig_1c,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1C/"),
                 file_name = "Figure_1C",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = fig_1c.raster,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1C/"),
                 file_name = "Figure_1C_raster",
                 height = 8,
                 width = 8,
                 dpi = 300)

# Figure 1D:
fig_1d <- Figure_1D_1E(type = "oligodendroglioma")
SCpubr::save_Plot(plot = fig_1d,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1D/"),
                 file_name = "Figure_1D",
                 height = 7,
                 width = 11,
                 dpi = 300)

# Figure 1E:
fig_1e <- Figure_1D_1E(type = "astrocytoma")
SCpubr::save_Plot(plot = fig_1e,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1E/"),
                 file_name = "Figure_1E",
                 height = 7,
                 width = 11,
                 dpi = 300)

# Figure 1F:
fig_1f <- SCpubr::do_DimPlot(sample = sample.oligo,
                             group.by = "New_NMF_labelling",
                             colors.use = cluster_cols,
                             plot.title = "Oligodendroglioma",
                             legend.position = "bottom",
                             idents.keep = c("OPC-like", "Astro-like", "Cycling", "RA", "Gradient"),
                             na.value = "grey50",
                             legend.ncol = 3,
                             fontsize = 22)
fig_1f.raster <- SCpubr::do_DimPlot(sample = sample.oligo,
                                    group.by = "New_NMF_labelling",
                                    colors.use = cluster_cols,
                                    plot.title = "Oligodendroglioma",
                                    pt.size = 3,
                                    raster = TRUE,
                                    legend.position = "bottom",
                                    idents.keep = c("OPC-like", "Astro-like", "Cycling", "RA", "Gradient"),
                                    na.value = "grey50",
                                    legend.ncol = 3,
                                    fontsize = 22)
SCpubr::save_Plot(plot = fig_1f,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1F/"),
                 file_name = "Figure_1F",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = fig_1f.raster,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1F/"),
                 file_name = "Figure_1F_raster",
                 height = 8,
                 width = 8,
                 dpi = 300)


# Figure 1G:
fig_1g <- SCpubr::do_DimPlot(sample = sample.astro,
                             group.by = "New_NMF_labelling",
                             colors.use = cluster_cols,
                             plot.title = "Astrocytoma",
                             legend.position = "bottom",
                             idents.keep = c("OPC-like", "Astro-like", "Cycling", "RA", "Gradient"),
                             na.value = "grey50",
                             legend.ncol = 3,
                             fontsize = 22)
fig_1g.raster <- SCpubr::do_DimPlot(sample = sample.astro,
                                    group.by = "New_NMF_labelling",
                                    colors.use = cluster_cols,
                                    plot.title = "Astrocytoma",
                                    pt.size = 3,
                                    raster = TRUE,
                                    legend.position = "bottom",
                                    idents.keep = c("OPC-like", "Astro-like", "Cycling", "RA", "Gradient"),
                                    na.value = "grey50",
                                    legend.ncol = 3,
                                    fontsize = 22)
SCpubr::save_Plot(plot = fig_1g,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1G/"),
                 file_name = "Figure_1G",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = fig_1g.raster,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1G/"),
                 file_name = "Figure_1G_raster",
                 height = 8,
                 width = 8,
                 dpi = 300)
# Figure 1H:

# Manually picked.
features_use <- c("NRG3", "ADGRV1", "SLC4A4", "SPARCL1", "ADCY2",
                  "MKI67", "CENPK", "EZH2", "POLQ", "EGFR",
                  "OPCML", "DSCAM", "FGF12", "SOX6", "DLGAP1",
                  "EEF2", "EEF1A1", "OLIG1", "ETV1", "RPL13"
)
features_use.list <- list("Astro-like" = c("NRG3", "ADGRV1", "SLC4A4", "SPARCL1", "ADCY2"),
                  "Cycling" = c("MKI67", "CENPK", "EZH2", "POLQ", "EGFR"),
                   "OPC-like" = c("OPCML", "DSCAM", "FGF12", "SOX6", "DLGAP1"),
                  "RA" = c("EEF2", "EEF1A1", "OLIG1", "ETV1", "RPL13")
)
scale.subtype <- c("AS" = "#b38b14",
                   "OD" = "#3c5b8b")


sample.oligo.tumor <- sample.oligo[, sample.oligo$New_NMF_labelling %in% c("Astro-like", "OPC-like", "RA", "Cycling", "Gradient")]
Seurat::Idents(sample.oligo.tumor) <- factor(sample.oligo.tumor$New_NMF_labelling, levels =rev(c("Astro-like", "Cycling", "Gradient", "OPC-like", "RA")))
fig_1h_v <- SCpubr::do_DotPlot(sample = sample.oligo.tumor,
                               features = features_use.list,
                               colors.use = c("grey75", "#3c5b8b"),
                               legend = TRUE,
                               legend.position = "right",
                               dot.scale = 10,
                               ylab = "OD",
                               fontsize = 22) +
             ggplot2::theme(strip.text = ggplot2::element_text(size = 21)) +
             ggpubr::rremove("x.text") +
             ggpubr::rremove("x.ticks") +
             ggpubr::rremove("x.axis")
fig_1h_v$guides$colour$title <- "Average Expression\nOligodendroglioma"
fig_1h_v$guides$size$order <- 1
fig_1h_v$guides$colour$order <- 2
fig_1h_v$guides$colour$label.hjust <- 1
fig_1h_v$guides$colour$label.hjust <- 0


sample.astro.tumor <- sample.astro[, sample.astro$New_NMF_labelling %in% c("Astro-like", "OPC-like", "RA", "Cycling", "Gradient")]
Seurat::Idents(sample.astro.tumor) <- factor(sample.astro.tumor$New_NMF_labelling, levels =rev(c("Astro-like", "Cycling", "Gradient", "OPC-like", "RA")))
fig_1i_v <- SCpubr::do_DotPlot(sample = sample.astro.tumor,
                               features = features_use.list,
                               colors.use = c("grey75", "#b38b14"),
                               ylab = "AS",
                               legend = TRUE,
                               legend.position = "right",
                               dot.scale = 10,
                               fontsize = 22)  +
             ggplot2::guides(size = "none") +
             ggplot2::theme(strip.text = ggplot2::element_blank(),
                            strip.background = ggplot2::element_blank())

fig_1i_v$guides$colour$title <- "Average Expression\nAstrocytoma"
fig_1i_v$guides$colour$label.hjust <- 1
fig_1i_v$guides$colour$label.hjust <- 0


layout <- "AAAC
BBBC"
patch_v <- patchwork::wrap_plots(A = fig_1h_v,
                                 B = fig_1i_v,
                                 C = patchwork::guide_area(),
                                 design = layout,
                                 guides = "collect")

SCpubr::save_Plot(plot = patch_v,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1H/"),
                 file_name = "Figure_1H",
                 height = 8,
                 width = 13,
                 dpi = 300)

# Figure 1I:
sample.lgg <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_integrated_LGGs/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_integrated_LGGs_after_normalization_to_use")
sample.lgg <- sample.lgg[, sample.lgg$New_NMF_labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA", "Cycling")]



sample.lgg$orig.ident <- stringr::str_remove_all(sample.lgg$orig.ident, "OE0145-")

# Order of bars.
labels.use <- c("IDH_ACB_AD_809",
                "IDH_NCH6341",
                "IDH_NCH536",
                "IDH_ACB_AD_540",
                "IDH_NCH6702",
                "IDH_NCH2111",
                "IDH_NCH781",
                "IDH_ACB_AD_883",
                "IDH_ACB_AD_785",
                "IDH_ACB_AD_832",
                "IDH_NCH2018",
                "IDH_ACB_AD_865",
                "IDH_NCH2157",
                "IDH_NCH2164")

# Bar plot cell identities.
out <- SCpubr::do_BarPlot(sample.lgg,
                          features = "orig.ident",
                          group.by = "New_NMF_labelling",
                          labels.order = labels.use,
                          order.by = "RA",
                          position = "fill",
                          colors.use = cluster_cols,
                          xlab = "",
                          ylab = "Proportion of each cell identity in sample",
                          horizontal = T,
                          return_data_matrix = T,
                          legend.position = "bottom",
                          legend.ncol = 2)
p.orig <- out$plot +
          ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2))
data.long <- out$data$orig.ident$long
data.wide <- out$data$orig.ident$wide

# Bar plot grade.
p.grade <- SCpubr::do_BarPlot(sample.lgg,
                      features = "orig.ident",
                      labels.order = labels.use,
                      group.by = "grade",
                      position = "fill",
                      colors.use = scale.grade,
                      xlab = "",
                      ylab = "",
                      horizontal = T,
                      legend.title = F,
                      legend.position = "bottom",
                      legend.ncol = 1) +
            ggpubr::rremove("axis") +
            ggpubr::rremove("axis.text") +
            ggpubr::rremove("ticks")

# Bar plot subtypes.
astro_samples <- c("IDH_ACB_AD_785",
                   "IDH_ACB_AD_832",
                   "IDH_NCH2018",
                   "IDH_ACB_AD_865",
                   "IDH_NCH2157",
                   "IDH_NCH2164")
sample.lgg$subtype <- ifelse(sample.lgg$orig.ident %in% astro_samples, "AS", "OD")
p.subtype <- SCpubr::do_BarPlot(sample.lgg,
                        feature = "orig.ident",
                        labels.order = labels.use,
                        group.by = "subtype",
                        position = "fill",
                        colors.use = scale.subtype.short,
                        xlab = "",
                        ylab = "",
                        horizontal = T,
                        legend.title = FALSE,
                        legend.ncol = 1,
                        legend.position = "bottom") +
              ggpubr::rremove("axis") +
              ggpubr::rremove("axis.text") +
              ggpubr::rremove("ticks")



layout <- "AAAAAAAAAAAAAABC"
patch <- patchwork::wrap_plots(A = p.orig,
                               B = p.grade,
                               C = p.subtype,
                               design = layout)

SCpubr::save_Plot(plot = patch,
                 figure_path = paste0(figure_path, "/Figure_1/Figure_1I/"),
                 file_name = "Figure_1I",
                 height = 8,
                 width = 11,
                 dpi = 300)

write.table(data.long, file = paste0(figure_path, "/Figure_1/Figure_1I/Figure_1I_data_long.tsv"),
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

write.table(data.wide, file = paste0(figure_path, "/Figure_1/Figure_1I/Figure_1I_data_wide.tsv"),
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")



# Figure 2X: Progeny.
out <- progeny_analysis()
h <- out$h
h.t <- out$h.t

SCpubr::save_Plot(plot = h,
                 figure_path = paste0(figure_path, "/Figure_2/progeny/"),
                 file_name = "progeny_tumor_bulk",
                 height = 7,
                 width = 11,
                 dpi = 300)
SCpubr::save_Plot(plot = h.t,
                 figure_path = paste0(figure_path, "/Figure_2/progeny/"),
                 file_name = "progeny_tumor_bulk_t",
                 height = 11,
                 width = 7,
                 dpi = 300)


# Figure 2X: Dorothea.

out <- dorothea_analysis()
h <- out$h
h.t <- out$h.t


SCpubr::save_Plot(plot = h,
                 figure_path = paste0(figure_path, "/Figure_2/dorothea/"),
                 file_name = "dorothea_tumor_bulk",
                 height = 7,
                 width = 11,
                 dpi = 300)
SCpubr::save_Plot(plot = h.t,
                 figure_path = paste0(figure_path, "/Figure_2/dorothea/"),
                 file_name = "dorothea_tumor_bulk_t",
                 height = 11,
                 width = 7,
                 dpi = 300)



# Figure 3X: CellphoneDB
p <- cellphoneDB_output()
SCpubr::save_Plot(plot = p,
                 figure_path = paste0(figure_path, "/Figure_3/cellphoneDB/"),
                 file_name = "dorothea_tumor_bulk",
                 height = 10,
                 width = 18,
                 dpi = 300)

# FIgure 2X: Diffusion maps only stem markers.
diffusion_maps_analysis()

# Figure 1X: Barplots.


savefig(plot = patch, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/barplots_celltype_composition/")), file_name = paste0("barplots_by_sample_annotated_ordered_by_LGG_subtype_and_grade_and_RA"), width = 15, heigth = 12)


patch_v <- p.oligo | p.astro


plot_h <- do_BarPlot(sample.lgg,
                     var.to.plot = "subtype",
                     group.by = "New_NMF_labelling",
                     position = "fill",
                     colors.use = cluster_cols[names(cluster_cols) %in% unique(sample.lgg$New_NMF_labelling)],
                     ylab = "Cell type composition",
                     xlab = "LGG subtype",
                     horizontal = T)
plot_v <- do_BarPlot(sample.lgg,
                     var.to.plot = "subtype",
                     group.by = "New_NMF_labelling",
                     position = "fill",
                     colors.use = cluster_cols[names(cluster_cols) %in% unique(sample.lgg$New_NMF_labelling)],
                     ylab = "Cell type composition",
                     xlab = "LGG subtype",
                     horizontal = F,
                     legend.position = "right")
savefig(plot = plot_h, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/barplots_celltype_composition/")), file_name = paste0("barplots_celltype_composition_horizontal"), width = 10, heigth = 5)
savefig(plot = plot_v, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/barplots_celltype_composition/")), file_name = paste0("barplots_celltype_composition_vertical"), width = 8, heigth = 10)
savefig(plot = plot_h, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/barplots_celltype_composition/")), file_name = paste0("barplots_celltype_composition_horizontal_squared"), width = 12, heigth = 12)
savefig(plot = plot_v, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/barplots_celltype_composition/")), file_name = paste0("barplots_celltype_composition_vertical_squared"), width = 12, heigth = 12)

# Figure SX: DimPlot integrated LGG.
plot <- do_DimPlot(sample.lgg, group.by = "New_NMF_labelling", label = F, colors.use = cluster_cols[names(cluster_cols) %in% unique(sample.lgg$New_NMF_labelling)])
savefig(plot = plot, figure_path = paste0(paste0("/omics/odcf/analysis/OE0145_projects/idh_gliomas/Figures_Science/second_iteration/supplementary/integrated_LGG_DimPlot/")), file_name = paste0("integrated_LGG_DimPlot"), width = 7, heigth = 8)




# Fig 2D modified.

min_value <- min(c(min(sample.oligo$Suva_stemness1), min(sample.astro$Suva_stemness1)))
max_value <- max(c(max(sample.oligo$Suva_stemness1), max(sample.astro$Suva_stemness1)))

p1 <- SCpubr::do_BeeSwarmPlot(sample.oligo,
                              feature_to_rank = "Suva_stemness1",
                              group.by = "New_NMF_labelling",
                              continuous_feature = T) +
  ggplot2::xlab("Stemness enrichment") +
  ggplot2::ggtitle("Oligodendroglioma") +
  ggplot2::scale_color_viridis_c(limits = c(min_value, max_value)) +
  ggpubr::rremove("legend")

p2 <- SCpubr::do_BeeSwarmPlot(sample.astro,
                              feature_to_rank = "Suva_stemness1",
                              group.by = "New_NMF_labelling",
                              continuous_feature = T) +
  ggplot2::xlab("Stemness enrichment") +
  ggplot2::ggtitle("Astrocytoma") +
  ggplot2::scale_color_viridis_c(limits = c(min_value, max_value)) +
  ggpubr::rremove("y.text") +
  ggpubr::rremove("y.ticks")

p1 | p2


p1 <- SCpubr::do_VlnPlot(sample.oligo, features = "Suva_stemness1", colors.use = cluster_cols, plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_VlnPlot(sample.astro, features = "Suva_stemness1", colors.use = cluster_cols, plot.title = "Astrocytoma")
p1 /p2










# FIX FOR INMA FIGURES

# FIGURE 2
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/OE0145-IDH_integrated_astrocytoma_peaks_activity_chromvar_IHL.rds")
sample.astro$predicted.id[sample.astro$predicted.id == "Undefined"] <- "Mixed"
sample.astro$labelling <- factor(sample.astro$predicted.id, levels = sort(unique(sample.astro$predicted.id)))
sample.astro$subtype <- "Astrocytoma"
grade_2_samples <- c("IDH_ACB_AD_809",
                     "IDH_NCH536",
                     "IDH_NCH6341",
                     "IDH_ACB_AD_832",
                     "IDH_ACB_AD_785")
sample.astro$grade <- ifelse(sample.astro$orig.ident %in% grade_2_samples, "2", "3")
Seurat::Idents(sample.astro) <- sample.astro$labelling
sample.astro$orig.ident <- stringr::str_remove_all(sample.astro$orig.ident, "OE0145-")



sample.oligo<- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/OE0145-IDH_integrated_oligodendroglioma_peaks_activity_chromvar_IHL.rds")
sample.oligo$predicted.id[sample.oligo$predicted.id == "Undefined"] <- "Mixed"
sample.oligo$labelling <- factor(sample.oligo$predicted.id, levels = sort(unique(sample.oligo$predicted.id)))
sample.oligo$subtype <- "Oligodendroglioma"
grade_2_samples <- c("IDH_ACB_AD_809",
                     "IDH_NCH536",
                     "IDH_NCH6341",
                     "IDH_ACB_AD_832",
                     "IDH_ACB_AD_785")
sample.oligo$grade <- ifelse(sample.oligo$orig.ident %in% grade_2_samples, "2", "3")
Seurat::Idents(sample.oligo) <- sample.oligo$labelling
sample.oligo$orig.ident <- stringr::str_remove_all(sample.oligo$orig.ident, "OE0145-")

saveRDS(sample.oligo, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/oligodendroglioma_atac_ready_to_use.rds", compress = F)
saveRDS(sample.astro, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use.rds", compress = F)

sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/oligodendroglioma_atac_ready_to_use.rds")
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use.rds")

# UMAP overview
layout <- "ABBBBB"
p.a <- SCpubr::do_DimPlot(sample.astro,
                          group.by = "labelling",
                          idents.keep = c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like"),
                          colors.use = cluster_cols,
                          na.value = "grey50",
                          legend.position = "bottom",
                          legend.ncol = 3)
p.b <- SCpubr::do_DimPlot(sample.astro,
                          split.by = "labelling",
                          idents.keep = c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like"),
                          colors.use = cluster_cols,
                          na.value = "grey50",
                          legend = F,
                          ncol = 5)
p <- patchwork::wrap_plots(A = p.a, B = p.b, design = layout)

layout <- "ABBBB"
p.a <- SCpubr::do_DimPlot(sample.oligo,
                          group.by = "labelling",
                          idents.keep = c("OPC-like", "Astro-like", "Gradient", "RA"),
                          colors.use = cluster_cols,
                          na.value = "grey50",
                          legend.position = "bottom",
                          legend.ncol = 3)
p.b <- SCpubr::do_DimPlot(sample.oligo,
                          split.by = "labelling",
                          idents.keep = c("OPC-like", "Astro-like", "Gradient", "RA"),
                          colors.use = cluster_cols,
                          na.value = "grey50",
                          legend = F,
                          ncol = 4)
p <- patchwork::wrap_plots(A = p.a, B = p.b, design = layout)

# Number of cells per sample.
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         features = "orig.ident",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = scale.ident,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per sample",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         features = "orig.ident",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = scale.ident,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per sample",
                         plot.title = "Astrocytoma")

p <- p1 | p2

# Number of cells per identity - fill
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         group.by = "orig.ident",
                         features = "subtype",
                         horizontal = T,
                         colors.use = scale.ident,
                         legend = T,
                         legend.ncol = 1,
                         legend.position = "right",
                         position = "fill") +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.text") +
  ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         group.by = "orig.ident",
                         features = "subtype",
                         horizontal = T,
                         colors.use = scale.ident,
                         legend = T,
                         legend.position = "right",
                         position = "fill",
                         legend.ncol = 1,
                         ylab = "Proportion of cells per sample")

p <- p1 / p2


# Number of cells per identity
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         features = "labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         features = "labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Astrocytoma")

p <- p1 | p2


# Number of cells per identity - fill
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         group.by = "labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = F,
                         legend.position = "right",
                         position = "fill") +
      ggpubr::rremove("x.axis") +
      ggpubr::rremove("x.text") +
      ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         group.by = "labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = T,
                         legend.position = "bottom",
                         position = "fill",
                         legend.nrow = 2,
                         ylab = "Proportion of cells per cell type")

p <- p1 / p2

# Number of cells per identity in the tumor bulk
p1 <- SCpubr::do_BarPlot(sample.oligo[, sample.oligo$labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         features = "labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro[, sample.astro$labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         features = "labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Astrocytoma")

p <- p1 | p2

# Number of cells per identity in the tumor bulk - fill
p1 <- SCpubr::do_BarPlot(sample.oligo[, sample.oligo$labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         group.by = "labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = F,
                         legend.position = "right",
                         position = "fill") +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.text") +
  ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro[, sample.astro$labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         group.by = "labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = T,
                         legend.position = "bottom",
                         position = "fill",
                         legend.nrow = 1,
                         ylab = "Proportion of cells per cell type")

p <- p1 / p2

# FIGURE 3

sample <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_3/rds_figure3/Clusters_microglia_onlyPaired.rds")




cluster_cols <- c("Homeostatic"= "#ECA809",                          # Marigold.
                  "Inflammatory" = "#043362",                    # Prussian Blue.
                  "Activated"= "#9A031E",                  # Ruby Red
                  "Resident-like"= "#009FF5",                     # Carolina Blue.
                  "Anti-inflammatory"= "#BC5210",            # Burnt Orange.
                  "Phagocytic"= "#279185",                  # Celadon Green.
                  "IFNg"= "#7EB356",                   # Bud Green.
                  "Stressed"= "#AC70FF",                   # Medium Purple.
                  "Inflammatory ICAM1+"= "#63412C",                   # Van Dyke Brown.
                  "Infiltrating"= "#D6D6D6",
                  "Hypoxic" = "#5F0F40")

astro_samples <- c("OE0145-IDH_ACB_AD_785",
                   "OE0145-IDH_ACB_AD_832",
                   "OE0145-IDH_ACB_AD_865",
                   "OE0145-IDH_NCH2018",
                   "OE0145-IDH_NCH2157",
                   "OE0145-IDH_NCH2164")

grade_2_samples <- c("IDH_ACB_AD_809",
                     "IDH_NCH536",
                     "IDH_NCH6341",
                     "IDH_ACB_AD_832",
                     "IDH_ACB_AD_785")
sample$orig.ident <- sample$gem_id
sample$subtype <- "Astrocytoma"
sample$grade <- ifelse(sample$orig.ident %in% grade_2_samples, "2", "3")

sample$labelling <- as.character(sample$Subclusters)
sample$Subclusters <- NULL
sample$gem_id <- NULL

sample$labelling[sample$labelling %in% "Mg Activated"] <- "Activated"
sample$labelling[sample$labelling %in% "Mg-IFNg TAMs"] <- "IFNg"
sample$labelling[sample$labelling %in% "Mg Inflammatory ICAM+"] <- "Inflammatory ICAM+"
sample$labelling[sample$labelling %in% "Mg homeostatic"] <- "Homeostatic"
sample$labelling[sample$labelling %in% "Mg inflammatory TAMs"] <- "Inflammatory"
sample$labelling[sample$labelling %in% "Mg phagocytic"] <- "Phagocytic"
sample$labelling[sample$labelling %in% "Mg resident-like TAMs"] <- "Resident-like"
sample$labelling[sample$labelling %in% "Mg stressed TAMs"] <- "Stressed"
sample$labelling[sample$labelling %in% "Mo-TAMs Infiltrating"] <- "Infiltrating"
sample$labelling[sample$labelling %in% "Mo-TAMs anti-inflammatory"] <- "Anti-inflammatory"





# FIGURE 4
sample <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_4/rds_figure4/Clusters_microglia_primary_annotated.rds")

astro_samples <- c("OE0145-IDH_ACB_AD_785",
                   "OE0145-IDH_ACB_AD_832",
                   "OE0145-IDH_ACB_AD_865",
                   "OE0145-IDH_NCH2018",
                   "OE0145-IDH_NCH2157",
                   "OE0145-IDH_NCH2164")

grade_2_samples <- c("IDH_ACB_AD_809",
                     "IDH_NCH536",
                     "IDH_NCH6341",
                     "IDH_ACB_AD_832",
                     "IDH_ACB_AD_785")
sample$orig.ident <- sample$gem_id
sample$subtype <- ifelse(sample$orig.ident %in% astro_samples, "Astrocytoma", "Oligodendroglioma")
sample$grade <- ifelse(sample$orig.ident %in% grade_2_samples, "2", "3")

sample$labelling <- as.character(sample$Subclusters)
sample$Subclusters <- NULL
sample$gem_id <- NULL

sample$labelling[sample$labelling %in% "Mg Activated"] <- "Activated"
sample$labelling[sample$labelling %in% "Mg-IFNg TAMs"] <- "IFNg"
sample$labelling[sample$labelling %in% "Mg Inflammatory ICAM+"] <- "Inflammatory ICAM+"
sample$labelling[sample$labelling %in% "Mg homeostatic"] <- "Homeostatic"
sample$labelling[sample$labelling %in% "Mg inflammatory TAMs"] <- "Inflammatory"
sample$labelling[sample$labelling %in% "Mg phagocytic"] <- "Phagocytic"
sample$labelling[sample$labelling %in% "Mg resident-like TAMs"] <- "Resident-like"
sample$labelling[sample$labelling %in% "Mg stressed TAMs"] <- "Stressed"
sample$labelling[sample$labelling %in% "Mo-TAMs Infiltrating"] <- "Infiltrating"
sample$labelling[sample$labelling %in% "Mo-TAMs anti-inflammatory"] <- "Anti-inflammatory"
interactions$clusters_plot <- factor(interactions$clusters_plot, levels = sort(unique(interactions$clusters_plot)))

Seurat::Idents(sample) <- sample$labelling
saveRDS(sample, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_4/rds_figure4/microglia_combined_to_use.rds", compress = F)
sample <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_4/rds_figure4/microglia_combined_to_use.rds")


# Composition of microglia per sample.

p <- SCpubr::do_BarPlot(sample,
                        features = "orig.ident",
                        horizontal = T,
                        add.summary_labels = T,
                        colors.use = scale.ident,
                        size.labels = 5,
                        legend = F,
                        plot.title = "Microglia composition by sample")

p.subtype <- SCpubr::do_BarPlot(sample,
                                features = "orig.ident",
                                group.by = "subtype_short",
                                horizontal = T,
                                position = "fill",
                                colors.use = scale.subtype.short,
                                labels.order = rev(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()),
                                legend.position = "bottom",
                                legend.ncol = 1,
                                legend.title = TRUE,
                                legend.title.position = "top") +
             Seurat::NoAxes()+
             ggplot2::labs(fill = "Subtype")

p.grade <- SCpubr::do_BarPlot(sample,
                              features = "orig.ident",
                              group.by = "grade",
                              horizontal = T,
                              position = "fill",
                              colors.use = scale.grade,
                              labels.order = rev(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()),
                              legend.position = "bottom",
                              legend.ncol = 1,
                              legend.title = TRUE,
                              legend.title.position = "top") +
  Seurat::NoAxes() +
  ggplot2::labs(fill = "Grade")

layout <- "AAAAAAAAAAAAAAAAAABC"
patch <- patchwork::wrap_plots(A = p,
                               B = p.grade,
                               C = p.subtype,
                               design = layout)




order <- c("Homeostatic",
           "Resident-like",
           "Activated",
           "Phagocytic",
           "Stressed",
           "IFNg",
           "Inflammatory",
           "Inflammatory ICAM1+",
           "Anti-inflammatory",
           "Infiltrating")
cluster_cols <- c("Homeostatic"= "#ECA809",                          # Marigold.
                  "Inflammatory" = "#043362",                    # Prussian Blue.
                  "Activated"= "#9A031E",                  # Ruby Red
                  "Resident-like"= "#009FF5",                     # Carolina Blue.
                  "Anti-inflammatory"= "#BC5210",            # Burnt Orange.
                  "Phagocytic"= "#279185",                  # Celadon Green.
                  "IFNg"= "#7EB356",                   # Bud Green.
                  "Stressed"= "#AC70FF",                   # Medium Purple.
                  "Inflammatory ICAM1+"= "#63412C",                   # Van Dyke Brown.
                  "Infiltrating"= "#D6D6D6")
cluster_cols <- cluster_cols[order]

# Microglia UMAP
p <- SCpubr::do_DimPlot(sample = sample,
                        colors.use = cluster_cols,
                        legend.position = "bottom",
                        legend.ncol = 4)

# Bar plots proportions.

p <- SCpubr::do_BarPlot(sample = sample,
                        feature = "subtype",
                        group.by = "labelling",
                        colors.use = cluster_cols,
                        position = "fill",
                        horizontal = TRUE)

SCpubr::save_plot(p = p,
                  figure_path = "~/test/",
                  create_path = TRUE,
                  file_name = "test",
                  width = 10,
                  height = 5)





# Markers assessment.
homeostatic <- c("PTPRC", "ITGAM", "P2RY12")
resident_like <- c("PLCL1", "CEBPD", "VIM")
activated <- c("CX3CR1", "GPR34", "FOXP2")
phagocytic <- c("GPNMB", "RGCC", "PPARG")
stressed <- c("HSPA1A", "HSPA1B", "SORCS2")
infg <- c("IFIT2", "IFIT3", "STAT1")
inflammatory <- c("CCL4", "IL1B", "CCL3")
inflammatory_icam <- c("RELB", "ICAM1", "TNFAIP3")
anti_inflammatory <- c("CD163", "SELENOP", "MRC1")
infiltrating <- c("F13A1", "LYVE1", "CD36")

#hypoxic <- c("HIF1A-AS2", "SLC16A10", "ANGPTL4", "VEGFA")

list.microglia.markers <- list(
  "Homeostatic" = homeostatic,
  "Resident-like" = resident_like,
  "Activated" = activated,
  "Phagocytic" = phagocytic,
  "Stressed" = stressed,
  "IFNg" = infg,
  "Inflammatory" = inflammatory,
  "Inflammatory ICAM1+" = inflammatory_icam,
  "Anti-inflammatory" = anti_inflammatory,
  "Infiltrating" = infiltrating
)
color_low <-"#001FA9"
color_high <- "#00A6A9"
p <- SCpubr::do_DotPlot(sample, features = list.microglia.markers) + ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
SCpubr::save_Plot(p = p,
                  figure_path = "/b06x-isilon/b06x-g/G703/eblanco/projects/test_figures",
                  create_path = TRUE,
                  file_name = "Microglia_marker_assignment",
                  width = 18,
                  height = 5,
                  dpi = 300)

p1 <- SCpubr::do_DotPlot(sample[, sample$subtype == "Oligodendroglioma"], features = list.microglia.markers, ylab = "Oligodendroglioma", colors.use = c("grey75", "#3c5b8b")) +
  ggplot2::theme(strip.text = ggplot2::element_text(size = 8)) +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.ticks") +
  ggpubr::rremove("x.text")

p1$guides$colour$title <- "Average Expression OD"
p1$guides$size$order <- 1
p1$guides$colour$order <- 2
p1$guides$colour$label.hjust <- 1
p1$guides$colour$label.hjust <- 0

p2 <- SCpubr::do_DotPlot(sample[, sample$subtype == "Astrocytoma"], features = list.microglia.markers, ylab = "Astrocytoma", colors.use = c("grey75", "#b38b14")) +
  ggplot2::theme(strip.text = ggplot2::element_blank(),
                 strip.background = ggplot2::element_blank())

p2$guides$colour$title <- "Average Expression AS"
p2$guides$size$order <- 1
p2$guides$colour$order <- 2
p2$guides$colour$label.hjust <- 1
p2$guides$colour$label.hjust <- 0

p <- p1 / p2
SCpubr::save_Plot(p = p,
                 figure_path = "/b06x-isilon/b06x-g/G703/eblanco/projects/test_figures",
                 create_path = TRUE,
                 file_name = "Microglia_marker_assignment_by_subtype",
                 width = 22,
                 height = 10,
                 dpi = 300)

# Unlisted mmm
p <- SCpubr::do_DotPlot(sample, features = unlist(list.microglia.markers), split.by = "subtype")


# ATAC extra:


#-------------------------------------------------------------------------------
#B. Astrocytoma
#B4. Motif Analysis
#------------------------------------------------------------------------------


sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use.rds")
Seurat::DefaultAssay(sample.astro) <- 'chromvar'
sample.astro <- sample.astro[, sample.astro$labelling %in% c("Astro-like", "OPC-like", "Gradient", "RA")]
sample.astro$labelling <- factor(sample.astro$labelling, levels = c("Astro-like", "Gradient", "OPC-like", "RA"))

# Find DA peaks for the tumor bulk.
da_regions <- Seurat::FindAllMarkers(object = sample.astro,
                                     only.pos = TRUE,
                                     min.pct = 0.1,
                                     test.use = 'LR',
                                     latent.vars = 'nCount_ATAC')

# Return top 50 DA peaks.
top50 <- da_regions %>%
         dplyr::group_by(.data$cluster) %>%
         dplyr::top_n(n = 50, wt = .data$avg_log2FC)

# Return the averaged data.
avg_data <- as.data.frame(Seurat::AverageExpression(sample.astro,
                                                    assay = "chromvar",
                                                    slot = "data",
                                                    return.seurat = F)[["chromvar"]])
avg_data$motifs <- rownames(avg_data)

# Subset the matrix for only the DA peaks.
motifs_use <- avg_data$motifs[avg_data$motifs %in% top50$gene]
avg_data <- avg_data[motifs_use, ]
avg_data$motifs <- NULL
avg_data <- as.matrix(avg_data)

# Compute correlation.
cc <- stats::cor(avg_data)

corrplot::corrplot(cc,
                   col = rev(corrplot::COL2('RdBu', 200)),
                   method='color',
                   addCoef.col = NULL,
                   tl.col = 'black',
                   bg = "white",
                   diag = T,
                   outline = "grey75",
                   addgrid.col = "black",
                   type = "lower")

motifs_names <- data.frame(BiocGenerics::sapply(rownames(avg_data), function(x) strsplit(x, "-", fixed=TRUE)[[1]][3]))[,1]
rownames(avg_data) <- motifs_names

pdf(file="top50motifs_variable_oligo.pdf")
pheatmap::pheatmap(as.matrix(avg_data), scale='row', cluster_cols=FALSE, col = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize = 6)
dev.off()


# Motif analysis.
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
pfm <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE, matrixtype = "PFM")
)

sample.astro <- Signac::AddMotifs(
  object = sample.astro,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "ATAC"
)
saveRDS(sample.astro, "/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use_TB_PFM.rds", compress = F)
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use_TB_PFM.rds")

# Sample with chromvar.
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/inma_obj/figure_2/rds_figure2/astrocytoma_atac_ready_to_use_TB_PFM_CV.rds")

pfm <- Signac::GetMotifData(object = sample.astro, assay = "ATAC", slot = "pwm")
ggseqlogo::ggseqlogo(data = pfm[["ENSG00000111046_LINE137_MYF6_D_N1"]],
                     method = "prob")

Signac::MotifPlot(object = sample.astro,
                  motifs = "ENSG00000111046_LINE137_MYF6_D_N1",
                  assay = 'ATAC',
                  stack_width = 0.95)

FeaturePlot(
  object = seurat,
  features = "ENSG00000111046-LINE137-MYF6-D-N1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1)

FeaturePlot(
  object = seurat,
  features = "ENSG00000130522-LINE411-JUND-D-N6",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1)

FeaturePlot(
  object = seurat,
  features = "ENSG00000130816-LINE1722-DNMT1-D",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1)


# GO? KEGG terms analysis.

do_invididual_Cluster_report <- function(list_name, genes){
  library(enrichR)
  # Set the search to Human genes.
  setEnrichrSite("Enrichr")
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  dbs_use <- c("Azimuth_Cell_Types_2021",
               "CellMarker_Augmented_2021",
               "PanglaoDB_Augmented_2021",
               "Descartes_Cell_Types_and_Tissue_2021",
               "MSigDB_Hallmark_2020",
               "GO_Biological_Process_2021",
               "GO_Molecular_Function_2021",
               "KEGG_2021_Human")

  enriched <- enrichr(genes, dbs_use)

  list_enrichr <- list()
  for (database in names(enriched)){
    data <- enriched[[database]] %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Count = {length(unlist(stringr::str_split(.data$Genes, ";")))}) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Adjusted.P.value) %>%
      dplyr::select(c("Adjusted.P.value", "Term", "Count")) %>%
      dplyr::slice_head(n = 10) %>%
      dplyr::mutate(Term = ifelse(nchar(Term) >= 60, modify_string(Term), Term)) %>%
      dplyr::mutate(Term = factor(Term))
    p <- ggplot2::ggplot(data, mapping = ggplot2::aes(x = Count, y = Term, color = Adjusted.P.value))  +
      ggplot2::scale_color_continuous(trans = 'reverse') +
      ggpubr::theme_pubr(legend = "right") +
      ggplot2::geom_point(size = 4) +
      geom_point(shape = 1,size = 4,colour = "black") +
      ggplot2::xlab("Genes supporting the Term") +
      ggplot2::ylab("Enriched Term") +
      ggplot2::ggtitle(database) +
      ggplot2::theme(axis.text = ggplot2::element_text(face = "bold"),
                     axis.title = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text(face = "bold"),
                     legend.title =  ggplot2::element_text(face = "bold"))

    list_enrichr[[database]] <- p
  }

  fig <- patchwork::wrap_plots(list_enrichr, ncol = 4, byrow = TRUE)  &
    ggplot2::theme(legend.text = ggplot2::element_text(size = 14)) & patchwork::plot_annotation(title = list_name)
  return(fig)
}



# SUPPLEMENTARY

# Supp bar plots primary RNA samples.
sample.oligo$orig.ident <- stringr::str_remove_all(sample.oligo$orig.ident, "OE0145-")
sample.astro$orig.ident <- stringr::str_remove_all(sample.astro$orig.ident, "OE0145-")
sample.oligo$subtype <- "Oligodendroglioma"
sample.astro$subtype <- "Astrocytoma"
# Number of cells per sample.
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         features = "orig.ident",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = scale.ident,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per sample",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         features = "orig.ident",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = scale.ident,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per sample",
                         plot.title = "Astrocytoma")

p <- p1 | p2

# Number of cells per identity - fill
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         group.by = "orig.ident",
                         features = "subtype",
                         horizontal = T,
                         colors.use = scale.ident,
                         legend = T,
                         legend.ncol = 1,
                         legend.position = "right",
                         position = "fill") +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.text") +
  ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         group.by = "orig.ident",
                         features = "subtype",
                         horizontal = T,
                         colors.use = scale.ident,
                         legend = T,
                         legend.position = "right",
                         position = "fill",
                         legend.ncol = 1,
                         ylab = "Proportion of cells per sample")

p <- p1 / p2


# Number of cells per identity
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         features = "New_NMF_labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         features = "New_NMF_labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Astrocytoma")

p <- p1 | p2


# Number of cells per identity - fill
p1 <- SCpubr::do_BarPlot(sample.oligo,
                         group.by = "New_NMF_labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = F,
                         legend.position = "right",
                         position = "fill") +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.text") +
  ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro,
                         group.by = "New_NMF_labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = T,
                         legend.position = "bottom",
                         position = "fill",
                         legend.nrow = 2,
                         ylab = "Proportion of cells per cell type")

p <- p1 / p2

# Number of cells per identity in the tumor bulk
p1 <- SCpubr::do_BarPlot(sample.oligo[, sample.oligo$New_NMF_labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         features = "New_NMF_labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Oligodendroglioma")
p2 <- SCpubr::do_BarPlot(sample.astro[, sample.astro$New_NMF_labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         features = "New_NMF_labelling",
                         horizontal = T,
                         add.summary_labels = T,
                         colors.use = cluster_cols,
                         legend = F,
                         size.labels = 5,
                         ylab = "Number of cells per cell type",
                         plot.title = "Astrocytoma")

p <- p1 | p2

# Number of cells per identity in the tumor bulk - fill
p1 <- SCpubr::do_BarPlot(sample.oligo[, sample.oligo$New_NMF_labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         group.by = "New_NMF_labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = F,
                         legend.position = "right",
                         position = "fill") +
  ggpubr::rremove("x.axis") +
  ggpubr::rremove("x.text") +
  ggpubr::rremove("x.ticks")
p2 <- SCpubr::do_BarPlot(sample.astro[, sample.astro$New_NMF_labelling %in% c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling-like")],
                         group.by = "New_NMF_labelling",
                         features = "subtype",
                         horizontal = T,
                         colors.use = cluster_cols,
                         legend = T,
                         legend.position = "bottom",
                         position = "fill",
                         legend.nrow = 1,
                         ylab = "Proportion of cells per cell type")

p <- p1 / p2

# Suva Markers
markers <- readxl::read_excel("/omics/odcf/analysis/OE0145_projects/idh_gliomas/scripts/main/06_scana/marker_genes/suva_markers.xlsx")
markers.suva <- list("Astro program" = markers$Astro_program[!is.na(markers$Astro_program)],
                     "Oligo program" = markers$Oligo_program[!is.na(markers$Oligo_program)],
                     "Stemness program" = markers$Stemness_program[!is.na(markers$Stemness_program)])

p.oligo <- SCpubr::do_EnrichmentHeatmap(sample = sample.oligo, list_genes = markers.suva, group.by = "seurat_clusters", row_title = "", column_title = "OD", transpose = T)
p.astro <- SCpubr::do_EnrichmentHeatmap(sample = sample.astro, list_genes = markers.suva, group.by = "seurat_clusters", row_title = "", column_title = "AS", transpose = T)

# Enrichment scores.
sample.oligo <- Seurat::AddModuleScore(sample.oligo, features = list(markers.suva$`Oligo program`), name = "Oligo_program")
sample.astro <- Seurat::AddModuleScore(sample.astro, features = list(markers.suva$`Oligo program`), name = "Oligo_program")

sample.oligo <- Seurat::AddModuleScore(sample.oligo, features = list(markers.suva$`Astro program`), name = "Astro_program")
sample.astro <- Seurat::AddModuleScore(sample.astro, features = list(markers.suva$`Astro program`), name = "Astro_program")

sample.oligo <- Seurat::AddModuleScore(sample.oligo, features = list(markers.suva$`Stemness program`), name = "Stemness_program")
sample.astro <- Seurat::AddModuleScore(sample.astro, features = list(markers.suva$`Stemness program`), name = "Stemness_program")


p.oligo.1 <- SCpubr::do_FeaturePlot(sample.oligo, features = "Oligo_program1", plot.title = "OD - Oligo program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)
p.oligo.2 <- SCpubr::do_FeaturePlot(sample.oligo, features = "Astro_program1", plot.title = "OD - Astro program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)
p.oligo.3 <- SCpubr::do_FeaturePlot(sample.oligo, features = "Stemness_program1", plot.title = "OD - Stemness program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)

p.astro.1 <- SCpubr::do_FeaturePlot(sample.astro, features = "Oligo_program1", plot.title = "AS - Oligo program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)
p.astro.2 <- SCpubr::do_FeaturePlot(sample.astro, features = "Astro_program1", plot.title = "AS - Astro program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)
p.astro.3 <- SCpubr::do_FeaturePlot(sample.astro, features = "Stemness_program1", plot.title = "OD - Stemness program", fontsize = 18, raster = T, pt.size = 1.5, raster.dpi = 1024)

patch <- (p.oligo.1 | p.oligo.2 | p.oligo.3) /(p.astro.1 | p.astro.2 | p.astro.3)

SCpubr::save_Plot(plot = p.oligo,
                 figure_path = paste0(figure_path, "/Supplementary/Suva_markers/"),
                 file_name = "Suva_markers_OD_heatmap",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = p.astro,
                 figure_path = paste0(figure_path, "/Supplementary/Suva_markers/"),
                 file_name = "Suva_markers_AS_heatmap",
                 height = 8,
                 width = 8,
                 dpi = 300)

SCpubr::save_Plot(plot = patch,
                 figure_path = paste0(figure_path, "/Supplementary/Suva_markers/"),
                 file_name = "Suva_markers_OD_and_AS_enrichment_scores",
                 height = 16,
                 width = 24,
                 dpi = 300)

# NMF metasignatures.
nmf <- as.list(read.table("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/final_marker_genes/IDH_gliomas_NMF_metasignatures_oligodendroglioma.tsv", header = T))
names(nmf) <- stringr::str_remove_all(names(nmf), "O_")


for (number_comparisons in c(1, 4)){
  message(number_comparisons)
  for (list.name in names(nmf)){
    message(list.name)
    out.oligo <- SCpubr:::do_PTEA(sample = sample.oligo, markers = nmf, list.name = list.name, idents.use = c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling"), number_comparisons = number_comparisons)
    out.astro <- SCpubr:::do_PTEA(sample = sample.astro, markers = nmf, list.name = list.name, idents.use = c("OPC-like", "Astro-like", "Gradient", "RA", "Cycling"), number_comparisons = number_comparisons)


    if (number_comparisons == 1){
      method <- "single_comparisons"
    } else if (number_comparisons == 4){
      method <- "multiple_comparisons"
    }
    SCpubr::save_Plot(plot = out.oligo$p.umap,
                     figure_path = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name),
                     file_name = paste0("OD_", list.name, "_umap"),
                     height = 5,
                     width = 10,
                     dpi = 300)
    SCpubr::save_Plot(plot = out.oligo$p.dist,
                     figure_path = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name),
                     file_name = paste0("OD_", list.name, "_distributions"),
                     height = 7,
                     width = 10,
                     dpi = 300)


    saveRDS(out.oligo$surpassed_cells, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_surpassing_cells.rds"))
    write.table(out.oligo$surpassed_cells, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_surpassing_cells.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
    saveRDS(out.oligo$empirical_data, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_empirical_data.rds"))
    write.table(out.oligo$empirical_data, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_empirical_data.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
    saveRDS(out.oligo$null_data, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_null_data.rds"))
    write.table(out.oligo$null_data, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/OD_", list.name, "_null_data.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")


    SCpubr::save_Plot(plot = out.astro$p.umap,
                     figure_path = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name),
                     file_name = paste0("AS_", list.name, "_umap"),
                     height = 5,
                     width = 10,
                     dpi = 300)
    SCpubr::save_Plot(plot = out.astro$p.dist,
                     figure_path = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name),
                     file_name = paste0("AS_", list.name, "_distributions"),
                     height = 7,
                     width = 10,
                     dpi = 300)


    saveRDS(out.astro$surpassed_cells, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_surpassing_cells.rds"))
    write.table(out.astro$surpassed_cells, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_surpassing_cells.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
    saveRDS(out.astro$empirical_data, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_empirical_data.rds"))
    write.table(out.astro$empirical_data, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_empirical_data.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")
    saveRDS(out.astro$null_data, paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_null_data.rds"))
    write.table(out.astro$null_data, file = paste0(figure_path, "/Supplementary/permutation_approach/", method, "/", list.name, "/AS_", list.name, "_null_data.tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

  }
}

# Supplementary - Correlation of clusters and final identities regarding HVG


p1 <- SCpubr::do_CorrelationPlot(sample = sample.oligo, group.by = "seurat_clusters")
p2 <- SCpubr::do_CorrelationPlot(sample = sample.astro, group.by = "seurat_clusters")

p3 <- SCpubr::do_CorrelationPlot(sample = sample.oligo, group.by = "New_NMF_labelling")
p4 <- SCpubr::do_CorrelationPlot(sample = sample.astro, group.by = "New_NMF_labelling")

SCpubr::save_Plot(plot = p1,
                 figure_path = paste0(figure_path, "/Supplementary/HVG_correlation/"),
                 file_name = "OD_HVG_correlation_seurat_clusters",
                 height = 10,
                 width = 10,
                 dpi = 300)

SCpubr::save_Plot(plot = p2,
                 figure_path = paste0(figure_path, "/Supplementary/HVG_correlation/"),
                 file_name = "AS_HVG_correlation_seurat_clusters",
                 height = 10,
                 width = 10,
                 dpi = 300)

SCpubr::save_Plot(plot = p3,
                 figure_path = paste0(figure_path, "/Supplementary/HVG_correlation/"),
                 file_name = "OD_HVG_correlation_final_clusters",
                 height = 10,
                 width = 10,
                 dpi = 300)

SCpubr::save_Plot(plot = p4,
                 figure_path = paste0(figure_path, "/Supplementary/HVG_correlation/"),
                 file_name = "AS_HVG_correlation_final_clusters",
                 height = 10,
                 width = 10,
                 dpi = 300)


# DimPlots merged samples.
sample.oligo <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_merged_oligodendroglioma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_merged_oligodendroglioma_labelled")
sample.astro <- readRDS("/omics/odcf/analysis/OE0145_projects/idh_gliomas/results/analysis/OE0145-IDH_merged_astrocytoma/snRNAseq/10x_3_prime_v3/saved_objects/RNBR/regressed/OE0145-IDH_merged_astrocytoma_labelled")


p.oligo <- SCpubr::do_DimPlot(sample = sample.oligo,
                              group.by = "first_labelling",
                              colors.use = cluster_cols,
                              legend.position = "bottom",
                              legend.ncol =  3,
                              raster = T,
                              pt.size = 1,
                              raster.dpi = 1024,
                              plot.title = "Oligodendroglioma") + ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5))

SCpubr::save_Plot(plot = p.oligo,
                  figure_path = paste0(figure_path, "/Supplementary/merged_samples_DimPlot/"),
                  file_name = "OD",
                  height = 8,
                  width = 7,
                  dpi = 300)

p.astro <- SCpubr::do_DimPlot(sample = sample.astro,
                              group.by = "first_labelling",
                              colors.use = cluster_cols,
                              legend.position = "bottom",
                              legend.ncol =  3,
                              raster = T,
                              pt.size = 1,
                              raster.dpi = 1024,
                              plot.title = "Astrocytoma") + ggplot2::theme(plot.title = ggtext::element_markdown(hjust = 0.5))

SCpubr::save_Plot(plot = p.astro,
                  figure_path = paste0(figure_path, "/Supplementary/merged_samples_DimPlot/"),
                  file_name = "AS",
                  height = 8,
                  width = 7,
                  dpi = 300)
