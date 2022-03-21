# Enrique Blanco Carmona
# e.blancocarmona@kitz-heidelberg.de
# PhD Student – Clinical Bioinformatics
# Division of Pediatric Neurooncology (B062)




#--------------------------------------------------------------------
# 1 - READ IN COUNT MATRICES AND INDIVIDUAL QC
#--------------------------------------------------------------------

# ------------------------------------------------------------------
# All the following code including doublet detection has to be run
# for each sample individually.
# ------------------------------------------------------------------

# Global parameters. Used throughout the whole analysis.
sample_name <- "" # Name of the sample.

# Load samples.
mt_pattern <- "^MT-"
sample_path <- "" # Path to where the count matrix from cellranger (either filtered or raw, your choice) is located.

# Read the count matrix and generate a Seurat object.
sample <- Seurat::Read10X(data.dir = sample_path) %>%
      		Seurat::CreateSeuratObject(project = sample_name, min.cells = 3, min.features = 200)

# Compute percentage of mitochondrial RNA in the cells.
sample[["percent.mt"]] <- Seurat::PercentageFeatureSet(object = sample,
                                                       pattern = mt_pattern)

# Perform QC on the sample.
# Get cutoffs.
counts_lower_cutoff <- 1000 # Minimum amount of UMIs per cell.
genes_lower_cutoff <- 500 # Minimum amount of genes per cell.
mito_higher_cutoff <- 5 # Maximum amount of mitochondrial RNA per cell.

# Generate the first subset.
count_mask <- sample$nCount_RNA > counts_lower_cutoff
gene_mask <- sample$nFeature_RNA > genes_lower_cutoff
mito_mask <- sample$percent.mt < mito_higher_cutoff
mask <- count_mask & gene_mask & mito_mask
sample <- sample[, mask] # Subset the sample.

# Upper cutoffs are determined by the mean and standard deviation of the remaining cells.
counts_higher_cutoff <- mean(sample$nCount_RNA) + 3 * stats::sd(sample$nCount_RNA)
genes_higher_cutoff <- mean(sample$nFeature_RNA) + 3 * stats::sd(sample$nFeature_RNA)


# Second subset based on the mean and standard deviation of the remaining cells.
count_mask <- sample$nCount_RNA < counts_higher_cutoff
gene_mask <- sample$nFeature_RNA < genes_higher_cutoff
mask <- count_mask & gene_mask
sample <- sample[, mask] # Subset the sample.


#--------------------------------------------------------------------
# 2 - DOUBLET DETECTION WITH SCRUBBLET
#--------------------------------------------------------------------

# For this part we have scrublet installed in a conda environment.
conda_env <- "" #Path to your conda environment.
reticulate::use_condaenv(conda_env)

# Import scrublet.
scrublet <- reticulate::import("scrublet")

# Transpose the count matrix.
counts_transposed <- Matrix::t(sample@assays$RNA@counts)

# Run scrublet.
# Code adapted from: https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb
scrub = scrublet$Scrublet(counts_transposed, expected_doublet_rate = 0.06)

# Compute the doublets.
return_list = scrub$scrub_doublets() # List with the output from scrublet.
scrublet_score <- return_list[[1]] # Scrublet scores per cell.
scrublet_binary <- return_list[[2]] # Scrublet assignment for each cell.

# Add cell names to the output, so it can be integrated in the Seurat object.
row.names(scrublet_score) <- colnames(sample)
row.names(scrublet_binary) <- colnames(sample)

# Add the output as metadata.
sample$scrublet_score <- scrublet_score
sample$scrublet_binary <- scrublet_binary

# Visualize the doublet scores and the assignment.
h <- graphics::hist(sample$scrublet_score, breaks = "FD") # To compute the breaks according to Freedman–Diaconis rule. https://stats.stackexchange.com/a/383145
p <- ggplot2::ggplot(sample[[]], ggplot2::aes(x = sample$scrublet_score, fill = sample$scrublet_binary)) +
     ggplot2::geom_histogram(breaks = h$breaks) +
     ggplot2::scale_fill_manual(values = colortools::opposite("steelblue")) +
     ggplot2::geom_vline(ggplot2::aes(xintercept = mean_binary), colour = "grey", linetype = "dashed") +
     ggpubr::theme_pubclean()
p$labels$fill <- "Doublet assignment"
p$labels$y <- "Number of nuclei"
p$labels$x <- "scrublet score"
p$labels$subtitle <- sprintf("Binary prediction of doublets failed: %s", was_null)
p$theme$legend.position <- "bottom"
p$labels$subtitle <- sprintf("Cutoff: %s\t Predicted doublets: %s\t Number of singlets: %s\t Percentage of doublets: %s",
                                 mean_binary,
                                 sum(sample$scrublet_binary),
                                 sum(!sample$scrublet_binary),
                                 round(sum(sample$scrublet_binary) / length(colnames(sample)), 3))

# At this point, it might be the case that the cutoff decided by scrublet is suboptimal.
# Therefore, it might make more sense to decide it yourself based on the histogram.

doublet_cutoff <- 0.2 # Put your own value.

# Modify the binary assignment accordingly.
sample$scrublet_binary <- ifelse(sample$scrublet_score > doublet_cutoff, TRUE, FALSE)
# You can re-run the histogram above to get the new visualization and metrics.

# Generate a reporting df.
report_df <- data.frame(number_of_cells_after_qc = length(colnames(sample)),
						number_of_predicted_doublets = sum(sample$scrublet_binary),
        				cutoff = mean_binary,
        				percentage_of_predicted_doublets_in_sample = sum(sample$scrublet_binary) * 100 / length(colnames(sample)),
        				prediction_doublet_failed = was_null)

# Save individual samples.
output_path <- "" # Path where the samples will be stored.

saveRDS(sample, paste0(output_path, "/", sample_name, "_after_QC.rds"))

#--------------------------------------------------------------------
# 3 - MERGE INDIVIDUAL SAMPLES TOGETHER
#--------------------------------------------------------------------

# Generate a list containing all samples.
list_samples <- list()
sample_names_vector <- c() # Vector containing all sample names.

for (sample_name_use in sample_names_vector){
	sample.individual <- readRDS(paste0(output_path, "/", sample_name_use, "_after_QC.rds"))
	list.samples[[sample_name_use]] <- sample.individual
}

# Generate the merged sample.
sample_name <- "" # Name for the merged sample.
sample <- merge(list_samples[[1]],
		        y = list_samples[2: length(list_samples)],
		        add.cell.ids = names(list_samples),
		        project = sample_name)

merged_sample@meta.data$orig.ident <- stringr::str_replace_all(merged_sample@meta.data$orig.ident, "\\.", "-")

# Exclude doublets.
message(print0("Total number of doublets in the sample: ", sum(sample$scrublet_binary)))
sample <- sample[, sample$scrublet_binary == FALSE]

#--------------------------------------------------------------------
# 4 - NORMALIZATION USING SCTRANSFORM AND CELL CYCLE SCORING
#--------------------------------------------------------------------

# In this process, we want to regress out the effect of 5 variables:
# - nCount_RNA: UMIs. We observe different numbers across samples.
# - nFeature_RNA: genes. We observe different numbers across samples.
# - percent.mt: % of mitochondrial RNA in the cells.
# - S.Score and G2M.Score: Cell cycle effect. We want to be sure this is not a bias in the data.

# Get S phase genes.
s.genes <- Seurat::cc.genes.updated.2019$s.genes
# Get G2-M phase genes.
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes


# We normalize first using SCTransform and removing the effect of UMIs, genes and percent.mt.
# Process followed as in https://github.com/satijalab/seurat/issues/1679#issuecomment-557781838

sample <- Seurat::SCTransform(sample,
                  			  assay = "RNA",
                              new.assay.name = "SCT",
                              vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"))


# With the data normalized, we asses the cell cycle effect.
sample <- Seurat::CellCycleScoring(object = sample,
                                   s.features = s.genes,
                                   g2m.features = g2m.genes,
                                   set.ident = FALSE)


# And then, since now S.Score and G2M.Score are stored as metadata variables, we can normalize again using also these variables in the regress out.

sample <- Seurat::SCTransform(sample,
                  			  assay = "RNA",
                              new.assay.name = "SCT",
                              vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"))
# Please note that this second normalization uses the data in the "RNA" assay, this is, the unnormalized data.
# Therefore, the data we will work on is only normalized once.

#--------------------------------------------------------------------
# 5 - DIMENSIONAL REDUCTION
#--------------------------------------------------------------------

sample <- Seurat::RunPCA(sample, assay = "SCT")
elbow <- Seurat::ElbowPlot(sample) # Visualize the number of PCs and select an appropriate number.
npcs <- 15 # Selected number of PCs.

# Normalize using Harmony. The effect to remove is the orig.ident metadata column (from where each sample comes from).
sample <- harmony::RunHarmony(sample, dims = 1:npcs, group.by.vars = "orig.ident", assay.use = "SCT")

# Dimensional reduction: UMAP over harmony reduction.
sample <- Seurat::RunUMAP(sample, dims = 1:npcs, reduction = "harmony")

sample <- Seurat::FindNeighbors(sample, dims = 1:npcs, reduction = "harmony")

sample <- Seurat::FindClusters(sample,
							   random.seed = 42)



#--------------------------------------------------------------------
# 6 - CLUSTER LABELLING
#--------------------------------------------------------------------

# Samples are labelled according to external lists of marker genes:
# TME: PanglaoDB: https://panglaodb.se/index.html
# Tumor bulk: https://pubmed.ncbi.nlm.nih.gov/27806376/
# A table with the annotation for each cell will be provided as supplementary data.



#----------------------------------------------------------------------------
# 7 - CNV PROFILING FOR 10X SAMPLES GENERATING METACELLS FOR INCREASED SIGNAL
#----------------------------------------------------------------------------

# CNV profiling is done using inferCNV package.
# In addition, here we provide the code for the generation of metacells, to further increase the resolution of the output, which can be quite noisy for 10X datasets.


# First of all, we want to get rid of any clusters with less than 10 cells.
# Also, we assume that the current labelling is stored in the metadata variable "New_NMF_labelling" and this is also use as Idents in the seurat object.
Seurat::Idents(sample) <- sample$New_NMF_labelling


clusters_keep <- c()
for (cluster in levels(sample)){
    num <- sum(sample$New_NMF_labelling == cluster)
    if (num >= 10){
        clusters_keep <- append(clusters_keep, cluster)
    }
}

# Subset the sample.
sample <- sample[, sample$New_NMF_labelling %in% clusters_keep]


# We need to retrieve the count matrix. For this, we generate metacells as follows:

# Generate a new metadata column storing the mapping cell-metacell.
sample[["metacell_mapping"]] <- "not_mapped"

# Will store the complete annotation for the metacells.
whole_annotation <- data.frame(cluster_names = "test", row.names = "test")

meta_counter <- 0 # To keep a count of the metacells that are created.


for (cluster_id in levels(sample)){
	print(sprintf("Computing metacells for cluster %s.", cluster_id))
	# Will store the metacells per cluster.
	metacells <- data.frame(test = rownames(sample), row.names = rownames(sample))

	# Subset the sample by each cluster ID.
	chunksample <- sample[, sample$New_NMF_labelling == cluster_id]

	# Get the count data as a data frame and transpose it so columns are GENES and rows are CELLS.
	countdata <- t(as.data.frame(Seurat::GetAssayData(chunksample, "counts")))

	# Get the possible amount of metacells.
	times <- trunc(dim(countdata)[1] / metacell_content)

	for (i in seq(1,times)){
		meta_counter <- meta_counter + 1
		# Generate slice points for each metacell.
		start <- ((i -1) * metacell_content + 1)
		end <- i * metacell_content


		# Compute the slice as a data frame containing the sum of the subsetted cells. dims = 1 row (metacell), X columns (genes)
		slice <- as.data.frame(colSums(countdata[start:end, ]))

		# Get the name of the cells merged.
		cell_names <- rownames(countdata[start:end, ])

		# Add the metacell.
		col_name <- sprintf("metacell_%s", meta_counter)
		metacells[[col_name]] <- slice[,1]

		# Add the mapping.
		sample$metacell_mapping[colnames(sample) %in% cell_names] <- col_name
	}

	# Delete the test column as we already have more than 1 column in our data frame.
	metacells[["test"]] <- NULL

	# Will contain the annotation of the generated metacells. Columns: cluster identities. Rows: each metacell.
	annotation <- data.frame(cluster_names = colnames(metacells), row.names = colnames(metacells))
	# Replace the dummy cluster_names column's values for the actual label for the cluster.
	annotation$cluster_names <- cluster_id

	# Add the annotation data and the metacell data to the global containers. In the end: # Columns for metacell object = # rows for annotation object.
	whole_metacells <- cbind(whole_metacells, metacells)
	whole_annotation <- rbind(whole_annotation, annotation)
}

# Turn the names into characters for the sake of avoiding errors when subsetting.
whole_annotation$cluster_names <- as.character(whole_annotation$cluster_names)

# Delete the test row from the global annotation data.
whole_annotation <- whole_annotation[!rownames(whole_annotation) %in% c("test"), , drop = FALSE]

# Delete the test column from the global metacell data.
whole_metacells$test <- NULL

cnv_analysis_folder <- "" # Path to store the output of inferCNV
annotation_file <- sprintf("%s/annotation_metacells.tsv", cnv_analysis_folder)
# Save the annotation object.
utils::write.table(whole_annotation,
	               file = annotation_file,
	               sep = "\t",
	               row.names = TRUE,
	               col.names = FALSE,
	               quote = FALSE)

# Return the metacell object as a matrix (required for running inferCNV).
count_matrix <- as.matrix(whole_metacells)

# You might also want to save the Seurat object with the metacell mapping annotation.
saveRDS(sample, "") # Modify accordingly.

# Run inferCNV:
gene_ordering_file <- "" # Path to where the file with the order of the genes is stored. It can also be downloaded here: https://data.broadinstitute.org/Trinity/CTAT/cnv/
ref_clusters <- "" # Which clusters to use as a reference.

# Create the inferCNV object.
infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = count_matrix,
                                               annotations_file = annotation_file,
                                               delim = "\t",
                                               gene_order_file = gene_ordering_file,
                                               ref_group_names = ref_clusters)

# Run inferCNV.
cnv_analysis_folder_output <- "" # This path needs to not exist in your filesystem otherwise inferCNV will stop complaining that the path exists.
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              min_cells_per_gene = 3, # Default.
                              out_dir = cnv_analysis_folder,  # dir is auto-created for storing outputs
                              cluster_by_groups = TRUE,   # Cluster by groups.
                              denoise = TRUE,
                              HMM = TRUE,
                              HMM_type = "i6",
                              window_length = 201,
                              num_threads = 8,
                              resume_mode=FALSE)

# For further denoising, you can apply median filtering, but this might remove true signal from the plot.
infercnv_obj_median_filtered = infercnv::apply_median_filtering(infercnv_obj)

infercnv::plot_cnv(infercnv_obj_median_filtered,
                   out_dir = cnv_analysis_folder_output,
                   output_filename = 'infercnv.median_filtered',
                   x.range = "auto",
                   x.center = 1,
                   title = "infercnv",
                   color_safe_pal = FALSE)


# ------------------------------------------------------------------
# 8 - NMF ANALYSIS ON THE TUMOR BULK
# ------------------------------------------------------------------

# Notes:
# The first half of this code, up to individual NMF signature generation and also the function for
# scoring metasignatures has been adapted from Volker Hovestadt's code.
# Methodology explained in: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6754173/
#

# The following code is strongly suggested to be run in a cluster, for each sample individually as a job.
# It can take easily 100+ GB of RAM for each unique sample in sample$orig.ident.

nmf_output_folder <- "" # Folder to store the results.
sample_index <- 1 # This will select the sample out of unique(sample$orig.ident). Modify this parameter as you run it on the cluster.


tumor_clusters <- c() # Vector with the clusters that form the tumor bulk.
sample.tumor <- sample[, sample$New_NMF_labelling %in% tumor_clusters]

# Center the genes over malignant cells.
sample.tumor <- Seurat::ScaleData(sample.tumor)

nmf.list <- BiocGenerics::lapply(unique(sample.tumor$orig.ident)[sample_index],
                                   function(samp){
                                       message(samp)

                                       # Subset the tumor population pertaining to the sample.
                                       tum <- sample.tumor[, sample.tumor$orig.ident == samp]

                                       # Re-center the genes for this sample.
                                       tum <- Seurat::ScaleData(tum)

                                       # Extract the count matrix.
                                       cts <- tum@assays$SCT@scale.data

                                       # Set negative values to 0.
                                       cts[cts < 0] <- 0

                                       # Remove rows without expression in this sample.
                                       cts <- cts[rowSums(cts) > 0, ]

                                       # Compute NMF for different K (i.e: signatures).
                                       # From 2 to 10 NMF signatures.
                                       k <- 2:10
                                       nrun <- 10
                                       seed <- 777 # Let's have some luck :P.

                                       # Run NMF.
                                       message(sprintf("Computing NMF for sample: %s", samp))

                                           nmf <- NMF::nmf(x = cts,
                                                           rank = k,
                                                           nrun = nrun,
                                                           seed = seed,
                                                           method = "snmf/r",
                                                           .options = "p4v")
                                           # Save the results.
                                           save(list = c("cts", "nmf"),
                                                file = paste0(nmf_output_folder, "_", samp, "_nmf2to10", ".RData"))

                                   })



# Generate further directories for the results.
nmf_heatmap_summary <- sprintf("%s/signature_summary/", nmf_output_folder)
nmf_feature_plots <- sprintf("%s/feature_plots/", nmf_output_folder)
nmf_signature_heatmaps <- sprintf("%s/signature_heatmaps/", nmf_output_folder)
nmf_go_kegg <- sprintf("%s/signature_GO_KEGG/", nmf_output_folder)
nmf_top30 <- sprintf("%s/top30_genes_per_signature/", nmf_output_folder)
nmf_metaclusters <- sprintf("%s/metaclusters/", nmf_output_folder)

dir.create(nmf_heatmap_summary, recursive = T)
dir.create(nmf_feature_plots, recursive = T)
dir.create(nmf_signature_heatmaps, recursive = T)
dir.create(nmf_go_kegg, recursive = T)
dir.create(nmf_top30, recursive = T)
dir.create(nmf_metaclusters, recursive = T)


# Remove MT genes, as being a snRNAseq experiment they are prone to appear later on as meta-signatures.
# Retrieve MT genes.
mitochondrial_genes <- rownames(tumor_scaled)[grepl("^MT", rownames(tumor_scaled))]
# Filter out ribosomal genes.
tumor_scaled <- tumor_scaled[!(rownames(tumor_scaled) %in% mitochondrial_genes), ]

# Read in the individual NMF signatures.
lf <- list.files(nmf_output_folder,
                         pattern = "nmf2to10.R*",
                         full.names = TRUE)
nmf.list <- BiocGenerics::lapply(lf, function(x){
    load(x)
    nmf
})

fi.names <- lf <- list.files(nmf_output_folder,
                                     pattern = "nmf2to10.R*",
                                     full.names = FALSE)
fi.names <- BiocGenerics::sapply(fi.names, function(x) strsplit(x, "_nmf")[[1]][1])
fi.names <- BiocGenerics::sapply(fi.names, function(x) strsplit(x, "^_")[[1]][2])
names(nmf.list) <- fi.names


message("Analyzing NMF signatures from k = 2 until k = 10.")
times <- 0
only_one_k <- FALSE # Check if only NMF was computed for only one value of K.
for (k in seq(2, 10)){
    message(sprintf("Using k = %s", k))
    message("...")
    times <- times + 1
    # Start a report.
    nmf_heatmap_summary_out <- sprintf("%s/k_%s/", nmf_heatmap_summary, k)
    dir.create(nmf_heatmap_summary_out, recursive = T)

    grDevices::pdf(sprintf("%s/%s_nmf%s.heatmap.pdf", nmf_heatmap_summary, sample_name, k),
                   width = 6,
                   height = 6)
    graphics::par(mar = c(4, 4, 4, 4))

    # Initialize empty lists to store data.
    list.cor <- list()
    list.top30 <- list()
    for (i in names(nmf.list)){
        if (only_one_k == FALSE){
            X.nmf <- nmf.list[[i]][["fit"]][[times]]
        } else if (only_one_k == TRUE){
            X.nmf <- nmf.list[[i]]@fit
        }

        # Basis (metagenes).
        X.nmf.basis <- NMF::basis(X.nmf)

        # Exclude mitochondrial genes.
        X.nmf.basis <- X.nmf.basis[!(rownames(X.nmf.basis) %in% rownames(X.nmf.basis)[grepl("^MT", rownames(X.nmf.basis))]), ]

        colnames(X.nmf.basis) <- paste0("k", seq(ncol(X.nmf.basis)))
        X.nmf.basis.hc1 <- stats::hclust(stats::dist(X.nmf.basis), method = "average")

        # Plot heatmap.
        heatmap(X.nmf.basis,
                Rowv = stats::as.dendrogram(X.nmf.basis.hc1),
                Colv = NA,
                zlim = c(0, max(X.nmf.basis)),
                scale = "n",
                col = grDevices::grey.colors(1000, 1, 0),
                useRaster = TRUE,
                labRow = FALSE,
                main = paste0(i, ", base"))

        # Coef.
        X.nmf.coef <- NMF::coef(X.nmf)
        rownames(X.nmf.coef) <- paste0("k", seq(nrow(X.nmf.coef)))
        X.nmf.coef.hc1 <- stats::hclust(stats::dist(X.nmf.coef), method = "average")
        X.nmf.coef.hc2 <- stats::hclust(stats::dist(t(X.nmf.coef)), method = "average")
        heatmap(X.nmf.coef[rev(seq(nrow(X.nmf.coef))), ],
                Rowv = NA,
                Colv = stats::as.dendrogram(X.nmf.coef.hc2),
                zlim=c(0, max(X.nmf.coef)),
                scale="n",
                col = grDevices::grey.colors(1000, 1, 0),
                useRaster=FALSE,
                labCol=FALSE,
                main=paste0(i, ", coef"))

        # Top genes, just use basis.
        X.nmf.basis.cor <- t(X.nmf.basis)
        list.cor[[i]] <- X.nmf.basis.cor

        X.nmf.basis.top30 <- apply(X.nmf.basis.cor, 1, function(x) colnames(X.nmf.basis.cor)[order(x, decreasing = TRUE)[1:30]])
        list.top30[[i]] <- X.nmf.basis.top30

        Y <- tumor_scaled[as.vector(X.nmf.basis.top30), intersect(colnames(X.nmf.coef), colnames(tumor_scaled))] #
        Y <- Y - rowMeans(Y)  # re-center
        Y.hc <- stats::as.dendrogram(stats::hclust(stats::as.dist(1 - stats::cor(Y)), method = "ward.D"))
        #Y.hc <- stats::reorder(Y.hc, wts = colMeans(tumor_scaled[ribosomal_genes, ])-colMeans(tumor_scaled[c(x1, x2), ]), agglo.FUN = mean)
        #Z <- as.matrix(scaleMinMax(Y , -6, 6))
        Z <- as.matrix(scale(Y))
        heatmap(Z[rev(seq(nrow(Z))), ],
                Rowv = NA,
                Colv = as.dendrogram(Y.hc),
                scale = "n",
                zlim = c(0, 1),
                useRaster = TRUE,
                cexRow = 0.5,
                cexCol = 0.1,
                add.expr = graphics::abline(h = head(cumsum(rep(nrow(X.nmf.basis.top30),
                                                                        ncol(X.nmf.basis.top30))), -1) + 0.5), main=paste0(i, ", ", paste0(dim(Z), collapse="x")))
    }
    grDevices::dev.off()
    save(list = c("list.top30", "list.cor", "nmf.list"),
         file = paste0(nmf_top30, "all_samples_nmf", k, "_signatures.RData"))
    list.top30.cbind <- do.call(cbind, list.top30)
    colnames(list.top30.cbind) <- unlist(BiocGenerics::lapply(names(list.top30), function(i){paste0(i, ".", seq(ncol(list.top30[[i]])))}))

    write.table(list.top30.cbind,
                file = paste0(nmf_top30, "all_samples_nmf", k, "_signatures.txt"),
                sep = "\t",
                quote = FALSE)
}

# Define 2 custom functions.

# Define 2 functions.
## (supposed) input for "scoreSignature":
## X.center - centered expression matrix for all tumor cells
## X.mean - mean expression for each gene across all tumor cells
## s - signature genes

# Function from Volker Hovestadt.
scoreSignature <- function(X.center, X.mean, signature, n = 100, verbose = F) {
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

# Adaptation from Volker Hovestadt function.
# Pretty much the same contept as scoreSignature but with a twist in the end. Instead of doing colMeans to get the average score for each cell for all of the genes in the signature.
# rowMeans to get the averaged score for each gene for all of the cells in the dataset provided.
scoreGene <- function(X.center, X.mean, signature, n = 100, verbose = F) {
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

# The following code creates a list of heatmaps that will help you decide for which value of N (number of NMF signatures) is best.
# This is decided based on the number and quality of the groups in the correlation plot.

# List placeholder for plot outputs.
list_plots <- list()
list_nmf_signatures_heatmaps <- list()

for (k in seq(2,10)){
    # Read in signatures.
    print("Reading signatures.")
    message(paste0("Using k = ", k))
    nmf.sign <- read.table(paste0(nmf_top30, "all_samples_nmf", k, "_signatures.txt"), sep = "\t")

    # Calculate the scores for each signature in each cell.
    X.mean <- rowMeans(tumor_scaled) # Mean expression level of the scale.data dataset of the tumor bulk.
    sign.scores <- apply(nmf.sign, 2, function(signature){scoreSignature(X.center = tumor_scaled,
                                                                         X.mean = X.mean,
                                                                         signature)})

    # Correlate/Cluster signatures over all scores.
    corr <- stats::cor(sign.scores)
    range <- max(abs(corr))

    # Compute the heatmap.
    h <- pheatmap::pheatmap(corr,
                            color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                            breaks = seq(-range, range, length.out = 100),
                            show_colnames = F,
                            treeheight_col = 0)

    # Store the heatmap.
    list_nmf_signatures_heatmaps[[paste0("nmf", k)]] <- h
}

# Modify this values accordingly
ind_k <- 10 # Final k for individual signatures to use.
meta_k <- 4 # Number of meta-clusters to derive from the signatures generated with ind_k.

nmf.sign <- read.table(paste0(nmf_top30, "all_samples_nmf", ind_k, "_signatures.txt"), sep = "\t")
X.mean <- rowMeans(tumor_scaled)
sign.scores <- apply(nmf.sign, 2, function(signature){scoreSignature(X.center = tumor_scaled,
                                                                     X.mean = X.mean,
                                                                     s = signature)})

# Correlate/Cluster signatures over all scores.
corr <- stats::cor(sign.scores)
range <- max(abs(corr))



h <- pheatmap::pheatmap(corr,
                        color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                        breaks = seq(-range, range, length.out = 100),
                        show_colnames = F,
                        treeheight_col = 0,
                        cutree_rows = meta_k,
                        cutree_cols = meta_k)

# Colors for the metasignature annotation.
ann_colors <- colortools::wheel("steelblue", meta_k)
names(ann_colors) <- BiocGenerics::sapply(as.character(seq(1, meta_k)), function(x){paste0("Metasignature_", x)})

# Colors for the sample annotation.
samp_colors <- colortools::wheel("navyblue", length(unique(sample.tumor$orig.ident)))
names(samp_colors) <- unique(sample.tumor$orig.ident)

# Colors for the grade annotation.
subtype_colors <- c("#195A70", "#702F19")
names(subtype_colors) <- c("2", "3")

# List containing all colors.
colors <- list("Metasignature" = ann_colors,
               "Origin" = samp_colors,
               "Grade" = subtype_colors)

# Workaround to modify the names.
metaclust <- stats::cutree(h$tree_row, k = meta_k)
metaclust.s <- as.character(BiocGenerics::sapply(metaclust, function(x){paste0("Metasignature_", x)}))
names(metaclust.s) <- names(metaclust)
metaclust <- metaclust.s

# Annotation values for the sample of origin.
origin <- BiocGenerics::sapply(colnames(corr), function(x){substring(x, 1, nchar(x)-2)})
origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "[.]", "-")})
origin <- BiocGenerics::sapply(origin, function(x){stringr::str_replace_all(x, "-$", "")})

# Annotation vlaues for the grade.
subgroup <- origin
grade_2 <- unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident)
grade_3 <- unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident)

origin <- factor(origin)

subgroup[subgroup %in% grade_2] <- 2
subgroup[subgroup %in% grade_3] <- 3

cellwidth <- 10
cellheight <- 10

# First visualization of the correlation heatmap with annotations.
h.anno <- pheatmap::pheatmap(corr,
                      		 color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                      		 breaks = seq(-range, range, length.out = 100),
                      		 show_colnames = F,
                      		 treeheight_col = 0,
                      		 annotation_col = data.frame("Metasignature" = metaclust,
                      		                             "Origin" = origin,
                      		                             "Subtype" = subgroup),
                      		 annotation_colors = colors,
                      		 cellwidth = cellwidth,
                      		 cellheight = cellheight)

# Copy  to modify the names. This can also be done to subset only given metasignatures if needed (hence, the .small naming).
metaclust_small <- metaclust

# Recompute the correlation.
corr.small <- stats::cor(sign.scores[, names(metaclust_small)])
range.small <- max(abs(corr.small))

# Rename groups if needed.
names_sign <- names(metaclust_small)
assign_sign <- as.character(metaclust_small)
# Rename the signatures accordingly.
assign_sign[assign_sign %in% c("Metasignature_3", "Metasignature_5")] <- "Metasignature_1" # Rename several metasignatures at once.
assign_sign[assign_sign == "Metasignature_6"] <- "Metasignature_3" # Rename individual metasignatures at one.
# Regenerate the metaclust_small object.
names(assign_sign) <- names_sign
metaclust_small <- as.factor(assign_sign)

# Recompute annotation. Same process as before.
ann_colors_small <- colortools::wheel("steelblue", length(levels(metaclust_small)))
names(ann_colors_small) <- levels(metaclust_small)


origin_small <- BiocGenerics::sapply(colnames(corr.small), function(x){substring(x, 1, nchar(x)-2)})
origin_small  <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "[.]", "-")})
origin_small <- BiocGenerics::sapply(origin_small, function(x){stringr::str_replace_all(x, "-$", "")})

subgroup_small <- origin_small
grade_2 <- unique(sample.tumor[, sample.tumor$grade == "2"]$orig.ident)
grade_3 <- unique(sample.tumor[, sample.tumor$grade == "3"]$orig.ident)

origin_small <- factor(origin_small)

subgroup_small[subgroup_small %in% grade_2] <- 2
subgroup_small[subgroup_small %in% grade_3] <- 3


samp_colors <- colortools::wheel("navyblue", length(unique(sample.tumor$orig.ident)))
names(samp_colors) <- unique(sample.tumor$orig.ident)

colors_small <- list("Metasignature" = ann_colors_small,
                     "Origin" = samp_colors,
                     "Subtype" = subtype_colors)

# Final visualization with the annotation fixed. It should look pretty neat.
h.small <- pheatmap::pheatmap(corr.small,
                              color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                              breaks = seq(-range.small, range.small, length.out = 100),
                              show_colnames = F,
                              treeheight_col = 0,annotation_col = data.frame("Metasignature" = metaclust_small,
                                                          "Origin" = origin_small,
                                                          "Subtype" = subgroup_small),
                              annotation_colors = colors_small,
                              cellwidth = 10,
                              cellheight = 10)


# If the previous visualization is good, we can proceed with the metasignature assignment.
metaclust <- metaclust_small # Assign the good one.

# Compute a metasignature score for each of the metasignatures for each of the cells.
avg.sc.cell <- as.data.frame(BiocGenerics::sapply(as.character(levels(metaclust)), function(c) {rowMeans(sign.scores[,which(metaclust == c), drop = FALSE])}))

# Add the scores as metadata.
for (mc in colnames(avg.sc.cell)){
    scores <- avg.sc.cell[, mc]
    names(scores) <- rownames(avg.sc.cell)
    sample.tumor[[mc]] <- scores
}

# Assign to each cell the maximum score out of the different computed metasignatures.
# PLEASE NOTE: while this might give you a hint towards the different populations in the tumor bulk, the assignment might not be perfect.
# What if a cells is almost equally scored for several metasignatures? Use this as a first look, but not as a definitive result.
# For this, we recommend the use of the top 30 genes per metasignature paired with enrichment scores (shown afterwards).
sample.tumor$metaclust <- BiocGenerics::sapply(apply(avg.sc.cell, 1, which.max), function(x){paste0("Metasignature_", x)})
sample.tumor$metaclust <- as.character(sample.tumor$metaclust)


 # Compute the top 30 scoring genes per each metasignature
 final.sign <- list()
 for (mc in colnames(avg.sc.cell)){
     # c here is the FUNCTION to combine VALUES into a VECTOR.
     ## nmf.sign[which(metaclust == mc) => retrieve the signatures (top 30 genes) from nmf.sign that belong to the meta-signature with the scores per cell for those signatures.
     ## do.call(c, nmf.sign[which(metaclust == mc)])) => Generate a vector with the genes in each signature.
     ## unique(do.call(c, nmf.sign[which(metaclust == mc)])) => Avoid repeting genes.
     ## Therefore, sign is a vector of the unique genes present in all of the signatures that form the meta-signature.
     sign <- unique(do.call(c, nmf.sign[which(metaclust == mc)]))
     # Compute the gene scores for each gene for all the cells that were assigned to a given metacluster.
     gene.score <- scoreGene(tumor_scaled[, which(sample.tumor$metaclust == mc)], X.mean, sign)
     #expr <- rowMeans(Ec[sign, which(tumor.subs$metaclust == x)])
     top.genes <- sign[order(gene.score, decreasing = T)[1:30]]
     # Include the list of top.genes in the final signature list.
     final.sign[[mc]] <- top.genes
 }

# Make sure to save both the metaclust and final.sign object.
saveRDS(metaclust, "") # Modify the path accordingly.
saveRDS(final.sign, "") # Modify the path accordingly.

# With the top 30 scoring genes per MS, one can do multiple downstream analysis to try to infer the properties of the metasignatures.


# ----------------------------------------------------------------------------------------------------
# 9 - PERMUTATION ANALYSIS TO DEFINE ENRICHED CELLS ON GRADIENT-LIKE ENRICHMENT CASES FOR 10X DATASETS
# ----------------------------------------------------------------------------------------------------


# The whole code is designed in a function fashion.

#' Assess enrichment of a Seurat object in a list of genes by permutation testing.
#'
#' The following analysis aims to find a way to statistically select the cells that are truly enriched in a list of marker genes.
#' For this, a permutation analysis approach is followed, in which we compute a null distribution by shuffling the expression values
#' of the genes queried for enrichment. Therefore, this disrupts the equation for the enrichment scores, since they are, in short,
#' the result of the difference in means between the expression values of the genes queried (which we disrupted) and the control
#' selected randomly by the function. This is done until we have one million values, which would act as the null distribution.
#' Once we have a null distribution and an empirical one, we can assess how extreme our values in the empirical distribution are with
#' respect to the values in the null distribution, being the p-value the proportion of values in the null distribution that are higher
#' than the given value from the empirical distribution you are querying. These p-values are, then, adjusted for multiple testing, selecting
#' a FDR cutoff of 0.05/n, being n the total number of lists of marker genes we are going to query to the same cells and use altogether
#' to define a labelling (for instance, if we query two lists for the same tumor bulk, it would be 0.05/2). This is doing to avoid
#' overinflation of the alpha error.
#'
#' @param sample Seurat object to test.
#' @param markers Named list of marker genes. Can contain multiple lists, the important point here is that each list has to be named.
#' @param list.name Name of the list of markers to use.
#' @param outpath Output path of the method.
#' @param FDR_cutoff # FDR cutoff to apply.
#' @param number_comparisons # Number of different lists that are going to be queried to the same cells. FDR value will be divided by this.
#' @return
#' @export
#'
#' @examples
Null_Distribution_approach <- function(sample,
                                       markers,
                                       list.name,
                                       outpath,
                                       FDR_cutoff = 0.05,
                                       number_comparisons = 1){


    message(paste("Running list:", list.name))

    # Retrieve empirical distribution.
    enrichment_name <- paste0(stringr::str_replace_all(list.name, "-", "."))
    scores_name <- paste0(enrichment_name, "1") # Seurat::AddModuleScore adds a "1" to the given name because yes.
    test.dist <- sample.tumor[[]][, scores_name]
    names(test.dist) <- colnames(sample.tumor)

    # Compute null distribution.
    # What happens inside the replicate seems to work on a different environment level.
    set.seed(777) # Reproducibility.

    message("Computing permutations.")
    # We want a null distribution with at least 1.000.000 permutations.
    wanted_permutations <- 1000000
    number_cells <- ncol(sample.tumor)
    nruns <- trunc(wanted_permutations / number_cells) + 1
    null.dist <- pbapply::pbreplicate(nruns, {
        genes.query <- markers[[list.name]][markers[[list.name]] %in% rownames(sample.tumor)]
        # First, create a replacement object.
        sample.null <- sample.tumor
        # Get the normalized data assay (sparse matrix).
        data.use <- sample.null@assays$SCT@data
        row.order <- rownames(data.use) # We want to preserve the original order of the matrix.

        # Subset out the matrix we do not want to reshuffle.
        data.keep <- data.use[-which(rownames(data.use) %in% genes.query), ] # Remove the input genes from the matrix.

        # Get the subset for reshuffle.
        data.to.shuffle <- data.use[which(rownames(data.use) %in% genes.query), ]
        # For each gene in the list of markers.
        df.new <- list()
        for (gene in genes.query){
            # Get the scores.
            expression.scores <- data.to.shuffle[gene, ]
            # Permute the scores for all cells for that given gene.
            shuffled.scores <- sample(expression.scores, length(expression.scores))
            # As the cell names get shuffled as well, we have to change them back to the original order.
            names(shuffled.scores) <- names(expression.scores)
            df.new[[gene]] <- shuffled.scores
        }
        data.to.shuffle <- Matrix::as.matrix(t(sapply(df.new, unlist)))
        data.use <- Matrix::rBind(data.keep, data.to.shuffle) # Perform a rowbind of the two matrices (instead of directly modifying the first one, which takes forever to do)
        data.use <- data.use[row.order, ] # Reorder back the matrix.
        # Set the new matrix as the Assay data from which the enrichment scores will be computed on.
        sample.null <- Seurat::SetAssayData(sample.null, assay = "SCT", slot = "data", new.data = data.use)
        # Compute enrichment scores, which will be the null distribution of the iteration.
        sample.null <- AddModuleScore(sample.null, list(markers[[list.name]]), name = enrichment_name)
        # Retrieve teh null distribution.
        null.dist <- sample.null[[]][, scores_name]
        # Return it.
        return(null.dist)
    })

    # Assign the cell names back to the output matrix.
    rownames(null.dist) <- colnames(sample.tumor)

    # Prepare the data for plotting.
    data.plot <- as.data.frame(null.dist) %>% # Make the matrix into tidiverse accepted object.
        dplyr::mutate(Empirical = test.dist) %>% # Add the empirical distribution to the matrix of AddModuleScore() iterations.
        pivot_longer(dplyr::everything(), names_to = "Distribution", values_to = "Enrichment") %>%
        dplyr::mutate(Distribution = ifelse(Distribution == "Empirical", "Empirical", "Null")) # Assign any name in Distribution that is not "Empirical" into "Null". By default, non-labelled df columns are named V1, V2...

    # Visualize the density plot of both distributions.
    p.dist <- data.plot %>% # Transform from wide to long format.
        ggplot2::ggplot() +
        ggplot2::geom_density(mapping = aes(x = Enrichment, color = Distribution)) +
        ggplot2::scale_color_manual(values = colortools::opposite("steelblue")) +
        ggplot2::ggtitle(paste0(list.name)) +
        ggpubr::theme_pubr(legend = "right")


    # Generate the p-values for each enrichment score in the empirical distribution.
    num_permutations <- sum(data.plot$Distribution == "Null")
    null_dist_values <- data.plot %>% # Gather the null distribution.
        dplyr::filter(Distribution == "Null") %>% # Filter only the values for the NULL.
        dplyr::select(Enrichment)
    null_dist_values <- null_dist_values$Enrichment

    message("Computing p-values.")
    p.value.vector <- c() # Will store all the p-values and will become a column of dist.data.
    # This might also take 10 minutes, since the vector is of 1 million data points to be suuuuuuper exact on the null distribution side.
    p.value.vector <- unlist(pbapply::pblapply(test.dist, function(x){
        greater_values <- sum(null_dist_values > x)
        p.value <- (greater_values + 1) / (num_permutations + 1) #https://pubmed.ncbi.nlm.nih.gov/21044043/
        names(p.value) <- names(x) # Assign the cell name to the p-value.
        p.value.vector <- c(p.value.vector, p.value) # Add the p-value to the output vector.
    }))

    # Generate a reporting matrix for the given permutation test.
    dist.data <- tidyr::tibble(Cell = names(test.dist),
                               Empirical = test.dist,
                               p.value = p.value.vector)

    # FDR correction.
    FDR <- FDR_cutoff / number_comparisons
    message(paste0("Using the FDR cutoff of: ", FDR))
    fdr <- FDR  # FDR to use. Divided by 4 which is the total number of lists of markers that we are gonna test and compare to the same subset of cells (tumor bulk).
    dist.data <- dist.data %>%
                    dplyr::arrange(p.value) %>% # Order by ascending p-value.
                    dplyr::mutate(q.value = stats::p.adjust(p.value, method = "BH")) %>% # Adjust for multiple testing and produce q-values.
                    dplyr::mutate(significant = ifelse(q.value < fdr, TRUE, FALSE)) %>%  # Assign significance.
                    dplyr::mutate(significant_corrected = ifelse(q.value < fdr & Empirical > 0, TRUE, FALSE)) # Assess the outliers with negative enrichment scores that surpass the cutoffs.

    # Check if we had weird cases of significant cells with negative enrichments.
    message(paste("Number of significant results:",
                  sum(dist.data$significant),
                  "\nNumber of significant results with enrichment scores higher than 0:",
                  sum(dist.data$significant_corrected),
                  "\nNumber of outliers:", sum(dist.data$significant) - sum(dist.data$significant_corrected)))

    # Visualizations.
    # UMAP coloring the cells that surpassed the FDR correction.
    surpassing_cells <- dist.data %>%
                            dplyr::filter(significant_corrected == TRUE) %>%
                            dplyr::pull(Cell)
    p.umap <- Seurat::DimPlot(sample.tumor,
                              reduction = "umap",
                              cells.highlight = surpassing_cells,
                              pt.size = 0.5) +
        ggpubr::theme_pubr() +
        ggpubr::rremove("legend.title") +
        ggplot2::ggtitle(paste("Cells that surpassed FDR correction",
                               sum(dist.data$q.value < fdr),
                               "\nFDR applied: ",
                               fdr,
                               "\nFalse positives according to FDR:",
                               trunc(sum(dist.data$q.value < fdr) * fdr),
                               "\nEnrichment score cutoff:",
                               round(min(dist.data$Empirical[sum(dist.data$q.value < fdr)]), 2))) +
        ggpubr::rremove("legend") +
        Seurat::NoAxes() +
        ggplot2::theme(title = ggplot2::element_text(face = "bold", hjust = 0.5))
    # Feature plot of the enrichment scores.
    p.feature <- Seurat::FeaturePlot(sample.tumor, reduction = "umap", scores_name, order = T, pt.size = 0.5, combine = T, ncol = 1) + ggtitle(paste0("Enrichment scores for list ", list.name))  & Seurat::NoAxes() & viridis::scale_color_viridis() &
        ggplot2::theme(title = ggplot2::element_text(face = "bold", hjust = 0.5),
                       legend.text = ggplot2::element_text(size = 10, face = "bold"))
    # Join both plots.
    p.final <- p.umap | p.feature

    # Store the plots.
    ggplot2::ggsave(filename = sprintf("%s_umap_selected_cells_and_enrichment_scores.png", list.name),
                    plot = p.final,
                    path = outpath,
                    device = "png",
                    dpi = 300,
                    width = 14,
                    height = 7)

    ggplot2::ggsave(filename = sprintf("%s_umap_selected_cells_and_enrichment_scores.pdf", list.name),
                    plot = p.final,
                    path = outpath,
                    device = "pdf",
                    dpi = 300,
                    width = 14,
                    height = 7)

    ggplot2::ggsave(filename = sprintf("%s_umap_selected_cells_and_enrichment_scores.svg", list.name),
                    plot = p.final,
                    path = outpath,
                    device = "svg",
                    dpi = 300,
                    width = 14,
                    height = 7)

    # Prepare the output of the function.
    output <- list(empirical_data = dist.data,
                   null_data = data.plot,
                   surpassed_cells = surpassing_cells,
                   p.dist = p.dist)
    # Also, save the output for future use.
    saveRDS(output, paste0(permutation_test_folder, "/results_", list.name), compress = F)
    return(output)
}

# So, to get the output:
output.enrichment <- Null_Distribution_approach(sample = sample.tumor,
                                         		list.name = "your_name",
                                         		markers = final.sign,
                                         		outpath = permutation_test_folder)



# ------------------------------------------------------------------
# 10 - PROGENY ANALYSIS
# ------------------------------------------------------------------

# Method adapted from: https://github.com/saezlab/progeny/blob/master/vignettes/ProgenySingleCell.Rmd

outfolder_progeny <- "" # Path to the folder where the output will be stored.

output_heatmap <- sprintf("%s/heatmaps", outfolder_progeny)
dir.create(path = output_heatmap,
           showWarnings = FALSE,
           recursive = TRUE)

assay.use <- "SCT" # If normalized with SCTransform. Use "RNA" if not.
ref_use <- "Human"
model_human_full <<- progeny::model_human_full # Required object by progeny.

# If the sample is too big, we need to split it in chunks.
number_cells <- ncol(sample)
if (number_cells > 50000){
    if (sample_name != "OE0145-IDH_integrated_LGGs"){
        midpoint <- number_cells / 2
        sample.1 <- sample[, colnames(sample)[1:midpoint]]
        sample.2 <- sample[, colnames(sample)[midpoint + 1 : number_cells]]
        sample.1 <- progeny::progeny(sample.1,
                                     scale = FALSE,
                                     organism = ref_use,
                                     top = 500,
                                     perm = 1,
                                     assay = assay.use,
                                     return_assay = TRUE)
        sample.2 <- progeny::progeny(sample.2,
                                     scale = FALSE,
                                     organism = ref_use,
                                     top = 500,
                                     perm = 1,
                                     assay = assay.use,
                                     return_assay = TRUE)
        sample <- merge(sample.1, sample.2)
    } else {
        number_cells <- length(colnames(sample))
        sample.1 <- sample[, colnames(sample)[1:50000]]
        sample.2 <- sample[, colnames(sample)[50001 : 100000]]
        sample.3 <- sample[, colnames(sample)[100001:number_cells]]

        sample.1 <- progeny::progeny(sample.1,
                                     scale = FALSE,
                                     organism = ref_use,
                                     top = 500,
                                     perm = 1,
                                     assay = assay.use,
                                     return_assay = TRUE)
        sample.2 <- progeny::progeny(sample.2,
                                     scale = FALSE,
                                     organism = ref_use,
                                     top = 500,
                                     perm = 1,
                                     assay = assay.use,
                                     return_assay = TRUE)
        sample.3 <- progeny::progeny(sample.3,
                                     scale = FALSE,
                                     organism = ref_use,
                                     top = 500,
                                     perm = 1,
                                     assay = assay.use,
                                     return_assay = TRUE)
        sample <- merge(sample.1, sample.2)
        sample <- merge(sample, sample.3)
    }
} else {
    sample <- progeny::progeny(sample,
                               scale = FALSE,
                               organism = ref_use,
                               top = 500,
                               perm = 1,
                               assay = assay.use,
                               return_assay = TRUE)
}


Seurat::DefaultAssay(sample) <- "progeny"
# Center the scores.
sample <- Seurat::ScaleData(sample, assay = "progeny")

# Obtain the scores from the Seurat object.
progeny_scores <- as.data.frame(t(Seurat::GetAssayData(sample,
                                                       slot = "scale.data",
                                                       assay = "progeny"))) %>%
                  tibble::rownames_to_column("Cell") %>%
                  tidyr::gather(Pathway, Activity, -Cell)

saveRDS(progeny_scores, sprintf("%s/%s_progeny_scores", output_heatmap, sample_name))
# Generate a mapping of the cells and their cluster labels.
cluster_mapping <- data.frame(Cell = names(Seurat::Idents(sample)),
                              CellType = as.character(Seurat::Idents(sample)),
                              stringsAsFactors = FALSE)
saveRDS(cluster_mapping, sprintf("%s/%s_cluster_mapping", output_heatmap, sample_name))
# Match the scores to the cluster names.
progeny_scores <- dplyr::inner_join(progeny_scores, cluster_mapping)
saveRDS(progeny_scores, sprintf("%s/%s_progeny_scores_inner_join", output_heatmap, sample_name))
# Summarize the scores per cell population.
summarized_progeny_scores <- progeny_scores %>%
    dplyr::group_by(Pathway, CellType) %>%
    dplyr::summarise(avg = mean(Activity), std = stats::sd(Activity))
saveRDS(summarized_progeny_scores, sprintf("%s/%s_summarized_progeny_scores", output_heatmap, sample_name))

summarized_progeny_scores <- summarized_progeny_scores %>%
                             dplyr::select(-std) %>%
                             tidyr::spread(Pathway, avg) %>%
                             data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)


# Generate a heatmap visualization.
data.use <- summarized_progeny_scores[, unique(sample.tumor$New_NMF_labelling)]
range <- max(abs(data.use))

# Visualize it.
h <- pheatmap::pheatmap(as.matrix(data.use),
                        color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                        breaks = seq(-range, range, length.out = 100),
                        cellwidth = 20,
                        cellheight = 20)

# ------------------------------------------------------------------
# 11 - DOROTHEA ANALYSIS
# ------------------------------------------------------------------

# Generate output structure.
outfolder_dorothea <- "" # Will store the output.
outfolder_dorothea.individual <- sprintf("%s/dorothea/individual_plots", outfolder)
outfolder_dorothea.summary <- sprintf("%s/dorothea/summary", outfolder)
outfolder_dorothea.dimensional_reduction <- sprintf("%s/dorothea/dimensional_reduction", outfolder)
outfolder_dorothea.markers <- sprintf("%s/dorothea/tf_markers", outfolder)
outfolder_dorothea.feature_plots <- sprintf("%s/dorothea/feature_plots", outfolder)

dir.create(path = outfolder_dorothea,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(path = outfolder_dorothea.individual,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(path = outfolder_dorothea.summary,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(path = outfolder_dorothea.dimensional_reduction,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(path = outfolder_dorothea.markers,
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(path = outfolder_dorothea.feature_plots,
           showWarnings = FALSE,
           recursive = TRUE)

# Retrieve dorothea regulons.
dorothea_regulon <- get(utils::data("dorothea_hs", package = "dorothea"))
assay.use <- "SCT" # Select "RNA" if not normalized with SCTransform.

# Modify run_viper.Seurat() function to allow the use of SCT assay.
# This should be already fixed in recent versions.
run_viper <- function(input,
                      regulons,
                      options = list(),
                      tidy = FALSE,
                      assay.use = assay.use) {
    UseMethod("run_viper")
}
run_viper.Seurat <- function(input,
                             regulons,
                             options = list(),
                             tidy = FALSE,
                             assay.use = "RNA") {
    if (tidy) {
        tidy <- FALSE
        warning("The argument 'tidy' cannot be TRUE for 'Seurat' objects. ",
                "'tidy' is set to FALSE")
    }
    print(sprintf("Using %s assay.", assay.use))
    mat <- as.matrix(Seurat::GetAssayData(input, assay = assay.use, slot = "data"))

    tf_activities <- dorothea::run_viper(mat, regulons = regulons, options = options,
                                         tidy = FALSE)

    # include TF activities in Seurat object
    dorothea_assay <- Seurat::CreateAssayObject(data = tf_activities)
    Seurat::Key(dorothea_assay) <- "dorothea_"
    input[["dorothea"]] <- dorothea_assay

    return(input)
}

# Obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon %>%
    dplyr::filter(confidence %in% c("A", "B", "C"))


# Same scenario in which the sample is too big to be analyzed as a whole (> 50000 cells).
if (ncol(sample) > 50000){
    if (sample_name != "OE0145-IDH_integrated_LGGs"){
        number_cells <- length(colnames(sample))
        midpoint <- number_cells / 2
        sample.1 <- sample[, colnames(sample)[1:midpoint]]
        sample.2 <- sample[, colnames(sample)[midpoint + 1 : number_cells]]
        sample.1 <- run_viper(sample.1,
                              regulon,
                              assay.use = "SCT",
                              options = list(method = "scale",
                                             minsize = 4,
                                             eset.filter = FALSE,
                                             cores = 1,
                                             verbose = FALSE))
        sample.2 <- run_viper(sample.2,
                              regulon,
                              assay.use = "SCT",
                              options = list(method = "scale",
                                             minsize = 4,
                                             eset.filter = FALSE,
                                             cores = 1,
                                             verbose = FALSE))
        sample <- merge(sample.1, sample.2)
    } else {
        number_cells <- length(colnames(sample))
        sample.1 <- sample[, colnames(sample)[1:50000]]
        sample.2 <- sample[, colnames(sample)[50001 : 100000]]
        sample.3 <- sample[, colnames(sample)[100001:number_cells]]
        sample.1 <- run_viper(sample.1,
                              regulon,
                              assay.use = "SCT",
                              options = list(method = "scale",
                                             minsize = 4,
                                             eset.filter = FALSE,
                                             cores = 1,
                                             verbose = FALSE))
        sample.2 <- run_viper(sample.2,
                              regulon,
                              assay.use = "SCT",
                              options = list(method = "scale",
                                             minsize = 4,
                                             eset.filter = FALSE,
                                             cores = 1,
                                             verbose = FALSE))
        sample.3 <- run_viper(sample.3,
                              regulon,
                              assay.use = "SCT",
                              options = list(method = "scale",
                                             minsize = 4,
                                             eset.filter = FALSE,
                                             cores = 1,
                                             verbose = FALSE))
        sample <- merge(sample.1, sample.2)
        sample <- merge(sample, sample.3)
    }
} else {
    sample <- run_viper(sample,
                        regulon,
                        assay.use = "SCT",
                        options = list(method = "scale",
                                       minsize = 4,
                                       eset.filter = FALSE,
                                       cores = 1,
                                       verbose = FALSE))
}

Seurat::DefaultAssay(sample) <- "dorothea"
sample <- Seurat::ScaleData(sample, assay = "dorothea")

# Find differential enriched regulons.
dorothea.markers <- Seurat::FindAllMarkers(sample,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25,
                                           verbose = FALSE,
                                           base = exp(1)) # Seurat v3, if v4 use "2".
saveRDS(dorothea.markers, sprintf("%s/%s_dorothea_markers", outfolder_dorothea.individual, sample_name))

# Retrieve viper scores.
viper_scores_df <- Seurat::GetAssayData(sample,
                                        slot = "scale.data",
                                        assay = "dorothea") %>%
        		   data.frame() %>%
        		   t()

# Bugfix: change . into - in the cell names.
rownames(viper_scores_df) <- gsub("[.]", "-", rownames(viper_scores_df))
saveRDS(viper_scores_df, sprintf("%s/%s_viper_scores_df", outfolder_dorothea.individual, sample_name))

# Obtain a df with the mapping cell-identity.
cell_clusters <- data.frame(cell = names(Seurat::Idents(sample)),
                            cell_type = as.character(Seurat::Idents(sample)),
                            stringsAsFactors = FALSE)
saveRDS(cell_clusters, sprintf("%s/%s_cell_clusters", outfolder_dorothea.individual, sample_name))

# Generate a df with the viper score per cell and its clusters.
viper_scores_clusters <- viper_scores_df  %>%
        				 data.frame() %>%
        				 tibble::rownames_to_column("cell") %>%
        				 tidyr::gather(tf, activity, -cell) %>%
        				 dplyr::inner_join(cell_clusters)
saveRDS(viper_scores_clusters, sprintf("%s/%s_viper_scores_clusters", outfolder_dorothea.individual, sample_name))

# Summarise by cell population.
summarized_viper_scores <- viper_scores_clusters %>%
        				   dplyr::group_by(tf, cell_type) %>%
        				   dplyr::summarise(avg = mean(activity),
        				                    std = stats::sd(activity))
saveRDS(summarized_viper_scores, sprintf("%s/%s_summarized_viper_scores", outfolder_dorothea.individual, sample_name))

# Run for 20, 30 and 40 top regulons.
for (num_tf in c(20, 30, 40)){
    highly_variable_tfs <- summarized_viper_scores %>%
        dplyr::group_by(tf) %>%
        dplyr::mutate(var = stats::var(avg))  %>%
        dplyr::ungroup() %>%
        dplyr::top_n((num_tf * length(levels(sample))), var) %>%
        dplyr::distinct(tf)

    # Prepare the data for the heatmap plot.
    summarized_viper_scores_df <- summarized_viper_scores %>%
        dplyr::semi_join(highly_variable_tfs, by = "tf") %>%
        dplyr::select(-std) %>%
        tidyr::spread(tf, avg) %>%
        data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

    data.use <- summarized_viper_scores_df
    range <- max(abs(data.use))

    h <- pheatmap::pheatmap(as.matrix(data.use),
                            color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
                            breaks = seq(-range, range, length.out = 100),
                            cellwidth = 20,
                            cellheight = 20)
    saveRDS(summarized_viper_scores_df, sprintf("%s/%s_heatmap_using_top%s_TF", outfolder_dorothea.individual, sample_name, num_tf))

    ggplot2::ggsave(filename = sprintf("%s_Dorothea_scores_per_cluster_heatmap_using_top%s_TF.png", sample_name, num_tf),
                    plot = h,
                    path = outfolder_dorothea.individual,
                    device = "png",
                    dpi = 300,
                    height = 10,
                    width = 14)

    ggplot2::ggsave(filename = sprintf("%s_Dorothea_scores_per_cluster_heatmap_using_top%s_TF.svg", sample_name, num_tf),
                    plot = h,
                    path = outfolder_dorothea.individual,
                    device = "svg",
                    dpi = 300,
                    height = 10,
                    width = 14)
}


# ------------------------------------------------------------------
# 12 - DIFFUSION MAPS
# ------------------------------------------------------------------

sample.subset <- sample # Or the subset of the sample you want to compute it over. (i.e: Tumor Bulk)
sample.use <- Seurat::as.SingleCellExperiment(sample.subset)

diffusion_folder <- ""  # Output folder.
dir.create(diffusion_folder, recursive = T)

# Debugging version:
diff.map <- try({destiny::DiffusionMap(sample.use, verbose = T)})

if (class(diff.map) == "try-error"){
    message("The previous run failed for any given reason.")
    message("Retrying with k = nrow(sample) - 1")
    diff.map <- destiny::DiffusionMap(sample.use, verbose = T, k = nrow(sample.use) - 1)
}
saveRDS(diff.map, paste0(diffusion_folder, "/", sample_name, "_difussion_map_object"))

# Add the diffusion map as a reduction in the Seurat object.
sample.subset@reductions$diffusion <- Seurat::CreateDimReducObject(embeddings = diff.map.oligo@eigenvectors, key = "DC", assay = "SCT")

# Select a diffusion component. There is an error in which, for any reason, a low amount of cells are not present in the diffusion map object.
# Therefore, we will subset the sample to only the cells present in the diffusion map object (this subset will only be used here).
# The number of cells not present is really low, maybe 10 out of 65000.
dc_choice <- 1
dc_use <- paste0("DC", dc_choice)
selected<- names(diff.map@eigenvectors[, dc_use])
sample.subset <- sample.subset[, selected]
# The sample is ready for visualization.


