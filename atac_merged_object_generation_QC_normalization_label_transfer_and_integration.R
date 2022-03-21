# Enrique Blanco Carmona
# e.blancocarmona@kitz-heidelberg.de
# PhD Student – Clinical Bioinformatics
# Division of Pediatric Neurooncology (B062)


# Code heavily based on Signac's vignette: https://satijalab.org/signac/articles/merging.html


# Generate the following objects.

# Vector containig all the sample IDs.
atac_samples <- c("sample_a",
                  "sample_b",
                  "sample_c")


# List containing as names the IDs in atac_samples and as values the full path to the peaks file.
peak.list <- list("sample_a" = "full_path_to_peaks_file",
                  "sample_b" = "full_path_to_peaks_file",
                  "sample_c" = "full_path_to_peaks_file")

# Generate a list with the peak objects.
peak.out <- list()
for (s.name in atac_samples){
    db <- generate_database(sprintf("/icgc/dkfzlsdf/analysis/OE0145_projects/idh_gliomas/scripts/main/06_scana/config_files/%s_config.yml", s.name))
    peaks <- read.table(file = peak.list[[s.name]],
                       col.names = c("chr", "start", "end"))
    peak.out[[s.name]] <- peaks
}

# Convert peaks into GRanges (Adapt for the number of samples you have).
gr.1 <- GenomicRanges::makeGRangesFromDataFrame(peak.out$`sample_a`)
gr.2 <- GenomicRanges::makeGRangesFromDataFrame(peak.out$`sample_b`)
gr.3 <- GenomicRanges::makeGRangesFromDataFrame(peak.out$`sample_c`)

# Put them back to a vector.
reduce.vector <- c(gr.1, gr.2, gr.3, gr.4, gr.5, gr.6, gr.7, gr.8, gr.9, gr.10, gr.11)

# Create an unified set of peaks to quantify in each dataset.
combined.peaks <- GenomicRanges::reduce(x = reduce.vector)


# Generate a list for the fragment files location.
fragment.paths <- list("sample_a" = "full_path_to_fragments_file",
                       "sample_b" = "full_path_to_fragments_file",
                       "sample_c" = "full_path_to_fragments_file")

# Generate a list for the metadata files location.
metadata.paths <- list("sample_a" = "full_path_to_metadata_file",
                       "sample_b" = "full_path_to_metadata_file",
                       "sample_c" = "full_path_to_metadata_file")

# Generate a list to store the fragment objects.
fragment.list <- list() # Will contain the framgent objects.
md.list <- list() # Will contain the metadata objects.
for(s.name in atac_samples){
    md <- read.table(file = metadata.paths[[s.name]],
                     stringsAsFactors = FALSE,
                     sep = ",",
                     header = TRUE,
                     row.names = 1,)[-1, ]
    # Filter low count cells.
    md <- md[md$passed_filters > 500, ]
    md.list[[s.name]] <- md

    # Create fragment object.
    frag.md <- Signac::CreateFragmentObject(path = fragment.paths[[s.name]],
                                            cells = rownames(md))
    fragment.list[[s.name]] <- frag.md
}


# Generate count matrices and store them in a list.
count.list <- list()
for(s.name in atac_samples){
    ct <- Signac::FeatureMatrix(fragments = fragment.list[[s.name]],
                                features = combined.peaks,
                                cells = rownames(md.list[[s.name]]))
    count.list[[s.name]] <- ct
}


# Generate a list that contains the signac objects.
sample.list <- list()

# Generate an annotation object to include to the Signac object. 
GenomeInfoDb::seqlevelsStyle(EnsDb.Hsapiens.v86) <- "UCSC"
annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
Signac::genome(annotation) <- "hg38"

for(s.name in samples.use){
    s.assay <- Signac::CreateChromatinAssay(count.list[[s.name]],
                                            fragments = fragment.list[[s.name]],
                                            genome = "hg38",
                                            annotation = annotation)
    s <- Seurat::CreateSeuratObject(s.assay,
                                    assay = "ATAC",
                                    meta.data = md.list[[s.name]],
                                    project = s.name)

    sample.list[[s.name]] <- s
}

# Generate a merged sample.
# To merge, you need to select a "reference" object upon which you merge the others. 
ref_id <- "sample_a"
sample <- merge(x = sample.list[[ref_id]],
                y = sample.list[!(names(sample.list) %in% ref_id)],
                add.cell.ids = c(ref_id, names(sample.list)[!(names(sample.list) %in% ref_id)]))
# Reassign the annotation just in case. 
Signac::Annotation(sample) <- annotation

# Define QC cutoffs.
reference_organism <- "h.sapiens"
# QC cutoff: lower peak region fragments.
lower_peak_region_fragments_cutoff = 3000
# QC cutoff: higher peak region fragments.
higher_peak_region_fragments_cutoff = 20000
# QC cutoff: percentage of reads in peaks (higher than).
percentage_reads_in_peaks_cutoff = 15
# QC cutoff: blacklist ratio (lower than).
blacklist_ratio_cutoff = 0.05
# QC cutoff: nucleosome signal cutoff (lower than).
nucleosome_signal_cutoff = 10
# QC cutoff: TSS enrichment score cutoff (higher than).
TSS_enrichment_score_cutoff = 2



# Calculate the chromosome binding pattern. Stored as "nucleosome_signal" column. # Computes it for chr1.
sample <- Signac::NucleosomeSignal(object = sample)

# Compute TSS enrichment.
sample <- Signac::TSSEnrichment(sample, fast = FALSE)

# Compute the percentage of reads in peaks.
sample$percentage_reads_in_peaks <- sample$peak_region_fragments / sample$passed_filters * 100

# Compute ratio of reads in blacklist sites.
sample$blacklist_ratio <- sample$blacklist_region_fragments / sample$peak_region_fragments

# Classify in different nucleosome groups.
sample$nucleosome_group <- ifelse(sample$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# Assign cells to high or low TSS enrichment.
sample$high_tss <- ifelse(sample$TSS.enrichment > 2, "High", "Low")


# Compute the QC logical vectors.
peak_lower_cutoff <- sample$peak_region_fragments > lower_peak_region_fragments_cutoff
peak_higher_cutoff <- sample$peak_region_fragments < higher_peak_region_fragments_cutoff
percentage_reads_in_peaks_cutoff <- sample$percentage_reads_in_peaks > percentage_reads_in_peaks_cutoff
blacklist_ratio_cutoff <- sample$blacklist_ratio < blacklist_ratio_cutoff
nucleosome_signal_cutoff <- sample$nucleosome_signal < nucleosome_signal_cutoff
tss_cutoff <- sample$TSS.enrichment > TSS_enrichment_score_cutoff

# Combine them.
mask <- peak_lower_cutoff & peak_higher_cutoff & percentage_reads_in_peaks_cutoff & blacklist_ratio_cutoff & nucleosome_signal_cutoff & tss_cutoff

# Subset the sample. 
sample <- sample[, which(mask)]

# Normalize.
sample <- Signac::RunTFIDF(sample)

# Find top peaks (in our case, we use q0 to select all of them).
sample <- Signac::FindTopFeatures(sample, min.cutoff = "q0")

# Run SVD.
sample <- Signac::RunSVD(object = sample,
                         assay = "ATAC",
                         reduction.key = "LSI_",
                         reduction.name = "lsi")

# Check the LSI components.
Signac::DepthCor(sample) # If the component 1 is at -1, exclude it.
first_comp <- 2 # Set to 1 if component 1 is not at -1.


# Compute UMAP.
sample <- Seurat::RunUMAP(object = sample,
                          reduction = "lsi",
                          dims = first_comp:30)

# Find neighbors.
sample <- Seurat::FindNeighbors(object = sample,
                                reduction = "lsi",
                                dims = first_comp:30)

# Find clusters.
sample <- Seurat::FindClusters(object = sample,
                               verbose = FALSE,
                               algorithm = 1)

# Compute activity matrix.
gene.activities <- Signac::GeneActivity(object = sample)

# Add the activity matrix to the Seurat object as a new assay, and normalize it.
sample[["RNA"]] <- Seurat::CreateAssayObject(counts = gene.activities)

# Normalize new assay data - treating the Activity assay as if it was a scRNAseq experiment.
sample <- Seurat::SCTransform(object = sample,
                              assay = "RNA",
                              new.assay.name = "SCT")



# Transfer labels from true scRNAseq object to the activity assay from scATACseq object. 
# This assumes you have a corresponding scRNAseq experiment.

# Load snRNAseq dataset.
sample.rna <- readRDS("full_path_to_scRNAseq_object")

# Rename gene activity assay to ACTIVITY.
sample <- SeuratObject::RenameAssays(sample, SCT = "ACTIVITY")

# Find anchors between datasets.
transfer.anchors <- Seurat::FindTransferAnchors(reference = sample.rna, # Seurat object to use as reference.
                                                query = sample,    # Seurat object to use as querey.
                                                features = Seurat::VariableFeatures(object = sample.rna), # Features to use for Dimensional Reduction.
                                                reference.assay = "SCT", # Assay to use as reference.
                                                query.assay = "ACTIVITY", # Assay to use as query.
                                                reduction = "cca") # Dimensional reduction to use when finding anchors.

# Transfer cluster labels from snRNAseq dataset.
celltype.predictions <- Seurat::TransferData(anchorset = transfer.anchors,
                                             refdata = sample.rna$cluster_names, # Data to transfer (cluster names).
                                             weight.reduction = sample[["lsi"]],
                                             dims = 1:15) # Dimensional reduction to use.

# Add the predictions as metadata.
sample <- Seurat::AddMetaData(sample, metadata = celltype.predictions) # Predictions stored in "predicted.id" metadata column.

# Generate an integrated object removing sample-specific effects..
Seurat::DefaultAssay(sample) <- "ATAC"
Seurat::Idents(sample) <- sample$predicted.id # Set the new identities.

# Run Harmony. 
sample.integrated <- harmony::RunHarmony(object = sample,
                                        group.by.vars = 'orig.ident',
                                        reduction = 'lsi',
                                        assay.use = 'ATAC',
                                        project.dim = FALSE)
# Recompute UMAP.
sample.integrated <- Seurat::RunUMAP(sample.integrated, reduction = "harmony", dims = first_comp:30)

