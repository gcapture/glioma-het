#!/bin/bash

#
# CellPhoneDB
#


# Running statistical analysis

## Define parameters
#metadata_path="/scratch/devel/inmah/GLIOBLASTOMA/submission/CellPhoneDB/metadata_oligo_microgliaAll.txt"
#counts_path="/scratch/devel/inmah/GLIOBLASTOMA/submission/CellPhoneDB/counts_oligo_microgliaAll.txt"

## Run CellPhoneDB v.2.0 in statistical analysis mode

cellphonedb method statistical_analysis --pvalue=0.001 /scratch/devel/inmah/GLIOBLASTOMA/submission/CellPhoneDB/metadata_oligo_microgliaAll.txt  /scratch/devel/inmah/GLIOBLASTOMA/submission/CellPhoneDB/counts_oligo_microgliaAll.txt --counts-data=gene_name --iterations=1000 --threads=8 # because we use gene symbol 


#cellphonedb method statistical_analysis \
#    $metadata_path \
#    $counts_path \
#    --counts-data=gene_name \
#    --iterations=1000 \
#    --threads=8

## Run fixOrder.py script
python /scratch/groups/hheyn/software/anaconda3/envs/cellphonedb_ljimenez/lib/python3.7/site-packages/cellphonedb/utils/fixOrder.py out/means.txt > out/means_fixorder.txt

