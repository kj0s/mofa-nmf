# This file is mostly just to convert RDS file to something else that is more
# python friendly like anndata.

setwd("/vast/projects/Sisseq/givanna")

# From a bit of digging, it seems the following RDS files are the most "raw"
# as one can get in this messy code structure:
# 1. ST223.annotated_2024.rds
# 2. ST223.annotated.rds
# 3. GEX_ADT_TW040723.rds
# 
# For 1 and 3, no idea where they came from, but src/main.R script, which is
# the preparation script that was used for fused nmf, made use of these rds files.
# File 2, seem to be similar to 1. 
# 
# Let's start by checking what is the difference between 1 and 2.

library(Seurat)
library(dplyr)
library(tibble)
library(Matrix)

data_dir <- "/vast/projects/Sisseq/givanna/human-haematopoiesis-sis-seq/data"

# ==== Question 1 ====
# What is the difference between ST223.annotated.rds and ST223.annotated_2024.rds?
# ====================

st223_annotated <- readRDS(file.path(data_dir, "ST223.annotated.rds"))
st223_annotated_2024 <- readRDS(file.path(data_dir, "ST223.annotated_2024.rds"))

# both have the same number of cells i think?
sprintf("N cells in annotated: %d", ncol(st223_annotated))
sprintf("N cells in annotated 2024: %d", ncol(st223_annotated_2024))

# are the cell barcode the same in both?
cell_bcode_diff <- union(
    setdiff(rownames(st223_annotated_2024@meta.data), rownames(st223_annotated@meta.data)),
    setdiff(rownames(st223_annotated@meta.data), rownames(st223_annotated_2024@meta.data))
)

sprintf("How many cell barcodes that are in annotated but not in annotated 2024, and vice versa? %d", length(cell_bcode_diff))

# Are the layers the same?
layers_diff <- union(
    setdiff(Layers(st223_annotated), Layers(st223_annotated_2024)),
    setdiff(Layers(st223_annotated_2024), Layers(st223_annotated))
)
sprintf("Do they have the same layers? %s", length(layers_diff) == 0)

# what about the assays?
assays_diff <- union(
    setdiff(Assays(st223_annotated), Assays(st223_annotated_2024)),
    setdiff(Assays(st223_annotated_2024), Assays(st223_annotated))
)
sprintf("Do they have the same assays? %s", length(assays_diff) == 0)

# ok layers and assays are all the same.
# what about the metadata columns? are they the same?
metadata_col_diff <- union(
    setdiff(colnames(st223_annotated@meta.data), colnames(st223_annotated_2024@meta.data)),
    setdiff(colnames(st223_annotated_2024@meta.data), colnames(st223_annotated@meta.data))
)
sprintf("Are the metadata column the same? %s", length(metadata_col_diff) == 0)

# No they are not. So what is the difference?
meta_col_in_annotated_not_in_2024 <- setdiff(
    colnames(st223_annotated@meta.data), 
    colnames(st223_annotated_2024@meta.data)
)
sprintf("Metadata column in annotated but not in 2024: %d", length(meta_col_in_annotated_not_in_2024))

meta_col_in_2024_not_in_annotated <- setdiff(
    colnames(st223_annotated_2024@meta.data),
    colnames(st223_annotated@meta.data)
)
sprintf("Metadata column in 2024 but not in annotated: %d", length(meta_col_in_2024_not_in_annotated))
sprintf("What is the extra metadata column in 2024? %s", paste(meta_col_in_2024_not_in_annotated, sep = ", "))

# ==== Finding 1 ====
# It seems populations column is what makes 2024 useful for main.R script 
# which preprocess data for fused NMF.
# Hence the 2024 is the one to use
# ===================


# ==== Question 2 ====
# What is in GEX_ADT_TW040723.rds?
# ====================

gex_adt <- readRDS(file.path(data_dir, "GEX_ADT_TW040723.rds"))
sprintf(
    "N cells: %d, N columns: %d, N Layers: %d, N Assays: %d in GEX_ADT?",
    ncol(gex_adt), nrow(gex_adt), length(Layers(gex_adt)), length(Assays(gex_adt))
)

cat("Assays in GEX_ADT?", paste(Assays(gex_adt), sep = ", "))
for (assay in Assays(gex_adt)) {
    layer_in_assay <- Layers(gex_adt[[assay]])
    cat("Layers in", assay, ": ", paste(layer_in_assay), "\n")
}

# what assays are in 2024 annotated?
cat("Assays in 2024: ", paste(Assays(st223_annotated_2024), sep = ", "))

# Seems there are SCT, RNA, and ADT layers in GEX.
# All of which has count layer.
# The 2024 version has RNA and ADT as well. Are the content the same?
# Easiest way to find out is to just spot check some cells
# and see if the values are the same.
cat("N cells in GEX:", ncol(gex_adt), "\nN cells in 2024: ", ncol(st223_annotated_2024))

# So GEX has more cells. Are all 2024 cells in GEX?
cat(
    "Are all cells in 2024 also in GEX?",
    length(intersect(colnames(st223_annotated_2024), colnames(gex_adt))) == length(colnames(st223_annotated_2024))
)

# What about the row counts? Are they the same?
# Just do sum of all. They should be the same.
# Will take some time to run
for (assay in c("ADT", "RNA")) {
    # subset gex adt to just cells in 2024, as otherwise the sum won't be the same
    # as gex adt has more cells
    sum_gex_adt <- sum(GetAssayData(subset(gex_adt, cells = colnames(st223_annotated)), layer = "count"))
    sum_st223_2024 <- sum(GetAssayData(st223_annotated_2024, layer = "count"))
    print(paste(
        "Are the sum of Assay", assay, "between GEX and 2024 the same?",
        sum_gex_adt == sum_st223_2024
    ))
}


# ==== Finding 2 =====
# It seems the GEX has an extra few cells,
# But all annotated 2024 cells are all in GEX and
# The two ADT and RNA assays seem to be the same.
# ====================


# ==== Question 3 ====
# How did the clone in the fate readout mapped to cons_bc in these RDS files?
# There are 2 fate readout files: ST223_all_barcodes.<cell or cpm>.norm.txt
# Both version values have been normalised in that the raw readout is first
# divided by total readout in the sample, similar to library size normalisation.
# cell then multiply the propotion by number of cells sorted in a given sample.
# cpm then multiply by just one million.
# One sample here i think is one day one cell type.
# But why in main.R (and some other scripts), the barcode is trimmed to 15bp?
# ====================

# Could this 15bp is just the first 15 bp of cons_bc?
# Is it enough to differentiate the barcodes?
clone_bc <- st223_annotated_2024@meta.data$cons_bc
# remove na
clone_bc <- clone_bc[!is.na(clone_bc)]

clone_bc_15bp <- substr(clone_bc, 0, 15)

cat(
    "Is 15bp truncation yielded the same number of unique clone barcodes compared to no trunction?",
    length(unique(clone_bc)) == length(unique(clone_bc_15bp))
)

# Nope. But what if we concatenated the donor?
clone_bc_15bp_with_donorid <- paste(
    clone_bc_15bp, st223_annotated_2024@meta.data$donor_id,
    sep = "_"
)
cat(
    "Is 15bp truncation with donor id yielded the same number of unique clone barcodes compared to no trunction?",
    length(unique(clone_bc)) == length(unique(clone_bc_15bp_with_donorid))
)

cat(
    "N clone barcodes untouched: ", length(unique(clone_bc)), ".\n",
    "N clone barcodes 15bp: ", length(unique(clone_bc_15bp)), ".\n",
    "N clone barcodes 15bp with donor id: ", length(unique(clone_bc_15bp_with_donorid)), ".\n"
)

# There is nothing matching above.
# So taking just the first 15bp will not yield the same barcode variety as if it is untouched.
# So why the first 15? why not just use the entire 21bp?

# let's load up the clone barcodes
# note, the columns are not just days + cell type.
# there ar also plates, PA, PB, PC.
fate_readout <- read.table(
    file.path(data_dir, "ST223_all_barcodes.cell.norm.txt")
)
# assign clone bc as either the 1st 15bp or 4th to 18th bp
# why? because ST223_barcode_clustering.R use the latter....
fate_readout$bcode_first_15bp <- paste(
    substr(
        sapply(rownames(fate_readout), function(x) strsplit(x, " ")[[1]][[2]]),
        0, 15
    ),
    fate_readout$patient, sep = "_"
)
fate_readout$bcode_4_to_18bp <- paste(
    substr(
        sapply(rownames(fate_readout), function(x) strsplit(x, " ")[[1]][[2]]),
        4, 18
    ),
    fate_readout$patient, sep = "_"
)

# also add a column where we grabbed just the first however many characters
# there are in the clone barcode of thd RDS file untouched.
# i think this is 20bp..
fate_readout$bcode_first_20bp <- paste(
    substr(
        sapply(rownames(fate_readout), function(x) strsplit(x, " ")[[1]][[2]]),
        0, 20
    ),
    fate_readout$patient, sep = "_"
)

# create a data.frame for clone bc from 2024 which has just the first 15bp
# of the barcode and the full
clone_bc_with_donor <- st223_annotated_2024@meta.data[!is.na(st223_annotated_2024@meta.data$cons_bc), c("cons_bc", "donor_id")]
clone_bc_with_donor$bcode_with_donor <- paste(
    clone_bc_with_donor$cons_bc,
    clone_bc_with_donor$donor_id,
    sep = "_"
)
clone_bc_with_donor$bcode_first_15bp <- paste(
    substr(clone_bc_with_donor$cons_bc, 0, 15), clone_bc_with_donor$donor_id,
    sep = "_"
)
clone_bc_with_donor$bcode_4_to_18bp <- paste(
    substr(clone_bc_with_donor$cons_bc, 4, 18), clone_bc_with_donor$donor_id,
    sep = "_"
)
# clone_bc_with_donor gives one row per cell, but cells can share clone barcode,
# and we only care about the clone barcode here. so just remove the duplicated
# clone barcodes.
clone_bc_with_donor <- clone_bc_with_donor[!duplicated(clone_bc_with_donor$bcode_with_donor), ]
# remove redundcant cell id used as rownames
rownames(clone_bc_with_donor) <- NULL


# do some merge
merged_by_15bp <- inner_join(
    fate_readout,
    clone_bc_with_donor,
    by = "bcode_first_15bp"
)
# if we merged just by 15bp, how many barcodes do we retain?
sprintf(
    "N barcode retained if merging by first 15bp: %d out of %d in 2024 RDS file", 
    nrow(merged_by_15bp),
    length(unique(clone_bc_with_donor$bcode_with_donor))
)

# what if we join by 4-18bp?
merged_by_4_18bp <- inner_join(
    fate_readout,
    clone_bc_with_donor,
    by = "bcode_4_to_18bp"
)
# no, the 4-18 merge ended up as many to many join. it won't work.

# join by the full bcode in the RDS file?
merged_by_20bp <- inner_join(
    fate_readout,
    clone_bc_with_donor,
    by = join_by(bcode_first_20bp == bcode_with_donor)
)
sprintf(
    "N barcode retained if merging by first 20bp: %d out of %d in 2024 RDS file", 
    nrow(merged_by_20bp),
    length(unique(clone_bc_with_donor$bcode_with_donor))
)

# ==== Finding 3 ====
# So if we took just the first 15 bp for the barcode but concatenated it with the patient id,
# we retained only 843 out of 3131 barcode in the RDS file.
# Also, the fused nmf, used the clone barcode concatenated by patients when
# doing the claculations..
# But if you join by the whole barcode length in the RDS file, it ended up with
# only 806 clones, that is a lot less than using the first 15bp.
# So let's just use the first 15bp for now..
# ====================

# Let's export files for KJ so it is easier for her.
# First, fate read out with just the first 15bp.
# Simple, just save the data frame as csv.

out_dir <- "/vast/projects/Sisseq/givanna/code"
# modification
out_dir <- "/vast/projects/Sisseq/human-haematopoiesis-sis-seq/mofa_kj/files"
write.csv(fate_readout, file.path(out_dir, "fate_readout.csv"))

# attach the metadata to the assay export
cell_metadata <- st223_annotated_2024@meta.data

# save the raw counts for ADT and RNA in a csv file 
# transpose so features = columns
adt_raw <- as.data.frame(
    as.matrix(
        t(GetAssayData(st223_annotated_2024, assay = "ADT", layer = "count"))
    )
)
adt_raw <- as_tibble(inner_join(
    rownames_to_column(adt_raw, var = "cell"), 
    cell_metadata,
    by="cell"
))
write.csv(adt_raw, file.path(out_dir, "adt_raw.csv"))

# repeat for rna, but store as mtx as it is big..
writeMM(
    t(GetAssayData(st223_annotated_2024, assay = "RNA", layer = "counts")), 
    file.path(out_dir, "X_rna_counts.mtx")
)
writeLines(rownames(st223_annotated_2024), "X_rna_counts_gene_names.txt")
# modification
writeLines(rownames(st223_annotated_2024[["RNA"]]), file.path(out_dir, "X_rna_counts_gene_names.txt"))
writeLines(colnames(st223_annotated_2024), "X_rna_counts_cell_ids.txt")
# modification
writeLines(colnames(st223_annotated_2024), file.path(out_dir, "X_rna_counts_cell_ids.txt"))
