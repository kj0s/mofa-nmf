library(Seurat)

# For each file:
obj <- readRDS("/vast/projects/Sisseq/human-haematopoiesis-sis-seq/data/ST223.annotated_2024.rds")

# Is it even a Seurat object, or a matrix/df?
class(obj)

# Row count vs unique barcode count — this is the key diagnostic.
# If ncol(obj) >> number of unique clonal barcodes, it's single-cell.
# If ncol(obj) == number of unique clonal barcodes, it's pseudobulked/averaged.
dim(obj)
head(colnames(obj))          # cell barcodes (10x-style, e.g. AAACCTGAGCTAGTGG-1)
head(obj@meta.data)
colnames(obj@meta.data)      # look for a clonal-barcode column, cluster/fate columns

# Check for a lineage/clonal barcode metadata column specifically
grep("barcode|clone|lineage", colnames(obj@meta.data), ignore.case = TRUE, value = TRUE)

# Count cells per unique clonal barcode — if this shows real distributions
# (e.g. 1-50 cells per barcode) rather than all 1s, you've got single-cell data
# with clone info attached
table(obj@meta.data$clonal_barcode_column_name) |> table()  # distribution of clone sizes

obj <- readRDS(".../ST223.annotated_fate.rds")

# Is cons_bc unique per cell, or repeated (i.e. multiple cells sharing a clone)?
length(unique(obj$cons_bc))
ncol(obj)

# If unique < ncol, we have single-cell data nested under clones — exactly what we want.
# Distribution of clone sizes (cells per cons_bc):
table(table(obj$cons_bc))

# Sanity check: does cons_bc (minus the _P/_1 patient suffix) match the "barcode"
# column format in ST223.rna.singlets.rds?
head(sort(table(obj$cons_bc), decreasing = TRUE), 10)

library(Seurat)
obj <- readRDS("/vast/projects/Sisseq/human-haematopoiesis-sis-seq/data/ST223.annotated_fate.rds")

# Is cons_bc unique per cell, or repeated (i.e. multiple cells sharing a clone)?
length(unique(obj$cons_bc))
ncol(obj)

# If unique < ncol, we have single-cell data nested under clones — exactly what we want.
# Distribution of clone sizes (cells per cons_bc):
table(table(obj$cons_bc))

# Sanity check: does cons_bc (minus the _P/_1 patient suffix) match the "barcode"
# column format in ST223.rna.singlets.rds?
head(sort(table(obj$cons_bc), decreasing = TRUE), 10)

## new qs asked
head(unique(obj$cons_bc), 10)
nchar(unique(obj$cons_bc))
