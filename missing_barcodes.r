library(dplyr)
library(readr)
library(stringr)
library(Seurat)

# which of 112 bc absent from single-cell RNA data
# at which stage r they missing? 
# (never captured vs dropped during RNA+ADT merge/QC).


FATE_CSV <- "/vast/projects/Sisseq/human-haematopoiesis-sis-seq/mofa_kj/files/barcodes_with_fate.csv"
RAW_RDS  <- "/vast/projects/Sisseq/human-haematopoiesis-sis-seq/data/ST223.rna.singlets.rds"
FATE_RDS <- "/vast/projects/Sisseq/human-haematopoiesis-sis-seq/data/ST223.annotated_fate.rds"

fate <- read_csv(FATE_CSV, show_col_types = FALSE)
fate$barcode_core <- sub("_[^_]+$", "", fate$barcode)     # strip donor suffix
stopifnot(all(nchar(fate$barcode_core) == 15))            
stopifnot(!any(duplicated(fate$barcode)))                 # confirm 112 unique

n_fate_total <- nrow(fate)
cat(sprintf("Total fate-sorted barcodes to trace: %d\n", n_fate_total))


# there is alignment offset in data. avg file made in st Figure3_barcode_clustering.R
# where align?


raw <- readRDS(RAW_RDS)
raw_barcodes <- unique(na.omit(raw$barcode))   
stopifnot(all(nchar(raw_barcodes) == 15))

# issues here!! 2 of 3,040 unique raw barcodes are malformed "GATTATGCA" at 9 chars, "ATTATGCA" at 8 chars . 

sc <- readRDS(FATE_RDS)
sc$cons_bc_core <- sub("_[^_]+$", "", sc$cons_bc)
core_len <- unique(nchar(sc$cons_bc_core))
stopifnot(length(core_len) == 1)   # all cons_bc_core same length

offset_results <- data.frame()
for (start in 1:(core_len - 15 + 1)) {
  window <- unique(substr(sc$cons_bc_core, start, start + 14))
  match_count <- sum(fate$barcode_core %in% window)
  offset_results <- rbind(offset_results, data.frame(start = start, n_exact_matches = match_count))
}
print(offset_results)

best_offset <- offset_results$start[which.max(offset_results$n_exact_matches)]
cat(sprintf("\nBest alignment offset: start = %d (exact matches = %d / %d)\n",
            best_offset, max(offset_results$n_exact_matches), n_fate_total))
stopifnot(max(offset_results$n_exact_matches) > offset_results$n_exact_matches[offset_results$start != best_offset] %>% max())

# one clear winner offset=1 but still lots of barcodes unaccounted for. 

sc$barcode15 <- substr(sc$cons_bc_core, best_offset, best_offset + 14)

# trace where each goes missing ?? 
report <- fate %>%
  mutate(
    in_raw_pool        = barcode_core %in% raw_barcodes,
    in_fate_rds        = barcode_core %in% sc$barcode15,
    n_cells_in_raw_pool = sapply(barcode_core, function(b) sum(raw$barcode == b, na.rm = TRUE)),
    n_cells_in_fate_rds = sapply(barcode_core, function(b) sum(sc$barcode15 == b, na.rm = TRUE))
  ) %>%
  mutate(
    status = case_when(
      in_raw_pool  & in_fate_rds  ~ "present_end_to_end",
      in_raw_pool  & !in_fate_rds ~ "dropped_during_RNA_ADT_merge_or_QC",
      !in_raw_pool                ~ "never_captured_in_single_cell_pool"
    )
  )


                                 
summary_table <- report %>% count(status) %>% mutate(pct = round(100 * n / n_fate_total, 1))
print(summary_table)

stopifnot(sum(summary_table$n) == n_fate_total)  # are all bc present in table 

cat(sprintf(
  "\nFINAL: %d / %d (%.1f%%) fate-sorted barcodes are absent from the single-cell RNA object (%s).\n",
  sum(report$status != "present_end_to_end"), n_fate_total,
  100 * sum(report$status != "present_end_to_end") / n_fate_total,
  basename(FATE_RDS)
))

