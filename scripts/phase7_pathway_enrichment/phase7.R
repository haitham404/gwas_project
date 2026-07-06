# ==========================================
# Phase 7 - Pathway Enrichment Analysis
# Human - using gprofiler2
# ==========================================

library(gprofiler2)
library(ggplot2)

# ---------------------------------------------------------
# STEP 1: Load Phase 6 annotated SNPs (FIXED PATH)
# ---------------------------------------------------------
snps <- read.csv("../../results/phase6_annotation/phase6_annotated_snps.csv")

cat("Rows in Phase6:", nrow(snps), "\n")

# ---------------------------------------------------------
# STEP 2: Extract clean gene list
# ---------------------------------------------------------
genes <- unique(trimws(snps$GENE))
genes <- genes[genes != ""]
genes <- genes[!is.na(genes)]

cat("Unique genes:", length(genes), "\n")

# Safety check
if (length(genes) == 0) {
  stop("No valid genes found for enrichment analysis")
}

# ---------------------------------------------------------
# STEP 3: GO/Pathway enrichment
# ---------------------------------------------------------
res <- gost(
  query = genes,
  organism = "hsapiens",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  correction_method = "fdr"
)

# ---------------------------------------------------------
# STEP 4: Handle empty results safely
# ---------------------------------------------------------
if (is.null(res$result) || nrow(res$result) == 0) {
  cat("No enriched pathways found.\n")
  quit()
}

cat("Enriched pathways:", nrow(res$result), "\n")

# ---------------------------------------------------------
# STEP 5: Clean results table
# ---------------------------------------------------------
result <- as.data.frame(res$result)

# remove list columns (important fix)
result <- result[, !sapply(result, is.list)]

# ---------------------------------------------------------
# STEP 6: Save results (FIXED PATH)
# ---------------------------------------------------------
write.csv(
  result,
  "../results/phase7_pathway_enrichment.csv",
  row.names = FALSE
)

cat("Results saved successfully\n")

# ---------------------------------------------------------
# STEP 7: Dot plot (FIXED safe save)
# ---------------------------------------------------------
png("../plots/phase7_dotplot.png", width = 1200, height = 800)

gostplot(res, interactive = FALSE)

dev.off()

cat("Dot plot saved successfully\n")

# ---------------------------------------------------------
# STEP 8: Top pathways
# ---------------------------------------------------------
cat("\nTop Enriched Pathways:\n")

print(
  res$result[
    1:min(10, nrow(res$result)),
    c("source", "term_name", "p_value", "intersection_size")
  ]
)

# ---------------------------------------------------------
# DEBUG INFO
# ---------------------------------------------------------
cat("\nUnique genes used:\n")
print(length(unique(na.omit(snps$GENE))))

cat("\n========== Phase 7 Finished ==========\n")