# =========================================================
# Phase 6 - GWAS Annotation + TRUE Manhattan (ggplot)
# =========================================================

library(biomaRt)
library(dplyr)
library(ggplot2)
library(ggrepel)

options(scipen = 999)

# ---------------------------------------------------------
# STEP 0: Correct paths (FIXED for your tree)
# ---------------------------------------------------------
pc1 <- read.table("../phase4_population_gwas/assoc_pc1.PC1.glm.linear",
                   header = TRUE, comment.char = "")

pc2 <- read.table("../phase4_population_gwas/assoc_pc2.PC2.glm.linear",
                   header = TRUE, comment.char = "")

sex <- read.table("../../results/phase5_sex_gwas/sex_assoc.SEX.glm.logistic.hybrid",
                  header = TRUE, comment.char = "")

# ---------------------------------------------------------
# STEP 1: Standardize column names (robust fix)
# ---------------------------------------------------------
rename_cols <- function(df) {

  if ("#CHROM" %in% colnames(df)) {
    names(df)[names(df) == "#CHROM"] <- "CHR"
  }

  if ("X.CHROM" %in% colnames(df)) {
    names(df)[names(df) == "X.CHROM"] <- "CHR"
  }

  if ("BP" %in% colnames(df)) {
    names(df)[names(df) == "BP"] <- "POS"
  }

  if ("POS" %in% colnames(df)) {
    names(df)[names(df) == "POS"] <- "BP"
  }

  if ("ID" %in% colnames(df)) {
    names(df)[names(df) == "ID"] <- "SNP"
  }

  return(df)
}

pc1 <- rename_cols(pc1)
pc2 <- rename_cols(pc2)
sex <- rename_cols(sex)

# ---------------------------------------------------------
# STEP 2: Fix chromosome encoding
# ---------------------------------------------------------
fix_chr <- function(df) {

  df$CHR <- as.character(df$CHR)

  df$CHR[df$CHR == "X"] <- "23"
  df$CHR[df$CHR == "Y"] <- "24"
  df$CHR[df$CHR %in% c("MT", "M")] <- "25"

  df$CHR <- as.numeric(df$CHR)
  df <- df[!is.na(df$CHR), ]

  return(df)
}

pc1 <- fix_chr(pc1)
pc2 <- fix_chr(pc2)
sex <- fix_chr(sex)

# ---------------------------------------------------------
# STEP 3: Clean GWAS data
# ---------------------------------------------------------
clean_gwas <- function(df) {
  df <- df[!is.na(df$P), ]
  df <- df[is.finite(df$P), ]
  return(df)
}

pc1 <- clean_gwas(pc1)
pc2 <- clean_gwas(pc2)
sex <- clean_gwas(sex)

# ---------------------------------------------------------
# STEP 4: Select significant SNPs
# ---------------------------------------------------------
threshold <- 1e-5

pc1_top <- subset(pc1, P < threshold)
pc2_top <- subset(pc2, P < threshold)
sex_top <- subset(sex, P < threshold)

top_snps <- unique(c(
  pc1_top$SNP,
  pc2_top$SNP,
  sex_top$SNP
))

cat("Selected SNPs:", length(top_snps), "\n")

# ---------------------------------------------------------
# STEP 5: Gene annotation (biomaRt FIXED)
# ---------------------------------------------------------
mart <- useEnsembl(
  biomart = "snp",
  dataset = "hsapiens_snp",
  mirror  = "useast"
)

gene_info <- getBM(
  attributes = c(
    "refsnp_id",
    "chr_name",
    "chrom_start",
    "ensembl_gene_stable_id",
    "consequence_type_tv"
  ),
  filters = "snp_filter",
  values  = top_snps,
  mart    = mart
)

colnames(gene_info) <- c("SNP", "CHR", "BP", "GENE", "CONSEQUENCE")

# ---------------------------------------------------------
# STEP 6: Save results (FIXED PATH)
# ---------------------------------------------------------
write.csv(
  gene_info,
  "../../results/phase6_annotation/phase6_annotated_snps.csv",
  row.names = FALSE
)

# ---------------------------------------------------------
# DONE
# ---------------------------------------------------------
cat("Phase 6 completed successfully\n")