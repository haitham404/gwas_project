# ============================================================
# Phase 11 - Metabolite Quality Control
# ============================================================

# -----------------------------
# Paths
# -----------------------------

input_file <- "data/metabolites/Qatari_metabolomics(in).csv"
output_dir <- "results/phase11"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# -----------------------------
# Load data
# -----------------------------

df <- read.csv(
  input_file,
  check.names = FALSE
)

cat("Original shape:", nrow(df), "x", ncol(df), "\n")

id_col <- "mapped_id"
phenotype_col <- "Diabetes"

metabolite_cols <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

cat("Samples:", nrow(df), "\n")
cat("Metabolites:", length(metabolite_cols), "\n")


# ============================================================
# 1. Metabolite missingness QC
# ============================================================

met_missing <- colMeans(
  is.na(df[, metabolite_cols, drop = FALSE])
)

met_missing_df <- data.frame(
  metabolite = names(met_missing),
  missing_rate = met_missing
)

write.csv(
  met_missing_df,
  file.path(output_dir, "metabolite_missingness.csv"),
  row.names = FALSE
)

remove_metabolites <- names(
  met_missing[met_missing > 0.20]
)

cat(
  "\nMetabolites removed (>20% missing):",
  length(remove_metabolites),
  "\n"
)

df <- df[, !colnames(df) %in% remove_metabolites, drop = FALSE]


# ============================================================
# 2. Sample missingness QC
# ============================================================

remaining_metabolites <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

sample_missing <- rowMeans(
  is.na(df[, remaining_metabolites, drop = FALSE])
)

sample_missing_df <- data.frame(
  mapped_id = df[[id_col]],
  missing_rate = sample_missing
)

write.csv(
  sample_missing_df,
  file.path(output_dir, "sample_missingness.csv"),
  row.names = FALSE
)

remove_samples <- sample_missing > 0.20

cat(
  "Samples removed (>20% missing):",
  sum(remove_samples),
  "\n"
)

df <- df[!remove_samples, , drop = FALSE]


# ============================================================
# 3. Zero / near-zero variance QC
# ============================================================

remaining_metabolites <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

variances <- sapply(
  df[, remaining_metabolites, drop = FALSE],
  var,
  na.rm = TRUE
)

zero_variance <- names(
  variances[is.na(variances) | variances <= 1e-8]
)

cat(
  "Zero/near-zero variance metabolites:",
  length(zero_variance),
  "\n"
)

df <- df[, !colnames(df) %in% zero_variance, drop = FALSE]


# ============================================================
# 4. Export removed metabolites
# ============================================================

removed_df <- data.frame(
  metabolite = c(
    remove_metabolites,
    zero_variance
  ),
  reason = c(
    rep(
      "missingness_gt_20_percent",
      length(remove_metabolites)
    ),
    rep(
      "zero_or_near_zero_variance",
      length(zero_variance)
    )
  )
)

write.csv(
  removed_df,
  file.path(output_dir, "removed_metabolites.csv"),
  row.names = FALSE
)


# ============================================================
# 5. Export QC-passed dataset
# ============================================================

write.csv(
  df,
  file.path(output_dir, "metabolites_qc_passed.csv"),
  row.names = FALSE
)

final_metabolites <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

qc_passed_df <- data.frame(
  metabolite = final_metabolites
)

write.csv(
  qc_passed_df,
  file.path(output_dir, "qc_passed_metabolites.csv"),
  row.names = FALSE
)


# ============================================================
# QC Summary
# ============================================================

cat("\n")
cat("========== QC SUMMARY ==========\n")
cat("Final samples:", nrow(df), "\n")
cat("Final metabolites:", length(final_metabolites), "\n")
cat("Final shape:", nrow(df), "x", ncol(df), "\n")
cat("================================\n")


# ============================================================
# 6. Distribution and skewness analysis
# ============================================================

final_metabolites <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

calculate_skewness <- function(x) {

  x <- x[!is.na(x)]

  n <- length(x)

  if (n < 3) {
    return(NA)
  }

  mean_x <- mean(x)
  sd_x <- sd(x)

  if (sd_x == 0) {
    return(NA)
  }

  sum((x - mean_x)^3) / n / sd_x^3
}


skewness_values <- sapply(
  df[, final_metabolites, drop = FALSE],
  calculate_skewness
)


skewness_df <- data.frame(
  metabolite = names(skewness_values),
  skewness = skewness_values,
  absolute_skewness = abs(skewness_values)
)


skewness_df <- skewness_df[
  order(
    skewness_df$absolute_skewness,
    decreasing = TRUE
  ),
]


write.csv(
  skewness_df,
  file.path(output_dir, "metabolite_skewness.csv"),
  row.names = FALSE
)


cat("\n========== DISTRIBUTION SUMMARY ==========\n")

cat(
  "Metabolites with |skewness| > 1:",
  sum(abs(skewness_values) > 1, na.rm = TRUE),
  "\n"
)

cat(
  "Metabolites with |skewness| > 2:",
  sum(abs(skewness_values) > 2, na.rm = TRUE),
  "\n"
)

cat("\nTop 10 most skewed metabolites:\n")

print(
  head(skewness_df, 10)
)



# ============================================================
# 7. Z-score standardization
# ============================================================

final_metabolites <- setdiff(
  colnames(df),
  c(id_col, phenotype_col)
)

df_standardized <- df

df_standardized[, final_metabolites] <- scale(
  df[, final_metabolites, drop = FALSE]
)

write.csv(
  df_standardized,
  file.path(output_dir, "metabolites_standardized.csv"),
  row.names = FALSE
)

cat("\n========== STANDARDIZATION ==========\n")
cat("Standardized metabolites:", length(final_metabolites), "\n")

cat(
  "Mean range:",
  range(
    sapply(
      df_standardized[, final_metabolites, drop = FALSE],
      mean,
      na.rm = TRUE
    )
  ),
  "\n"
)

cat(
  "SD range:",
  range(
    sapply(
      df_standardized[, final_metabolites, drop = FALSE],
      sd,
      na.rm = TRUE
    )
  ),
  "\n"
)