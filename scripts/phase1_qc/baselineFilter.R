# Phase 1: QC Report - Baseline


freq <- read.table(
  "../../results/phase1_qc/qc/baseline/res/qc_baseline.afreq",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)

snp_miss <- read.table(
  "../../results/phase1_qc/qc/baseline/res/qc_baseline.vmiss",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)

sample_miss <- read.table(
  "../../results/phase1_qc/qc/baseline/res/qc_baseline.smiss",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)


freq$MAF <- pmin(freq$ALT_FREQS, 1 - freq$ALT_FREQS)





cat("Minimum MAF:", min(freq$MAF), "\n")
cat("Maximum MAF:", max(freq$MAF), "\n\n")

cat("Minimum SNP Missingness:", min(snp_miss$F_MISS), "\n")
cat("Maximum SNP Missingness:", max(snp_miss$F_MISS), "\n\n")

cat("Minimum Sample Missingness:", min(sample_miss$F_MISS), "\n")
cat("Maximum Sample Missingness:", max(sample_miss$F_MISS), "\n")


hist(freq$MAF,breaks = 100,main = "Minor Allele Frequency Distribution",xlab = "Minor Allele Frequency",ylab = "Number of SNPs")

hist(snp_miss$F_MISS,breaks = 20,main = "Distribution of SNP Missingness",xlab = "Missingness Rate",ylab = "Number of SNPs")

hist(sample_miss$F_MISS,breaks = 20,main = "Distribution of Sample Missingness",xlab = "Missingness Rate",ylab = "Number of Samples")