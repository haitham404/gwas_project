freq <- read.table(
  "../../results/phase1_qc/freq.afreq",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)

snp_miss <- read.table(
  "../../results/phase1_qc/missing.vmiss",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)

sample_miss <- read.table(
  "../../results/phase1_qc/missing.smiss",
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)




freq$MAF <- pmin(freq$ALT_FREQS, 1 - freq$ALT_FREQS)




min(freq$MAF)
max(freq$MAF)
min(snp_miss$F_MISS)
max(snp_miss$F_MISS)
min(sample_miss$F_MISS)
max(sample_miss$F_MISS)
hist(freq$MAF, breaks = 100, main = "Minor Allele Frequency Distribution", xlab = "Minor Allele Frequency", ylab = "Count" )
hist(snp_miss$F_MISS,breaks = 20, main = "Distribution of SNP Missingness", xlab = "Missingness Rate", ylab = "Number of SNPs" )
hist(sample_miss$F_MISS,breaks = 20, main = "Distribution of Sapmle Missingness", xlab = "Missingness Rate", ylab = "Number of Samples" )
