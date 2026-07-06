assoc <- read.table(
    "../../results/phase5_sex_gwas/sex_assoc.SEX.glm.logistic.hybrid",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

# ----------------------------
# Fix chromosome format (X -> 23)
# ----------------------------
assoc$`#CHROM`[assoc$`#CHROM` == "X"] <- "23"
assoc$`#CHROM` <- as.numeric(assoc$`#CHROM`)

# ----------------------------
# Load qqman
# ----------------------------
library(qqman)

# ----------------------------
# Manhattan Plot (Sex GWAS)
# ----------------------------
png(
    "../../plots/phase5_sex_gwas/phase5_manhattan.png",
    width = 800,
    height = 600
)

manhattan(
    assoc,
    chr = "#CHROM",
    bp = "POS",
    snp = "ID",
    p = "P",
    main = "Phase 5 - Logistic GWAS (Sex)"
)

dev.off()

# ----------------------------
# QQ Plot
# ----------------------------
png(
    "../../plots/phase5_sex_gwas/phase5_qq.png",
    width = 800,
    height = 800
)

qq(assoc$P)

dev.off()

# ----------------------------
# Top hits export
# ----------------------------
assoc_sorted <- assoc[order(assoc$P), ]

write.table(
    assoc_sorted[1:20, ],
    "../../results/phase5_sex_gwas/top_hits.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)