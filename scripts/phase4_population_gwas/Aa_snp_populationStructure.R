eigenvec <- read.table(
    "../../results/phase2_pca/pca/raw_pca.eigenvec",
    header = FALSE
)

colnames(eigenvec) <- c("FID", "IID", paste0("PC", 1:(ncol(eigenvec)-2)))

write.table(
    eigenvec[, c("FID", "IID", "PC1")],
    "pc1.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

write.table(
    eigenvec[, c("FID", "IID", "PC2")],
    "pc2.txt",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

library(qqman)

assoc <- read.table(
    "assoc_pc2.PC2.glm.linear",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

assoc$`#CHROM`[assoc$`#CHROM` == "X"] <- "23"
assoc$`#CHROM` <- as.numeric(assoc$`#CHROM`)

png(
    "../../plots/phase4_population_gwas/phase4_manhattan.png",
    width = 800,
    height = 600
)

manhattan(
    assoc,
    chr = "#CHROM",
    bp = "POS",
    snp = "ID",
    p = "P",
    main = "Phase 4 - Linear GWAS (PC2)"
)

dev.off()

assoc <- read.table(
    "assoc_pc1.PC1.glm.linear",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

assoc$`#CHROM`[assoc$`#CHROM` == "X"] <- "23"
assoc$`#CHROM` <- as.numeric(assoc$`#CHROM`)

png(
    "../../plots/phase4_population_gwas/phase4_manhattan_pc1.png",
    width = 800,
    height = 600
)

manhattan(
    assoc,
    chr = "#CHROM",
    bp = "POS",
    snp = "ID",
    p = "P",
    main = "Phase 4 - Linear GWAS (PC1)"
)

dev.off()


assoc <- read.table(
    "assoc_pc1_with_pc2.PC1.glm.linear",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

assoc$`#CHROM`[assoc$`#CHROM` == "X"] <- "23"
assoc$`#CHROM` <- as.numeric(assoc$`#CHROM`)

png(
    "../../plots/phase4_population_gwas/phase4_manhattan_pc1_with_pc2.png",
    width = 800,
    height = 600
)

manhattan(
    assoc,
    chr = "#CHROM",
    bp = "POS",
    snp = "ID",
    p = "P",
    main = "Phase 4 - Linear GWAS (PC1 + PC2 Covariate)"
)

dev.off()


assoc <- read.table(
    "assoc_pc2_with_pc1.PC2.glm.linear",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

assoc$`#CHROM`[assoc$`#CHROM` == "X"] <- "23"
assoc$`#CHROM` <- as.numeric(assoc$`#CHROM`)

png(
    "../../plots/phase4_population_gwas/phase4_manhattan_pc2_with_pc1.png",
    width = 800,
    height = 600
)

manhattan(
    assoc,
    chr = "#CHROM",
    bp = "POS",
    snp = "ID",
    p = "P",
    main = "Phase 4 - Linear GWAS (PC2 + PC1 Covariate)"
)

dev.off()