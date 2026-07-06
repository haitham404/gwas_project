pcs <- read.table(
    "../../results/phase2_pca/pca/raw_pca.eigenvec",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

library(ggplot2)

p <- ggplot(
    pcs,
    aes(PC1, PC2)
) +
geom_point(size = 2) +
theme_classic()

ggsave(
    "../../plots/phase2_pca/pca_scatter.png",
    plot = p
)

eig <- scan("../../results/phase2_pca/pca/raw_pca.eigenval")
variance <- eig / sum(eig) * 100

df <- data.frame(
    PC = 1:10,
    Variance = variance
)

write.csv(
    df,
    "../../results/phase2_pca/pca/variance_explained.csv",
    row.names = FALSE
)

ggplot(df, aes(x = PC, y = Variance)) +
    geom_col() +
    geom_point() +
    geom_line() +
    theme_classic()

ggsave(
    "../../plots/phase2_pca/pca_scree.png",
    width = 7,
    height = 5
)