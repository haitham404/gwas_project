library(ggplot2)
library(plotly)
library(htmlwidgets)

pca <- read.table(
    "../../results/phase2_pca/pca/raw_pca.eigenvec",
    header = TRUE,
    comment.char = "",
    check.names = FALSE
)

pcs <- pca[, c("PC1", "PC2", "PC3")]
wss <- c()

for(i in 1:10){
    result <- kmeans(pcs, centers = i, nstart = 25)
    wss[i] <- result$tot.withinss
}

pdf("../../plots/phase3_clustering/elbow_plot.pdf")
plot(
    1:10,
    wss,
    type = "b",
    xlab = "Number of Clusters",
    ylab = "WSS",
    main = "Elbow Method"
)
dev.off()

pcs1 <- pca[, c("PC1", "PC2")]
pcs2 <- pca[, c("PC1", "PC2", "PC3")]

Kmeans1 <- kmeans(pcs1, centers = 4, nstart = 25)
Kmeans2 <- kmeans(pcs2, centers = 4, nstart = 25)

pca$Cluster1 <- as.factor(Kmeans1$cluster)
pca$Cluster2 <- as.factor(Kmeans2$cluster)

ggplot(pca, aes(x = PC1, y = PC2, color = Cluster1)) +
    geom_point(size = 2) +
    labs(title = "PCA Clustering (2D)", x = "PC1", y = "PC2") +
    theme_minimal()

ggsave(
    "../../plots/phase3_clustering/pca_2d.png",
    width = 7,
    height = 5
)

p <- plot_ly(
    data = pca,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~Cluster2,
    type = "scatter3d",
    mode = "markers"
)

saveWidget(
    p,
    "../../plots/phase3_clustering/pca_3d.html",
    selfcontained = FALSE
)