library(ggplot2)
library(reshape2)

# =====================================================
# PATHS
# =====================================================

ld_file <- "results/phase8_ld/rs4131079_ld.unphased.vcor2"

vars_file <- "results/phase8_ld/rs4131079_ld.unphased.vcor2.vars"

output_file <- "plots/phase8_ld/ld_heatmap.png"


# =====================================================
# LOAD DATA
# =====================================================

ld <- as.matrix(
  read.table(
    ld_file,
    header = FALSE
  )
)

snps <- readLines(vars_file)


# =====================================================
# ADD SNP NAMES
# =====================================================

rownames(ld) <- snps
colnames(ld) <- snps


# =====================================================
# LONG FORMAT
# =====================================================

ld_long <- melt(ld)

colnames(ld_long) <- c(
  "SNP1",
  "SNP2",
  "R2"
)


# =====================================================
# PLOT
# =====================================================

p <- ggplot(
  ld_long,
  aes(
    x = SNP1,
    y = SNP2,
    fill = R2
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 1)
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    )
  ) +
  labs(
    title = "LD Heatmap Around rs4131079",
    subtitle = "Chr13: 43.86–44.86 Mb",
    x = "SNP",
    y = "SNP",
    fill = expression(R^2)
  )


# =====================================================
# SAVE
# =====================================================

ggsave(
  output_file,
  p,
  width = 12,
  height = 10,
  dpi = 300
)

print(p)

cat("\nLD heatmap saved to:\n")
cat(output_file, "\n")