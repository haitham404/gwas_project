library(tidyverse)

# ============================================================
# Phase 12: Metabolite-Diabetes Phenotype Association
# Step 0: Load and merge data
# ============================================================

metabolites <- read.csv(
  "results/phase11/metabolites_qc_passed.csv",
  check.names = FALSE
)

mapping <- read.csv(
  "data/metabolites/mapping(Sheet1).csv"
)

pca <- read.table(
  "results/phase2_pca/pca/raw_pca.eigenvec",
  header = TRUE,
  sep = "\t",
  comment.char = "",
  check.names = FALSE
)

colnames(pca)[colnames(pca) == "#FID"] <- "FID"

fam <- read.table(
  "data/Qatari156_filtered_pruned.fam",
  header = FALSE
)

colnames(fam) <- c(
  "FID",
  "IID",
  "Father",
  "Mother",
  "Sex",
  "Phenotype"
)


# ============================================================
# Map metabolite IDs to genotype IDs
# ============================================================

data <- metabolites %>%
  left_join(mapping, by = "mapped_id")


# ============================================================
# Add PCA covariates
# ============================================================

data <- data %>%
  left_join(
    pca,
    by = c("main_id" = "IID")
  )


# ============================================================
# Add sex
# ============================================================

data <- data %>%
  left_join(
    fam %>% select(IID, Sex),
    by = c("main_id" = "IID")
  )


# ============================================================
# Basic merge QC
# ============================================================

cat("Samples:", nrow(data), "\n")

cat(
  "Missing mapped IDs:",
  sum(is.na(data$main_id)),
  "\n"
)

cat(
  "Missing PCA:",
  sum(is.na(data$PC1)),
  "\n"
)

cat(
  "Missing Sex:",
  sum(is.na(data$Sex)),
  "\n"
)

cat("\nDiabetes distribution:\n")
print(table(data$Diabetes))

cat("\nSex distribution:\n")
print(table(data$Sex))


# ============================================================
# Save merged analysis data
# ============================================================

write.csv(
  data,
  "results/phase12/phase12_analysis_data.csv",
  row.names = FALSE
)




# ============================================================
# Step 1: Exploratory bar plot
# Metabolite levels by diabetes status
# ============================================================

metadata_cols <- c(
  "mapped_id",
  "Diabetes",
  "main_id",
  "FID",
  "PC1", "PC2", "PC3", "PC4", "PC5",
  "PC6", "PC7", "PC8", "PC9", "PC10",
  "Sex"
)

metabolite_cols <- setdiff(
  colnames(data),
  metadata_cols
)

cat(
  "\nNumber of metabolites:",
  length(metabolite_cols),
  "\n"
)


# Convert metabolite data from wide to long format
metabolite_long <- data %>%
  select(Diabetes, all_of(metabolite_cols)) %>%
  pivot_longer(
    cols = all_of(metabolite_cols),
    names_to = "Metabolite",
    values_to = "Level"
  )


# Calculate mean and standard error
bar_summary <- metabolite_long %>%
  group_by(Metabolite, Diabetes) %>%
  summarise(
    Mean = mean(Level, na.rm = TRUE),
    SD = sd(Level, na.rm = TRUE),
    N = sum(!is.na(Level)),
    SE = SD / sqrt(N),
    .groups = "drop"
  )


# Label diabetes groups
bar_summary$Diabetes <- factor(
  bar_summary$Diabetes,
  levels = c(0, 1),
  labels = c("Non-diabetic", "Diabetic")
)


# Save summary table
write.csv(
  bar_summary,
  "results/phase12/metabolite_diabetes_bar_summary.csv",
  row.names = FALSE
)


# Create exploratory bar plot
# Create exploratory bar plot
p <- ggplot(
  bar_summary,
  aes(
    x = Diabetes,
    y = Mean,
    fill = Diabetes
  )
) +
  geom_col(
    width = 0.7
  ) +

  # Mean ± SE
  geom_errorbar(
    aes(
      ymin = Mean - SE,
      ymax = Mean + SE
    ),
    width = 0.2,
    linewidth = 0.4
  ) +

  # Show exact mean above each bar
  geom_text(
    aes(
      label = sprintf("%.2f", Mean)
    ),
    vjust = ifelse(bar_summary$Mean >= 0, -0.5, 1.5),
    size = 2.5,
    fontface = "bold"
  ) +

  facet_wrap(
    ~ Metabolite,
    scales = "free_y",
    ncol = 4
  ) +

  labs(
    title = "Metabolite Levels by Diabetes Status",
    subtitle = "Bars show mean standardized metabolite level; error bars show ± SE",
    x = "Diabetes Status",
    y = "Mean Standardized Metabolite Level",
    fill = "Diabetes Status"
  ) +

  theme_bw() +

  theme(
    legend.position = "bottom",

    strip.text = element_text(
      size = 7,
      face = "bold"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 6
    ),

    plot.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "plots/phase12/metabolite_diabetes_barplot.pdf",
  plot = p,
  width = 20,
  height = 45
  )

cat(
  "Exploratory bar plot saved successfully.\n"
)





# ============================================================
# Step 2: Formal metabolite-diabetes association analysis
# Model:
# Metabolite ~ Diabetes + Sex + PC1 + PC2 + PC3
# ============================================================

association_results <- map_dfr(
  metabolite_cols,
  function(metabolite) {

    model_data <- data %>%
      select(
        Diabetes,
        Sex,
        PC1,
        PC2,
        PC3,
        all_of(metabolite)
      ) %>%
      drop_na()

    colnames(model_data)[
      colnames(model_data) == metabolite
    ] <- "MetaboliteLevel"

    model <- lm(
      MetaboliteLevel ~ Diabetes + Sex + PC1 + PC2 + PC3,
      data = model_data
    )

    model_summary <- summary(model)$coefficients

    tibble(
      Metabolite = metabolite,
      Beta = model_summary["Diabetes", "Estimate"],
      SE = model_summary["Diabetes", "Std. Error"],
      P_value = model_summary["Diabetes", "Pr(>|t|)"]
    )
  }
)


# ============================================================
# Step 3: Bonferroni correction
# ============================================================

number_of_metabolites <- nrow(association_results)

bonferroni_threshold <- 0.05 / number_of_metabolites

association_results <- association_results %>%
  mutate(
    Bonferroni_significant = ifelse(
      P_value < bonferroni_threshold,
      "Y",
      "N"
    )
  ) %>%
  arrange(P_value)


# ============================================================
# Save full association table
# ============================================================

write.csv(
  association_results,
  "results/phase12/full_metabolite_diabetes_association.csv",
  row.names = FALSE
)


# ============================================================
# Print summary
# ============================================================

cat(
  "\nNumber of metabolites tested:",
  number_of_metabolites,
  "\n"
)

cat(
  "Bonferroni threshold:",
  bonferroni_threshold,
  "\n"
)

cat(
  "Bonferroni-significant metabolites:",
  sum(
    association_results$Bonferroni_significant == "Y"
  ),
  "\n"
)

cat("\nTop associations:\n")

print(
  association_results %>%
    head(10)
)




# ============================================================
# Step 4: Volcano plot
# ============================================================

association_results <- association_results %>%
  mutate(
    NegLog10P = -log10(P_value),
    Significance = ifelse(
      Bonferroni_significant == "Y",
      "Bonferroni significant",
      "Not significant"
    )
  )

volcano_plot <- ggplot(
  association_results,
  aes(
    x = Beta,
    y = NegLog10P,
    color = Significance
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  geom_hline(
    yintercept = -log10(bonferroni_threshold),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Metabolite-Diabetes Association Volcano Plot",
    subtitle = paste(
      "Bonferroni threshold =",
      signif(bonferroni_threshold, 3)
    ),
    x = "Diabetes Effect Size (Beta)",
    y = "-log10(P-value)",
    color = "Association"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )


ggsave(
  "plots/phase12/metabolite_diabetes_volcano_plot.pdf",
  plot = volcano_plot,
  width = 10,
  height = 8
)


# ============================================================
# Save significant metabolites
# ============================================================

significant_metabolites <- association_results %>%
  filter(Bonferroni_significant == "Y") %>%
  arrange(P_value)


write.csv(
  significant_metabolites,
  "results/phase12/significant_metabolites.csv",
  row.names = FALSE
)


cat(
  "\nVolcano plot saved successfully.\n"
)

cat(
  "Significant metabolites saved:",
  nrow(significant_metabolites),
  "\n"
)