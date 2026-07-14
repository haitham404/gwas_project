library(dplyr)
library(ggplot2)

# ============================================================
# Paths
# ============================================================

KINSHIP_FILE <- "results/phase9_kinship/king.kin0"

RESULT_DIR <- "results/phase9_kinship"
PLOT_DIR <- "plots/phase9_kinship"

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Load KING kinship table
# ============================================================

kinship <- read.table(
  KINSHIP_FILE,
  header = TRUE,
  comment.char = "",
  check.names = FALSE
)


names(kinship) <- sub("^#", "", names(kinship))

print(head(kinship))
print(names(kinship))


# ============================================================
# Classify relatedness
# ============================================================

kinship <- kinship %>%
  mutate(
    Relationship = case_when(
      KINSHIP > 0.354  ~ "Duplicate/MZ twin",
      KINSHIP >= 0.177 ~ "1st-degree",
      KINSHIP >= 0.0884 ~ "2nd-degree",
      KINSHIP >= 0.0442 ~ "3rd-degree",
      TRUE ~ "Unrelated"
    )
  )


# ============================================================
# Count pairs
# ============================================================

category_order <- c(
  "Duplicate/MZ twin",
  "1st-degree",
  "2nd-degree",
  "3rd-degree",
  "Unrelated"
)

pair_counts <- kinship %>%
  count(Relationship, name = "Number_of_pairs") %>%
  right_join(
    data.frame(Relationship = category_order),
    by = "Relationship"
  ) %>%
  mutate(
    Number_of_pairs = ifelse(
      is.na(Number_of_pairs),
      0,
      Number_of_pairs
    ),
    Relationship = factor(
      Relationship,
      levels = category_order
    )
  ) %>%
  arrange(Relationship)

print(pair_counts)

write.csv(
  pair_counts,
  file.path(RESULT_DIR, "relatedness_pair_counts.csv"),
  row.names = FALSE
)


# ============================================================
# Relative pairs
# ============================================================

relative_pairs <- kinship %>%
  filter(Relationship != "Unrelated")

write.csv(
  relative_pairs,
  file.path(RESULT_DIR, "possible_relative_pairs.csv"),
  row.names = FALSE
)

n_relative_pairs <- nrow(relative_pairs)


# ============================================================
# Individuals involved in relative pairs
# ============================================================

relative_individuals <- unique(
  c(
    as.character(relative_pairs$IID1),
    as.character(relative_pairs$IID2)
  )
)

n_relative_individuals <- length(relative_individuals)

writeLines(
  as.character(relative_individuals),
  file.path(RESULT_DIR, "relative_individuals.txt")
)

# ============================================================
# Build kinship matrix
# ============================================================

ids <- unique(
  c(
    as.character(kinship$IID1),
    as.character(kinship$IID2)
  )
)

kinship_matrix <- matrix(
  0,
  nrow = length(ids),
  ncol = length(ids),
  dimnames = list(ids, ids)
)

diag(kinship_matrix) <- 0.5

for (i in seq_len(nrow(kinship))) {

  id1 <- as.character(kinship$IID1[i])
  id2 <- as.character(kinship$IID2[i])
  k <- kinship$KINSHIP[i]

  kinship_matrix[id1, id2] <- k
  kinship_matrix[id2, id1] <- k
}


# ============================================================
# Convert matrix to long format
# ============================================================

heatmap_data <- expand.grid(
  Individual1 = rownames(kinship_matrix),
  Individual2 = colnames(kinship_matrix),
  stringsAsFactors = FALSE
)

heatmap_data$Kinship <- as.vector(kinship_matrix)

print(dim(kinship_matrix))
print(head(heatmap_data))


# ============================================================
# Kinship heatmap
# ============================================================

p <- ggplot(
  heatmap_data,
  aes(
    x = Individual1,
    y = Individual2,
    fill = Kinship
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  labs(
    title = "KING-Robust Kinship Matrix",
    subtitle = "Qatari QC'd Cohort",
    x = "Individual",
    y = "Individual",
    fill = "Kinship"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  file.path(PLOT_DIR, "kinship_heatmap.png"),
  p,
  width = 9,
  height = 8,
  dpi = 300
)


# ============================================================
# Final report
# ============================================================

cat("\n")
cat("====================================================\n")
cat("PHASE 9 KINSHIP SUMMARY\n")
cat("====================================================\n")

print(pair_counts)

cat("\n")
cat(
  n_relative_pairs,
  "possible relative pairs identified in the cohort.\n"
)

cat(
  n_relative_individuals,
  "individuals are involved in at least one non-unrelated pair.\n"
)

cat("\n")

if (n_relative_pairs > 0) {

  cat(
    "Close relatives were detected. ",
    "These relationships should be considered before Phase 11. ",
    "A mixed-model analysis can account for genetic relatedness ",
    "through the GRM; duplicate samples should be investigated ",
    "and potentially excluded.\n"
  )

} else {

  cat(
    "No possible relative pairs were identified. ",
    "No relatedness-based sample exclusion is indicated before Phase 11.\n"
  )
}