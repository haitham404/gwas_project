library(tidyverse)
library(randomForest)
library(pROC)

set.seed(42)

# ============================================================
# Directories
# ============================================================

dir.create(
  "results/phase13",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "plots/phase13",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Load data
# ============================================================

metabolites <- read.csv(
  "results/phase11/metabolites_qc_passed.csv",
  check.names = FALSE
)

mapping <- read.csv(
  "data/metabolites/mapping(Sheet1).csv"
)

pca <- read.delim(
  "results/phase2_pca/pca/raw_pca.eigenvec",
  check.names = FALSE
)

fam <- read.table(
  "data/Qatari156_filtered_pruned.fam"
)


# ============================================================
# Prepare metadata
# ============================================================

colnames(fam) <- c(
  "FID",
  "IID",
  "Father",
  "Mother",
  "Sex",
  "Phenotype"
)

colnames(pca)[1:2] <- c(
  "FID",
  "IID"
)


# ============================================================
# Merge data
# ============================================================

data <- metabolites %>%
  left_join(mapping, by = "mapped_id") %>%
  left_join(
    pca,
    by = c("main_id" = "IID")
  ) %>%
  left_join(
    fam %>% select(IID, Sex),
    by = c("main_id" = "IID")
  )


# ============================================================
# Define metabolite predictors
# ============================================================

metadata_cols <- c(
  "mapped_id",
  "main_id",
  "Diabetes",
  "FID",
  "Sex",
  paste0("PC", 1:10)
)

metabolite_cols <- setdiff(
  colnames(data),
  metadata_cols
)


cat(
  "Samples:",
  nrow(data),
  "\n"
)

cat(
  "Metabolites:",
  length(metabolite_cols),
  "\n"
)

cat("\nDiabetes distribution:\n")

print(
  table(data$Diabetes)
)


# ============================================================
# Predictor dataset
# ============================================================

model_data <- data %>%
  select(
    Diabetes,
    Sex,
    PC1,
    PC2,
    PC3,
    all_of(metabolite_cols)
  ) %>%
  drop_na()


cat(
  "\nSamples available for classification:",
  nrow(model_data),
  "\n"
)



# ============================================================
# Step 2: Stratified 80/20 train-test split
# ============================================================

set.seed(42)

train_data <- model_data %>%
  group_by(Diabetes) %>%
  slice_sample(prop = 0.80) %>%
  ungroup()

test_data <- model_data %>%
  anti_join(
    train_data,
    by = colnames(model_data)
  )


cat(
  "\nTraining samples:",
  nrow(train_data),
  "\n"
)

cat(
  "Testing samples:",
  nrow(test_data),
  "\n"
)


cat("\nTraining diabetes distribution:\n")

print(
  table(train_data$Diabetes)
)


cat("\nTesting diabetes distribution:\n")

print(
  table(test_data$Diabetes)
)


# ============================================================
# Save split
# ============================================================

write.csv(
  train_data,
  "results/phase13/train_data.csv",
  row.names = FALSE
)

write.csv(
  test_data,
  "results/phase13/test_data.csv",
  row.names = FALSE
)



# ============================================================
# Step 3A: Logistic Regression
# ============================================================

logistic_model <- glm(
  Diabetes ~ .,
  data = train_data,
  family = binomial(link = "logit")
)


# ============================================================
# Predict on test set
# ============================================================

logistic_probability <- predict(
  logistic_model,
  newdata = test_data,
  type = "response"
)

logistic_prediction <- ifelse(
  logistic_probability >= 0.5,
  1,
  0
)


# ============================================================
# Save test predictions
# ============================================================

logistic_predictions <- tibble(
  Actual = test_data$Diabetes,
  Predicted = logistic_prediction,
  Probability = logistic_probability
)

write.csv(
  logistic_predictions,
  "results/phase13/logistic_test_predictions.csv",
  row.names = FALSE
)


cat("\nLogistic Regression predictions:\n")

print(logistic_predictions)




# ============================================================
# Step 4: Evaluate Logistic Regression
# ============================================================

actual <- test_data$Diabetes
predicted <- logistic_prediction


# Confusion matrix values

TP <- sum(actual == 1 & predicted == 1)
TN <- sum(actual == 0 & predicted == 0)
FP <- sum(actual == 0 & predicted == 1)
FN <- sum(actual == 1 & predicted == 0)


cat("\nConfusion Matrix Values:\n")

cat("TP:", TP, "\n")
cat("TN:", TN, "\n")
cat("FP:", FP, "\n")
cat("FN:", FN, "\n")


# ============================================================
# Classification metrics
# ============================================================

accuracy <- (TP + TN) / length(actual)

sensitivity <- TP / (TP + FN)

specificity <- TN / (TN + FP)

balanced_accuracy <- (
  sensitivity + specificity
) / 2


# ============================================================
# AUC
# ============================================================

roc_object <- roc(
  actual,
  logistic_probability,
  quiet = TRUE
)

auc_value <- as.numeric(
  auc(roc_object)
)


# ============================================================
# Metrics table
# ============================================================

logistic_metrics <- tibble(
  Model = "Logistic Regression",
  Accuracy = accuracy,
  Sensitivity = sensitivity,
  Specificity = specificity,
  Balanced_Accuracy = balanced_accuracy,
  AUC = auc_value
)


print(logistic_metrics)


# ============================================================
# Save metrics
# ============================================================

write.csv(
  logistic_metrics,
  "results/phase13/classification_metrics.csv",
  row.names = FALSE
)


# ============================================================
# Step 3B: Random Forest
# ============================================================

train_rf <- train_data
test_rf <- test_data

train_rf$Diabetes <- factor(
  train_rf$Diabetes,
  levels = c(0, 1)
)

test_rf$Diabetes <- factor(
  test_rf$Diabetes,
  levels = c(0, 1)
)


# ============================================================
# Train Random Forest
# ============================================================

set.seed(42)

rf_model <- randomForest(
  Diabetes ~ .,
  data = train_rf,
  ntree = 500,
  importance = TRUE
)

cat("\nRandom Forest Model:\n")

print(rf_model)


# ============================================================
# Test predictions
# ============================================================

rf_probability <- predict(
  rf_model,
  newdata = test_rf,
  type = "prob"
)[, "1"]

rf_prediction <- predict(
  rf_model,
  newdata = test_rf,
  type = "response"
)


rf_actual <- as.numeric(
  as.character(test_rf$Diabetes)
)

rf_predicted <- as.numeric(
  as.character(rf_prediction)
)


# ============================================================
# Confusion matrix
# ============================================================

rf_TP <- sum(
  rf_actual == 1 &
  rf_predicted == 1
)

rf_TN <- sum(
  rf_actual == 0 &
  rf_predicted == 0
)

rf_FP <- sum(
  rf_actual == 0 &
  rf_predicted == 1
)

rf_FN <- sum(
  rf_actual == 1 &
  rf_predicted == 0
)


cat("\nRandom Forest Confusion Matrix Values:\n")

cat("TP:", rf_TP, "\n")
cat("TN:", rf_TN, "\n")
cat("FP:", rf_FP, "\n")
cat("FN:", rf_FN, "\n")


# ============================================================
# Evaluation metrics
# ============================================================

rf_accuracy <- (
  rf_TP + rf_TN
) / length(rf_actual)

rf_sensitivity <- rf_TP / (
  rf_TP + rf_FN
)

rf_specificity <- rf_TN / (
  rf_TN + rf_FP
)

rf_balanced_accuracy <- (
  rf_sensitivity +
  rf_specificity
) / 2


rf_roc <- roc(
  rf_actual,
  rf_probability,
  quiet = TRUE
)

rf_auc <- as.numeric(
  auc(rf_roc)
)


# ============================================================
# Random Forest metrics
# ============================================================

rf_metrics <- tibble(
  Model = "Random Forest",
  Accuracy = rf_accuracy,
  Sensitivity = rf_sensitivity,
  Specificity = rf_specificity,
  Balanced_Accuracy = rf_balanced_accuracy,
  AUC = rf_auc
)


cat("\nRandom Forest Performance:\n")

print(rf_metrics)


# ============================================================
# Combine model metrics
# ============================================================

classification_metrics <- bind_rows(
  logistic_metrics,
  rf_metrics
)

write.csv(
  classification_metrics,
  "results/phase13/classification_metrics.csv",
  row.names = FALSE
)


cat("\nModel Comparison:\n")

print(classification_metrics)


# ============================================================
# Save Random Forest predictions
# ============================================================

rf_predictions <- tibble(
  Actual = rf_actual,
  Predicted = rf_predicted,
  Probability = rf_probability
)

write.csv(
  rf_predictions,
  "results/phase13/random_forest_test_predictions.csv",
  row.names = FALSE
)



# ============================================================
# Step 5: Random Forest Variable Importance
# ============================================================

rf_importance <- importance(rf_model)

importance_table <- data.frame(
  Feature = rownames(rf_importance),
  MeanDecreaseAccuracy = rf_importance[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = rf_importance[, "MeanDecreaseGini"]
) %>%
  arrange(desc(MeanDecreaseGini))


cat("\nTop 20 Important Features:\n")

print(
  head(importance_table, 20)
)


write.csv(
  importance_table,
  "results/phase13/random_forest_variable_importance.csv",
  row.names = FALSE
)

# ============================================================
# Variable Importance Plot
# ============================================================

top_importance <- importance_table %>%
  filter(
    !Feature %in% c("Sex", "PC1", "PC2", "PC3")
  ) %>%
  slice_head(n = 20)


importance_plot <- ggplot(
  top_importance,
  aes(
    x = reorder(Feature, MeanDecreaseGini),
    y = MeanDecreaseGini
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 20 Random Forest Metabolite Predictors",
    x = "Metabolite",
    y = "Mean Decrease Gini"
  ) +
  theme_minimal()


ggsave(
  "plots/phase13/random_forest_variable_importance.png",
  importance_plot,
  width = 10,
  height = 7,
  dpi = 300
)
