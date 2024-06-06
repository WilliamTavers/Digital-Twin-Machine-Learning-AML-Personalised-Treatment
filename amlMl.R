library(readr)
library(caret)
library(randomForest)
setwd("~/Desktop/amlDrug")
      
data <- read_csv("~/Desktop/amlDrug/AML_Drug_Response_Dataset2.csv")

data <- data[, !colnames(data) %in% "FLT3-ITD"]

# Scaling the features
preproc <- preProcess(data[, -ncol(data)], method = c("center", "scale"))
data_scaled <- predict(preproc, data[, -ncol(data)])
data_scaled$Drug_Response <- data$Drug_Response
data_scaled$Drug <- data$Drug


# Predicting drug response
set.seed(123)  
training_samples <- createDataPartition(data_scaled$Drug_Response, p = 0.8, list = FALSE)
train_data <- data_scaled[training_samples, ]
test_data <- data_scaled[-training_samples, ]

train_data$Drug_Response <- as.factor(train_data$Drug_Response)
test_data$Drug_Response <- as.factor(test_data$Drug_Response)

model <- randomForest(Drug_Response ~ ., data = train_data, ntree = 100)

predictions <- predict(model, test_data[, -which(colnames(test_data) == "Drug_Response")])

confusionMatrix(predictions, test_data$Drug_Response)

library(pROC)

probabilities <- predict(model, test_data[, -which(colnames(test_data) == "Drug_Response")], type = "prob")

roc_result <- roc(as.numeric(test_data$Drug_Response) == 1, probabilities[,2])
auc_result <- auc(roc_result)
print(auc_result)


roc_result <- roc(test_data$Drug_Response, probabilities)
auc(roc_result)
plot(roc_result, main="ROC Curve", col="#1c61b6")

legend("bottomright", legend=paste("AUC =", round(auc(roc_result), 2)), col="#1c61b6", lwd=10)
varImpPlot(model) 

# Extract feature importance
importance_data <- importance(model)

importance_df <- data.frame(
  Feature = rownames(importance_data),
  Importance = importance_data[, "%IncMSE"]  
)

ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_col(fill = "steelblue") +
  labs(title = "Feature Importance", x = "Features", y = "Increase in MSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  



# Convert Drug_Response to numeric
numeric_responses <- as.numeric(test_data$Drug_Response)

# Calculate Spearman Correlation
spearman_correlation <- cor(numeric_responses, probabilities[,2], method = "spearman")
print(paste("Spearman Correlation: ", spearman_correlation))

conf_matrix <- confusionMatrix(as.factor(predictions), test_data$Drug_Response)
precision <- conf_matrix$byClass['Pos Pred Value']
recall <- conf_matrix$byClass['Sensitivity']
f1_score <- 2 * precision * recall / (precision + recall)

print(paste("Precision: ", precision))
print(paste("Recall: ", recall))
print(paste("F1 Score: ", f1_score))

library(boot)

boot_data <- data.frame(Drug_Response = as.numeric(test_data$Drug_Response), predicted = probabilities[,2])

# Function to obtain Spearman correlation
spearman_function <- function(data, indices) {
  d <- data[indices,] 
  return(cor(d$Drug_Response, d$predicted, method = "spearman"))
}

# Running the bootstrap
results <- boot(data = boot_data, statistic = spearman_function, R = 1000)
boot_ci <- boot.ci(results, type = "perc")

print(boot_ci)



############### Drugs



data_scaled$Drug <- as.factor(data_scaled$Drug)  

# Drug Type
training_samples <- createDataPartition(data_scaled$Drug, p = 0.8, list = FALSE)
train_data <- data_scaled[training_samples, ]
test_data <- data_scaled[-training_samples, ]

train_data$Drug <- as.factor(train_data$Drug)
test_data$Drug <- as.factor(test_data$Drug)

model <- randomForest(Drug ~ ., data = train_data, ntree = 100)

importance <- varImp(model, scale = FALSE)
plot(importance)

# Predict the Drug
predictions <- predict(model, test_data[, -which(colnames(test_data) == "Drug")])
conf_matrix <- confusionMatrix(predictions, test_data$Drug)
print(conf_matrix)


probabilities <- predict(model, test_data[, -which(colnames(test_data) == "Drug")], type = "prob")

confidence_scores <- apply(probabilities, 1, max)

results <- data.frame(
  Actual = test_data$Drug,
  Predicted = predictions,
  Confidence = confidence_scores
)


print(head(results))

evaluate_at_threshold <- function(threshold) {
  filtered_results <- results[results$Confidence >= threshold, ]
  if (nrow(filtered_results) == 0) {
    return(c(Precision = NA, Recall = NA, F1_Score = NA, Count = 0))
  }
  
  conf_matrix <- confusionMatrix(filtered_results$Predicted, filtered_results$Actual)
  precision <- conf_matrix$byClass['Pos Pred Value']
  recall <- conf_matrix$byClass['Sensitivity']
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(c(Precision = precision, Recall = recall, F1_Score = f1_score, Count = nrow(filtered_results)))
}

thresholds <- seq(0, 1, by = 0.05)
evaluation_results <- t(sapply(thresholds, evaluate_at_threshold))
colnames(evaluation_results) <- c("Precision", "Recall", "F1_Score", "Count")
evaluation_results_df <- as.data.frame(evaluation_results)
evaluation_results_df$Threshold <- thresholds

print(evaluation_results_df)


optimal_threshold <- 0.75
optimal_results <- evaluate_at_threshold(optimal_threshold)

cat("Performance at Confidence Threshold 0.75:\n")
cat("Precision: ", optimal_results["Precision"], "\n")
cat("Recall: ", optimal_results["Recall"], "\n")
cat("F1-Score: ", optimal_results["F1_Score"], "\n")
cat("Number of Predictions: ", optimal_results["Count"], "\n")





# Predict probabilities for multiclass classification (which drug)
probabilities <- predict(model, test_data[, -which(colnames(test_data) == "Drug")], type = "prob")

accuracy <- sum(predict(model, test_data) == test_data$Drug) / nrow(test_data)
print(accuracy)

predictions <- predict(model, test_data)
confusionMatrix(predictions, test_data$Drug)
fourfoldplot(confusionMatrix$table, color=c("#6baed6", "#eff3ff", "#bdd7e7", "#3182bd"), conf.level=0)

varImpPlot(model)

roc_response <- roc(response = as.factor(test_data$Drug), predictor = as.numeric(predictions))
plot(roc_response)

accuracy <- sum(predictions == test_data$Drug) / nrow(test_data)
precision <- posPredValue(predictions, test_data$Drug, positive = "specific class")
recall <- sensitivity(predictions, test_data$Drug, positive = "specific class")
f1_score <- (2 * precision * recall) / (precision + recall)

data <- data.frame(
  metric = c("Accuracy", "Precision", "Recall", "F1-Score"),
  value = c(accuracy, precision, recall, f1_score)
)

ggplot(data, aes(x = metric, y = value, group = 1)) +
  geom_polygon(fill = "skyblue", colour = "blue") +
  geom_line(colour = "blue") +
  theme_minimal()


install.packages("pheatmap")
library(pheatmap)

conf_matrix <- table(predictions, test_data$Drug)

pheatmap(conf_matrix, color = colorRampPalette(c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd"))(n = 99))


library(pROC)
results <- lapply(levels(test_data$Drug), function(class) {
  roc(as.numeric(test_data$Drug == class), probabilities[, class])
})
# Plot ROC curve for each class
lapply(results, plot)


library(caret)
stats <- lapply(levels(test_data$Drug), function(class) {
  actual <- test_data$Drug == class
  predicted <- predictions == class
  data.frame(
    Accuracy = sum(predicted == actual) / length(actual),
    Precision = posPredValue(predicted, actual, positive = TRUE),
    Recall = sensitivity(predicted, actual, positive = TRUE),
    F1_Score = f1Score(predicted, actual)
  )
})

stats_df <- do.call(rbind, stats)
rownames(stats_df) <- levels(test_data$Drug)
print(stats_df)


library(ggplot2)
data <- melt(stats_df, varnames = c("Metric", "Drug"))
ggplot(data, aes(x = Metric, y = value, group = Drug, color = Drug)) +
  geom_line() +
  geom_point() +
  coord_polar() +
  theme_minimal()



#xgboost
library(caret)
library(xgboost)

control <- trainControl(method="cv", number=5, search="random")

tunegrid <- expand.grid(
  nrounds = c(50, 100, 200),
  max_depth = c(4, 6, 8),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0, 1, 5),
  colsample_bytree = c(0.5, 0.7, 1.0),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.5, 0.7, 1.0)
)

xgb_caret_model <- train(
  Drug ~ ., data = train_data,
  method = "xgbTree",
  trControl = control,
  tuneGrid = tunegrid,
  verbose = FALSE
)

best_model <- xgb_caret_model$bestTune

xgb_model <- xgboost(
  data = train_matrix,
  max_depth = best_model$max_depth,
  eta = best_model$eta,
  nrounds = best_model$nrounds,
  gamma = best_model$gamma,
  colsample_bytree = best_model$colsample_bytree,
  min_child_weight = best_model$min_child_weight,
  subsample = best_model$subsample,
  objective = "multi:softprob",
  num_class = length(unique(train_data$Drug)),
  verbose = 0
)

xgb_probabilities <- predict(xgb_model, test_matrix)
xgb_probabilities <- matrix(xgb_probabilities, nrow = num_classes, ncol = length(test_data$Drug), byrow = TRUE)
xgb_predictions <- max.col(t(xgb_probabilities)) - 1

xgb_predictions <- factor(xgb_predictions, levels = 0:(num_classes-1), labels = levels(train_data$Drug))

if (length(xgb_predictions) == length(test_data$Drug)) {
  xgb_conf_matrix <- confusionMatrix(xgb_predictions, test_data$Drug)
  print(xgb_conf_matrix)
} else {
  print("Error: Lengths of predictions and actual labels do not match.")
}

xgb_confidence_scores <- apply(xgb_probabilities, 2, max)

xgb_results <- data.frame(
  Actual = test_data$Drug,
  Predicted = xgb_predictions,
  Confidence = xgb_confidence_scores
)

evaluate_at_threshold <- function(threshold) {
  filtered_results <- xgb_results[xgb_results$Confidence >= threshold, ]
  if (nrow(filtered_results) == 0) {
    return(c(Precision = NA, Recall = NA, F1_Score = NA, Count = 0))
  }
  
  conf_matrix <- confusionMatrix(filtered_results$Predicted, filtered_results$Actual)
  precision <- conf_matrix$byClass['Pos Pred Value']
  recall <- conf_matrix$byClass['Sensitivity']
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(c(Precision = precision, Recall = recall, F1_Score = f1_score, Count = nrow(filtered_results)))
}

thresholds <- seq(0, 1, by = 0.05)
evaluation_results <- t(sapply(thresholds, evaluate_at_threshold))
colnames(evaluation_results) <- c("Precision", "Recall", "F1_Score", "Count")
evaluation_results_df <- as.data.frame(evaluation_results)
evaluation_results_df$Threshold <- thresholds

print(evaluation_results_df)

optimal_threshold <- 0.75
optimal_results <- evaluate_at_threshold(optimal_threshold)

cat("Performance at Confidence Threshold 0.75:\n")
cat("Precision: ", optimal_results["Precision"], "\n")
cat("Recall: ", optimal_results["Recall"], "\n")
cat("F1-Score: ", optimal_results["F1_Score"], "\n")
cat("Number of Predictions: ", optimal_results["Count"], "\n")
