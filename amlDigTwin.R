# Set working directory
setwd("~/Desktop/newTestR")

# Loading necessary libraries
install.packages("caret")
install.packages("randomForest")
install.packages("e1071")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")
install.packages("xgboost")
library(caret)
library(randomForest)
library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(xgboost)
library(pROC)

# Load data in R
gene_expression <- read.table("~/Desktop/amlFolder/aml/exp", header = TRUE, sep = " ")
dna_methylation <- read.table("~/Desktop/amlFolder/aml/methy", header = TRUE, sep = " ")
mirna_expression <- read.table("~/Desktop/amlFolder/aml/mirna", header = TRUE, sep = " ")
patient_survival <- read.table("~/Desktop/amlFolder/aml/survival", header = TRUE, sep = "\t")

# fixing as table 
gene_dt <- as.data.table(gene_expression)
methy_dt <- as.data.table(dna_methylation)
mirna_dt <- as.data.table(mirna_expression)

convert_to_long <- function(df, measure_variable_name) {
  # Create a variable for Gene (row names) and add it to the dataframe
  df$Gene <- rownames(df)
  # Convert the dataframe to long format
  df_long <- melt(df, id.vars = "Gene", variable.name = "PatientID", value.name = measure_variable_name)
  # Adjust the PatientID to ensure it matches the standard TCGA ID format (TCGA-XXXX-XXXX)
  df_long$PatientID <- toupper(gsub("^(TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}).*", "\\1", df_long$PatientID))
  return(df_long)
}

# Reapply the function to each dataset
gene_long <- convert_to_long(gene_expression, "Expression")
methy_long <- convert_to_long(dna_methylation, "Methylation")
mirna_long <- convert_to_long(mirna_expression, "miRNA_Expression")

# Standardize PatientID in patient_survival
patient_survival$PatientID <- gsub("-", ".", toupper(patient_survival$PatientID))


#merge gene expression data with survival data
common_ids <- intersect(gene_long$PatientID, patient_survival$PatientID)
gene_survival <-merge(gene_long, patient_survival, by  = "PatientID")

#merge methy data with survival data
common_ids <- intersect(methy_long$PatientID, patient_survival$PatientID)
methy_survival <-merge(methy_long, patient_survival, by  = "PatientID")

#merge mirna data with survival data
common_ids <- intersect(mirna_long$PatientID, patient_survival$PatientID)
mirna_survival <-merge(mirna_long, patient_survival, by  = "PatientID")



# Check for missing values
sum(is.na(gene_survival$Survival))
sum(is.na(gene_survival$Death))
sum(is.na(methy_survival$Survival))
sum(is.na(methy_survival$Death))
sum(is.na(mirna_survival$Survival))
sum(is.na(mirna_survival$Death))


# Remove rows with missing Survival or Death values
gene_survival_clean <- gene_survival[!is.na(gene_survival$Survival) & !is.na(gene_survival$Death), ]
mirna_survival_clean <- mirna_survival[!is.na(mirna_survival$Survival) & !is.na(mirna_survival$Death), ]
methy_survival_clean <- methy_survival[!is.na(methy_survival$Survival) & !is.na(methy_survival$Death), ]
head(gene_survival_clean)
tail(gene_survival_clean)
# Verify the dimensions after removing missing values
print(dim(gene_survival_clean))
print(dim(mirna_survival_clean))
print(dim(methy_survival_clean))


# need to download aml clinical patient
clinical_data <- read.table("~/Desktop/amlFolder/clinical/aml", header = TRUE, sep = "\t")

# Replace dashes with dots in sampleID
clinical_data$sampleID <- gsub("-", ".", clinical_data$sampleID)
gene_clinical <- merge(clinical_data, gene_survival_clean, by.x = "sampleID", by.y = "PatientID", all = TRUE)
head(gene_clinical)
#
# Merge the datasets on PatientID, retaining all gene columns
# omics_data <- merge(gene_survival_clean, methy_survival_clean, by = "PatientID", all = TRUE, suffixes = c("_gene", "_meth"))

# Combine Survival columns
# omics_data$Survival <- ifelse(is.na(omics_data$Survival_gene), omics_data$Survival_meth, omics_data$Survival_gene)

# Combine Death columns
# omics_data$Death <- ifelse(is.na(omics_data$Death_gene), omics_data$Death_meth, omics_data$Death_gene)

# Drop the old Survival and Death columns
# omics_data <- omics_data[, !colnames(omics_data) %in% c("Survival_gene", "Survival_meth", "Death_gene", "Death_meth")]

# Check the combined result
# head(omics_data)
# str(omics_data)



# agg
# gene_survival_agg <- gene_survival_clean %>%
  # group_by(PatientID) %>%
  # summarize(Expression = mean(Expression, na.rm = TRUE), 
            # Gene = Gene,
            # Survival = mean(Survival, na.rm = TRUE), 
            # Death = mean(Death, na.rm = TRUE))

# head(gene_survival_agg)

# methy_survival_agg <- methy_survival_clean %>%
 # group_by(PatientID) %>%
 # summarize(Methylation = mean(Methylation, na.rm = TRUE), 
    #        Gene = Gene,
 #           Survival = mean(Survival, na.rm = TRUE), 
#            Death = mean(Death, na.rm = TRUE))

# head(methy_survival_agg)

# omics_data <- merge(gene_survival_agg, methy_survival_agg, by = "PatientID", all = TRUE)
# head(omics_data)



# Merge the omics datasets
# omics_data <- merge(gene_survival_clean, mirna_survival_clean, by = "PatientID", all = TRUE)
# omics_data <- merge(omics_data, mirna_survival_clean, by = "PatientID", all = TRUE)

# Check the structure
# head(omics_data)
# str(omics_data)



# head(gene_clinical$gender)


# Convert data frames to data tables for efficient handling
dt_gene <- as.data.table(gene_survival_clean)
dt_methy <- as.data.table(methy_survival_clean)
dt_mirna <- as.data.table(mirna_survival_clean)

# Remove exact duplicates across all columns
dt_gene <- unique(dt_gene)
dt_methy <- unique(dt_methy)
dt_mirna <- unique(dt_mirna)

head(dt_gene)
head(dt_methy)
# Aggregating by taking the mean (if necessary)
# If values are the same for duplicates, the mean will not change the result
dt_gene <-  dt_gene[, lapply(.SD, mean), by = PatientID, .SDcols = -"PatientID"]
dt_methy <- dt_methy[, .(Methylation = Methylation), by = .(PatientID, Gene, Survival, Death)]
dt_mirna <- dt_mirna[, .(miRNA_Expression = miRNA_Expression), by = .(PatientID, Gene, Survival, Death)]
head(dt_gene)
head(dt_methy)
head(dt_mirna)

# Merge on PatientID, align by 'Death' if it's included in all tables
# Using .EACHI in data.table (only works if you're using data.table syntax)
omics_data <- dt_gene[dt_methy, on = "PatientID", allow.cartesian = TRUE, by = .EACHI]
omics_data <- omics_data[dt_mirna, on = "PatientID", allow.cartesian = TRUE, by = .EACHI]


# Optionally rename columns to clarify which omics data they represent
setnames(omics_data, old = c("Gene", "Gene_gene", "Gene_meth"), 
         new = c("Gene_mirna", "Gene_expression", "Gene_methylation"))

# Check the structure and the head of the final data table
str(omics_data)
head(omics_data)
tail(omics_data)

# (1) predicting survival

sapply(list(gene_survival_clean, mirna_survival_clean, methy_survival_clean), function(x) format(object.size(x), units = "MB"))

# Merge datasets on PatientID
omics_data <- Reduce(function(x, y) merge(x, y, by = "PatientID", all = TRUE), list(gene_survival_clean, mirna_survival_clean, methy_survival_clean))

# Check for missing values and impute or remove
omics_data <- na.omit(omics_data)  # Simple way to handle missing values

# Select predictors and the outcome variable
predictors <- setdiff(names(omics_data), c("PatientID", "Survival", "Death"))
outcome <- "Death"

# Split data into training and testing sets
set.seed(123)  # for reproducibility
training_rows <- createDataPartition(omics_data[[outcome]], p = 0.8, list = FALSE)
train_data <- omics_data[training_rows, ]
test_data <- omics_data[-training_rows, ]









library(Matrix) 

library(digest)
str(dt_gene)
# Correcting and verifying the hashing function
hash_function <- function(x, num_features) {
  # Convert string to integer vector
  int_vector <- utf8ToInt(x)
  
  # Sum the integer vector to get a single integer
  total <- sum(int_vector, na.rm = TRUE)
  
  # Apply modulo operation to reduce dimensionality
  return((total %% num_features) + 1)
}

# Number of features (buckets for hashing)
num_features <- 100000000  # You can adjust this number based on your needs
head(dt_methy)
# Apply hashing to the 'Gene' column
dt_gene[, hashed_gene := sapply(Gene, hash_function, num_features = num_features)]
head(dt_mirna)
head(dt_gene)
tail(dt_gene)

set.seed(123)  # for reproducibility

# Prepare the data (select only necessary columns)
data_model <- dt_gene[, .(Expression, Gene, Survival, Death)]
data_mirna_model <- dt_mirna[, .(Gene, Death, miRNA_Expression, Death)]

#
# Split data into training and testing sets
trainIndex <- createDataPartition(data_model$Death, p = 0.8, list = FALSE)
train_data <- data_model[trainIndex, ]
test_data <- data_model[-trainIndex, ]

  # mirna
mirna_trainIndex <- createDataPartition(data_mirna_model$Death, p = 0.8, list = FALSE)
mirna_train_data <- data_mirna_model[mirna_trainIndex, ]
mirna_test_data <- data_mirna_model[-mirna_trainIndex, ]
tail(dt_mirna)

# Convert Death to a factor if it's not already
train_data$Death <- as.factor(train_data$Death)
test_data$Death <- as.factor(test_data$Death)
mirna_train_data$Death <- as.factor(mirna_train_data$Death)
mirna_test_data$Death <- as.factor(mirna_test_data$Death)

# Try running the model again with fewer trees to check for memory issues
rf_model <- randomForest(Death ~ ., data = train_data, ntree = 25, importance = TRUE)
mirna_rf_model <- randomForest(Death ~ ., data = mirna_train_data, ntree = 110, importance = TRUE)
# Check the model
print(rf_model)
print(mirna_rf_model)
# Predict on the test set
predictions <- predict(rf_model, test_data)
mirna_predictions <- predict(mirna_rf_model, mirna_test_data)

# Calculate confusion matrix to see the accuracy
conf_mat <- confusionMatrix(predictions, test_data$Death)
print(conf_mat)

mirna_conf_mat <- confusionMatrix(mirna_predictions, mirna_test_data$Death)
print(mirna_conf_mat)


# Calculate ROC AUC

roc_result <- roc(test_data$Death, as.numeric(predictions))
mirna_roc_result <- roc(mirna_test_data$Death, as.numeric(mirna_predictions))

auc_value <- auc(roc_result)
mirna_auc_value <-auc(mirna_roc_result)

print(paste("AUC value:", auc_value))
print(paste("mirna AUC value:", mirna_auc_value))

varImpPlot(rf_model)
varImpPlot(mirna_rf_model)

# Plot ROC curves
plot(roc_result, main="ROC Curve for Gene Expression Model", col="#1c61b6", print.auc=TRUE, print.auc.x=0.6)
plot(mirna_roc_result, main="ROC Curve for miRNA Expression Model", col="#ff9900", print.auc=TRUE, print.auc.x=0.6)


# (2) predicting days till death


# List of relevant columns based on the provided column names
gene_columns <- grep("Gene|Expression", colnames(gene_clinical), value = TRUE)
demographics_columns <- c("age_at_initial_pathologic_diagnosis", "gender")
clinical_columns <- c("days_to_death", "disease_detection_molecular_analysis_method_type", 
                      "acute_myeloid_leukemia_calgb_cytogenetics_risk_category", "vital_status")

# Combine all relevant columns
relevant_columns <- c("sampleID", gene_columns, demographics_columns, clinical_columns)
print(relevant_columns)

# Check if all relevant columns exist in the dataset
missing_columns <- setdiff(relevant_columns, colnames(gene_clinical))
if (length(missing_columns) > 0) {
  cat("Missing columns:", missing_columns, "\n")
} else {
  # Select relevant columns from the dataset
  selected_data <- gene_clinical[, relevant_columns]
}

# Print the structure of the selected data
if (exists("selected_data")) {
  str(selected_data)
}


if (exists("selected_data")) {
  # Remove rows with missing target values
  selected_data <- na.omit(selected_data)
  
  # Split the dataset into training and testing sets
  set.seed(123) # For reproducibility
  trainIndex <- createDataPartition(selected_data$days_to_death, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  train_data <- selected_data[trainIndex, ]
  test_data <- selected_data[-trainIndex, ]
}

if (exists("train_data") && exists("test_data")) {
  # Train a Random Forest model
  rf_model <- randomForest(days_to_death ~ ., data = train_data, ntree = 64)
  
  # Predict on the test set
  rf_predictions <- predict(rf_model, test_data)
  
  # Evaluate the model
  evaluation <- postResample(rf_predictions, test_data$days_to_death)
  print(evaluation)
}



# Calculate evaluation metrics
mae <- mean(abs(rf_predictions - test_data$days_to_death))
mse <- mean((rf_predictions - test_data$days_to_death)^2)
r2 <- cor(rf_predictions, test_data$days_to_death)^2

cat("Mean Absolute Error:", mae, "\n")
cat("Mean Squared Error:", mse, "\n")
cat("R-squared:", r2, "\n")

# Plot variable importance
varImpPlot(rf_model)
importance(rf_model)

# Plot ROC curves
plot(roc_result, main="ROC Curve for Gene Expression Model", col="#1c61b6", print.auc=TRUE, print.auc.x=0.6)
plot(mirna_roc_result, main="ROC Curve for miRNA Expression Model", col="#ff9900", print.auc=TRUE, print.auc.x=0.6)







# treatment response
all_data_mirna <- merge(dt_mirna, clinical_data, by.x="PatientID", by.y="sampleID", all.x=TRUE)

head(all_data_mirna)
all_data_mirna <- unique(all_data_mirna)



# Example of creating a treatment response feature
all_data_mirna$high_risk_treatment_response <- ifelse(all_data_mirna$FISH_test_component_percentage_value > 50 &
                                                       all_data_mirna$acute_myeloid_leukemia_calgb_cytogenetics_risk_category == "High",
                                                     "High Risk", "Low Risk")

# Load necessary library
install.packages("missForest")
library(missForest)
# Convert FISH_test_component_percentage_value to numeric
final_data$FISH_test_component_percentage_value <- as.numeric(as.character(final_data$FISH_test_component_percentage_value))

# Check for any conversion failures which result in NAs
sum(is.na(final_data$FISH_test_component_percentage_value))

# If there are many NAs after conversion, you may need to decide how to handle these
# For instance, you can fill them with the median or mean
median_value <- median(final_data$FISH_test_component_percentage_value, na.rm = TRUE)
final_data$FISH_test_component_percentage_value[is.na(final_data$FISH_test_component_percentage_value)] <- median_value

# Impute missing values using missForest, which can handle mixed types (numerical and categorical)
final_data_imputed <- missForest(final_data)

# Train the model on the imputed data
rf_treatment_response <- randomForest(high_risk_treatment_response ~ ., data=final_data_imputed$ximp, ntree=10)
-
head(dt_mirna)
# Prepare data for modeling
final_data <- all_data_mirna[, .(miRNA_Expression, hashed_gene, FISH_test_component_percentage_value, high_risk_treatment_response, Survival, Death)]

# Train models
rf_treatment_response <- randomForest(high_risk_treatment_response ~ ., data=final_data, ntree=10)
rf_prognosis <- randomForest(Death ~ ., data=final_data, ntree=50)











# Convert 'high_risk_treatment_response' to factor explicitly and verify the conversion
final_data$high_risk_treatment_response <- as.factor(final_data$high_risk_treatment_response)

# Check the structure immediately after conversion to confirm the type
str(final_data$high_risk_treatment_response)

# If still showing as character, investigate if any steps prior might be changing its type or check the environment issues

final_data_imputed <- missForest(final_data)



# Proceed with your model training
rf_treatment_response <- randomForest(high_risk_treatment_response ~ ., data=final_data_imputed$ximp, ntree=10)

# Evaluate the model
print(summary(rf_treatment_response))

# Check the levels and their counts in the high_risk_treatment_response variable
table(final_data_imputed$ximp$high_risk_treatment_response)

# If the dataset has been imputed, check the original data as well
table(final_data$high_risk_treatment_response)



# Example: redefine "high_risk_treatment_response" with a different threshold or criteria
final_data$high_risk_treatment_response <- ifelse(
  final_data$FISH_test_component_percentage_value > 50,  # Adjusted from a previous strict criteria
  "High Risk",
  "Low Risk"
)
head(final_data)
# Split data into training and testing sets
set.seed(123)  # for reproducibility
trainIndex <- createDataPartition(final_data$high_risk_treatment_response, p = 0.8, list = FALSE)
train_data <- final_data[trainIndex, ]
test_data <- final_data[-trainIndex, ]

## Ensure the response variable is a factor
train_data$high_risk_treatment_response <- as.factor(train_data$high_risk_treatment_response)
test_data$high_risk_treatment_response <- as.factor(test_data$high_risk_treatment_response)

# Now, try training the Random Forest again
rf_treatment_response <- randomForest(high_risk_treatment_response ~ ., data = train_data, ntree = 10)
print(rf_treatment_response)

# Proceed with prediction and evaluation
predictions <- predict(rf_treatment_response, test_data)
conf_mat <- confusionMatrix(predictions, test_data$high_risk_treatment_response)
print(conf_mat)

# Calculate ROC and AUC
prob_predictions <- predict(rf_treatment_response, test_data, type = "prob")
roc_result <- roc(test_data$high_risk_treatment_response, prob_predictions[, "High Risk"])
auc_value <- auc(roc_result)
print(paste("AUC value:", auc_value))


# Variable importance plot
png("variable_importance_plot.png", width=800, height=600)
varImpPlot(rf_treatment_response)
dev.off()  # Save and close the plotting device


set.seed(123)  # for reproducibility

# Define training control
train_control <- trainControl(
  method = "cv",  # k-fold cross-validation
  number = 10,    # number of folds
  savePredictions = "final",
  summaryFunction = twoClassSummary,
  classProbs = TRUE  # to calculate AUC later
)

# Train the model
model <- train(
  high_risk_treatment_response ~ ., 
  data = final_data, 
  method = "rf",  # using random forest
  trControl = train_control,
  metric = "ROC"
)

# Check results
print(model)




# Downsample the data
set.seed(123)
train_indices_downsampled <- sample(1:nrow(train_data), size = 0.1 * nrow(train_data))
train_data_downsampled <- train_data[train_indices_downsampled, ]

# Prepare data for xgboost
train_matrix <- xgb.DMatrix(data = as.matrix(train_data_downsampled$Expression), label = train_data_downsampled$Survival)
test_matrix <- xgb.DMatrix(data = as.matrix(test_data$Expression), label = test_data$Survival)

# Define parameters
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse"
)

# Train the model
xgb_model <- xgb.train(params = params, data = train_matrix, nrounds = 500)

# Predict on the test set
predictions <- predict(xgb_model, test_matrix)

# Evaluate the model
actual <- test_data$Survival
results <- data.frame(Actual = actual, Predicted = predictions)
head(results)

# Calculate evaluation metrics
mse <- mean((results$Actual - results$Predicted)^2)
rmse <- sqrt(mse)
mae <- mean(abs(results$Actual - results$Predicted))

cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")




# Perform cross-validation
cv <- xgb.cv(
  params = params,
  data = train_matrix,
  nrounds = 500,
  nfold = 5,
  showsd = TRUE,
  stratified = TRUE,
  print.every.n = 10,
  early_stopping_rounds = 20,
  maximize = FALSE
)

# Get the best number of rounds from cross-validation
best_nrounds <- cv$best_iteration

# Train the model with the best number of rounds
xgb_model <- xgb.train(
  params = params,
  data = train_matrix,
  nrounds = best_nrounds
)

# Predict on the test set
predictions <- predict(xgb_model, test_matrix)

# Evaluate the model
actual <- test_data$Survival
results <- data.frame(Actual = actual, Predicted = predictions)
mse <- mean((results$Actual - results$Predicted)^2)
rmse <- sqrt(mse)
mae <- mean(abs(results$Actual - results$Predicted))

cat("Mean Squared Error (MSE):", mse, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")















