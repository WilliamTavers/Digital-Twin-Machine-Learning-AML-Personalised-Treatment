setwd("~/Desktop/newTestR")

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

gene_expression <- read.table("~/Desktop/amlFolder/aml/exp", header = TRUE, sep = " ")
dna_methylation <- read.table("~/Desktop/amlFolder/aml/methy", header = TRUE, sep = " ")
mirna_expression <- read.table("~/Desktop/amlFolder/aml/mirna", header = TRUE, sep = " ")
patient_survival <- read.table("~/Desktop/amlFolder/aml/survival", header = TRUE, sep = "\t")

gene_dt <- as.data.table(gene_expression)
methy_dt <- as.data.table(dna_methylation)
mirna_dt <- as.data.table(mirna_expression)

convert_to_long <- function(df, measure_variable_name) {
  df$Gene <- rownames(df)
  df_long <- melt(df, id.vars = "Gene", variable.name = "PatientID", value.name = measure_variable_name)
  df_long$PatientID <- toupper(gsub("^(TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}).*", "\\1", df_long$PatientID))
  return(df_long)
}

gene_long <- convert_to_long(gene_expression, "Expression")
methy_long <- convert_to_long(dna_methylation, "Methylation")
mirna_long <- convert_to_long(mirna_expression, "miRNA_Expression")

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

print(dim(gene_survival_clean))
print(dim(mirna_survival_clean))
print(dim(methy_survival_clean))

clinical_data <- read.table("~/Desktop/amlFolder/clinical/aml", header = TRUE, sep = "\t")

clinical_data$sampleID <- gsub("-", ".", clinical_data$sampleID)
gene_clinical <- merge(clinical_data, gene_survival_clean, by.x = "sampleID", by.y = "PatientID", all = TRUE)
head(gene_clinical)

# Convert data frames to data tables for efficient handling
dt_gene <- as.data.table(gene_survival_clean)
dt_methy <- as.data.table(methy_survival_clean)
dt_mirna <- as.data.table(mirna_survival_clean)

dt_gene <- unique(dt_gene)
dt_methy <- unique(dt_methy)
dt_mirna <- unique(dt_mirna)

head(dt_gene)
head(dt_methy)

# Aggregating by taking the mean (if necessary)
dt_gene <-  dt_gene[, lapply(.SD, mean), by = PatientID, .SDcols = -"PatientID"]
dt_methy <- dt_methy[, .(Methylation = Methylation), by = .(PatientID, Gene, Survival, Death)]
dt_mirna <- dt_mirna[, .(miRNA_Expression = miRNA_Expression), by = .(PatientID, Gene, Survival, Death)]
head(dt_gene)
head(dt_methy)
head(dt_mirna)

omics_data <- dt_gene[dt_methy, on = "PatientID", allow.cartesian = TRUE, by = .EACHI]
omics_data <- omics_data[dt_mirna, on = "PatientID", allow.cartesian = TRUE, by = .EACHI]


setnames(omics_data, old = c("Gene", "Gene_gene", "Gene_meth"), 
         new = c("Gene_mirna", "Gene_expression", "Gene_methylation"))

str(omics_data)
head(omics_data)
tail(omics_data)

# (1) predicting survival

sapply(list(gene_survival_clean, mirna_survival_clean, methy_survival_clean), function(x) format(object.size(x), units = "MB"))

omics_data <- Reduce(function(x, y) merge(x, y, by = "PatientID", all = TRUE), list(gene_survival_clean, mirna_survival_clean, methy_survival_clean))

omics_data <- na.omit(omics_data)  

predictors <- setdiff(names(omics_data), c("PatientID", "Survival", "Death"))
outcome <- "Death"

set.seed(123) 
training_rows <- createDataPartition(omics_data[[outcome]], p = 0.8, list = FALSE)
train_data <- omics_data[training_rows, ]
test_data <- omics_data[-training_rows, ]


library(Matrix) 

library(digest)
str(dt_gene)
hash_function <- function(x, num_features) {
  int_vector <- utf8ToInt(x)
  total <- sum(int_vector, na.rm = TRUE)
  return((total %% num_features) + 1)
}

num_features <- 100000000  
head(dt_methy)
dt_gene[, hashed_gene := sapply(Gene, hash_function, num_features = num_features)]
head(dt_mirna)
head(dt_gene)
tail(dt_gene)

data_model <- dt_gene[, .(Expression, Gene, Survival, Death)]
data_mirna_model <- dt_mirna[, .(Gene, Death, miRNA_Expression, Death)]


trainIndex <- createDataPartition(data_model$Death, p = 0.8, list = FALSE)
train_data <- data_model[trainIndex, ]
test_data <- data_model[-trainIndex, ]

mirna_trainIndex <- createDataPartition(data_mirna_model$Death, p = 0.8, list = FALSE)
mirna_train_data <- data_mirna_model[mirna_trainIndex, ]
mirna_test_data <- data_mirna_model[-mirna_trainIndex, ]
tail(dt_mirna)

train_data$Death <- as.factor(train_data$Death)
test_data$Death <- as.factor(test_data$Death)
mirna_train_data$Death <- as.factor(mirna_train_data$Death)
mirna_test_data$Death <- as.factor(mirna_test_data$Death)

rf_model <- randomForest(Death ~ ., data = train_data, ntree = 25, importance = TRUE)
mirna_rf_model <- randomForest(Death ~ ., data = mirna_train_data, ntree = 110, importance = TRUE)

print(rf_model)
print(mirna_rf_model)

predictions <- predict(rf_model, test_data)
mirna_predictions <- predict(mirna_rf_model, mirna_test_data)

conf_mat <- confusionMatrix(predictions, test_data$Death)
print(conf_mat)

mirna_conf_mat <- confusionMatrix(mirna_predictions, mirna_test_data$Death)
print(mirna_conf_mat)


roc_result <- roc(test_data$Death, as.numeric(predictions))
mirna_roc_result <- roc(mirna_test_data$Death, as.numeric(mirna_predictions))

auc_value <- auc(roc_result)
mirna_auc_value <-auc(mirna_roc_result)

print(paste("AUC value:", auc_value))
print(paste("mirna AUC value:", mirna_auc_value))

varImpPlot(rf_model)
varImpPlot(mirna_rf_model)

plot(roc_result, main="ROC Curve for Gene Expression Model", col="#1c61b6", print.auc=TRUE, print.auc.x=0.6)
plot(mirna_roc_result, main="ROC Curve for miRNA Expression Model", col="#ff9900", print.auc=TRUE, print.auc.x=0.6)


# (2) predicting days till death


gene_columns <- grep("Gene|Expression", colnames(gene_clinical), value = TRUE)
demographics_columns <- c("age_at_initial_pathologic_diagnosis", "gender")
clinical_columns <- c("days_to_death", "disease_detection_molecular_analysis_method_type", 
                      "acute_myeloid_leukemia_calgb_cytogenetics_risk_category", "vital_status")

relevant_columns <- c("sampleID", gene_columns, demographics_columns, clinical_columns)
print(relevant_columns)

missing_columns <- setdiff(relevant_columns, colnames(gene_clinical))
if (length(missing_columns) > 0) {
  cat("Missing columns:", missing_columns, "\n")
} else {
  selected_data <- gene_clinical[, relevant_columns]
}

if (exists("selected_data")) {
  str(selected_data)
}


if (exists("selected_data")) {
  selected_data <- na.omit(selected_data)
  
  set.seed(123) 
  trainIndex <- createDataPartition(selected_data$days_to_death, p = .8, 
                                    list = FALSE, 
                                    times = 1)
  train_data <- selected_data[trainIndex, ]
  test_data <- selected_data[-trainIndex, ]
}

if (exists("train_data") && exists("test_data")) {
  rf_model <- randomForest(days_to_death ~ ., data = train_data, ntree = 64)
  
  rf_predictions <- predict(rf_model, test_data)
  
  evaluation <- postResample(rf_predictions, test_data$days_to_death)
  print(evaluation)
}



mae <- mean(abs(rf_predictions - test_data$days_to_death))
mse <- mean((rf_predictions - test_data$days_to_death)^2)
r2 <- cor(rf_predictions, test_data$days_to_death)^2

cat("Mean Absolute Error:", mae, "\n")
cat("Mean Squared Error:", mse, "\n")
cat("R-squared:", r2, "\n")

varImpPlot(rf_model)
importance(rf_model)

plot(roc_result, main="ROC Curve for Gene Expression Model", col="#1c61b6", print.auc=TRUE, print.auc.x=0.6)
plot(mirna_roc_result, main="ROC Curve for miRNA Expression Model", col="#ff9900", print.auc=TRUE, print.auc.x=0.6)
