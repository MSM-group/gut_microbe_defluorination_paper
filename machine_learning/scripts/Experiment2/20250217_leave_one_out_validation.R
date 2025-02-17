# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "caret", 
               "tidymodels", "ranger", "tree", "seqinr",
               "rsample", "randomForest", "readxl", "ggpubr")

# Read in the dataset
dat <- read_csv("data/20250209_full_set_452_training_seqs.csv")

# Set random seed
set.seed(1234) 

# Split into test and training data - random option
dat_split <- rsample::initial_split(dat[1:448,], prop = 0.8, strata = "truth")
dat_train <- rsample::training(dat_split)
dat_test0 <- rsample::testing(dat_split)
dat_test <- bind_rows(dat_test0, dat[448:452,])

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("nams", "truth")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "truth")]

# Dependent variable
y_train <- dat_train$truth
y_test <- dat_test$truth

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nams)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nams)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$nams)

# Optional tuning of random forest parameters
mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0), ncol(df_train))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("gini", "extratrees"),
                       min.node.size = 1)

# Weight the classes by their importance
wp118 <- 1/length(grep("WP_118709078", rownames(x_train)))
wp118 * length(grep("WP_118709078", rownames(x_train))) # weight corresponds to 1
wp178 <- 1/length(grep("WP_178618037", rownames(x_train)))
wp178 * length(grep("WP_178618037", rownames(x_train))) # weight corresponds to 1
case_wts <- case_when(grepl("WP_118709078", rownames(x_train)) ~ wp118,
                      grepl("WP_178618037", rownames(x_train)) ~ wp178,
                      TRUE ~ 1)
bind_cols(case_wts, rownames(x_train))

# Train a model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, 
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")
rf$bestTune

# Best model
rf_full <- ranger(as.factor(y_train) ~., data = form_train, num.trees = 1000, 
                  splitrule = "gini",
                  case.weights = case_wts,
                  mtry = 15, min.node.size = 1,
                  importance = "permutation", probability = TRUE)

# Training set accuracy
rf_full$predictions
rf_full$prediction.error # out-of-bag error
mod_train_pred <- as.factor(ifelse(rf_full$predictions[,1] >= 0.5, "defluor", "nondefluor"))
confusionMatrix(mod_train_pred, as.factor(dat_train$truth))
saveRDS(rf_full, "data/20250209_classification_random_forest.rds")

rf_pred <- predict(rf_full, data = form_test)
mod_test_pred <- as.factor(ifelse(rf_pred$predictions[,1] >= 0.5, "defluor", "nondefluor"))
mod_test_pred
cm_rf <- confusionMatrix(mod_test_pred, as.factor(dat_test$truth))
cm_rf

see_preds <- bind_cols(rownames(dat_test), y_test, mod_test_pred, rf_pred$predictions)
see_preds

# Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp
rf_imp

rf1 <- ggplot(rf_imp, top = 10) + 
  xlab("") +
  theme_classic()
rf1

table(as.numeric(gsub("residue", "", word(rownames(rf_imp$importance), sep = "_", 1))) > 196)
80/(154 + 80) # 34%
table(as.numeric(gsub("residue", "", word(rownames(rf_imp$importance), sep = "_", 1))) > 225) # PTLY
55/(154 + 80) # 24%
rf_imp

rf2 <- ggplot(rf_imp, top = 10) + 
  xlab("") +
  theme_classic() 
rf2

rf_roc_train <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                          predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                          ci = T)

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

see_preds <- bind_cols(rownames(dat_test), form_test$y_test, rf_pred)
write_csv(see_preds, "output/classification_model_predictions.csv")
cm_rf <- confusionMatrix(rf_pred, as.factor(y_test))
cm_rf

# ROC curve
rf$pred$obs
rf$pred$pred

rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                    predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                    plot = TRUE)
rf_roc

# AUC
plot(rf_roc, type = "s", 
     col = "#529DCF", xaxs = "i", yaxs="i",
     print.auc = TRUE, print.auc.x = 0.8, print.auc.y = 0.6)

rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                    predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                    ci = T)

rf_roc$sensitivities
rf_roc$specificities
pl1 <- plot(rf_roc, type = "s", 
            col = "#529DCF", xaxs = "i", yaxs="i",
            print.auc = TRUE, print.auc.x = 0.95 
            , print.auc.y = 0.8,
            xlim = c(1.1,-0.1), ylim = c( 0, 1.1))
pl1

rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                    predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                    ci = F)


