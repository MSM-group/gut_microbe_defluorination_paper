# Install packages
library("tidyverse")
library("readxl")
library("ggpubr")
library("data.table")
library("ranger")
library("caret")

# Read in the template
temp1 <- read_excel("data/machine_learning/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res <- read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                  col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res <- res[!is.na(res)]

# Pair
dat <- bind_cols(label = temp1, activity = res) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::filter(aa != "O")
table(dat$aa)
dat_list <- as.list(dat$aa)
names(dat_list) <- dat$label

source("scripts/convert_seq_5aap.R")
extract_feat_list <- lapply(1:length(dat_list), function(x) { convert_seq_5aap(dat_list[[x]]) })
n <- 10

extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)
extract_feat_df
colnames(extract_feat_df) <- c("polarity", "secondary_structure", "size", "codon_diversity", "charge")
extract_feat_df
feat_df <- dat %>%
  bind_cols(extract_feat_df) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(position = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::filter(!is.na(position))

# Try predicting the activity of 20%
set.seed(20240801)
# Split into test and training data 
# 80% training
# 20% test
dat <- feat_df 
dat_split <- rsample::initial_split(dat, prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)
nrow(dat_test)
nrow(dat_train)/nrow(dat) 

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("label", "activity")]
x_test <- dat_test[,!colnames(dat_test) %in% c("label", "activity")]

# Dependent variable
y_train <- dat_train$activity
y_test <- dat_test$activity
y_test # check there is a mix of your two variables

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$label)

# Optional tuning of random forest parameters
mtrys <- c(round(log2(ncol(df_train)), 0), round(ncol(df_train), 0))
mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = mtrys,
                       splitrule = c("variance", "extratrees"),
                       min.node.size = 1)

# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 3, # increase this to 3 when you run the code 
                           verboseIter = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
#saveRDS(rf, "data/Machine_learning/20240701_regression_random_forest.rds")
rf$finalModel$prediction.error # out-of-bag error

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp

#pdf("data/machine_learning/variable_importance_plot.pdf", width = 10)
ggplot(rf_imp, top = 7) + 
  xlab("") +
  theme_pubr(base_size = 14)
#dev.off()

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

# Predict for new data
Metrics::rmse(y_test, rf_pred) # 0.21
rf_df <- data.frame(cbind(rf_pred, y_test))

my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))

rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)

ggplot(rf_df_resid, aes(x = rf_pred, y = resid)) +
  geom_point(alpha = .3) + 
  geom_line(col = "red", y = 0.0,
            lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Predicted enzyme activity") +
  ylab("Residuals") 

ggplot(rf_df, aes(x = y_test, y = rf_pred)) +
  geom_point(alpha = .3) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed") +
  ylab("Predicted") 

