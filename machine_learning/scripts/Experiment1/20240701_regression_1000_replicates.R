# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Template
temp1 <- read_excel("data/Experiment1/template_exp1.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res <-  read_excel('data/Experiment1/20240427_rep1_2_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res <- as.numeric(res[!is.na(res)])

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
dat <- feat_df

# Try predicting the activity of 20%
set.seed(20240711)
rf_models <- vector(mode = "list",
                    length = 1000)
dat$label
# Split randomly into 10 test/training sets
for(i in 1:1000) {
  # Split into test and training data
  dat_split <- initial_split(dat)
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-label, -activity)
  x_test <- dat_test %>%
    dplyr::select(-label, -activity)
  
  y_train <- dat_train$activity
  y_test <- dat_test$activity
  
  df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$id)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)
  
  ## Train model using set parameters (determined by cross-validation)
  rf_1 <- ranger(y_train ~., data = form_train, 
                 num.trees = 1000, 
                 splitrule = "extratrees",
                 mtry = 7, 
                 min.node.size = 1,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$train_rmse <- sqrt(rf_1$prediction.error) # Classifcation accuracy 1 - OOB error
  
  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$test_rmse <- Metrics::rmse(rf_1_pred$predictions, y_test)
}

mean_training_vec <- unlist(sapply(1:length(rf_models), function(x) { rf_models[[x]]$train_rmse })) 
traindf <- tibble(mean_training_vec)
traindf


mean_testing_vec <- unlist(sapply(1:length(rf_models), function(x) { rf_models[[x]]$test_rmse}))
testdf <- tibble(mean_testing_vec)


p1 <- ggplot(traindf) +
  geom_histogram(aes(x = mean_training_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean training set accuracy") +
  ylab("Frequency")
p1

p2 <- ggplot(testdf) +
  geom_histogram(aes( x = mean_testing_vec, y = ..density..), alpha = 0.7,
                 binwidth = 0.005, fill = "gray80", color = "black") +
  theme_pubr(base_size = 17) +
  xlab("Mean testing set accuracy") +
  ylab("Frequency")
p2

# Combine into a single histogram
traindf$split <- "training"
colnames(traindf) <- c("score", "split")
trainsumm <- data.frame(tibble(t(summary(traindf$score))), stringsAsFactors = F)
trainsumm

testdf$split <- "testing"
colnames(testdf) <- c("score", "split")
summary(testdf$score)
testsumm <- data.frame(tibble(t(summary(testdf$score))), stringsAsFactors = F)
testsumm

allsumm <- data.frame(c(trainsumm,testsumm)) 

alldf <- traindf %>%
  bind_rows(testdf)
head(alldf)
alldf$split <- as.factor(alldf$split)

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 
colMeans(dat_df)

sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df



