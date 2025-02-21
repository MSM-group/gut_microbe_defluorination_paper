# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", 
               "colorspace", "cowplot", "tidymodels", "ranger", "tree", 
               "rsample", "randomForest", "readxl", "ggpubr")

# Read in the dataset
reg_df <- read_csv('data/Experiment2/20250220_defluorinases_C_term_for_classification.csv')
truth <- data.frame(nams = reg_df$nams, 
                    truth = reg_df$truth)
rawdat <- reg_df %>%
  select(-2, -truth) %>%
  column_to_rownames("nams") %>% 
  filter(!duplicated(.)) %>%
  mutate(nams = rownames(.))
rawdat <- rawdat %>%
  left_join(., truth, by = c("nams"))

# Remove variables with nonzero variance (should already be removed)
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, 
                             uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
length(which_rem) # cut none

# Check for duplicates
dat <- rawdat  %>%
  select(-which_rem)
dat <- dat[!duplicated(dat),] 
nrow(dat)

colnames(rawdat) <- paste0(colnames(rawdat), "_", as.character(rawdat[1,]))
colnames(rawdat)[grepl("truth", colnames(rawdat))] <- "truth"
colnames(rawdat)[grepl("nams", colnames(rawdat))] <- "nams"

# Set random seed 
set.seed(20250220)

rf_models <- vector(mode = "list",
                    length = 1000) # 100 different splits

# Split randomly into 10 test/training sets
for(i in 1:1000) { # 100
  # Split into test and training data
  dat_split <- initial_split(dat, strata = "truth", prop = 0.8)
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  
  # Define variables of interest
  x_train <- dat_train %>%
    dplyr::select(-nams, -truth)
  x_test <- dat_test %>%
    dplyr::select(-nams, -truth)
  
  y_train <- as.factor(dat_train$truth)
  y_test <- as.factor(dat_test$truth)
  
  df_train <- data.frame(x_train, stringsAsFactors = F)
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F) # row.names = dat_train$nams)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F) # row.names = dat_test$nams)
  
  # Weight the classes by their importance
  wp118 <- 1/length(grep("WP_118709078", dat_train$nams))
  wp118 * length(grep("WP_118709078", dat_train$nams)) # weight corresponds to 1
  wp178 <- 1/length(grep("WP_178618037", dat_train$nams))
  wp178 * length(grep("WP_178618037", dat_train$nams)) # weight corresponds to 1
  case_wts <- case_when(grepl("WP_118709078", dat_train$nams) ~ wp118,
                        grepl("WP_178618037", dat_train$nams) ~ wp178,
                        TRUE ~ 1)
  
  ## Train model
  rf_1 <- ranger(y_train ~., data = form_train, 
                 num.trees = 1000, 
                 splitrule = "extratrees", 
                 mtry = 5, 
                 min.node.size = 1, 
                 case.weights = case_wts,
                 importance = "permutation")
  
  rf_models[[i]]$models <- rf_1
  rf_models[[i]]$oob <- 1 - rf_1$prediction.error # Classification accuracy 1 - OOB error

  rf_1_pred <- predict(rf_1, form_test)
  
  rf_models[[i]]$confmat <- confusionMatrix(rf_1_pred$predictions, as.factor(dat_test$truth))
}

mean_training_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$oob }) 
mean_training_vec
traindf <- tibble(mean_training_vec)
mean(mean_training_vec)
sd(mean_training_vec) 

mean_testing_vec <- sapply(1:length(rf_models), function(x) { rf_models[[x]]$confmat$overall[1]})
testdf <- tibble(mean_testing_vec)
mean(mean_testing_vec)
sd(mean_testing_vec) 

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

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) 

dat_df
colMeans(dat_df, na.rm = T)

sd_df <- dat_df %>%
  dplyr::summarise_each(sd)
sd_df


