# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr")

# Read in the dataset
dat <- read_csv("data/20250123_regression_feature_table.csv") %>%
  dplyr::select(-fd_uname) # remove label
dat$activity <- as.numeric(dat$activity)

# Set random seed 
seeds <- 1:1000
rf_models <- vector(mode = "list",
                    length = length(seeds))

# Define the function to train the model and extract top features
for(i in 1:length(seeds)) {
  set.seed(i)
  # Split into test and training data
  dat_split <- rsample::initial_split(dat, prop = 0.8)
  dat_train <- rsample::training(dat_split)
  dat_test  <- rsample::testing(dat_split)
  
  # Independent variables
  x_train <- dat_train[,!colnames(dat_train) %in% c("label", "activity")]
  x_test <- dat_test[,!colnames(dat_test) %in% c("label", "activity")]
  
  # Dependent variable
  y_train <- dat_train$activity
  y_test <- dat_test$activity
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)
  
  # Make a data frame for prediction
  df_train <- data.frame(x_train, stringsAsFactors = F, 
                         row.names = dat_train$label)

  # Train the Ranger model
  rf_full <- ranger(y_train ~., data = form_train, num.trees = 1000, 
                    mtry = 3, min.node.size = 1,
                    importance = "permutation")
  
  
  rf_imp <- rf_full$variable.importance
  names(rf_imp)
  
  # Extract top 6 features and their ranks
  importance_df <- data.frame(
    feature = names(rf_imp),
    importance = rf_imp,
    aa = word(names(rf_imp), sep = "_", -1),
    stringsAsFactors = FALSE
  )
  
  importance_df <- importance_df %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number()) %>%
    mutate(iteration = i)
    head(7)

  rf_pred <- predict(rf_full, data = form_test)
  rf_models[[i]]$importance_df <- importance_df
  rf_models[[i]]$train_rmse <- sqrt(rf_full$prediction.error) # Classifcation accuracy 1 - OOB error

  rf_models[[i]]$test_rmse <- Metrics::rmse(rf_pred$predictions, y_test)
  rf_models[[i]]$preds_df <- bind_cols(iteration = i, nams = dat_test$label, truth = y_test, 
                        predicted = rf_pred$predictions)

}

# Perform bootstrapping
mean_training_vec <- unlist(sapply(1:length(seeds), function(x) { rf_models[[x]]$train_rmse })) 
summ_training <- data.frame(mean = mean(mean_training_vec),
                            sd = sd(mean_training_vec))
summ_training
mean_testing_vec <- unlist(sapply(1:length(seeds), function(x) { rf_models[[x]]$test_rmse}))
summ_testing <- data.frame(mean = mean(mean_testing_vec),
                            sd = sd(mean_testing_vec))

combdf <- bind_cols(training_mean = summ_training$mean, training_sd = summ_training$sd, 
                    testing_mean = summ_testing$mean,  testing_sd = summ_testing$sd)
write_csv(combdf, "output/20250211_testing_results.csv")
  
# Look at performance 
combined_df <- bind_rows(lapply(rf_models, `[[`, 4))


# Combine the importance results into a single data frame
importance_df_all <- data.frame(bind_rows(lapply(rf_models, `[[`, 1)))
importance_df_summ <- importance_df_all %>%
  group_by(feature) %>%
  dplyr::summarise(mean_value = mean(importance),
                   sd_value = sd(importance))
importance_df_summ
ggplot(data = importance_df_summ) +
  geom_bar(aes(x = feature, y = mean_value), position = "stack", stat = "identity") +
  theme_pubr() +
  theme(element_text(angle = 45))



