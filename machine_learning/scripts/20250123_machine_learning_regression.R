# Install packages
library("tidyverse")
library("readxl")
library("ggpubr")
library("data.table")
library("ranger")
library("caret")
library("ggpmisc")

# Read in the protein normalized data
dat <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                  range = c("C135:N156"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
dat
attr(dat, "na.action") <- NULL
attr(dat, "class") <- NULL

# Read in the template
temp1 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() %>%
  na.omit() 
attr(temp1, "na.action") <- NULL
attr(temp1, "class") <- NULL
temp1

rawdf <-  bind_cols(label = temp1, 
                    value = dat)
wt <- rawdf$value[rawdf$label == "WT"]  


specdf <- bind_cols(label = temp1, activity = dat) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(truth = case_when(activity >= 1.0 ~ "defluor", 
                                  activity <= 0.3 ~ "nondefluor"))  %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(fd_uname = paste0("WP_178618037_1_", label)) %>%
  dplyr::mutate(fd_uname = gsub("_g", "_", fd_uname)) %>%
  dplyr::mutate(activity = log10(activity))# take the log10 of the activity

table(specdf$aa)
dat_list <- as.list(specdf$aa)
names(dat_list) <- specdf$label
dat_list

# Convert AAs to features
source("scripts/convert_seq_5aap.R")
extract_feat_list <- lapply(1:length(dat_list), function(x) { convert_seq_5aap(dat_list[[x]]) })


extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)
extract_feat_df
colnames(extract_feat_df) <- c("polarity", "secondary_structure", "size", "codon_diversity", "charge")
wtact <- specdf$activity[specdf$label == "WT"]
hist(specdf$activity) 
mean(specdf$activity) # 2.42

feat_df <- specdf %>%
  bind_cols(extract_feat_df) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::mutate(position = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::filter(!is.na(position)) %>%
  dplyr::filter(activity > mean(activity)) %>% 
  dplyr::filter(!label %in% c("P20", "P21", "P22", "P23", "P24"))
hist(feat_df$activity)

write_csv(feat_df, "data/20250123_regression_feature_table.csv")


# Try predicting the activity of 20%
set.seed(20250124)

# Split into test and training data 
# 80% training
# 20% test
dat <- feat_df %>%
  dplyr::select(-truth, -fd_uname)
colnames(dat)
dat_split <- rsample::initial_split(dat, prop = 0.8)
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("label", "activity")]
x_test <- dat_test[,!colnames(dat_test) %in% c("label", "activity")]

# Dependent variable
y_train <- dat_train$activity
y_test <- dat_test$activity
y_test # check there is a mix

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$label)
y_train
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
rf$finalModel$prediction.error # out-of-bag error 0.04

#Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp

pdf("data/variable_importance_plot.pdf", width = 10)
ggplot(rf_imp, top = 7) + 
  xlab("") +
  theme_pubr(base_size = 14)
dev.off()

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

# Predict for new data
Metrics::rmse(y_test, rf_pred) # 0.23

rf_df <- data.frame(cbind(rf_pred, y_test))
my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))


# Residual plot to check for skew
rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)

ggplot(rf_df_resid, aes(x = rf_pred, y = resid)) +
  geom_point(alpha = .3) + 
  geom_line(col = "red", y = 0.0,
            lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Predicted enzyme activity") +
  ylab("Residuals") 

my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))
ggplot(rfdat, aes(x = obs, y = pred)) +
  geom_point(aes(x = obs, y = pred, alpha = 0.8)) + 
  geom_smooth(se = FALSE, col = "red", method = "lm",
               lty = 2, lwd = 1, alpha = .5) +
  stat_poly_eq(formula = my.formula,
              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
              parse = TRUE) +
  theme_pubr() +
  xlab("Observed") +
  ylab("Predicted") 

