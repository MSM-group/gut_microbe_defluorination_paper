# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggpmisc", "caret", 
               "colorspace", "cowplot", "tidymodels", "ranger", "tree", 
               "rsample", "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Read in the template
temp1 <- read_excel("data/Machine_learning/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res <- read_excel('data/Machine_learning/20240427_alanine_scanning_replicate_average_fluoride_data.xlsx', 
                  col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 

res <- res[!is.na(res)]

# Delta calculation relative to WT
wt <- res[temp1 == "WT"]
wt

# Delta
delta <- wt - res # delta between WT -> Ala_mut

dat <- bind_cols(label = temp1, 
                 delta = delta) %>%
  dplyr::filter(label != "WT") %>%
  dplyr::mutate(position_index = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::mutate(position_from = substr(label, 2, 2)) %>%
  dplyr::mutate(position_to = substr(label, nchar(label), nchar(label))) %>%
  dplyr::mutate(unique_id = paste0(position_from, "_", position_index, "_", position_to, "_fwd"))
table(dat$position_from)

# Delta mirrored
delta_mirror <- res - wt # delta between Ala_mut -> WT
hist(delta_mirror)
dat_mirror <- bind_cols(label = temp1,
                        delta = delta_mirror) %>%
  dplyr::filter(label != "WT") %>%
  dplyr::mutate(position_index = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::mutate(position_from = substr(label, nchar(label), nchar(label))) %>%
  dplyr::mutate(position_to = substr(label, 2, 2)) %>%
  dplyr::mutate(unique_id = paste0(position_from, "_", position_index, "_", position_to, "_rev"))

# Make feature vector
# source("scripts/R/convert_seq_5aap.R")
# dat_list <- as.list(unique(dat$position_from))
# names(dat_list) <- unique(dat$position_from)
# extract_feat_list <- lapply(1:length(dat_list), function(x) { convert_seq_5aap(dat_list[[x]]) })
# extract_feat_df <- data.frame(matrix(unlist(extract_feat_list), nrow = length(extract_feat_list), byrow=T), stringsAsFactors=FALSE)
# extract_feat_df
# colnames(extract_feat_df) <- c("polarity", "secondary_structure", "size", "codon_diversity", "charge")
#extract_feat_df

# Read in 5 AA properties
props <- read_csv("data/Priestia_ML_analysis/5_aa_properties.csv")
props15 <- read_table("data/Priestia_ML_analysis/15_physicochemical_properties.txt")
colnames(props15)[1] <- "AA_ABREV"
colnames(props15)
feat_df <- dat %>%
 # bind_rows(dat_mirror) %>%
  dplyr::left_join(., props15, by = c("position_to" = "AA_ABREV")) %>%
 # dplyr::mutate(position_to = Factor_1_polarity) %>%
  #dplyr::mutate(paste0(c('to_', colnames(props15)) = props15))
  #dplyr::select(1:6) %>%
 # dplyr::left_join(., props15, by = c("position_from" = "AA_ABREV")) %>%
  #dplyr::mutate(position_from = Factor_1_polarity) %>%
  #dplyr::select(1:6) %>%
  #column_to_rownames(., var = "unique_id") %>%
  dplyr::mutate(label = rownames(.))
feat_df


# Try predicting the activity of 20%
set.seed(78910)
# Split into test and training data 
# 80% training
# 20% test
# TODO exclude the samples which HAD SEQUENCING ERRORS
errors <- read_excel("data/Priestia_ML_analysis/Batch353c_Robinson_Sanger_Samples_v2.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(sanger_status == "re-clone") %>%
  pull(gene_id)
errors
m8037 <- read_excel("data/Priestia_ML_analysis/Batch353c_Robinson_m8037_T7Express.xlsx") 

tosamp <- m8037 %>%
  janitor::clean_names() %>%
  dplyr::filter(top_construct_name %in% errors) %>%
  mutate(index1 = word(fd_uname, sep = "_", 4)) %>%
  dplyr::mutate(index2 = substr(index1, 2, nchar(index1))) %>%
  dplyr::mutate(index3 = gsub("A", "", index2)) %>%
  pull(index3) %>%
  as.numeric()
tosamp

dat <- feat_df
colnames(dat)
# dat_check <- dat %>%
#   dplyr::select(1, 246, 220:245)

# Split with known hold-out test set
dat_train <- dat %>% 
  dplyr::filter(!position_index %in% tosamp)
dat_test <- dat %>% 
  dplyr::filter(position_index %in% tosamp) # 18 test
nrow(dat_test) # 23 sequences held out
nrow(dat_train)/nrow(dat)  # 96% train/test split

# Independent variables
x_train <- dat_train[,!colnames(dat_train) %in% c("unique_id", "label", "delta")]
x_test <- dat_test[,!colnames(dat_test) %in% c("unique_id", "label", "delta")]
colnames(x_train)

# Dependent variable
y_train <- dat_train$delta
y_test <- dat_test$delta
y_test # check there is a mix of your two variables

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$label)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$label)

# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, 
                       row.names = dat_train$unique_id)

# Optional tuning of random forest parameters
#mtrys <- c(round(log2(ncol(df_train)), 0), round(sqrt(ncol(df_train)), 0), round(ncol(df_train)/2, 0))
#mtrys # number of variables available for splitting at each tree node

rf_grid <- expand.grid(mtry = 3,
                       splitrule = c("extratrees"),
                       min.node.size = 1)

# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 1, # increase this to 3 when you run the code 
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")
warnings()


# Training set accuracy
getTrainPerf(rf) # Training set accuracy 
rf$finalModel$prediction.error # out-of-bag error

# Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred

#pdf("data/machine_learning/variable_importance_plot.pdf", width = 10)
ggplot(rf_imp, top = 3) + 
  xlab("") +
  theme_pubr(base_size = 20)
#dev.off()

Metrics::rmse(y_test, rf_pred) # 0.22
rf_df <- data.frame(cbind(rf_pred, y_test))

my.formula <- y ~ x
summary(lm(rf_df$y_test ~ rf_pred))
summary(lm(rf_pred ~ rf_df$y_test))

rf_df_resid <- rf_df %>%
  mutate(resid = y_test - rf_pred)

ggplot(rf_df_resid, aes(x = rf_pred, y = resid)) +
  # geom_abline(col = "green", alpha = .5) + 
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
  ylab("Predicted") +
  xlim(-0.5, 0.5) +
  ylim(-0.5, 0.5)
# stat_poly_eq(formula = my.formula,
#              aes(label = paste(..rr.label.., sep = "~~~")),
#              parse = TRUE)

rfdf_train <- data.frame(cbind(rf$finalModel$predictions, y_train), stringsAsFactors = F)
colnames(rfdf_train) <- c("Pred", "Obs")


rfdf_test <- data.frame(cbind(rf$finalModel$predictions, y_train), stringsAsFactors = F)
colnames(rfdf_test) <- c("Pred", "Obs")

ggplot(rfdf_train, aes(x = Pred, y = Obs)) +
  geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  # geom_smooth(se = FALSE, col = "red", method = "lm",
  #              lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  xlim(-0.6, 0.6) +
  ylim(-0.6, 0.6) +
  stat_poly_eq(formula = my.formula,
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = TRUE)

ggplot(rfdf_test, aes(x = Pred, y = Obs)) +
  # geom_abline(col = "green", alpha = .5) + 
  geom_point(alpha = .3) + 
  # geom_smooth(se = FALSE, col = "red", method = "lm",
  #             lty = 2, lwd = 1, alpha = .5) +
  theme_pubr() +
  xlab("Observed enzyme activity") +
  ylab("Predicted enzyme activity") +
  xlim(-0.4, 0.4) +
  ylim(-0.4, 0.4)



