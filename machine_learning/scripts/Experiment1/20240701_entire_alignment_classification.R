# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "caret", 
                "tidymodels", "ranger", "tree", "seqinr",
               "rsample", "randomForest", "readxl", "ggpubr")


# Read in the template
temp1 <- read_excel("data/Experiment1/template_exp1.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
dat <-  read_excel('data/Experiment1/20240427_rep1_2_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
dat <- as.numeric(dat[!is.na(dat)])

specdf <- bind_cols(label = temp1, activity = dat) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(truth = case_when(activity >= 0.3 ~ "defluor", 
                                  activity <= 0.1 ~ "nondefluor"))  %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(fd_uname = paste0("WP_178618037_1_", label))

# Read in sequences
m8037 <- read_excel("data/Batch353c_Robinson_m8037_T7Express.xlsx") %>%
  dplyr::left_join(., specdf, by = c("fd_uname")) %>%
  dplyr::mutate(truth = case_when(activity >= 1.0 ~ "defluor", 
                                  activity <= 0.3 ~ "nondefluor"))  
m9078 <- read_excel("data/Batch353c_Robinson_m9078_T7Express.xlsx") %>%
  dplyr::mutate(truth = "nondefluor")
comball <- m8037 %>%
  bind_rows(m9078) %>%
  janitor::clean_names() %>%
  dplyr::mutate(nuc = as.character(toupper(dna_sequence))) %>%
  dplyr::mutate(protein = microseq::translate(nuc)) %>%
  dplyr::filter(!is.na(truth))

dat_list <- as.list(specdf$aa)
names(dat_list) <- specdf$label
mutlib <- AAStringSet(x = comball$protein)
names(mutlib) <- comball$fd_uname

# Align all sequences
pos <- readAAStringSet("data/pos_defluorinases_deduplicated.fasta")
neg <- readAAStringSet("data/neg_defluorinases_deduplicated.fasta")
comb <- AlignSeqs(AAStringSet(c(pos, neg, mutlib)))
names(comb) <- c(names(pos), names(neg), names(mutlib))
BrowseSeqs(comb)
writeXStringSet(comb, "data/mutlib_all_seqs_aligned.fasta")

# Load alignment file
seqs <- read.alignment("data/mutlib_all_seqs_aligned.fasta", format = "fasta") 
lb_rham<-seqs$seq[[1]] 
seqs$nam

# Identify the start of the enzyme
tri.pos<-words.pos("mdvsnv",lb_rham) 
tri.pos
nchar(lb_rham)

# Find it in all sequences
nuc<-lapply(seqs$seq,function(x) { substr(x,tri.pos,tri.pos+nchar(lb_rham)) })
nucr<-unlist(nuc)
names(nucr)<-seqs$nam
appendf <- data.frame(nams = names(nucr), motif = nucr)

# Merge with the defluorination data
merg <- appendf %>%
  dplyr::mutate(motif = toupper(motif))
merg_split <- merg %>%
  tidyr::separate(motif, into = paste0("residue", 1:width(nucr)[1]), sep = "") 
nrow(merg_split)

reg_df <- merg_split %>%
  dplyr::mutate(truth = c(rep("defluor", length(pos)), 
                          rep("nondefluor", length(neg)),
                          comball$truth[!is.na(comball$truth)]))


# Remove columns that the machine learning model should not see, e.g., nams, truth
rawdat <- reg_df %>%
  select(-2)

colnames(rawdat) <- paste0(colnames(rawdat), "_", as.character(rawdat[1,]))
colnames(rawdat)[grepl("truth", colnames(rawdat))] <- "truth"
colnames(rawdat)[grepl("nams", colnames(rawdat))] <- "nams"

# Remove variables with nonzero variance 
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
length(which_rem) # cut 42
which_rem

dat <- rawdat  %>%
  select(-which_rem) %>%
  dplyr::mutate(truth = rawdat$truth) %>%
  select(-contains("nams"))

# Check for duplicates
dat <- dat[!duplicated(dat),] 
colnames(dat)

# Set random seed 
set.seed(1234) 

# Split into test and training data - random option
dat_split <- rsample::initial_split(dat, prop = 0.8, strata = "truth")
dat_train <- rsample::training(dat_split)
dat_test  <- rsample::testing(dat_split)

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

# Train a machine learning model
rf <- train(
  x = df_train,
  y = y_train,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, # this is how many folds
                           repeats = 3, # increase this to 3 when you run the code 
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  verbose = TRUE,
  importance = "permutation")

# Training set accuracy
getTrainPerf(rf) # Training set accuracy 95% !
rf$finalModel$prediction.error # out-of-bag error 5%

# Plot of variable importance
rf_imp <- varImp(rf, scale = FALSE, 
                 surrogates = FALSE, 
                 competes = FALSE)
rf_imp
rf_imp

rf1 <- ggplot(rf_imp, top = 20) + 
  xlab("") +
  theme_classic()
rf1


rf_roc_train <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                          predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                          ci = T)

# Testing set
rf_pred <- predict(rf, newdata = form_test)
rf_pred
see_preds <- bind_cols(rownames(dat_test), form_test$y_test, rf_pred)
see_preds
cm_rf <- confusionMatrix(rf_pred, as.factor(y_test))
cm_rf

# ROC curve
rf$pred$obs
rf$pred$pred

rf_roc <- pROC::roc(response = ifelse(rf$pred$obs == "nondefluor", 0, 1),
                    predictor = ifelse(rf$pred$pred == "nondefluor", 0, 1),
                    plot = TRUE)
rf_roc
cm_rf$table

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

auc <- round(rf_roc$auc, 4)
auc

