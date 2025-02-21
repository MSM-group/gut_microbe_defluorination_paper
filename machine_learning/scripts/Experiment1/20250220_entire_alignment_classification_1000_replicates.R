# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", 
               "tidymodels", "ranger", "tree", 
               "rsample", "randomForest", "readxl", "ggpubr", "seqinr")

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
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(fd_uname = gsub("g", "", paste0("WP_178618037_1_", label))) %>%
  dplyr::mutate(truth = case_when(activity >= 0.3 ~ "defluor", 
                                  activity <= 0.1 ~ "nondefluor"))  
specdf$fd_uname
m8037$fd_uname
# Read in sequences
m8037 <- read_excel("data/Batch353c_Robinson_m8037_T7Express.xlsx") %>%
  dplyr::left_join(., specdf, by = c("fd_uname"))

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
writeXStringSet(comb, "data/mutlib_all_seqs_aligned.fasta") # optional write to file

# Load alignment file
seqs <- read.alignment("data/mutlib_all_seqs_aligned.fasta", format = "fasta") 
lb_rham<-seqs$seq[[1]] 

# Identify the start of the enzyme
tri.pos <- words.pos("mdvsnv",lb_rham) 
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
  select(-2) %>%
  filter(!duplicated(.))


# Remove variables with nonzero variance (should already be removed)
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, 
                             uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
which_rem <- which_rem[which_rem!="truth"] # cut none
which_rem
# Check for duplicates
dat <- rawdat  %>%
  select(-all_of(which_rem))
dat <- dat[!duplicated(dat),] 
dat$truth <- as.factor(dat$truth)

# Set random seed 
set.seed(20250220)

rf_models <- vector(mode = "list",
                    length = 1000) # 100 different splits

# Split randomly into 10 test/training sets
for(i in 1:1000) { # 100
  # Split into test and training data
  dat_split <- initial_split(dat, strata = "truth", prop = 0.8)
  dat_train <- training(dat_split)
  dat_train$nams
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
                 mtry = 7, 
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

mean_F1_vec <- sapply(1:length(rf_models), function(x) {  rf_models[[x]]$confmat$byClass[['F1']] }) 

mean_oob_vec <- 1-sapply(1:length(rf_models), function(x) {  rf_models[[x]][['oob']] }) 
mean(mean_oob_vec)
sd(mean_oob_vec)

mean_F1_vec[mean_F1_vec == "NaN"] <- NA

dat_df <- tibble(mean_training_vec) %>%
  bind_cols(tibble(mean_testing_vec)) %>%
  bind_cols(tibble(mean_F1_vec)) 
dat_df
colMeans(dat_df, na.rm = T)

sd_df <- dat_df %>%
  dplyr::filter(!is.na(mean_F1_vec)) %>%
  dplyr::summarise_each(sd)
sd_df


