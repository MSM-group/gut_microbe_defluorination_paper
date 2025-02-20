# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "caret", 
               "tidymodels", "ranger", "tree", "seqinr", "ggseqlogo",
               "rsample", "randomForest", "readxl", "ggpubr")

# Read in the protein normalized data
dat <- read_csv("data/Experiment2/20250220_rep1_2_3_linearized_fluoride_data.csv")
specdf <- dat %>%
  dplyr::filter(!label %in% c("WT", "P20", "P21", "P22", "P23", "P24")) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>%
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(activity = (activity - min(activity, na.rm=T))/(max(activity,na.rm=T) - min(activity,na.rm=T))) %>%
  dplyr::mutate(truth = case_when(activity >= 0.3 ~ "defluor",  
                                  activity <= 0.1 ~ "nondefluor"))  %>%
  dplyr::filter(complete.cases(.)) %>%
  dplyr::mutate(fd_uname = paste0("WP_178618037_1_", label)) %>%
  dplyr::mutate(fd_uname = gsub("_g", "_", fd_uname))
table(specdf$truth)
summary(specdf$activity)
hist(specdf$activity,breaks = 20)

# Read in sequences
m8037 <- read_excel("data/Batch353c_Robinson_m8037_T7Express.xlsx") %>%
  dplyr::left_join(., specdf, by = c("fd_uname")) %>%
  dplyr::mutate(truth = case_when(activity >= 0.3 ~ "defluor",   
                                  activity <= 0.1 ~ "nondefluor"))  

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

# Align positive and negative examples
pos <- readAAStringSet("data/pos_defluorinases_deduplicated.fasta")
neg <- readAAStringSet("data/neg_defluorinases_deduplicated.fasta")
comb <- AlignSeqs(AAStringSet(c(pos, neg, mutlib)))
BrowseSeqs(comb)
#writeXStringSet(comb, "data/mutlib_all_seqs_aligned.fasta")

# Load alignment file
seqs <- read.alignment("data/mutlib_all_seqs_aligned.fasta", format = "fasta") 
lb_rham<-seqs$seq[[1]] 
seqs$nam

# Identify the start of the enzyme
nchar("PTLYVPRPLEYGAVNNHVEPEAKYDHETVADFRELAARLGV") # full region is 41 amino acids, however alignment has a gap
tri.pos1 <- words.pos("ptl", lb_rham) # C-terminal motif starts with 'PTL'
end.pos1 <- words.pos("nnh", lb_rham)
tri.pos1
end.pos1

# Section 1
nuc1 <- lapply(seqs$seq, function(x) { substr(x,tri.pos1,end.pos1) })
nuc1
nucr1 <- unlist(nuc1)
nucr1

# Section 2
tri.pos2 <- words.pos("vepe",lb_rham) # ptl #mdvsnv
end.pos2 <- words.pos("gv-", lb_rham)
tri.pos2
end.pos2
nuc2<-lapply(seqs$seq, function(x) { substr(x,tri.pos2,end.pos2+1) })
nuc2
nucr2<-unlist(nuc2)

# Concatenate two sections
nucr <- paste0(nucr1, nucr2)
names(nucr) <- seqs$nam
appendf <- data.frame(nams = names(nucr), motif = nucr)

# Merge with the defluorination data
merg <- appendf %>%
  dplyr::mutate(motif = toupper(motif))
nchar(merg$motif)

merg_split <- merg %>%
  tidyr::separate(motif, into = paste0("residue", c((tri.pos1-12):(end.pos1-11), (tri.pos2-11):(290-11))), sep = "") 

# Write data frame to file
reg_df <- merg_split %>%
  dplyr::mutate(truth = c(rep("defluor", length(pos)), 
                          rep("nondefluor", length(neg)),
                          comball$truth[!is.na(comball$truth)]))

# write_csv(reg_df, "data/Experiment2/20250220_defluorinases_C_term_for_classification.csv")

# Remove columns that the machine learning model should not learn e.g., nams, truth
reg_df <- read_csv("data/Experiment2/20250220_defluorinases_C_term_for_classification.csv") 

rawdat <- reg_df %>%
  select(-(1:2)) %>%
  as.data.frame()
rownames(rawdat) <- reg_df$nams
colnames(rawdat) <- paste0(colnames(rawdat), "_", as.character(rawdat[1,]))
colnames(rawdat)[grepl("truth", colnames(rawdat))] <- "truth"
colnames(rawdat)[grepl("nams", colnames(rawdat))] <- "nams"

# Remove variables with nonzero variance 
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 

dat <- rawdat  %>%
  select(-all_of(which_rem)) 
# Check and remove duplicates
dat <- dat[!duplicated(dat),] %>%
  dplyr::filter(!is.na(truth))
rownames(dat)
colnames(dat)

# Set random seed 
set.seed(20250220) 

# Split into test and training data - random option
dat_split <- rsample::initial_split(dat, prop = 0.8, strata = "truth")
dat_train <- rsample::training(dat_split)
dat_test <- rsample::testing(dat_split)

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

# Weight the classes by their importance
wp118 <- 1/length(grep("WP_118709078", rownames(x_train)))
wp118 * length(grep("WP_118709078", rownames(x_train))) # weight corresponds to 1
wp178 <- 1/length(grep("WP_178618037", rownames(x_train)))
wp178 * length(grep("WP_178618037", rownames(x_train))) # weight corresponds to 1
case_wts <- case_when(grepl("WP_118709078", rownames(x_train)) ~ wp118,
                      grepl("WP_178618037", rownames(x_train)) ~ wp178,
                      TRUE ~ 1)
bind_cols(case_wts, rownames(x_train))
length(case_wts)

rf_full <- ranger(as.factor(y_train) ~., data = form_train, num.trees = 1000,
                  splitrule = "gini",
                  case.weights = case_wts,
                  mtry = 5, min.node.size = 1,
                  importance = "permutation", probability = TRUE)

# Training set accuracy
rf_full$predictions
rf_full$prediction.error # out-of-bag error
mod_train_pred <- as.factor(ifelse(rf_full$predictions[,1] >= 0.5, "defluor", "nondefluor"))
confusionMatrix(mod_train_pred, as.factor(dat_train$truth))
saveRDS(rf, "data/Experiment2/20250220_C_term_classification_random_forest.rds")

# Testing set
rf_pred <- predict(rf_full, data = form_test)
mod_test_pred <- as.factor(ifelse(rf_pred$predictions[,1] >= 0.5, "defluor", "nondefluor"))
mod_test_pred
cm_rf <- confusionMatrix(mod_test_pred, as.factor(dat_test$truth))
cm_rf
see_preds <- bind_cols(rownames(dat_test), y_test, mod_test_pred, rf_pred$predictions)





