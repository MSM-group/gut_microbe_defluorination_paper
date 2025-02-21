# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", 
               "colorspace", "cowplot", "tidymodels", "ranger", "tree", 
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

expdat <- data.frame(activity = dat) %>%
  bind_cols(label = temp1) %>%
  dplyr::mutate(aa = substr(label, 2, 2)) %>%
  arrange(desc(activity)) %>% # mean activity
  dplyr::filter(aa != "O") %>%
  dplyr::mutate(newlab = paste0("WP_178618037_1_", substr(label, 2, nchar(label))))

# Read in the sequences
m8037 <- read_excel("data/Batch353c_Robinson_m8037_T7Express.xlsx") %>%
  dplyr::left_join(., expdat, by = c("fd_uname" = "newlab")) %>%
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

mutlib <- AAStringSet(x = comball$protein)
names(mutlib) <- comball$fd_uname

# Align all sequences
pos <- readAAStringSet("data/pos_defluorinases_deduplicated.fasta")
neg <- readAAStringSet("data/neg_defluorinases_deduplicated.fasta")
comb <- AlignSeqs(AAStringSet(c(pos, neg, mutlib)))
names(comb) <- c(names(pos), names(neg), names(mutlib))

# BrowseSeqs(comb)
writeXStringSet(comb, "data/mutlib_all_seqs_aligned.fasta") # optional write to file

# Load alignment file
seqs <- read.alignment("data/mutlib_all_seqs_aligned.fasta", format = "fasta") 
lb_rham <- seqs$seq[[1]] 

# Identify the start and end of the C-terminal region
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
nuc2<-lapply(seqs$seq,function(x) { substr(x,tri.pos2,end.pos2+1) })
nuc2
nucr2<-unlist(nuc2)

# Concatenate two sections
nucr <- paste0(nucr1, nucr2)
names(nucr) <- seqs$nam
appendf <- data.frame(nams = names(nucr), motif = nucr)

# Make a logo of the C-terminal motif 
motif_seqs <- data.frame(toupper(substr(appendf$motif, 
                                        start = 23,
                                        stop = nchar(appendf$motif))))

fig1 <- ggplot() + 
  geom_logo(motif_seqs, method = "p", 
            col_scheme = "auto") +
  theme_logo() +
  theme(axis.text.x = element_text(angle = 45))
fig1
fig1$scales$scales[[1]] <- scale_x_continuous(breaks = c(1:18), labels=as.character(c(274:(275+16))))
fig1

# Merge with the defluorination data
merg <- appendf %>%
  dplyr::mutate(motif = toupper(motif))
nchar(merg$motif)

merg_split <- merg %>%
  tidyr::separate(motif, into = paste0("residue", c((tri.pos1-12):(end.pos1-11), (tri.pos2-11):(290-11), "_", (292-11):(end.pos2-10))), sep = "")
appendf <- data.frame(nams = names(nucr), motif = nucr)

# Write unencoded data frame to file
reg_df <- merg_split %>%
  dplyr::mutate(truth = c(rep("defluor", length(pos)), rep("nondefluor", length(neg)),
                          comball$truth))

# Remove columns that the machine learning model should not learn e.g., nams, truth
rawdat <- reg_df %>%
  select(-(2)) 

colnames(rawdat) <- paste0(colnames(rawdat), "_", as.character(rawdat[1,]))
colnames(rawdat)[grepl("truth", colnames(rawdat))] <- "truth"
colnames(rawdat)[grepl("nams", colnames(rawdat))] <- "nams"

# Remove variables with nonzero variance 
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 

dat <- rawdat  %>%
  select(-all_of(which_rem))

# Check and remove duplicates
dups <- dat %>%
  dplyr::select(-1) %>%
  dplyr::filter(duplicated(.))
dat <- dat[!duplicated(dat),] %>%
  dplyr::filter(!is.na(truth)) %>%
  dplyr::filter(!nams %in% rownames(dups))
rownames(dat)
colnames(dat)
table(dat$truth)

# Set random seed 
set.seed(20240521)

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


