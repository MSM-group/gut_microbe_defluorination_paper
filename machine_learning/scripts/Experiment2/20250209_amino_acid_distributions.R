# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", "ranger", 
               "doParallel", "permute", "coin")

# Read in the dataset
reg_df <- read_csv('data/Experiment2/20250220_defluorinases_entire_alignment_for_classification.csv')
truthdf <- data.frame(nams = reg_df$nams, 
                      truth = reg_df$truth)

rawdat <- reg_df %>%
  dplyr::select(-2, -truth) %>%
  column_to_rownames("nams") %>% 
  dplyr::filter(!duplicated(.)) %>%
  mutate(nams = rownames(.))

rawdat <- rawdat %>%
  left_join(., truthdf, by = c("nams"))

# Remove variables with nonzero variance (should already be removed)
nozdat <- caret::nearZeroVar(rawdat, saveMetrics = TRUE, 
                             uniqueCut = 1)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE] 
length(which_rem) # cut none

# Check for duplicates
dat <- rawdat  %>%
  dplyr::select(-which_rem)
dat <- dat[!duplicated(dat),] 
nrow(dat)

# Register parallel backend
registerDoParallel(cores = parallel::detectCores())

# Define the function to train the model and extract top features
train_and_extract_features <- function(seed) {
  set.seed(seed)
  
  # Split into test and training data - random option
  dat_split <- rsample::initial_split(dat, prop = 0.8, strata = "truth")
  dat_train <- rsample::training(dat_split)
  dat_test <-rsample::testing(dat_split)
  
  # Independent variables
  x_train <- data.frame(dat_train[,!colnames(dat_train) %in% c("nams", "truth")])
  x_test <- dat_test[,!colnames(dat_test) %in% c("nams", "truth")]

  
  # Dependent variable
  y_train <- dat_train$truth
  y_test <- dat_test$truth
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F)# row.names = dat_train$nams)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F)# row.names = dat_test$nams)
    
  # Make a data frame for prediction
  df_train <- data.frame(x_train, stringsAsFactors = F, 
                         row.names = dat_train$nams)
  
  # Weight the classes by their importance
  wp118 <- 1/length(grep("WP_118709078", dat_train$nams))
  wp118 * length(grep("WP_118709078", dat_train$nams)) # weight corresponds to 1
  wp178 <- 1/length(grep("WP_178618037", dat_train$nams))
  wp178 * length(grep("WP_178618037", dat_train$nams)) # weight corresponds to 1
  case_wts <- case_when(grepl("WP_118709078", dat_train$nams) ~ wp118,
                        grepl("WP_178618037", dat_train$nams) ~ wp178,
                        TRUE ~ 1)
  

  # Train the Ranger model
  rf_full <- ranger(as.factor(y_train) ~., data = form_train, num.trees = 1000, 
                    splitrule = "gini",
                    case.weights = case_wts,
                    mtry = 15, min.node.size = 1,
                    importance = "permutation", probability = TRUE)
  
  
  rf_imp <- rf_full$variable.importance
  names(rf_imp)
  
  # Extract top 50 features and their ranks
  importance_df <- data.frame(
    feature = names(rf_imp),
    importance = rf_imp,
    aa = word(names(rf_imp), sep = "_", -1),
    stringsAsFactors = FALSE
  )
  
  importance_df <- importance_df %>%
    arrange(desc(importance)) %>%
    mutate(rank = row_number()) %>%
    head(95)
  
  rf_pred <- predict(rf_full, data = form_test)
  
 
  preds_df <- bind_cols(nams = dat_test$nams, truth = y_test, 
                        predicted_defluor = rf_pred$predictions[,1],
                        predicted_nondefluor = rf_pred$predictions[,2])  
  return(list(importance_df))
}

# Set random seed 
seeds <- 1:1000
all_top_features_rf <- lapply(seeds, train_and_extract_features)
combined_df <- bind_rows(all_top_features_rf, .id = "iteration")

# Combine the results into a single data frame
combined_df_rf <- combined_df %>%
  dplyr::select(1:5) %>%
  dplyr::filter(!is.na(aa)) %>%
  dplyr::filter(aa!= ".") %>%
  dplyr::filter(rank %in% c(1:20)) %>%
  dplyr::mutate(position = as.numeric(str_extract(feature, "\\-?\\d+\\.?\\d*")))


checksumms <- combined_df_rf %>%
  group_by(iteration) %>%
  dplyr::mutate(summs = sum(position))

checksumms$iteration[which.max(checksumms$summs)] # 778 #531 # 910
combmax <- combined_df_rf %>%
  dplyr::filter(iteration == 984) %>%
  arrange(importance)
combmax$feature <- gsub("residue", "motif", combmax$feature)

combmax$feature <- factor(combmax$feature, levels = as.character(combmax$feature))

pdf("output/revision_figures_final/Experiment2_importance_plot.pdf", height = 6, width = 6)
ggplot(combmax) +
  geom_bar(aes(x = feature, y = importance), stat = "identity") +
  theme_pubr() +
  #theme(text = element_text(color = "gray40")) +
  coord_flip() +
  xlab("") +
  ylab("Importance")
dev.off()

write_csv(combined_df_rf, "output/20250220_1000_classification_predictions.csv")



# Barplot of important AAs
ggplot(combined_df_rf) +
  geom_bar(aes(aa)) +
  theme_pubr()

# How this compares to the distribution of residues in WP_178618037 
wp178 <- readAAStringSet("data/Experiment2/wp178618037.fasta")
wpdf <- data.frame(nam = names(wp178)[1], seq = wp178[1])
merg_split <- wpdf %>%
  tidyr::separate(seq, into = paste0("residue", 1:(width(wp178)[1] + 1)), sep = "") %>%
  reshape2::melt(., id.vars = "nam") %>%
  dplyr::filter(value != "")
nrow(merg_split)

ggplot(merg_split) +
  geom_bar(aes(value)) +
  theme_pubr()


# Calculate proportions
combined_df_rf

amino_acids <- sort(c(unique(combined_df_rf$aa), "C"))[unique(combined_df_rf$aa) != "."] # Example amino acid names
amino_acids
observed_prop_raw <- data.frame(table(combined_df_rf$aa)/length(combined_df_rf$aa))
observed_prop_raw
new_row1<- data.frame(Var1 = c("C"), Freq = c(0))
new_row2<- data.frame(Var1 = c("W"), Freq = c(0))

observed_proportions_long <- bind_rows(observed_prop_raw, new_row1, new_row2) 
observed_proportions <- observed_proportions_long$Freq
names(observed_proportions) <- observed_proportions_long$Var1
observed_proportions <- observed_proportions[order(names(observed_proportions))]

n_total <- table(merg_split$value)  # Total number of occurrences of each amino acid
p_background <- c(table(merg_split$value)/width(wp178))  # Expected proportion of important amino acids (background)
p_background
n_total <- table(merg_split$value)  # Total number of occurrences of each amino acid
p_background <- c(table(merg_split$value)/width(wp178))  # Expected proportion of important amino acids (background)
p_background_df <- data.frame(Var1 = names(p_background), p_background)

df2 <-  observed_proportions_long %>%
  dplyr::left_join(p_background_df, by = "Var1") %>%
  reshape2::melt() %>%
  dplyr::mutate(source = ifelse(variable == "Freq", "important (random forest)", "background"))
df2$variable                  

pdf("~/Downloads/20250209_underlying_amino_acid_distribution.pdf", height = 4, width = 6)
ggplot(df2, aes(x = Var1, y = value, fill = source)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_pubr() +
  xlab("Amino acid") +
  ylab("Proportion") +
  scale_fill_manual(values = c("gray80", "gray20")) +
  theme(legend.title = element_blank())
dev.off()

# Create dataframe
length(observed_proportions)
observed_proportions
length(p_background)
unique(df2$Var1)

data <- data.frame(AminoAcid = unique(df2$Var1), Observed_Prop = observed_proportions, Expected_Prop = c(p_background))

# Set seed for reproducibility
set.seed(123)
?kruskal_test
library("coin")
# Proper permutation test function
perm_test <- function(observed_value, expected_value) {
  # Ensure input is a properly structured data frame with two groups
  df <- data.frame(Group = as.factor(c("Observed", "Expected")),
                   Value = as.numeric(c(observed_value, expected_value)))  # Ensure numeric values
  
  # Run permutation test (using nresample instead of B)
  perm_test_result <- kruskal_test(Value ~ Group, data = df)#distribution = approximate(nresample = 10000))
  
  return(pvalue(perm_test_result))
}

# Apply the permutation test correctly using mapply
p_values <- mapply(perm_test, as.numeric(data$Observed_Prop), as.numeric(data$Expected_Prop))
p_values
# Assign the result to the p_value column
data$p_value <- p_values

# Display results
print(data)

# Identify significant amino acids (p < 0.05)
significant_aa <- data[data$p_value < 0.05, ]
print("Significantly different amino acids:")
print(significant_aa)
?independence_test
