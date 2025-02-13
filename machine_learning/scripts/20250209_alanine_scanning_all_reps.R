# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", "ggridges")

# Read in the protein normalized data
rep1 <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("C110:N132"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
rep1

rep2 <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("P110:AA132"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
rep2

rep3 <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("AC110:AN132"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
rep3

# Read in the template
temp1 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() %>%
  na.omit() 
temp1
attr(temp1, "na.action") <- NULL
attr(temp1, "class") <- NULL
temp1
raw_wide <-  bind_cols(label = temp1,
                       rep1 = rep1,
                       rep2 = rep2,
                       rep3 = rep3) %>%
  dplyr::filter(!label %in% c("P20", "P21", "P22", "P23", "P24")) 
rawdf <- reshape2::melt(raw_wide)

# Calculate the activity relative to wild-type
wt <- mean(rawdf$value[rawdf$label == "WT"])
attr(rawdf, "na.action") <- NULL
attr(rawdf, "class") <- NULL
attr(rawdf, "row.names") <- NULL

specdf <- data.frame(rawdf) %>%
  dplyr::mutate(delta = wt-value) %>%
  dplyr::mutate(position_index = as.numeric(gsub('[[:alpha:]]', "", substr(label, 3, 5)))) %>%
  dplyr::mutate(position_from = substr(label, 2, 2)) %>%
  dplyr::mutate(position_to = substr(label, nchar(label), nchar(label))) %>%
  dplyr::mutate(unique_id = paste0(position_from, "_", position_index, "_", position_to, "_fwd")) %>%
  dplyr::mutate(equal_segment = ntile(position_index, n = 5)) %>%
  dplyr::mutate(chimera_segment_six = ntile(position_index, n = 6)) %>%
  dplyr::mutate(chimera_segment = case_when(position_index %in% 1:35 ~ 1,
                                            position_index %in% 36:69 ~ 2,
                                            position_index %in% 70:109 ~ 3,
                                            position_index %in% 110:150 ~ 4,
                                            position_index %in% 150:196 ~ 5,
                                            position_index %in% 197:238 ~ 6,
                                            TRUE ~ NA)) %>%
  dplyr::filter(label != "WT")
attr(specdf$delta, "na.action") <- NULL
attr(specdf$delta, "class") <- NULL
attr(specdf$delta, "row.names") <- NULL
specdf <- data.frame(specdf)
specdf$delta

pdf("output/revision_figures_final/geom_violin_significance.pdf", width = 4, height = 6)
ggplot(specdf, aes(x = position_index, y = delta)) +
  geom_violin(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95))+
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Amino acid residue number") +
  ylab("Wild-type enzyme - alanine variant activity") 
dev.off()

pdf("output/revision_figures_final/geom_point_significance.pdf", width = 4, height = 6)
ggplot(specdf, aes(x = as.factor(chimera_segment), y = delta)) +
  #geom_boxplot(coef = 0, aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
   geom_jitter(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Chimera segment") +
  ylab("Wild-type enzyme - alanine variant activity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_signif(comparisons = list(c(6,5),
                               c(6,4),
                               c(6,3),
                               c(6,2),
                               c(6,1)),
            y_position = 900, tip_length = 0, vjust = 0.7,
            step_increase = 0.05,  # To stagger the bars
            textsize = 1,
            map_signif_level = T)
dev.off()
# 
# identify_outliers <- function(x) {
#   Q1 <- quantile(x, 0.25)
#   Q3 <- quantile(x, 0.75)
#   IQR <- Q3 - Q1
#   lower_bound <- Q1 - 1.5 * IQR
#   upper_bound <- Q3 + 1.5 * IQR
#   return(x < lower_bound | x > upper_bound)  # TRUE for outliers
# }
# 
# outlierdf <- specdf %>%
#   mutate(outlier = identify_outliers(delta))  # 'hwy' is the numeric variable of interest
# outlierdf
# 
# # Step 3: Create a contingency table for outliers by class
# outlier_table <- outlierdf %>%
#   count(chimera_segment, outlier) %>%
#   spread(outlier, n, fill = 0)  # Spread the outliers into separate columns
# outlier_table
# outlier_counts_pairwise <- outlier_table %>%
#   filter(chimera_segment %in% c("6", "5")) 
# 
# # Step 3: Perform the Chi-Square test
# chisq_test_pairwise <- chisq.test(outlier_counts_pairwise[, -1])  # Exclude 'class' column
# print(chisq_test_pairwise)
# # 6 and 3 are significant, 6 and 4 are significant