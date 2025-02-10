# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", 
               "doParallel", "permute", "coin")

# Make a table
maketab <- read_csv("output/20250209_1000_classification_predictions.csv") %>%
  dplyr::mutate(prediction_class = as.factor(ifelse(predicted_defluor >= 0.5, "defluor", "nondefluor")))

# Calculate overall statistics
maketab$truth
cm_rf <- confusionMatrix(as.factor(maketab$truth), maketab$prediction_class)
cm_rf


maketab_trim <- maketab %>%
  dplyr::filter(!grepl("WP_178|WP_118", nams))
cm_rf <- confusionMatrix(as.factor(maketab_trim$truth), maketab_trim$prediction_class)
cm_rf

maketab_P <- maketab_trim %>%
  dplyr::filter(grepl("P2", nams)) %>%
  dplyr::mutate(proteins = word(nams, sep = "_", 1))
maketab_P
cm_rf <- confusionMatrix(as.factor(maketab_P$truth), maketab_P$prediction_class)
cm_rf

summtab <- maketab_P %>%
  dplyr::group_by(nams) %>%
  summarise(
    mean_value = mean(predicted_defluor),
    sd_value = sd(predicted_defluor)
  )
maketab_P$predicted_defluor

pdf("output/20250209_ML_predictions_P2025.pdf", width = 3, height = 3)
ggplot(maketab_P, aes(x = proteins, y = predicted_defluor, fill = truth, color = truth)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  #geom_hline(aes(yintercept = 0.5), "dashed") +
 # geom_jitter(width = 0.20) +
  #geom_errorbar(data = summtab, aes(x = nams, ymin = mean_value - sd_value, ymax = mean_value + sd_value)) +
  theme_pubr() +
#  theme(axis.text.x = element_text(angle = 45, vjust = -0.0005)) +
  ylab("Defluorination probability") +
  xlab("Validation set HADs")
dev.off()



write_csv(summtab, "output/scores_for_5_predicted_HADs.csv")
