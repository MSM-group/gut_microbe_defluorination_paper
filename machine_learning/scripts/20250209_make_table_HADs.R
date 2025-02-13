# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", "caret",
               "doParallel", "permute", "coin", "viridis")

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

pdf("output/supplemental_figures_final/validation_HAD_boxplot.pdf", width = 3, height = 3)
ggplot(maketab_P, aes(x = proteins, y = predicted_defluor, fill = truth, color = truth)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  geom_hline(aes(yintercept = 0.5), linetype = "dashed") +
  scale_color_manual(values = c("goldenrod", "gray80", "gray80", "gray80", "gray80")) +
  scale_fill_manual(values = c("goldenrod", "gray80", "gray80", "gray80", "gray80")) +
  theme_pubr() +
  ylab("Defluorination probability") +
  xlab("Validation set HADs") +
  theme(legend.title = element_blank())
dev.off()

pdf("output/supplemental_figures_final/validation_HAD_densityplot.pdf", width = 5, height = 3)
ggplot(maketab_P, aes(x = predicted_defluor, fill = proteins)) +
  geom_density(alpha = 0.7) +
  theme_pubr() +
  scale_color_manual(values = c("goldenrod", "orchid1", "skyblue2", "palegreen2", "#FB9A99")) +
  scale_fill_manual(values = c("goldenrod", "orchid1", "skyblue2", "palegreen2", "#FB9A99")) +
  ylab("Density") +
  xlab("Defluorination probability") +
  theme(legend.title = element_blank())
dev.off()

write_csv(summtab, "output/scores_for_5_predicted_HADs.csv")

