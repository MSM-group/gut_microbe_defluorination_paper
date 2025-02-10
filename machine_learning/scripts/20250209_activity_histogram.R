# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", "diptest")

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

attr(dat, "na.action") <- NULL
attr(dat, "class") <- NULL

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


mean(wt)
summary(rawdf$value)
vtert <- quantile(rawdf$value, c(0:5/5))
vtert

pdf("output/log10_activity_histogram.pdf")
ggplot(rawdf, aes(x = log10(value))) +
  geom_histogram(binwidth = 0.1, fill = "gray80") +
  theme_pubr() +
  geom_vline(aes(xintercept = log10(830.7176)), size = 2, color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = 2.322878), size = 2, color = "blue", linetype = "dashed")
dev.off()

dip_result <- dip(log10(rawdf$value), full.result = TRUE)
dip_test <- dip.test(log10(rawdf$value))
kmeans_result <- kmeans(log10(rawdf$value), centers = 2)

# Find the centers (means) of the two clusters using k-means
cluster_centers <- sort(kmeans_result$centers)
# The valley would be roughly in the middle of the two cluster centers
valley_point <- mean(cluster_centers)
valley_point # 2.322878
 
ggplot(rawdf, aes(x = log10(value))) +
  geom_freqpoly() +
  theme_pubr() 

ggplot(rawdf, aes(x = log10(value))) +
  geom_density() +
  theme_pubr() 


# Calculate the activity
wt <- mean(rawdf$value[rawdf$label == "WT"])
attr(rawdf, "na.action") <- NULL
attr(rawdf, "class") <- NULL
attr(rawdf, "row.names") <- NULL

specdf <- data.frame(rawdf) %>%
  dplyr::mutate(delta = log10(value)) %>%
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

pdf("output/geom_violin_significance.pdf", width = 3, height = 5)
ggplot(specdf, aes(x = as.factor(chimera_segment), y = value)) +
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) + 
  geom_violin(aes(group = as.factor(chimera_segment), color = as.factor(chimera_segment), fill = as.factor(chimera_segment)), alpha = 0.2) +
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Chimera segment") +
  ylab("Enzyme activity") +
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) 
  # geom_signif(comparisons = list(c(6,5),
  #                                c(6,4),
  #                                c(6,3),
  #                                c(6,2),
  #                                c(6,1)),
  #             y_position = 3500, tip_length = 0, vjust = 0.7,
  #             step_increase = 0.05,  # To stagger the bars
  #             textsize = 0,
  #             map_signif_level = T)
dev.off()
