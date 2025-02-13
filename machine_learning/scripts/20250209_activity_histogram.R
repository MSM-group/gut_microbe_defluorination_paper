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

pdf("output/revision_figures_final/log10_activity_histogram.pdf", height = 5, width = 5)
ggplot(rawdf, aes(x = log10(value))) +
  geom_histogram(binwidth = 0.1, fill = "gray80") +
  theme_pubr() +
  geom_vline(aes(xintercept = log10(830.7176)), size = 2, color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = 2.322878), size = 2, color = "blue", linetype = "dashed") +
  ylab("count in alanine scanning library") +
  xlab("log10 of fluoride concentration / \n protein concentration in supernatant \n[µM*mg^(-1)*ml]")
dev.off()


dip_result <- dip(log10(rawdf$value), full.result = TRUE)
dip_test <- dip.test(log10(rawdf$value))
kmeans_result <- kmeans(log10(rawdf$value), centers = 2)

# Find the centers (means) of the two clusters using k-means
cluster_centers <- sort(kmeans_result$centers)
# The valley point is the middle of the two cluster centers
valley_point <- mean(cluster_centers)
valley_point # 2.322878
 
