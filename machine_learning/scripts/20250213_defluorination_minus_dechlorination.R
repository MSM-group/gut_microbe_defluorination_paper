# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "caret", 
               "tidymodels", "ranger", "tree", "seqinr", "ggseqlogo",
               "rsample", "randomForest", "readxl", "ggpubr")

# Read in the protein normalized data
dat <- read.table("data/Dechlorination-Defluorination.csv", sep = ";") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
dat
attr(dat, "na.action") <- NULL
attr(dat, "class") <- NULL

# Read in the template
temp1 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() %>%
  na.omit() 
  
temp1 <- temp1[!temp1 %in% c("P20", "P21", "P22", "P23", "P24", "WT")]
temp1
attr(temp1, "na.action") <- NULL
attr(temp1, "class") <- NULL
temp1
raw_wide <-  bind_cols(label = temp1,
                       dat = dat) 
rawdf <- reshape2::melt(raw_wide)
rawdf

# Calculate the activity relative to wild-type
# wt <- mean(rawdf$value[rawdf$label == "WT"])
attr(rawdf, "na.action") <- NULL
attr(rawdf, "class") <- NULL
attr(rawdf, "row.names") <- NULL


specdf <- data.frame(rawdf) %>%
  dplyr::mutate(delta = value) %>%
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

ggplot(specdf, aes(x = position_index, y = delta)) +
  geom_violin(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95))+
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Amino acid residue number") +
  ylab("Wild-type enzyme - alanine variant activity") 

ggplot(specdf, aes(x = as.factor(chimera_segment), y = abs(delta))) +
  geom_jitter(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment), alpha = 0.95)) +
  theme_pubr() +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  theme(legend.position = "none") +
  xlab("Chimera segment") +
  ylab("Dechlorination - defluorination- activity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_signif(comparisons = list(c(6,5),
                                 c(6,4),
                                 c(6,3),
                                 c(6,2),
                                 c(6,1)),
               tip_length = 0, vjust = 0.7,
              step_increase = 0.05,  # To stagger the bars
              textsize = 3,
              map_signif_level = T)
