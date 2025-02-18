# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "ggpmisc", "caret", "gplots",
               "colorspace", "cowplot", "tidymodels", "ranger", "tree",
               "rsample", "randomForest", "readxl", "ggpubr", "RColorBrewer")

# Read in the template
temp1 <- read_excel("data/machine_learning/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
temp1 <- temp1[!is.na(temp1)]

# Read in the CSV
res3 <- read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                  col_names = F, sheet = "rep3_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res3 <- res3[!is.na(res3)]

res4 <- read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "rep4_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 
res4 <- as.numeric(res4[!is.na(res4)])

res5 <-  read_excel('data/machine_learning/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                    col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() 

res5 <- as.numeric(res4[!is.na(res5)])

# Delta calculation relative to WT
wt3 <- res3[temp1 == "WT"]
wt4 <- as.numeric(res4[temp1 == "WT"])
wt5 <- as.numeric(res5[temp1 == "WT"])

# Delta
delta3 <- wt3 - res3 # delta between WT -> Ala_mut
delta4 <- wt4 - res4
delta5 <- wt5 - res5

delta <- c(delta3, delta4)
label <- c(temp1, temp1)
delta <- delta5
label <- temp1
# Delta
dat <- bind_cols(label = label, 
                 delta = delta) %>%
  dplyr::filter(label != "WT") %>%
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
                                            position_index %in% 151:196 ~ 5,
                                            position_index %in% 197:238 ~ 6,
                                            TRUE ~ NA))

dat$chimera_segment
whichkeep <- dat$delta[dat$chimera_segment == "6"]
whichkeep <- whichkeep[!is.na(whichkeep)]
length(whichkeep[whichkeep >= 0]) / length(whichkeep)

p1 <-  ggplot(dat, aes(x = position_index, y = delta)) +
  geom_point(aes(color = as.factor(chimera_segment), 
                 fill = as.factor(chimera_segment),
                 alpha = 0.9),
             size = 2)+
  theme_pubr() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  xlab("Protein position") +
  ylab("Enzyme activity: wild-type - variant")+
  geom_line(aes(y=0), color = "gray20", linetype = "dashed", linewidth = 0.5) 
p1  


ggsave(p1, file = "output/dot_plot.png", height = 4, width = 3)


