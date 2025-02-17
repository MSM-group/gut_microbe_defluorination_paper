# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr", "diptest")

# Experiment1
temp1 <- read_excel("data/Experiment1/template_linearized.xlsx", col_names = F) %>%
  dplyr::select(-31) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit()
temp1 <- temp1[!is.na(temp1)]
temp1

# Read in the CSV
exp1_rep1 <- read_excel('data/Experiment1/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "rep3_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit()
exp1_rep2 <- read_excel('data/Experiment1/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                        col_names = F, sheet = "rep4_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit()

exp1_wide <-  bind_cols(label = temp1,
                       rep1 = as.numeric(exp1_rep1),
                       rep2 = as.numeric(exp1_rep2))
sd(exp1_wide$rep1[1], exp1_wide$rep2[1])

exp1dat <- exp1_wide %>%
  rowwise() %>%
  mutate(
    mean1 = mean(c_across(2:3), na.rm =T),
    sd1 = sd(c_across(2:3), na.rm =T)
  ) %>%
  ungroup() %>%
  mutate(exp = "Experiment 1") %>%
  dplyr::mutate(activity = (mean1 - min(mean1, na.rm=T))/(max(mean1,na.rm=T) - min(mean1,na.rm=T))) 

#### EXPERIMENT 2 DILUTED ######
# Read in the protein normalized data
rep1 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("B2:M22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep1

rep2 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("P2:AA22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep2

rep3 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("AD2:AO22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep3

# Read in the template
temp2 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() 
temp2
temp2 <- temp2[!temp2 %in% c("P20", "P21", "P22", "P23", "P24")]
temp2<- temp2[1:240]

raw_wide <-  bind_cols(label = temp2,
                       rep1 = as.numeric(rep1),
                       rep2 = as.numeric(rep2),
                       rep3 = as.numeric(rep3)) %>%
  
  rowwise() %>%
  mutate(
    mean2 = mean(c_across(2:4), na.rm =T),
    sd2 = sd(c_across(2:4), na.rm =T)
  ) %>%
  ungroup() %>%
  mutate(exp = "Experiment 2") %>%
  dplyr::mutate(activity = (mean2 - min(mean2, na.rm=T))/(max(mean2,na.rm=T) - min(mean2,na.rm=T))) %>%
  na.omit()

combdf <- raw_wide %>%
  bind_rows(exp1dat) %>%
  dplyr::mutate(position = as.numeric(str_extract(label, "\\-?\\d+\\.?\\d*")))

meandf <- raw_wide %>%
  left_join(exp1dat, by = "label") %>%
  dplyr::select(-contains("rep")) %>%
  dplyr::select(label, mean_exp1 = mean1, sd_exp1 = sd1, mean_exp2 = mean2, sd_exp2 = sd2)
write_csv(meandf, "output/revision_figures_final/alanine_scanning_activity_data.csv")
combdf2 <- combdf %>%
  dplyr::filter(label!= "WT") 

pdf("output/revision_figures_final/experiment1_2_comparison.pdf", width =5, height =6)
ggplot(combdf2) +
  geom_bar(aes(x = position, y = activity, color = exp, fill = exp), stat = "identity") + 
  facet_grid(exp~.) +
  theme_pubr() +
  scale_color_manual(values = c("orchid4", "forestgreen", "orange3")) +
  scale_fill_manual(values = c("orchid4", "forestgreen", "orange3")) +
  xlab("Amino acid position in the P6 alanine scanning library") +
  ylab("Defluorination activity (min-max normalized)") +
  theme(legend.title = element_blank())
dev.off()
