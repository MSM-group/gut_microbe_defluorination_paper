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
exp1 <- read_excel('data/Experiment1/20240427_rep3_4_linearized_fluoride_data.xlsx', 
                   col_names = F, sheet = "average_linearized_fluoride") %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit()

exp1dat <- data.frame(activity = exp1) %>%
  bind_cols(label = temp1) %>%
  mutate(exp = "Experiment 1") %>%
  arrange(desc(activity)) %>%
  dplyr::filter(label!= "WT")

#### EXPERIMENT 2 DILUTED ######
# Read in the protein normalized data
rep1 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("B2:M22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() #%>%
  #na.omit()
rep1

rep2 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("P2:AA22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() #%>%
  #na.omit()
rep2

rep3 <- read_excel("data/Experiment2/Fluoride Library screening 1 to 3 diluted.xlsx",
                   range = c("AD2:AO22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() # %>%
  #na.omit()
rep3

# Read in the template
temp2 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() 
temp2
temp2 <- temp2[!temp2 %in% c("P20", "P21", "P22", "P23", "P24", "WT")]
temp2<- temp2[1:240]

raw_wide <-  bind_cols(label = temp2,
                       rep1 = as.numeric(rep1),
                       rep2 = as.numeric(rep2),
                       rep3 = as.numeric(rep3)) 

rawsum <- rowMeans(raw_wide[, 2:4]) 
rawdf_diluted <- data.frame(label = temp2, exp = "Experiment 2", mean1 = rawsum) %>%
  dplyr::mutate(activity = (mean1 - min(mean1, na.rm=T))/(max(mean1,na.rm=T) - min(mean1,na.rm=T))) %>%
  na.omit()
combdf <- rawdf_diluted %>%
  bind_rows(exp1dat) %>%
  dplyr::mutate(position = as.numeric(str_extract(label, "\\-?\\d+\\.?\\d*")))

##### EXPERIMENT 2 UNDILUTED ######
# Read in the protein normalized data
rep1 <- read_excel("data/Experiment2/Fluoride Library screening undiluted.xlsx",
                   range = c("B2:M22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() #%>%
#na.omit()
rep1

rep2 <- read_excel("data/Experiment2/Fluoride Library screening undiluted.xlsx",
                   range = c("O2:Z22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() #%>%
#na.omit()
rep2

rep3 <- read_excel("data/Experiment2/Fluoride Library screening undiluted.xlsx",
                   range = c("AB2:AM22"), sheet = 1) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() # %>%
#na.omit()
rep3

# Read in the template
temp2 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() 
temp2
temp2 <- temp2[!temp2 %in% c("P20", "P21", "P22", "P23", "P24", "WT")]
temp2<- temp2[1:240]
raw_wide <-  bind_cols(label = temp2,
                       rep1 = as.numeric(rep1),
                       rep2 = as.numeric(rep2),
                       rep3 = as.numeric(rep3)) 
rawsum <- rowMeans(raw_wide[, 2:4]) 
rawdf_undiluted <- data.frame(label = temp2, exp = "Exp2 undiluted", mean1 = rawsum) %>%
  dplyr::mutate(activity = (mean1 - min(mean1, na.rm=T))/(max(mean1,na.rm=T) - min(mean1,na.rm=T))) %>%
  na.omit()
combdf2 <- rawdf_undiluted %>%
  bind_rows(combdf) %>%
  dplyr::mutate(position = as.numeric(str_extract(label, "\\-?\\d+\\.?\\d*"))) %>%
  dplyr::filter(!exp %in% ("Exp2 undiluted")) 

pdf("output/revision_figures_final/experiment1_2_comparison.pdf", width =5, height =6)
ggplot(combdf) +
  geom_bar(aes(x = position, y = activity, color = exp, fill = exp), stat = "identity") + 
  facet_grid(exp~.) +
  theme_pubr() +
  scale_color_manual(values = c("orchid4", "forestgreen", "orange3")) +
  scale_fill_manual(values = c("orchid4", "forestgreen", "orange3")) +
  xlab("Amino acid position in the P6 alanine scanning library") +
  ylab("Defluorination activity (min-max normalized)") +
  theme(legend.title = element_blank())
dev.off()
