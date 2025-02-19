# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr")

### EXPERIMENT 1 #####
temp1 <- read_excel("data/Experiment1/template_exp1.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit()
temp1 <- temp1[!is.na(temp1)]
temp1

# Read in the CSV
exp1_rep1 <- read_excel('data/Experiment1/20240427_rep1_2_linearized_fluoride_data.xlsx', 
                        col_names = F, sheet = "rep1_linearized_fluoride", range = c("A1:AD8")) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit() 
exp1_rep2 <- read_excel('data/Experiment1/20240427_rep1_2_linearized_fluoride_data.xlsx', 
                        col_names = F, sheet = "rep2_linearized_fluoride", range = c("A1:AD8")) %>%
  janitor::clean_names() %>%
  as.matrix(bycol = T) %>%
  t %>%
  as.vector() %>%
  na.omit() 

exp1_wide <-  bind_cols(label = temp1,
                        rep1 = as.numeric(exp1_rep1),
                        rep2 = as.numeric(exp1_rep2)) %>%
  dplyr::filter(label != "WT")

exp1dat <- exp1_wide %>%
  rowwise() %>%
  mutate(
    mean1 = mean(c_across(2:3), na.rm =T)
  ) %>%
  ungroup() %>%
  mutate(exp = "Experiment 1") %>%
  dplyr::mutate(activity1 = (rep1 - min(mean1, na.rm=T))/(max(mean1,na.rm=T) - min(mean1,na.rm=T)),
                activity2 = (rep2 - min(mean1, na.rm=T))/(max(mean1,na.rm=T) - min(mean1,na.rm=T))) %>%
  rowwise() %>%
  dplyr::mutate(
    mean_activity1 = mean(c_across(activity1:activity2), na.rm =T),
    sd_activity1 = sd(c_across(activity1:activity2), na.rm =T))

#### EXPERIMENT 2 #####
# Averages to check
# Read in the protein normalized data
dat <- read_excel("data/Experiment2/Fluoride_concentrations+normalized_activities.xlsx",
                  range = c("C135:N156"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
temp2 <- read_excel("data/Experiment2/template_exp2.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() 
temp2 <- temp2[1:242]
checkdata <- bind_cols(dat, temp2)
attr(dat, "na.action") <- NULL
attr(dat, "class") <- NULL

#### EXPERIMENT 2 ######
# Read in the individual replicates
rep1 <- read_excel("data/Experiment2/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("C110:N130"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep1

rep2 <- read_excel("data/Experiment2/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("P110:AA130"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep2

rep3 <- read_excel("data/Experiment2/Fluoride_concentrations+normalized_activities.xlsx",
                   range = c("AC110:AN130"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() 
rep3

# Read in the template
temp2 <- read_excel("data/Experiment2/template_exp2.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() 
temp2 <- temp2[1:length(rep1)]
temp2


exp2dat <-  bind_cols(label = temp2,
                        rep1 = as.numeric(rep1),
                        rep2 = as.numeric(rep2),
                        rep3 = as.numeric(rep3)) %>%
  dplyr::filter(!label %in% c("WT", "P20", "P21", "P22", "P23", "P24")) %>%
  rowwise() %>%
  dplyr::mutate(mean2 = mean(c_across(2:4), na.rm =T)) %>%
  ungroup() %>%
  dplyr::mutate(activity1 = (rep1 - min(mean2, na.rm=T))/(max(mean2,na.rm=T) - min(mean2,na.rm=T)),
                  activity2 = (rep2 - min(mean2, na.rm=T))/(max(mean2,na.rm=T) - min(mean2,na.rm=T)),
                  activity3 = (rep3 - min(mean2, na.rm=T))/(max(mean2,na.rm=T) - min(mean2,na.rm=T))) %>%
  rowwise() %>%
  dplyr::mutate(
    mean_activity2 = mean(c_across(6:8), na.rm =T),
    sd_activity2 = sd(c_across(6:8), na.rm =T)) %>%
  ungroup() %>% 
  mutate(exp = "Experiment 2") %>%
  dplyr::filter(!duplicated(label)) %>%
  dplyr::filter(!is.na(mean2))

combdf <- exp2dat %>%
  left_join(exp1dat, by = "label") %>%
  dplyr::mutate(position = as.numeric(str_extract(label, "\\-?\\d+\\.?\\d*")))

meandf <- combdf %>%
  dplyr::select(-contains("rep")) %>%
  dplyr::select(label, Experiment1 = mean_activity1, sd_exp1 = sd_activity1, Experiment2 = mean_activity2, sd_exp2 = sd_activity2)

exp2_data_final <- meandf %>%
  dplyr::select(label, activity = Experiment2)
write_csv(exp2_data_final, "data/Experiment2/20250220_rep1_2_3_linearized_fluoride_data.csv")

# write_csv(meandf, "output/revision_figures_final/alanine_scanning_activity_diluted_data.csv")

combdf2 <- meandf %>%
  dplyr::select(-contains("sd")) %>%
  reshape2::melt(id.var = "label") %>%
  dplyr::mutate(Position = as.numeric(str_extract(label, "\\-?\\d+\\.?\\d*")))
colnames(combdf2) <- c("Label", "Experiment", "Activity", "Position")

pdf("output/revision_figures_final/experiment1_2_comparison.pdf", width =5, height =6)
ggplot(combdf2) +
  geom_bar(aes(x =Position, y = Activity, color = Experiment, fill = Experiment), stat = "identity") + 
  facet_grid(Experiment~.) +
  theme_pubr() +
  scale_color_manual(values = c("orchid4", "forestgreen")) +
  scale_fill_manual(values = c("orchid4", "forestgreen")) +
  xlab("Amino acid position in the P6 alanine scanning library") +
  ylab("Defluorination activity (min-max normalized)") +
  theme(legend.title = element_blank())
dev.off()
