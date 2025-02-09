# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr")

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
                    rep3 = rep3)
rawdf <- reshape2::melt(raw_wide)

wt <- rawdf$value[rawdf$label == "WT"] 
specdf <- rawdf %>%
  dplyr::mutate(delta = value - wt) %>%
  dplyr::filter(!label %in% c("P20", "P21", "P22", "P23", "P24")) %>%
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
                                            TRUE ~ NA))
