# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr")

# Read in the protein normalized data
dat <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                  range = c("C135:N156"), sheet = "Data normalized Defluorination") %>%
  as.matrix(byrow = T) %>%
  t() %>%
  c() %>%
  na.omit()
dat
attr(dat, "na.action") <- NULL
attr(dat, "class") <- NULL
