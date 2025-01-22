# Read in the packages
pacman::p_load("tidyverse", "readxl", "ggpubr")

# Read in the protein normalized data
dat <- read_excel("data/Fluoride_concentrations+normalized_activities.xlsx",
                  range = c("C131:N152")) %>%
      as.matrix(byrow = T) %>%
      as.vector() %>%
      na.omit()
dat

# Read in the template
temp1 <- read_excel("data/template_raw.xlsx", col_names = F) %>%
  janitor::clean_names() %>%
  as.matrix(byrow = T) %>%
  as.vector() %>%
  na.omit()
temp1

# Assign the values to the template
specdf <- bind_cols(label = temp1, 
                 delta = dat) %>%
  dplyr::filter(!label %in% c("WT", "P20", "P21", "P22", "P23", "P24")) %>%
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
  dplyr::arrange(desc(delta))
specdf$label <- factor(specdf$label, levels = as.character(specdf$label))

# Make a barplot
pdf("output/alanine_variant_activity_barplot.pdf", height= 40)
ggplot(data = specdf, aes(x = label, y = delta-1)) +
  geom_bar(stat = "identity") +
  theme_pubr() +
  coord_flip()
dev.off()

# Make a segment plot
p1 <-  ggplot(specdf, aes(x = position_index, y = delta)) +
  geom_point(aes(color = as.factor(chimera_segment), fill = as.factor(chimera_segment),alpha = 0.95)) +
  theme_pubr() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  scale_color_manual(values = c("yellow4","darkblue", "magenta", "limegreen","gray60", "maroon")) + 
  xlab("Protein position") +
  ylab("Enzyme activity normalized by wild-type")+
  geom_line(aes(y=1), color = "gray20", linetype = "dashed", linewidth = 0.5)
p1  
