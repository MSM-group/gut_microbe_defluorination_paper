# Install packages
# install.packages("pacman")
pacman::p_load(nlstools, 
               data.table,
               ggplot2,
               drc,
               rvg,
               officer,
               stringr,
               readxl,
               tidyverse,
               ggpubr)

# Read in the data and calculate rates
dat <- read_excel("kinetics/20250218_kinetic_raw_data.xlsx") %>%
        janitor::clean_names() %>%
        dplyr::mutate(FAc_mM = as.numeric(gsub("mM", "", word(ident, sep = "_", 3)))) %>%
        dplyr::slice(1:nrow(.)) %>%
        dplyr::mutate(timepoint_min = as.numeric(gsub("min", "", word(ident, sep = "_", 4)))) %>%
        dplyr::mutate(fluoride_mM = as.numeric(anion_fluoride_m_m)) %>%
        dplyr::mutate(id_rep = word(ident, sep = "_", -1)) %>%
        dplyr::filter(!grepl("Kinetik_Protein3_50mM_120min|Kinetik_Protein3_20mM_120min", ident))

# Plot the results for the different concentrations
ggplot(dat) +
  geom_point(aes(x = timepoint_min, y = fluoride_mM, 
                 group = FAc_mM, color = FAc_mM)) +
  theme_pubr()

result <- dat %>%
  group_by(replicate, FAc_mM) %>%
  summarize(
    model = list(lm(fluoride_mM ~  timepoint_min, data = cur_data())),  # Fit model within each group
    .groups = "drop"
  ) %>%
  mutate(
    slope = sapply(model, function(m) coef(m)[2])  # Extract slope for each model
  )  %>%
  dplyr::select(-model)
colnames(result)


# Plot the results for the different concentrations
ggplot(result)  +
  geom_point(aes(x = FAc_mM, y = slope, color = replicate)) +
  theme_pubr()

# Our reactions had a protein concentration of 0.5 mg / mL in a Tris-Buffer pH 8.5 
mich <- result %>%
    mutate(v = slope/60) # mM per second

# Plot the results for the different concentrations
pdf("kinetics/rate_per_conc.pdf", width = 4, height = 3)
ggplot(mich)  +
  geom_point(aes(x = FAc_mM, y = v, color = replicate)) +
  xlab("conc FAc (mM)") +
  ylab("Rate / mM [F-] per minute") + 
  theme_pubr()
dev.off()

trim_dat <- mich %>%
  dplyr::mutate(S = FAc_mM) %>% 
  dplyr::mutate(enzymes = "Protein3") %>%
  dplyr::select(v, S, enzymes)

# See Michaelis script
source("kinetics/michaelis.r")
raw <- michaelisRaw(trim_dat, filname = paste0(getwd(), "/kinetics/"))
raw

# Bootstrapping
mm <- trim_dat
boot1 <- michaelisBoot(mm = mm, filname = paste0(getwd(),  "/kinetics/"))
boot1
model.drm <- drm(mm$v ~ mm$S, data = mm, fct = MM.2())
mml <- data.frame(S = seq(0, max(mm$S), length.out = 100))
mml$v <- predict(model.drm, newdata = mml)
mml

michaelis <- formula(v~Vmax*S/(Km+S)) 
bestfit <- nls(michaelis, mm, start=list(Vmax=1,Km=50))
bestfit

row1 <- c(boot1$enzymes, boot1$Km, 
          boot1$Km.std, boot1$vmax, 
          boot1$vmax.std)

dtf <- data.frame(t(row1))
dtf
colnames(dtf) <- c("enzymes", "Km", "Km.std",
                   "Vmax", "Vmax.std")
dtf

dtf$Vmax <- as.numeric(dtf$Vmax) 
dtf$Km <- as.numeric(dtf$Km)
dtf$Vmax.std <- as.numeric(dtf$Vmax.std) 
dtf$Km.std <- as.numeric(dtf$Km.std) 
dtf$Vmax

dtfv <- dtf %>%
  mutate(kcat = (Vmax)/((0.5 / 24705.22) * 1000)) %>% #  24705.22 is the enzyme molecular weight in Da
  mutate(kcat.std = (Vmax.std)/((0.5 / 24705.22) * 1000)) %>% 
  mutate(Km.M = Km/1e3) %>% # convert from mM to M (divide by 1e3)
  mutate(Km.M.std = Km.std/1e3) %>% 
  mutate(cat.effic = (kcat)/Km.M)
# Standard error propagation when dividing two numbers with standard errors
dtfv$cat.effic.std <- (dtfv$cat.effic) * sqrt((dtfv$kcat.std^2) + (dtfv$Km.M.std^2)) 
dtfv

write.table(dtfv, "kinetics/Protein3_table_normalized.tsv",
            sep = "\t", quote = F, row.names = F)

