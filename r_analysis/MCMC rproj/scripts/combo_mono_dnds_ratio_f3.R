#load libraries
library(tidyverse)
library(coda)

#check working directory
getwd()

#run wrangling script to get dataframe
source("scripts/01_Data_Wrangling_dnds.R")

patients_to_use <- c("p1", "p100", "p106", "p126", "p132", "p1", "p60", "p21", "p26", "p28", "p55", "p69")

fig3_ratio <- ggplot(combine_data_hpd_therapy %>% filter(patient %in% patients_to_use), aes(x=patient, y=diffDNDS_mean)) +
  theme_bw()+
  geom_dotplot(binaxis="y", stackdir='center', dotsize=0.5)+
  geom_errorbar(aes(ymin = diffDNDS_lowerhpd, ymax = diffDNDS_upperhpd), width = 0.15, linewidth = 0.5)+
  labs(y = "Combo dN/dS / Mono dN/dS")
  
ggsave("Figures/Fig3_V9.1_hpd_ratio_20240913.pdf", fig3_ratio, height=11, width=8.5, unit="in")
