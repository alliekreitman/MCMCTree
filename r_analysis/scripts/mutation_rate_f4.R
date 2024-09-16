#load libraries
library(tidyverse)
library(plotrix)

#check working directory
getwd()

#turn on scientific notation
options(scipen = -100)

#run wrangling script to get dataframe
source("scripts/01_Data_Wrangling_mutation.R")

patients_to_use <- c("p1", "p100", "p106", "p126", "p132", "p1", "p60", "p21", "p26", "p28", "p55", "p69")

data <- combine_data_hpd %>% 
  filter(patient %in% patients_to_use)

fig4 <- ggplot(data, aes(x=patient, y=mutation_mean)) + 
  theme_bw()+
  geom_dotplot(binaxis="y", stackdir='center', dotsize=0.5)+
  geom_errorbar(aes(ymin = mutation_lowerhpd, ymax = mutation_upperhpd), width = 0.15, size = 0.5)+
  labs(y = "Mutation Rate", x="Patient")
fig4

# ggsave("Figures/Fig4_V7.0_hpd.pdf", fig4, height=11, width=8.5, unit="in")
