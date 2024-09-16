#load libraries
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(coda)

#check working directory
getwd()

#run wrangling script to get dataframe
source("scripts/01_Data_Wrangling_dnds.R")
rm(patient_dnds_dataframe)

# little bit of data wrangling to split by therapy 
mono_therapy_df <- combine_data_hpd_therapy %>% 
  select(patient, mono_therapy_mean, mono_therapy_lowerhpd, mono_therapy_upperhpd, Mono) %>% 
  rename(dnds_mean = mono_therapy_mean) %>% 
  rename(dnds_lowerhpd = mono_therapy_lowerhpd) %>% 
  rename(dnds_upperhpd = mono_therapy_upperhpd) %>% 
  rename(therapy = Mono) %>% 
  mutate(mono_combo = "mono")

combo_therapy_df <- combine_data_hpd_therapy %>% 
  select(patient, combo_therapy_mean, combo_therapy_lowerhpd, combo_therapy_upperhpd, Combo) %>% 
  rename(dnds_mean = combo_therapy_mean) %>% 
  rename(dnds_lowerhpd = combo_therapy_lowerhpd) %>% 
  rename(dnds_upperhpd = combo_therapy_upperhpd) %>% 
  rename(therapy = Combo) %>% 
  mutate(mono_combo = "combo")

data <- rbind(mono_therapy_df, combo_therapy_df)

rm(mono_therapy_df)
rm(combo_therapy_df)

patients <- unique(data$patient)

figure_list <- c()

for (i in 1:length(patients)){
  patient_current <- patients[i]
  
  
  data_patient_temp <- data %>% 
    filter(patient == patient_current)
  
  assign(paste(patient_current, sep="_", "figure"), 
         ggplot(data = data_patient_temp, aes(x = therapy, y = dnds_mean))+
           theme_bw()+
           geom_dotplot(binaxis="y", stackdir='center', dotsize=1)+
           geom_errorbar(aes(ymin = dnds_lowerhpd, ymax = dnds_upperhpd), width = 0.15, size = 0.5)+
           theme(axis.text = element_text(size = 8))+
           xlab("")+
           ylab("")+
           facet_wrap(~patient)
  )
  
  figure_list <- c(figure_list, paste(patient_current, sep="_", "figure"))
}

print(figure_list)

#combine figures
plot <- cowplot::plot_grid(p100_figure, p106_figure, p126_figure, p132_figure, p1_figure, p60_figure, 
                           p21_figure, p26_figure,p28_figure, p55_figure, p69_figure)
x.grob <- textGrob("Therapy", gp=gpar(fontsize=12))
y.grob <- textGrob("dN/dS", gp=gpar(fontsize=12))

plot <- plot_grid(arrangeGrob(plot, bottom = x.grob, left = y.grob))


# ggsave("Figures/Fig2_v7.1_hpd.pdf", plot, height = 8.5, width = 11, unit="in")

