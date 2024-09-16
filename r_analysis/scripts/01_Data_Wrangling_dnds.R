#load libraries
library(tidyverse)
library(plotrix)
library(coda)

options(scipen = 999)

raw_data <- read.csv("data/converged_data_from_gewek.csv")

combined_data_frame <- raw_data %>% 
  mutate(combo_mono_ratio = combo_therapy_dnds/mono_therapy_dnds) %>% 
  mutate(combo_mono_diff = combo_therapy_dnds - mono_therapy_dnds) %>% 
  mutate(diffDNDS = combo_therapy_dnds / mono_therapy_dnds)

# ------------ mcmc and hpd -------------


#list of patients
data_patient_list <- unique(combined_data_frame$patient)

combine_data_hpd <- data.frame(patient = as.character(),
                               mono_therapy_mean = as.numeric(), 
                               mono_therapy_lowerhpd = as.numeric(), 
                               mono_therapy_upperhpd = as.numeric(),
                               combo_therapy_mean = as.numeric(), 
                               combo_therapy_lowerhpd = as.numeric(), 
                               combo_therapy_upperhpd = as.numeric(),
                               comboMonoDiff_mean = as.numeric(), 
                               comboMonoDiff_lowerhpd = as.numeric(), 
                               comboMonoDiff_upperhpd = as.numeric(),
                               comboMonoRatio_mean = as.numeric(), 
                               comboMonoRatio_lowerhpd = as.numeric(), 
                               comboMonoRatio_upperhpd = as.numeric(),
                               diffDNDS_mean = as.numeric(), 
                               diffDNDS_lowerhpd = as.numeric(), 
                               diffDNDS_upperhpd = as.numeric())

#make dataframes for each patient 
for (i in 1:length(unique(combined_data_frame$patient))){
  #assign current patient
  current_patient <- data_patient_list[i]
  
  # make object for gillman statistic
  temp_df <- combined_data_frame %>% 
    filter(patient == current_patient)
  #list of unique seeds
  unique_seeds_temp <- unique(temp_df$seed)
  
  assign(paste0(current_patient, "_mcmc_list"),
         mcmc.list(
           mcmc(temp_df %>% 
                  filter(seed == unique_seeds_temp[1]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds)
                ),
           mcmc(temp_df %>% 
                  filter(seed == unique_seeds_temp[2]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds)
           ),
           mcmc(temp_df %>% 
                  filter(seed == unique_seeds_temp[3]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds)
           )
         )
         )
  
  # make dataframe for graphing
  
  #make dataframe of relevant information
  temp_df <- combined_data_frame %>% 
    filter(patient == current_patient) %>% 
    select(mono_therapy_dnds, combo_therapy_dnds, combo_mono_diff, combo_mono_ratio, diffDNDS)
  
  #make dataframe into mcmc object
  temp_mcmc <- mcmc(temp_df)
  # calculate hpd interval
  temp_hpdint <- HPDinterval(temp_mcmc, prob=0.95)
  
  #add values into master dataframe structure
  combine_data_hpd_temp = data.frame(patient = current_patient,
                                     mono_therapy_mean = mean(temp_df$mono_therapy_dnds), 
                                     mono_therapy_lowerhpd = temp_hpdint["mono_therapy_dnds", "lower"], 
                                     mono_therapy_upperhpd = temp_hpdint["mono_therapy_dnds", "upper"],
                                     combo_therapy_mean = mean(temp_df$combo_therapy_dnds), 
                                     combo_therapy_lowerhpd = temp_hpdint["combo_therapy_dnds", "lower"], 
                                     combo_therapy_upperhpd = temp_hpdint["combo_therapy_dnds", "upper"],
                                     comboMonoDiff_mean = mean(temp_df$combo_mono_diff), 
                                     comboMonoDiff_lowerhpd = temp_hpdint["combo_mono_diff", "lower"], 
                                     comboMonoDiff_upperhpd = temp_hpdint["combo_mono_diff", "upper"],
                                     comboMonoRatio_mean = mean(temp_df$combo_mono_ratio), 
                                     comboMonoRatio_lowerhpd = temp_hpdint["combo_mono_ratio", "lower"], 
                                     comboMonoRatio_upperhpd = temp_hpdint["combo_mono_ratio", "upper"],
                                     diffDNDS_mean = mean(temp_df$diffDNDS), 
                                     diffDNDS_lowerhpd = temp_hpdint["diffDNDS", "lower"], 
                                     diffDNDS_upperhpd = temp_hpdint["diffDNDS", "upper"])
  
  combine_data_hpd <- rbind(combine_data_hpd, combine_data_hpd_temp)
}

rm(combine_data_hpd_temp)

#group results by patient to calculate mean and standard error 
combined_data_frame_grouped <- combined_data_frame %>% 
  group_by(patient) %>% 
  summarise(count = n(), mean_combo_mono_diff = mean(combo_mono_diff), se_combo_mono_diff = std.error(combo_mono_diff), mean_diffDNDS = mean(diffDNDS), se_diffDNDS = std.error(diffDNDS))

#load in list of therapies
therapies_by_patient <- read.csv("Data_MCTree12/Therapies_by_patients.csv")

therapies_by_patient <- therapies_by_patient %>% 
  mutate(Patient = paste("p", Patient, sep=""))

#make monotherapy and combo therapy dataframes so dnds can become one column
combined_data_frame_mono <- combined_data_frame %>% 
  select(mono_therapy_dnds, patient) %>% 
  mutate(mono_or_combo = "mono")

combined_data_frame_mono <- left_join(combined_data_frame_mono, therapies_by_patient, by = c("patient" =  "Patient"))

combined_data_frame_mono <- combined_data_frame_mono %>% 
  select(-"Combo") %>% 
  rename("therapy" = "Mono") %>% 
  rename("dnds"  = "mono_therapy_dnds")

combined_data_frame_combo <- combined_data_frame %>% 
  select(combo_therapy_dnds, patient) %>% 
  mutate(mono_or_combo = "combo")

combined_data_frame_combo <- left_join(combined_data_frame_combo, therapies_by_patient, by = c("patient" =  "Patient"))

combined_data_frame_combo <- combined_data_frame_combo %>% 
  select(-"Mono") %>% 
  rename("therapy" = "Combo") %>% 
  rename("dnds"  = "combo_therapy_dnds")

patient_dnds_dataframe <- rbind(combined_data_frame_mono, combined_data_frame_combo)

#combine patient data with therapy data
combine_data_hpd_therapy <- left_join(combine_data_hpd, therapies_by_patient %>% select(-Explicit_switch_in_genebank) %>% rename("patient" = "Patient"), by="patient")
  
#save resulting dataframes ----------
# write.csv(combined_data_frame_grouped, "Data_MCTree12/Cleaned_Data/combo_mono_mean_SE.csv", row.names=FALSE)
# write.csv(combined_data_frame, "Data_MCTree12/Cleaned_Data/combo_mono_diff_all_data.csv", row.names=FALSE)
# write.csv(patient_dnds_dataframe, "Data_MCTree12/Cleaned_Data/patient_dnds.csv", row.names=FALSE)

# remove unecessary dataframes ----------
rm(combined_data_frame)
rm(combined_data_frame_combo)
rm(combined_data_frame_grouped)
rm(combined_data_frame_mono)
rm(combined_data_frame_raw)
rm(temp_df)
rm(temp_hpdint)
rm(therapies_by_patient)
rm(current_patient)
rm(data_patient_list)
rm(file)
rm(omega_txt_files)
rm(temp_mcmc)
rm(i)

