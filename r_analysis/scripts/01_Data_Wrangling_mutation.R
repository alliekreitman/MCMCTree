#load libraries
library(tidyverse)
library(plotrix)

#make into scientific notation
options(scipen = -100)

raw_data <- read.csv("data/converged_data_from_gewek.csv")

combined_data_frame <- raw_data %>% select(patient, mutation_rate)

# make a look to split each df up and make it a mcmc object, then put into a new dataframe ----
#list of patients
data_patient_list <- unique(combined_data_frame$patient)

combine_data_hpd <- data.frame(patient = as.character(),
                               mutation_mean = as.numeric(), 
                               mutation_lowerhpd = as.numeric(), 
                               mutation_upperhpd = as.numeric())

#make dataframes for each patient 
for (i in 1:length(unique(combined_data_frame$patient))){
  #assign current patient
  current_patient <- data_patient_list[i]
  
  #make dataframe of relevant information
  temp_df <- combined_data_frame %>% 
    filter(patient == current_patient) %>% 
    select(mutation_rate)
  
  #make dataframe into mcmc object
  temp_mcmc <- mcmc(temp_df)
  # calculate hpd interval
  temp_hpdint <- HPDinterval(temp_mcmc, prob=0.95)
  
  #add values into master dataframe structure
  combine_data_hpd_temp = data.frame(patient = current_patient,
                                     mutation_mean = mean(temp_df$mutation_rate), 
                                     mutation_lowerhpd = temp_hpdint["mutation_rate", "lower"], 
                                     mutation_upperhpd = temp_hpdint["mutation_rate", "upper"])
  
  combine_data_hpd <- rbind(combine_data_hpd, combine_data_hpd_temp)
  
  #remove temp data
  rm(temp_df)
  rm(temp_hpdint)
  rm(temp_mcmc)
  rm(current_patient)
  rm(combine_data_hpd_temp)
}

#save dataframe
# write.csv(combined_data_frame, "Data_MCTree12/Cleaned_Data/Clean_mutation_data.csv", row.names = FALSE)

rm(combined_data_frame_raw)

