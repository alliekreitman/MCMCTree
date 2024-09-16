###
# GOAL: Use geweke function (https://www.rdocumentation.org/packages/coda/versions/0.19-4/topics/geweke.diag) 
# to analyze convergence - T test of post-burn in data
# Author: Allie Kreitman
# Updated: 2024-06-20
####

# load libraries
library(tidyverse)
library(plotrix)
library(coda)

# user input
percent_burn_in <-  0.6 # percent of data to burn in 

#load data list
omega_txt_files <- list.files(path = "data/Data_MCTree12/", pattern = "*_omega.txt", recursive=TRUE)
mu_txt_files <- list.files(path = "data/Data_MCTree12/", pattern = "*_mu.txt", recursive=TRUE)

# ----------- CALCULATE SHORTEST LENGTH OF EACH PATIENT BY SEED ----------- 
#create an empty dataframe for the number of rows per file 
length_of_data_df <- data.frame(rows = as.numeric(), patient=as.character(), seed = as.character())

#load in data in a loop and compile into one large dataframe
for (file in omega_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = file, delim = "  ", col_names=FALSE)
  
  length_of_data_df <- rbind(length_of_data_df,
                             data.frame(rows = nrow(temp_data), file = file) %>% 
                               separate(file, into = c("data_folder", "parent_folder", "patient", "results_folder", "p", "seed")) %>% 
                               select(rows, patient, seed)
  )
  rm(temp_data)
}

#find min number of rows per patient
entries_per_patient <- length_of_data_df %>% 
  group_by(patient) %>% 
  summarise(min_rows = min(rows)) %>% #shortest seed per patient
  mutate(rows_after_burn = ceiling(min_rows*(1-percent_burn_in))) %>% #second half of data, rounding up if an odd number of rows
  mutate (nrow_for_convergence = ifelse(rows_after_burn %% 2 == 0, rows_after_burn, rows_after_burn + 1)) %>% # make convergence set even my rounding up odd number of rows after burn in
  ungroup()

# -------- LOAD AND COMBINE MU AND OMEGA FILES -------------

#create empty dataframe to load all data
combined_data_frame_omega_raw <- data.frame(mono_therapy_dnds = as.numeric(), 
                                            combo_therapy_dnds = as.numeric(), 
                                            file = as.character(), 
                                            id=as.character(),
                                            patient=as.character(),
                                            seed=as.numeric())


combined_data_frame_mu_raw <- data.frame(mutation_rate = as.numeric(), 
                                         no_idea = as.numeric(), 
                                         file = as.character(), 
                                         id=as.character(),
                                         patient=as.character(),
                                         seed=as.numeric())

# Compile Omega Text Files
for (file in omega_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = file, delim = "  ", col_names=FALSE)
  
  #clean up data so it fits with combined dataframe format
  temp_data <- temp_data %>% 
    rename(mono_therapy_dnds = X1) %>% 
    rename(combo_therapy_dnds = X2) %>% 
    mutate(file = file) %>% 
    mutate(file_copy=file) %>% 
    separate(file_copy, into = c("data_folder", "parent_folder", "patient", "results_folder", "p", "seed"))
  
  current_patient <- temp_data$patient[1]
  
  #number of data entries per patient
  current_patient_rows <- entries_per_patient %>% filter(patient==current_patient)
  
  # dumber of data for convergence 
  datapoints_to_include <- current_patient_rows$nrow_for_convergence[1]
  
  # filter dataframe for current patient, and assign burned and unburned data
  temp_data <- temp_data %>%
    mutate(id = row_number()) %>%
    mutate(burned_data = ifelse(id <= (nrow(temp_data) - datapoints_to_include), "Burned", "data_for_convergence")) %>% #use this line of code to have a burn in! (1/2)
    select(mono_therapy_dnds, combo_therapy_dnds, file, id, patient, seed, burned_data)
  
  #combine current dataset into combined_data_frame
  combined_data_frame_omega_raw <- rbind(combined_data_frame_omega_raw, temp_data)
}

# Compile Mu Text Files
for (file in mu_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = file, delim = "  ", col_names=FALSE)
  
  #clean up data so it fits with combined dataframe format
  temp_data <- temp_data %>% 
    rename(mutation_rate = X1) %>% 
    rename(no_idea = X2) %>% 
    mutate(file = file) %>% 
    mutate(file_copy=file) %>% 
    separate(file_copy, into = c("data_folder", "parent_folder", "patient", "results_folder", "p", "seed"))
  
  current_patient <- temp_data$patient[1]
  
  current_patient_rows <- entries_per_patient %>% filter(patient==current_patient)
  
  datapoints_to_include <- current_patient_rows$nrow_for_convergence[1]
  
  temp_data <- temp_data %>%
    mutate(id = row_number()) %>%
    mutate(burned_data = ifelse(id <= (nrow(temp_data) - datapoints_to_include), "Burned", "data_for_convergence")) %>% #use this line of code to have a burn in! (2/2)
    select(mutation_rate, no_idea, file, id, patient, seed, burned_data)
  
  #combine current dataset into combined_data_frame
  combined_data_frame_mu_raw <- rbind(combined_data_frame_mu_raw, temp_data)
}

data_omega_mu <- full_join(combined_data_frame_mu_raw, combined_data_frame_omega_raw, by = c("patient", "id", "seed", "burned_data"))

# remove burned
data_omega_mu <- data_omega_mu %>% filter(burned_data != "Burned") %>% 
# make a new ID that is only for the convergence dataset
  group_by(patient, seed) %>% 
  mutate(id_convergence = 1:n()) %>% 
  ungroup()

#count total number of datapoints per patient 
number_convergence_datapoints_by_patient <- data_omega_mu %>% #count total number of datapoints per patient 
  group_by(patient) %>% 
  summarise(total_convergence_datapoints = n()/3) %>% 
  ungroup()

#add total number of datapoints per patient to dataset
data_omega_mu <- left_join(data_omega_mu, number_convergence_datapoints_by_patient, by = "patient") 

# Save converged data for figure making
# write.csv(data_omega_mu, file = "Data_gewek/converged_data_from_gewek.csv")

# split the convergence data into a "A" and "B" set, which includes the first 50% of each seed then the second 50% of each seed
data_omega_mu <- data_omega_mu %>% 
  mutate(convergence_comparison_group = ifelse(id_convergence <= total_convergence_datapoints/2, "A", "B"))
  
  # confirm A and B groups are even (ues for confirmation as needed)
  # group_by(patient, seed, convergence_comparison_group) %>% 
  # summarise(n = n()) %>% 
  # ungroup()
  

#list of patients
data_patient_list <- unique(data_omega_mu$patient)

#define the A and B sets
data_groups <- c("A", "B")

for (j in 1:length(data_groups)){ # make a loop for A and B dataset
  
# set A or B group  
current_data_group <- data_groups[j]

# create DF for A gewek
gewek_stat_df <- data.frame(patient = as.character(),
                            seed = as.numeric(),
                            mono_therapy_dnds_gewek = as.numeric(),
                            combo_therapy_dnds_gewek = as.numeric(),
                            mutation_rate_gewek = as.numeric(),
                            var_4_gewek = as.numeric())

#make dataframes for each patient 
for (i in 1:length(unique(data_omega_mu$patient))){
  #assign current patient
  current_patient <- data_patient_list[i]
  
  # make object for gillman statistic
  temp_df <- data_omega_mu %>% 
    filter(patient == current_patient) %>% 
    filter(convergence_comparison_group == current_data_group)
  #list of unique seeds
  unique_seeds_temp <- unique(temp_df$seed)
  seed_temp_1 = unique_seeds_temp[1]
  seed_temp_2 = unique_seeds_temp[2]
  seed_temp_3 = unique_seeds_temp[3]
  
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[1]),
         mcmc(temp_df %>% 
                filter(seed == unique_seeds_temp[1]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea)
         )
  )
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[2]),
         mcmc(temp_df %>% 
                filter(seed == unique_seeds_temp[2]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea)
         )
  )
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[3]),
         mcmc(temp_df %>% 
                filter(seed == unique_seeds_temp[3]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea)
         )
  )
  
  
  temp_gewek_1 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[1])))
  temp_gewek_2 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[2])))
  temp_gewek_3 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[3])))
  
  gewek_stat_df <- rbind(gewek_stat_df, 
                         data.frame(patient = current_patient,
                                    seed = unique_seeds_temp[1],
                                    mono_therapy_dnds_gewek = temp_gewek_1[["z"]][["mono_therapy_dnds"]],
                                    combo_therapy_dnds_gewek = temp_gewek_1[["z"]][["combo_therapy_dnds"]],
                                    mutation_rate_gewek = temp_gewek_1[["z"]][["mutation_rate"]],
                                    var_4_gewek = temp_gewek_1[["z"]][["no_idea"]]),
                         data.frame(patient = current_patient,
                                    seed = unique_seeds_temp[2],
                                    mono_therapy_dnds_gewek = temp_gewek_2[["z"]][["mono_therapy_dnds"]],
                                    combo_therapy_dnds_gewek = temp_gewek_2[["z"]][["combo_therapy_dnds"]],
                                    mutation_rate_gewek = temp_gewek_2[["z"]][["mutation_rate"]],
                                    var_4_gewek = temp_gewek_2[["z"]][["no_idea"]]),
                         data.frame(patient = current_patient,
                                    seed = unique_seeds_temp[3],
                                    mono_therapy_dnds_gewek = temp_gewek_3[["z"]][["mono_therapy_dnds"]],
                                    combo_therapy_dnds_gewek = temp_gewek_3[["z"]][["combo_therapy_dnds"]],
                                    mutation_rate_gewek = temp_gewek_3[["z"]][["mutation_rate"]],
                                    var_4_gewek = temp_gewek_3[["z"]][["no_idea"]]))
  
  
# rename df based on current variable
assign(paste("gewek_stat_df", current_data_group, sep = "_"), gewek_stat_df)

}
}


# remove environment except for final object with gewek result
# all_objects <- ls()
# objects_to_remove <- setdiff(all_objects, c("gewek_stat_df_A", "gewek_stat_df_B", "data_omega_mu", "entries_per_patient", "number_convergence_datapoints_by_patient"))
# rm(list = objects_to_remove)
# rm(all_objects, objects_to_remove)

# --- T Test ---------------

options(scipen=999)

# T test by patient by variable --------
# make a table of p values
T_Test_Results_df <- data.frame(patient = as.character(), variable = as.character(), p_value = as.numeric(), t = as.numeric(), df = as.numeric())
patient_list <- unique(data_omega_mu$patient)
for (i in 1:length(patient_list)){
  current_patient <- patient_list[i] #current patient
  gewek_stat_df_A_patient <- gewek_stat_df_A %>% filter(patient == current_patient) #filter data for current patient in no burn dataset
  gewek_stat_df_B_patient <- gewek_stat_df_B %>% filter(patient == current_patient) #filter data for current patient in burned dataset
  
  
  # T tests
  # monotherapy
  monotherapy_ttest <- t.test(gewek_stat_df_A_patient$mono_therapy_dnds_gewek, gewek_stat_df_B_patient$mono_therapy_dnds_gewek)
  
  # combotherapy
  combotherapy_ttest <- t.test(gewek_stat_df_A_patient$combo_therapy_dnds_gewek, gewek_stat_df_B_patient$combo_therapy_dnds_gewek)
  
  # mutation
  mu_ttestt <- t.test(gewek_stat_df_A_patient$mutation_rate_gewek, gewek_stat_df_B_patient$mutation_rate_gewek)
  
  # put all the T test results into a dataframe
  T_Test_Results_df <- rbind(T_Test_Results_df, 
                             data.frame(patient = current_patient,
                                        variable = "monotherapy", 
                                        p_value = monotherapy_ttest[["p.value"]],
                                        t = monotherapy_ttest[["statistic"]][["t"]], 
                                        df = monotherapy_ttest[["parameter"]][["df"]]),
                             data.frame(patient = current_patient,
                                        variable = "combotherapy", 
                                        p_value = combotherapy_ttest[["p.value"]],
                                        t = combotherapy_ttest[["statistic"]][["t"]], 
                                        df = combotherapy_ttest[["parameter"]][["df"]]),
                             data.frame(patient = current_patient,
                                        variable = "mutation", 
                                        p_value = mu_ttestt[["p.value"]],
                                        t = mu_ttestt[["statistic"]][["t"]], 
                                        df = mu_ttestt[["parameter"]][["df"]])
  )
  
  rm(gewek_stat_df_A_patient, gewek_stat_df_B_patient)
}

#reorder cols so there is 1 row per patient
T_test_results_wider <- pivot_wider(T_Test_Results_df, names_from = variable, values_from = c(p_value, t, df)) %>% 
  select(patient, contains("mutation"), contains("monotherapy"), contains("combotherapy")) %>% 
  mutate(across(where(is.numeric), round, 3)) # shorten to 4 decimals
  
# remove the underscores in the col names
colnames(T_test_results_wider) <- gsub("_", " ", colnames(T_test_results_wider))

write.csv(T_test_results_wider, "Data_gewek/T-Test_0.6_burn_forpublication.csv", row.names=FALSE,quote = FALSE)

# write.csv(T_Test_Results_df, "Data_gewek/Gewek_Ttest_V5.3_0.54Burn.csv", row.names=FALSE)

# not_convergent_entries <- T_Test_Results_df %>% filter(p_value<0.05)
