###
# GOAL: Use geweke function (https://www.rdocumentation.org/packages/coda/versions/0.19-4/topics/geweke.diag) 
# to analyze convergence usign gewek.diag z score 
# Author: Allie Kreitman
# Updated: 2024-11-07
####

# load libraries
library(tidyverse)
library(plotrix)
library(coda)

# user input
percent_burn_in <- 0.6 # percent of data to burn in 

#load data list
omega_txt_files <- list.files(path = "data/Data_MCTree12/", pattern = "*_omega.txt", recursive=TRUE) # list of omega files
mu_txt_files <- list.files(path = "data/Data_MCTree12/", pattern = "*_mu.txt", recursive=TRUE) # list of mu files 

# ----------- CALCULATE SHORTEST LENGTH OF EACH PATIENT BY SEED ----------- 
#create an empty dataframe for the number of rows per file 
length_of_data_df <- data.frame(rows = as.numeric(), patient=as.character(), seed = as.character())

#load in data in a loop and compile into one large dataframe
for (file in omega_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = paste0("data/Data_MCTree12/", file, sep = ""), delim = "  ", col_names=FALSE, show_col_types = FALSE)
  
  #join the data to the empty data frame from above to count the number of rows for each patient and seed 
  length_of_data_df <- rbind(length_of_data_df, # join with empty/prior files
                             # format new files
                             data.frame(rows = nrow(temp_data), file = file) %>% # count the number of rows in the chain and list the file path 
                              separate(file, into = c("patient", "result", "file"), sep = "/") %>% # separate the file path to extract rows, patient, seed
                              separate(file, into = c("p", "seed", "omega"), sep = "_") %>% # separate the file path to extract rows, patient, seed
                               select(rows, patient, seed)) # select just the num rows, patient, and seed (chain)
  
  rm(temp_data) # clean up unecessary data frames
}

#find min number of rows per patient
entries_per_patient <- length_of_data_df %>% #define new dataframe
  group_by(patient) %>% # group by patient 
  summarise(min_rows = min(rows)) %>% #shortest seed per patient
  mutate(nrow_for_convergence = ceiling(min_rows*(1-percent_burn_in))) %>% #second half of data, rounding up if an odd number of rows
  # use following 2 lines (but not line above) for even number after convergence (not currently needed): 
  # mutate(rows_after_burn = ceiling(min_rows*(1-percent_burn_in))) %>% #second half of data, rounding up if an odd number of rows
  # mutate (nrow_for_convergence = ifelse(rows_after_burn %% 2 == 0, rows_after_burn, rows_after_burn + 1)) %>% # make convergence set even my rounding up odd number of rows after burn in
  ungroup()

# -------- LOAD AND COMBINE MU AND OMEGA FILES -------------

#create empty dataframe to load all data
#one for the omega files
combined_data_frame_omega_raw <- data.frame(mono_therapy_dnds = as.numeric(), 
                                            combo_therapy_dnds = as.numeric(), 
                                            file = as.character(), 
                                            id=as.character(),
                                            patient=as.character(),
                                            seed=as.numeric())

# and one for the mu files
combined_data_frame_mu_raw <- data.frame(mutation_rate = as.numeric(), 
                                         no_idea = as.numeric(), 
                                         file = as.character(), 
                                         id=as.character(),
                                         patient=as.character(),
                                         seed=as.numeric())

# Compile Omega Text Files
for (file in omega_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = paste0("data/Data_MCTree12/", file, sep = ""), delim = "  ", col_names=FALSE, show_col_types = FALSE)
  
  #clean up data so it fits with combined dataframe format
  temp_data <- temp_data %>% 
    rename(mono_therapy_dnds = X1) %>% 
    rename(combo_therapy_dnds = X2) %>% 
    mutate(file = file) %>% # add teh file path
    mutate(file_copy=file) %>% # duplicate the file path so i can break it apart
    separate(file_copy, into = c("data_folder", "parent_folder", "patient", "seed", "results_folder", "p")) # extract patient, seed from file path 
  
  current_patient <- temp_data$patient[1] # current patient can be pulled from first (or any row) of dataframe as its the same for all rows
  
  #number of data entries per patient
  current_patient_rows <- entries_per_patient %>% filter(patient==current_patient)
  
  # number of data for convergence 
  datapoints_to_include <- current_patient_rows$nrow_for_convergence[1]
  
  # filter dataframe for current patient, and assign burned and unburned data
  temp_data <- temp_data %>%
    mutate(id = row_number()) %>% # name row numbers
    # this next line does the burn in and marks early data for burn
      # by just marking the rows up to the count of 'datapoints to include' based on the burn in percent defined at top
    mutate(burned_data = ifelse(id <= (nrow(temp_data) - datapoints_to_include), "Burned", "data_for_convergence")) %>% 
    select(mono_therapy_dnds, combo_therapy_dnds, file, id, patient, seed, burned_data) # select variables 
  
  #combine current dataset into combined_data_frame
  combined_data_frame_omega_raw <- rbind(combined_data_frame_omega_raw, temp_data)
}

# Compile Mu Text Files
# this code is basically the same as the chunk above
for (file in mu_txt_files){
  #read in dataframe from list of results
  temp_data <- read_delim(file = paste0("data/Data_MCTree12/", file, sep = ""), delim = "  ", col_names=FALSE, show_col_types = FALSE)
  
  #clean up data so it fits with combined dataframe format
  temp_data <- temp_data %>% 
    rename(mutation_rate = X1) %>% 
    rename(no_idea = X2) %>% 
    mutate(file = file) %>% 
    mutate(file_copy=file) %>% 
    separate(file_copy, into = c("data_folder", "parent_folder", "patient", "seed", "results_folder", "p"))
  
  current_patient <- temp_data$patient[1]
  
  current_patient_rows <- entries_per_patient %>% filter(patient==current_patient)
  
  datapoints_to_include <- current_patient_rows$nrow_for_convergence[1]
  
  temp_data <- temp_data %>%
    mutate(id = row_number()) %>%
    mutate(burned_data = ifelse(id <= (nrow(temp_data) - datapoints_to_include), "Burned", "data_for_convergence")) %>% 
    select(mutation_rate, no_idea, file, id, patient, seed, burned_data)
  
  #combine current dataset into combined_data_frame
  combined_data_frame_mu_raw <- rbind(combined_data_frame_mu_raw, temp_data)
}
# join omega and mu data by patient/seed/burned
data_omega_mu <- full_join(combined_data_frame_mu_raw, combined_data_frame_omega_raw, by = c("patient", "id", "seed", "burned_data"))

# # make traceplot ---
# # reformat data for ggplot by pivoting longer
traceplot_data <- data_omega_mu %>%
  select(-no_idea) %>%
  pivot_longer(cols = c(mutation_rate, mono_therapy_dnds, combo_therapy_dnds), names_to = "statistic", values_to = "value") %>% # move all statistics to one column
  group_by(patient,statistic,burned_data, seed) %>% # group by patient, seed, statistic 
  mutate(separation_point = max(case_when(burned_data == "Burned" ~ id)))%>% # add separation point column that is when the burn in ends
  ungroup() %>% 
  filter(patient != "p156") %>%  # remove patient that we are not analyzing
  fill(separation_point, .direction = "down") # fill in NA values from burn in point
#
# #make ggplot for each patient
#list of patients
data_patient_list <- unique(data_omega_mu$patient)
data_patient_list <- data_patient_list[ !data_patient_list == 'p156']

for (i in 1:length(data_patient_list)){ # loop through each patient

  # make traceplot:
  curr_plot <- ggplot(data = traceplot_data %>% filter(patient == data_patient_list[i]))+ # plot from data filtered for just current patient
    theme_bw()+ # looks clean
    geom_line(aes(x = id, y = value), size = 1)+ # plot line of statistic values over time (id)
    geom_vline(aes(xintercept = separation_point), color = "blue")+ # add vertical line for each point where it goes to only covergence data
    facet_grid(statistic~seed, scales = "free")+ # facet by patient and chain (seed)
    labs(title = (paste0("Traceplot for ", data_patient_list[i])))

  # save traceplot for each patient
  assign(paste0("traceplot_", data_patient_list[i]), curr_plot)
}

# # save all the traceplots
traceplot_list <- mget(ls(pattern = "^traceplot_p")) # List all variables that start with "traceplot_p"
pdf("data/traceplots/all_traceplots.pdf", width = 8, height = 6) # Open a PDF device to save all plots into one file
lapply(traceplot_list, print) # Loop through each plot and print it to the PDF
dev.off() # Close the PDF device


# format data for gewek.diag --- 
# remove burned
data_omega_mu <- data_omega_mu %>% filter(burned_data != "Burned") %>% 
# make a new ID that is only for the convergence dataset
  group_by(patient, seed) %>% 
  mutate(id_convergence = 1:n()) %>% 
  filter(id_convergence %% 2 != 0) %>%  # thin mcmc chains by selecting only odd id (every other)
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

# ------------------------ big changes
# define fractions for gewek.diag
f1 <- 0.33 # fraction from beginning of chain
f2 <- 0.33 # fraction from the end of the chain

# make an empty dataframe for gewek outputs 
gewek_stat_df <- data.frame(patient = as.character(), # patient number
                            seed = as.numeric(), # seed of the 3 runs
                            mono_therapy_dnds_gewek = as.numeric(),
                            combo_therapy_dnds_gewek = as.numeric(),
                            mutation_rate_gewek = as.numeric(),
                            var_4_gewek = as.numeric())

# make a loop that runs through each patient to format data then run the gewek.diag, then add it to a growing df
for (i in 1:length(unique(data_omega_mu$patient))){ # go iteratively through each patient
  #assign current patient
  current_patient <- data_patient_list[i] # current patient is ith in the list of patients
  
  # format data 
  temp_df <- data_omega_mu %>% 
    filter(patient == current_patient)  # remove data from other patients
    
  #list of unique seeds
  unique_seeds_temp <- unique(temp_df$seed) # list 3 unique seeds, then define them as seed 1, 2, 3
  seed_temp_1 = unique_seeds_temp[1]
  seed_temp_2 = unique_seeds_temp[2]
  seed_temp_3 = unique_seeds_temp[3]
  
  # make mcmc object for seed # 1
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[1]), # this line makes the name
         mcmc(temp_df %>% # pull from dataframe with only current patient data
                filter(seed == unique_seeds_temp[1]) %>% # use only seed 1
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea) # select data rows
         )
  )
  # make mcmc object for seed # 2
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[2]),
         mcmc(temp_df %>% 
                filter(seed == unique_seeds_temp[2]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea)
         )
  )
  # make mcmc object for seed # 3
  assign(paste0(current_patient, "_mcmc", unique_seeds_temp[3]),
         mcmc(temp_df %>% 
                filter(seed == unique_seeds_temp[3]) %>% 
                select(mono_therapy_dnds, combo_therapy_dnds, mutation_rate, no_idea)
         )
  )
  
  # run the gewek.diag functions for each seed of the patient
  # use the fractions f1 and f2 defined above
  # input is the mcmc object defined in the chunk above 
  temp_gewek_1 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[1])), frac1 = f1, frac2 = f2) # seed 1
  temp_gewek_2 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[2])), frac1 = f1, frac2 = f2) # seed 2
  temp_gewek_3 <- geweke.diag(get(paste0(current_patient, "_mcmc", unique_seeds_temp[3])), frac1 = f1, frac2 = f2) # seed 3
  
  # combine output of gewek.diag into a big dataframe
  gewek_stat_df <- rbind(gewek_stat_df, 
                         data.frame(patient = current_patient, # patient
                                    seed = unique_seeds_temp[1], # chain/seed number
                                    mono_therapy_dnds_gewek = temp_gewek_1[["z"]][["mono_therapy_dnds"]], # z score for mono_therapy 
                                    combo_therapy_dnds_gewek = temp_gewek_1[["z"]][["combo_therapy_dnds"]], # z score for combo_therapy
                                    mutation_rate_gewek = temp_gewek_1[["z"]][["mutation_rate"]], # z score for mutation rate
                                    var_4_gewek = temp_gewek_1[["z"]][["no_idea"]]), # z score for IDK?
                         data.frame(patient = current_patient,# patient
                                    seed = unique_seeds_temp[2], # chain/seed number
                                    mono_therapy_dnds_gewek = temp_gewek_2[["z"]][["mono_therapy_dnds"]], # z score for mono_therapy 
                                    combo_therapy_dnds_gewek = temp_gewek_2[["z"]][["combo_therapy_dnds"]], # z score for combo_therapy
                                    mutation_rate_gewek = temp_gewek_2[["z"]][["mutation_rate"]], # z score for mutation rate
                                    var_4_gewek = temp_gewek_2[["z"]][["no_idea"]]), # z score for IDK?
                         data.frame(patient = current_patient,# patient
                                    seed = unique_seeds_temp[3], # chain/seed number
                                    mono_therapy_dnds_gewek = temp_gewek_3[["z"]][["mono_therapy_dnds"]], # z score for mono_therapy 
                                    combo_therapy_dnds_gewek = temp_gewek_3[["z"]][["combo_therapy_dnds"]], # z score for combo_therapy
                                    mutation_rate_gewek = temp_gewek_3[["z"]][["mutation_rate"]], # z score for mutation rate
                                    var_4_gewek = temp_gewek_3[["z"]][["no_idea"]])) # z score for IDK?
  
}
# reformat data
gewek_stat_df <- gewek_stat_df %>% 
  select(-var_4_gewek) %>% # remove that random variable we are not using
  pivot_longer(!c(patient, seed), names_to = "var", values_to = "z") %>% # pivot longer to one row per var 
  mutate(converged_ys = ifelse(z < -1.96 | z > 1.96, "no", "yes")) %>% # make a y/n column for whether it is outside the range of convergence 
  mutate(frac1 = f1, frac2 = f2, burn = percent_burn_in) %>%  # save variables of convergence analysis
  filter(patient != "p156") # remove patient not included in figures

# save output
write.csv(gewek_stat_df, "data/20241107_Convergence_outputs/convergence_f1_0.33_f2_0.33_burn_0.33_thinned.csv", row.names = FALSE)





