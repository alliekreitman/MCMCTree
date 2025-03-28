######################
# Author: Allie Kreitman
# Script Function: Use API to access genebank to collection information needed for init files
######################

# load packages ------
library(tidyverse)
# install.packages("rentrez")
library(rentrez) # package specifically designed to interact with NCBI Entrez databases, including GenBank

# --- use this code for making init files 

# list of patients in the study
patients_in_study <- c("patient100", "patient106", "patient126", "patient132", "patient1", "patient60", "patient21", "patient26", "patient28", "patient55", "patient69")

### read in file of accession IDs and parse ------

# Read in data file
text_data <- read_file("data/record.txt")

# Extract patient entries
lines <- str_split(text_data, "\n")[[1]]

# Initialize storage
data_list <- list()
current_patient <- NULL
current_therapy <- NULL

# Process each line
for (line in lines) {
  line <- str_trim(line)
  
  # Identify patient lines
  if (str_detect(line, "^patient\\d+:")) {
    current_patient <- str_extract(line, "patient\\d+")
    
    # Identify therapy type
  } else if (str_detect(line, "mono-therapy:|como-therapy:")) {
    current_therapy <- ifelse(str_detect(line, "mono-therapy"), "mono-therapy", "combo-therapy")
    ids <- str_extract_all(line, "AY\\d+")[[1]]
    
    # Store extracted data
    if (length(ids) > 0) {
      data_list <- append(data_list, list(tibble(patient = current_patient, therapy = current_therapy, id = ids)))
    }
    
    # Capture additional combo-therapy IDs
  } else if (str_detect(line, "\\{")) {
    ids <- str_extract_all(line, "AY\\d+")[[1]]
    
    if (length(ids) > 0) {
      data_list <- append(data_list, list(tibble(patient = current_patient, therapy = "combo-therapy", id = ids)))
    }
  }
}

# Combine all extracted data into a single table
genebank_data_df <- bind_rows(data_list) %>% 
  mutate(study_day_collected = NA) %>%  # make empty column for study day collected
  filter(patient %in% patients_in_study) # remove patients that were not included in the study 
  
### loop of accessing genebank data for each id -----

# list of patients
id_list <- unique(genebank_data_df$id)

# make a loop that cycles through each patient data
for (i in 1:length(id_list)){
  
  # Retrieve Full GenBank Record for ith id value
  sequence_gb <- entrez_fetch(db="nucleotide", id=id_list[i], rettype="gb", retmode="text")
  # cat(sequence_gb) # print all metadata for the sample
  
  # parse the study day from the comment 
  study_day <- str_extract(sequence_gb, "study day \\d+")
  
  # update study day collected into dataframe 
  genebank_data_df <- genebank_data_df %>% 
    mutate(study_day_collected = ifelse(id == id_list[i], study_day, study_day_collected))
}

# make study day a numeric column
genebank_data_df <- genebank_data_df %>% 
  mutate(study_day_collected = gsub("study day ", "", study_day_collected)) %>% # make study day just a number
  mutate(study_day_collected = as.numeric(study_day_collected))

### patient 26 ---- 
# patient 26 was not included in the list of samples, so the below script will use the fasta file to rebuild the list of samples and collectiond dates

# Read in data file
p26_fasta_file <- read_file("data/p26sequence.fasta")

# extract list of ids
p26_ids <- stringr::str_extract_all(p26_fasta_file, ">AY\\d+")
p26_ids <- gsub(">", "", unlist(p26_ids))  # Remove the ">" symbol

genebank_data_df_p26 <- data.frame("id" = p26_ids) %>% 
  mutate(patient = "patient26") %>% 
  mutate(therapy = NA, study_day_collected = NA)

# make a loop that cycles through each patient data
for (i in 1:length(p26_ids)){
  
  # Retrieve Full GenBank Record for ith id value
  sequence_gb <- entrez_fetch(db="nucleotide", id=p26_ids[i], rettype="gb", retmode="text")
  # cat(sequence_gb) # print all metadata for the sample
  
  # parse the study day from the comment 
  study_day <- str_extract(sequence_gb, "study day \\d+")
  
  
  # update study day collected into dataframe 
  genebank_data_df_p26 <- genebank_data_df_p26 %>% 
    mutate(study_day_collected = ifelse(id == p26_ids[i], study_day, study_day_collected))
}

# monotherapy duration was 14 days from genebank
genebank_data_df_p26 <- genebank_data_df_p26 %>% 
  mutate(study_day_collected = gsub("study day ", "", study_day_collected)) %>% # make study day just a number
  mutate(study_day_collected = as.numeric(study_day_collected)) %>% 
  # assign monotherapy if within first 2 weeks, then combotherapy for later than 2 weeks of study
  mutate(therapy = ifelse(study_day_collected > 14, "combo-therapy", "mono-therapy"))


# combine p26 with the rest of the data
genebank_data_df <- rbind(genebank_data_df, genebank_data_df_p26)

## format for init files -----
genebank_data_df_summarised <- genebank_data_df %>% 
  group_by(patient, study_day_collected) %>% 
  summarise(datanum = n()) %>% 
  ungroup() %>% 
  rename(datatime = study_day_collected)

# print the datanum and datatime lists per patient 
for (i in 1:length(patients_in_study)){
  curr_data <- genebank_data_df %>% 
    filter(patient == patients_in_study[i]) %>%  # make curr data frame with one patient
    arrange(study_day_collected) # sort data by patient and by study day 
  
  # list number of unique collection days
  collection_days_list <- unique(curr_data$study_day_collected)
  
  print(patients_in_study[i]) # print patient number
  
  # print list of datatime 
  print(genebank_data_df_summarised %>% 
          filter(patient == patients_in_study[i]) %>% 
          summarise(datatime_list = paste(datatime, collapse = ",")))
  
  # print list of datanum
  print(genebank_data_df_summarised %>% 
    filter(patient == patients_in_study[i]) %>% 
    summarise(datanum_list = paste(datanum, collapse = ",")))
  
  for (j in 1:length(collection_days_list)){ # loop through each of the collection days groups to print the data appropriately for the init file
    print(curr_data_collection_ids <- curr_data %>% 
      filter(study_day_collected == collection_days_list[j]) %>% 
        summarise(id = paste(dQuote(id), collapse = ",")))
    
    print(paste0("{", curr_data_collection_ids, "}")) # paste list of accession ids for that datatime
  }
}

# this code is to compile a table of the patient's notes from genebank
# two sets of patients - some included in study and others not, but still have fasta files compiled
# this calls one sample from each patient because comments are the same for all samples of the same patient 
patients_ID_in_study <- c(
  "AY001661.1", #patient 100
  "AY000682.1", #patient 106
  "AY001619.1", #patient 126
  "AY001563.1", #patient 132
  "AY000421.1", #patient 1
  "AY000828.1", #patient 21
  "AY000587.1", #patient 26
  "AY000540.1", #patient 28
  "AY000759.1", #patient 55
  "AY001157.1", #patient 60
  "AY000913.1" #patient 69
)

patients_ID_not_in_study <- c(
  "AY000221.1", #patient 11
  "AY001527.1", #patient 141
  "AY000157.1", #patient 14
  "AY001961.1", #patient 155
  "AY001381.1", #patient 156
  "AY000282.1", #patient 15
  "AY000363.1", #patient 20
  "AY000252.1", #patient 3 
  "AY001081.1", #patient 85
  "AY002087.1", #patient 86
  "AY001305.1", #patient 87
  "AY001497.1" #patient 93
)

patients_all <- c(patients_ID_in_study, patients_ID_not_in_study)
therapy_hx_df <- data.frame("ID" = as.character(), "therapy_hx" = as.character())

for (i in 1:length(patients_all)){
  # pull metadata from genebank using API 
  text <- entrez_fetch(db="nucleotide", id=patients_all[i], rettype="gb", retmode="text")
  
  # Extract the match
  extracted_text <- str_extract(text, "/note=\"([^\"]+)\"")
  
  # Remove the `/note="` part to keep only the relevant information
  extracted_text <- gsub('/note="', '', extracted_text)
  
  # put text of note into growing dataframe
  therapy_hx_df_tmp <- data.frame("ID" = patients_all[i], "therapy_hx" = extracted_text)
  therapy_hx_df <- rbind(therapy_hx_df, therapy_hx_df_tmp)
  
  # remove tmp df 
  rm(therapy_hx_df_tmp)
}

# add more metadata to df
therapy_hx_df <- therapy_hx_df %>% 
  # clean up double spaces in text
  mutate(therapy_hx = gsub("\\s+", " ", therapy_hx)) %>% 
  # add patient number
  mutate(patient = str_extract(therapy_hx, "P\\d+")) %>% 
  # add column for whether they were included in the OG analysis
  mutate(in_OG_analyis = case_when(ID %in% patients_ID_in_study ~ "y",
                                   ID %in% patients_ID_not_in_study ~ "n", 
                                   TRUE ~ "ruh-roh")
  ) %>% 
# therapy change day present? 
  mutate(therapy_change_date_present_yn = ifelse(
    str_detect(therapy_hx, "then later switched to"),  # Check if phrase exists
    "n",  # If phrase exists
    "y"   # If phrase does not exist
  ))
  
# save therapy_hx_df table
write.csv(x = therapy_hx_df, file = "data/20250328_genebank_note_df.csv", row.names=FALSE)
