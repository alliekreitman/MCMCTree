library(tidyverse)
library(Biostrings)
library(readxl)

### NEXT UP ###
# add list.files for all fasta
# then run loop fxn
# then summarise by patient for detected resistance mutations
# and join therapy they received

fasta_combotherapy_files <- list.files("data/Li_DnDs", pattern="comotherapy.fasta", full.names=T)
fasta_monotherapy_files <- list.files("data/Li_DnDs", pattern="monotherapy.fasta", full.names=T)
fasta_files <- c(fasta_combotherapy_files, fasta_monotherapy_files)
# table therapy per patient
therapy_per_patient <- read.csv("data/Data_MCTree12/Therapies_by_patients.csv") %>% 
  filter(Explicit_switch_in_genebank == "y") %>% 
  select(-Explicit_switch_in_genebank)
  # mutate(all_meds = paste(Mono, Combo, sep = " + ")) %>%
  # separate_rows(all_meds, sep = " \\+ ") %>%
  # mutate(medication = str_trim(all_meds)) %>%
  # distinct(Patient, medication)

therapy_per_patient_simplified <- therapy_per_patient %>%
  pivot_longer(!Patient, names_to = "combo_mono", values_to = "therapy") %>%
  mutate(combo_mono = ifelse(combo_mono == "Combo", "combotherapy", "monotherapy"))


patient_meds <- therapy_per_patient %>%
  mutate(all_meds = paste(Mono, Combo, sep = " + ")) %>%
  separate_rows(all_meds, sep = " \\+ ") %>%
  mutate(medication = str_trim(all_meds)) %>%
  distinct(Patient, medication) %>% 
  mutate(medication = case_when(medication == "indinavir" ~ "Indinavir",
                                medication == "efavirenz" ~ "Efavirenz",
                                .default=as.character(medication)))
  

# table mutations per therapy

resistance_mutations_table_raw <- read_excel("data/known_resistance_mutations_per_therapy.xlsx")
# now separate reistance mutation into its two AA and pos components 
resistance_mutations_table <- resistance_mutations_table_raw %>% 
  mutate(
    reference_aa = str_extract(mutation, "^[A-Z]"),
    position = str_extract(mutation, "[0-9]+") %>% as.integer(),
    resistance_aa = str_extract(mutation, "[A-Z]$")
  )
  

# for now, just use 1 fasta file (for debugging)
# LATER, use. list.files for all fasta files
# fasta_file <- "data/Data_MCTree12/p1/fasta/p60sequence.fasta"

results <- data.frame()

for (fasta_file in fasta_files){
  print(fasta_file)
  
  sequences <- readDNAStringSet(fasta_file)
  aa_sequences <- translate(sequences)
  
  
  for (i in seq_along(aa_sequences)) {
    seq_name <- names(aa_sequences)[i]
    seq <- as.character(aa_sequences[[i]])
    
    # Check each mutation
    for (j in 1:nrow(resistance_mutations_table)) {
      pos <- resistance_mutations_table$position[j]
      ref_aa <- resistance_mutations_table$reference_aa[j]
      res_aa <- resistance_mutations_table$resistance_aa[j]
      
      # Extract amino acid at position
      if (pos <= nchar(seq)) {
        aa_at_pos <- substr(seq, pos, pos)
        
          results <- rbind(results, data.frame(
            sequence_id = seq_name,
            fasta_file = basename(fasta_file),
            combo_mono = case_when(fasta_file %in% fasta_combotherapy_files ~ "combotherapy",
                                   fasta_file %in% fasta_monotherapy_files ~ "monotherapy",
                                   .default = fasta_file),
            therapy = resistance_mutations_table$therapy[j],
            mutation = resistance_mutations_table$mutation[j],
            position = pos,
            reference_aa = ref_aa,
            observed_aa = aa_at_pos,
            # Log if mutation is present
            resistance_detected = ifelse(aa_at_pos == res_aa, TRUE, FALSE),
            stringsAsFactors = FALSE
          ))
      }
    }
  }
}


results_edited <- results %>% 
  mutate(therapy_mut_associated_with = ifelse(therapy == "AZT (Zidovudine)" | therapy == "3TC (Lamivudine)", "ZDV/3TC", therapy)) %>%
  select(-therapy) %>% 
  group_by(fasta_file, combo_mono, therapy_mut_associated_with, mutation, position, reference_aa, observed_aa) %>% 
  summarise(mutation_proportion = mean(resistance_detected, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(mutation_proportion > 0) %>% 
  rename("combo_mono" = "isolate_time") %>% 
  mutate(Patient = as.numeric(str_extract(fasta_file, "(?<=p)\\d+"))) %>% 
  # keep only the resistance mutations associated with therapy they received
  inner_join(patient_meds, by = c("Patient", "therapy_mut_associated_with" = "medication")) %>% 
  select(-mutation_proportion, -position, -observed_aa, -reference_aa, -fasta_file) %>% 
  left_join(therapy_per_patient, by = "Patient", relationship = "many-to-many") %>% 
  rename("Mono" = "monotherapy",
         "Combo" = "combotherapy") %>% 
  select(Patient, monotherapy, combotherapy, mutation, therapy_mut_associated_with, isolate_time)


write.csv(results_edited, "out/resistance_mutations_identified_in_isolates.csv", row.names=F)

