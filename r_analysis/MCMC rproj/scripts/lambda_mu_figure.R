# load libraries
require(tidyverse)
require(cowplot)

# load data
converged_data <- read.csv("Data_gewek/converged_data_from_gewek.csv")

# convert data to long format
converged_data_long <- converged_data %>% 
  pivot_longer(c(mono_therapy_dnds, combo_therapy_dnds), names_to = "dnds_therapy", values_to = "dnds") %>% 
  filter(patient != ("p156"))

# facet by patient
# for monotherapy 
lambda_v_mu_by_patient_mono <- ggplot(data = converged_data_long %>% filter(dnds_therapy == "mono_therapy_dnds"), aes(x = mutation_rate, y = dnds))+
  theme_bw()+
  facet_wrap(.~patient, ncol = 1, scales = "free", strip.position="right")+
  geom_point(size = 0.5)+
  ylab("dnds - monotherapy")

# and combotherapy
lambda_v_mu_by_patient_combo <- ggplot(data = converged_data_long %>% filter(dnds_therapy == "combo_therapy_dnds"), aes(x = mutation_rate, y = dnds))+
  theme_bw()+
  facet_wrap(.~patient, ncol = 1, scales = "free", strip.position = "right")+
  geom_point(size = 0.5)+
  ylab("dnds - combotherapy")

# put mono and combo therapy side by side in single plot
lambda_v_mu_by_patient <- plot_grid(lambda_v_mu_by_patient_mono, lambda_v_mu_by_patient_combo, ncol = 2)


# combine patients
# monotherapy 
lambda_v_mu_combine_patients_mono <- ggplot(data = converged_data_long %>% filter(dnds_therapy == "mono_therapy_dnds"), aes(x = mutation_rate, y = dnds))+
  theme_bw()+
  facet_grid(.~factor(dnds_therapy, levels = c("mono_therapy_dnds", "combo_therapy_dnds")), scales = "free")+
  geom_point(size = 0.5)

# combo therapy 
lambda_v_mu_combine_patients_combo <- ggplot(data = converged_data_long %>% filter(dnds_therapy == "combo_therapy_dnds"), aes(x = mutation_rate, y = dnds))+
  theme_bw()+
  facet_grid(.~factor(dnds_therapy, levels = c("mono_therapy_dnds", "combo_therapy_dnds")), scales = "free")+
  geom_point(size = 0.5)

# combine into single graph
lambda_v_mu_combine_patients <- plot_grid(lambda_v_mu_combine_patients_mono, lambda_v_mu_combine_patients_combo, ncol = 2)

# place all four layouts into one figure
final_layout <- plot_grid(lambda_v_mu_by_patient, lambda_v_mu_combine_patients, ncol = 1, rel_heights = c(3,1), labels = c('A', 'B'))

ggsave("Figures/20240915_lambda_v_mu_v3.0.pdf", plot = final_layout, height = 11, width = 8.5, unit = "in")
ggsave("Figures/20240915_lambda_v_mu_combined_patients_v3.0.pdf", plot = lambda_v_mu_combine_patients, height = 11, width = 8.5, unit = "in")
ggsave("Figures/20240915_lambda_v_mu_by_patient_v3.0.pdf", plot = lambda_v_mu_by_patient, height = 11, width = 8.5, unit = "in")
