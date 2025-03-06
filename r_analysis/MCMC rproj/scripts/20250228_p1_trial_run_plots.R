#load libraries
library(tidyverse)
library(cowplot)

# load data
# Read the file
omega_lines <- readLines("data/20250228_p1_test_run/p1_76807_411274_omega.txt")
mu_lines <- readLines("data/20250228_p1_test_run/p1_76807_411274_mu.txt") 

# Process data
omega_results <- tibble(
  omega1 = as.numeric(str_extract(omega_lines[seq(1, length(omega_lines), by = 3)], "^[0-9\\.]+")),
  omega2 = as.numeric(str_extract(omega_lines[seq(1, length(omega_lines), by = 3)], "(?<=\\s)[0-9\\.]+$")),
  Unique_Tree_count = as.numeric(gsub("Unique_Tree_count=", "", omega_lines[seq(2, length(omega_lines), by = 3)])),
  Iter_count = as.numeric(gsub("Iter_count=", "", omega_lines[seq(3, length(omega_lines), by = 3)]))
)

mu_results <- tibble(
  lambda = as.numeric(str_extract(mu_lines[seq(1, length(mu_lines), by = 3)], "^[0-9\\.]+")),
  Nbegin = as.numeric(str_extract(mu_lines[seq(1, length(mu_lines), by = 3)], "(?<=\\s)[0-9\\.]+$")),
  Unique_Tree_count = as.numeric(gsub("Unique_Tree_count=", "", mu_lines[seq(2, length(mu_lines), by = 3)])),
  Iter_count = as.numeric(gsub("Iter_count=", "", mu_lines[seq(3, length(mu_lines), by = 3)]))
) %>% 
  mutate(log_Nbegin = log(Nbegin))

###### traceplots
# plot the entire chain trace plot -- i.e. make one figure with 4 subplots; 
# one subplot for each of the 4 parameters (omega1, omega2, lambda, Nbegin)
omega1_traceplot <- ggplot(omega_results, aes(x = Iter_count, y = omega1))+
  theme_bw()+
  geom_point()

omega2_traceplot <- ggplot(omega_results, aes(x = Iter_count, y = omega2))+
  theme_bw()+
  geom_point()

lambda_traceplot <- ggplot(mu_results, aes(x = Iter_count, y = lambda))+
  theme_bw()+
  geom_point()

log_Nbegin_traceplot <- ggplot(mu_results, aes(x = Iter_count, y = log_Nbegin))+
  theme_bw()+
  geom_point()

traceplots <- cowplot::plot_grid(omega1_traceplot, omega2_traceplot, lambda_traceplot, log_Nbegin_traceplot)
traceplots

# plot the ACFs (also 4 subplots as above) 

# Compute ACF using dplyr
acf_data <- full_join(omega_results, mu_results, by = "Iter_count", suffix = c(".omega", ".mu")) %>%
  select(omega1, omega2, lambda, log_Nbegin) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
  group_by(parameter) %>%
  summarise(acf_values = list(acf(value, plot = FALSE)$acf),
            lags = list(acf(value, plot = FALSE)$lag),
            .groups = "drop") %>%
  unnest(c(acf_values, lags))

# Plot the ACFs
ggplot(acf_data, aes(x = lags, y = acf_values)) +
  theme_bw()+
  geom_line() +
  geom_point() +
  facet_wrap(~parameter, scales = "free_y") +
  labs(title = "Autocorrelation Function (ACF) for Parameters",
       x = "Lag",
       y = "Autocorrelation") 


# # combine dataframes 
# omega_mu_results <- rbind(
#   omega_results %>%
#     pivot_longer(cols = c(omega1, omega2), names_to = "variable", values_to = "value"), 
#   mu_results %>% 
#     pivot_longer(cols = c(lambda, Nbegin), names_to = "variable", values_to = "value")
# )





