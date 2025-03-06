#load libraries
library(tidyverse)

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
)


# plot the entire chain trace plot -- i.e. make one figure with 4 subplots; 
  # one subplot for each of the 4 parameters (omega1, omega2, lambda, Nbegin)


# plot the ACFs (also 4 subplots as above) 