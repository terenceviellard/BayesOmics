
# Load useful packages
library(tidyverse)

# Load data
qData <- readxl::read_xlsx(path = "data/pr5c01206_si_002.xlsx", sheet = "Table S1")

# Data preprocessing
data_clean <- qData |> 
  # Long format
  pivot_longer(cols = starts_with(c("Relapse", "Cured")),
    names_to = "Sample_Name", values_to = "Output") |> 
  # Clean the sample names
  separate(col = Sample_Name, into = c("Outcome", "Patient_ID", "Timestamp"),sep = "\\.") |> 
  # Define Input, Sample and Group
  mutate(Input = as.numeric(Timestamp),
    Sample = paste(Outcome, Patient_ID, sep = "_"),
    Group = if_else(str_detect(Outcome, pattern = "Cured"),
      true = "Cured", false = "Relapsed"), .keep = "unused") |> 
  # Rename Feature identifier
  rename(ID = UniProt_Accession) |> 
  # Remove NA values
  filter(!is.na(Output))


pak::pak("terenceviellard/BayesOmics")
library(BayesOmics)
library(keRnel)
set.seed(17)

# Fit the kernel
kern <- variance_kernel(variance = 1) * se_kernel(length_scale = 50) +
  white_noise_kernel(noise = 1)

group1 <- data_clean[data_clean$Group == unique(data_clean$Group)[1],]
prior_mean <- mean(group1$Output)
hp_opt <- fit_kernel(kern, group1, prior_mean = prior_mean, prior_cov = 1)
kern_opt   <- kupdate(kern,
                       variance     = hp_opt[["variance"]],
                       length_scale = hp_opt[["length_scale"]],
                       noise        = hp_opt[["noise"]])