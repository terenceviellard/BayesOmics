
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
  # Define ID, Input, Sample and Group
  mutate(ID = paste0("Week_", Timestamp),
    Input = as.numeric(Timestamp),
    Sample = paste(Outcome, Patient_ID, sep = "_"),
    Group = if_else(str_detect(Outcome, pattern = "Cured"),
      true = "Cured", false = "Relapsed"), .keep = "unused") |> 
  # Rename Feature identifier
  rename(Protein = UniProt_Accession) |> 
  # Remove NA values
  filter(!is.na(Output))

data_clean |> 
  count(Protein, Sample, ID) |> 
  filter(n > 1)

pak::pak("terenceviellard/BayesOmics")
library(BayesOmics)
library(keRnel)
set.seed(17)
kern <- variance_kernel(variance = 1) * se_kernel(length_scale = 12) +
  white_noise_kernel(noise = 1)

# Only one protein
weeks <- c(0, 2, 4, 8, 17, 26, 52)

data_dpp4 <- data_clean |>
  filter(Protein == "P27487") |>
  group_by(Group, Sample) |>
  filter(
    n() == 7,
    n_distinct(Input) == 7,
    setequal(Input, weeks),
    all(is.finite(Output))
  ) |>
  ungroup()

data_dpp4 |>
  count(Group, ID) |>
  pivot_wider(
    names_from = Group,
    values_from = n
  )

group1 <- data_dpp4 |>
  filter(Group == "Cured")
prior_mean <- mean(group1$Output)

hp_opt <- fit_kernel(
  kern,
  group1,
  prior_mean = mean(group1$Output),
  prior_cov = 1
)

kern_opt   <- kupdate(kern,
                       variance     = hp_opt[["variance"]],
                       length_scale = hp_opt[["length_scale"]],
                       noise        = hp_opt[["noise"]])

posterior <- posterior_mean(data_dpp4, kern_opt, mu_0 = prior_mean, lambda_0 = 1)
group_diff(posterior)

# Tirages des distributions postérieures
set.seed(17)
samples <- sample_posterior(posterior, n = 2000)

# Comparaison à un temps particulier
plot_posterior_overlap(
  samples,
  group1 = "Cured",
  group2 = "Relapsed",
  id = "Week_0"  # Adapter à la casse de tes ID
)

# Différence conjointe sur les 7 temps
multi_diff <- compute_multi_diff(
  samples,
  results = posterior
)

plot_multi_diff(multi_diff)


# Fit the kernel on all data
group1 <- data_clean[data_clean$Group == unique(data_clean$Group)[1],]
prior_mean <- mean(group1$Output)
hp_opt <- fit_kernel(kern, group1, prior_mean = prior_mean, prior_cov = 1)
kern_opt   <- kupdate(kern,
                       variance     = hp_opt[["variance"]],
                       length_scale = hp_opt[["length_scale"]],
                       noise        = hp_opt[["noise"]])