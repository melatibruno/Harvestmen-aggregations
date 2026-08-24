# ==============================================================================
# 01_species_assortment_null_model.R
#
# Species assortment in Prionostemma aggregations
#
# Tests overall assortment by species identity at two levels:
#   1. Aggregation level: number of MSA and SSA.
#   2. Individual-record level: cryptic and conspicuous records in MSA and SSA.
#
# The null model reallocates the daily pool of individual records among the
# observed aggregation positions. Species and Color are shuffled jointly because
# each species has a fixed color phenotype. This preserves daily abundances,
# number of sites, and aggregation sizes while breaking the association between
# species identity and aggregation membership.
#
# This script contains ANALYSES ONLY. Figure code is kept in a separate script.
# Conditional color-assortativity analyses P1-P6 are also kept separately.
# ==============================================================================

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

# ==============================================================================
# 2. USER SETTINGS
# ==============================================================================

data_file <- "dados.xlsx"
n_sims <- 1000
random_seed <- 123

# ==============================================================================
# 3. IMPORT AND CLEAN DATA
# ==============================================================================

data_clean <- read_excel(data_file) %>%
  transmute(
    Site = trimws(as.character(Site)),
    Day = Day,
    Species = trimws(as.character(Species)),
    Color = trimws(toupper(as.character(Color)))
  )

stopifnot(
  all(!is.na(data_clean$Site)),
  all(!is.na(data_clean$Day)),
  all(!is.na(data_clean$Species)),
  all(!is.na(data_clean$Color))
)

# ==============================================================================
# 4. CLASSIFY AGGREGATION TYPE
#
# Species_composition:
#   0 = solitary record
#   1 = single-species aggregation (SSA)
#   2 = multispecies aggregation (MSA)
# ==============================================================================

classify_species_composition <- function(data) {
  data %>%
    group_by(Site, Day) %>%
    mutate(
      Species_composition = case_when(
        n() == 1 ~ 0L,
        n_distinct(Species) == 1 ~ 1L,
        n_distinct(Species) > 1 ~ 2L
      )
    ) %>%
    ungroup()
}

data_obs <- classify_species_composition(data_clean)

# ==============================================================================
# 5. CALCULATE SPECIES-ASSORTMENT STATISTICS
#
# Aggregation level:
#   n_MSA, n_SSA
#
# Individual-record level:
#   records_CRY_MSA, records_CRY_SSA
#   records_APO_MSA, records_APO_SSA
#
# CRY = cryptic; APO = conspicuous (dataset codes retained).
# ==============================================================================

calc_species_assortment <- function(data) {

  aggregation_data <- data %>%
    group_by(Site, Day) %>%
    summarise(
      Species_composition = first(Species_composition),
      .groups = "drop"
    )

  n_MSA <- sum(aggregation_data$Species_composition == 2, na.rm = TRUE)
  n_SSA <- sum(aggregation_data$Species_composition == 1, na.rm = TRUE)

  records_CRY_MSA <- sum(
    data$Color == "CRY" & data$Species_composition == 2, na.rm = TRUE
  )
  records_CRY_SSA <- sum(
    data$Color == "CRY" & data$Species_composition == 1, na.rm = TRUE
  )
  records_APO_MSA <- sum(
    data$Color == "APO" & data$Species_composition == 2, na.rm = TRUE
  )
  records_APO_SSA <- sum(
    data$Color == "APO" & data$Species_composition == 1, na.rm = TRUE
  )

  tibble(
    metric = c(
      "n_MSA", "n_SSA",
      "records_CRY_MSA", "records_CRY_SSA",
      "records_APO_MSA", "records_APO_SSA"
    ),
    scale = c(
      rep("Aggregation level", 2),
      rep("Individual-record level", 4)
    ),
    description = c(
      "Number of multispecies aggregations (MSA)",
      "Number of single-species aggregations (SSA)",
      "Cryptic records occurring in MSA",
      "Cryptic records occurring in SSA",
      "Conspicuous records occurring in MSA",
      "Conspicuous records occurring in SSA"
    ),
    value = c(
      n_MSA, n_SSA,
      records_CRY_MSA, records_CRY_SSA,
      records_APO_MSA, records_APO_SSA
    )
  )
}

observed_species_assortment <- calc_species_assortment(data_obs)

# ==============================================================================
# 6. TEMPORAL NULL MODEL
#
# Within each Day, Species and Color are shuffled jointly among records.
# Site remains fixed, preserving the observed number and size of aggregations.
# Species_composition is recalculated after each randomization.
# ==============================================================================

randomize_once <- function(data) {
  data %>%
    group_by(Day) %>%
    group_modify(~ {
      shuffled_rows <- sample.int(nrow(.x))
      .x %>%
        mutate(
          Species = Species[shuffled_rows],
          Color = Color[shuffled_rows]
        )
    }) %>%
    ungroup() %>%
    classify_species_composition()
}

# ==============================================================================
# 7. RUN NULL-MODEL SIMULATIONS
# ==============================================================================

set.seed(random_seed)

null_species_assortment <- map_dfr(
  seq_len(n_sims),
  function(i) {
    randomize_once(data_clean) %>%
      calc_species_assortment() %>%
      mutate(run = i, .before = 1)
  }
)

# ==============================================================================
# 8. CALCULATE SES AND TWO-TAILED MONTE CARLO P-VALUES
#
# SES = (observed - mean(null)) / SD(null)
#
# p_upper = (sum(null >= observed) + 1) / (N + 1)
# p_lower = (sum(null <= observed) + 1) / (N + 1)
# p_two_sided = min(2 * min(p_upper, p_lower), 1)
# ==============================================================================

calc_null_results <- function(null_data, observed_data) {

  metric_info <- observed_data %>%
    select(metric, scale, description)

  null_data %>%
    group_by(metric) %>%
    summarise(
      mean_null = mean(value, na.rm = TRUE),
      sd_null = sd(value, na.rm = TRUE),
      n_valid_simulations = sum(!is.na(value)),
      null_values = list(value[!is.na(value)]),
      .groups = "drop"
    ) %>%
    left_join(
      observed_data %>% select(metric, observed = value),
      by = "metric"
    ) %>%
    mutate(
      SES = (observed - mean_null) / sd_null,
      p_ge_obs = map2_dbl(
        null_values, observed,
        ~ (sum(.x >= .y) + 1) / (length(.x) + 1)
      ),
      p_le_obs = map2_dbl(
        null_values, observed,
        ~ (sum(.x <= .y) + 1) / (length(.x) + 1)
      ),
      p_two_sided = pmin(2 * pmin(p_ge_obs, p_le_obs), 1)
    ) %>%
    select(
      metric, observed, mean_null, sd_null, SES,
      p_ge_obs, p_le_obs, p_two_sided, n_valid_simulations
    ) %>%
    left_join(metric_info, by = "metric") %>%
    select(
      metric, scale, description,
      observed, mean_null, sd_null, SES,
      p_ge_obs, p_le_obs, p_two_sided, n_valid_simulations
    )
}

results_species_assortment <- calc_null_results(
  null_species_assortment,
  observed_species_assortment
)

# ==============================================================================
# 9. INSPECT RESULTS
# ==============================================================================

print(observed_species_assortment)
print(results_species_assortment)

# ==============================================================================
# 10. SAVE OUTPUTS
#
# The RDS contains the complete null distributions so figures can be regenerated
# without rerunning the 1,000 randomizations.
# ==============================================================================

saveRDS(
  null_species_assortment,
  "long_simulations_SPECIES_ASSORTMENT.rds"
)

write.csv(
  observed_species_assortment,
  "observed_SPECIES_ASSORTMENT.csv",
  row.names = FALSE
)

write.csv(
  results_species_assortment,
  "results_SPECIES_ASSORTMENT.csv",
  row.names = FALSE
)

# ==============================================================================
# END
# ==============================================================================
