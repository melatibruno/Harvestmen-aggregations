################################################################################
# COLOR ASSORTATIVITY NULL MODEL: PREDICTIONS P1-P6
################################################################################


# ==============================================================================
# 02_color_assortativity_P1_P6.R
#
# Color assortativity conditional on MSA formation
#
# Purpose
# -------
# This script tests color assortativity within multispecies aggregations (MSA)
# using six conditional predictions at two hierarchical levels:
#
#   1. Aggregation level:
#      - P1: proportion of all MSA that are mixed-color
#      - P2: among MSA containing cryptic species, proportion also containing
#            conspicuous species
#      - P3: among MSA containing conspicuous species, proportion also containing
#            cryptic species
#
#   2. Individual-record level:
#      - P4: proportion of all records within MSA occurring in mixed-color MSA
#      - P5: among cryptic records within MSA, proportion occurring in
#            mixed-color MSA
#      - P6: among conspicuous records within MSA, proportion occurring in
#            mixed-color MSA
#
# Null model
# ----------
# Within each sampling day, Species and Color are shuffled jointly among the
# observed records while Site remains fixed. This preserves daily species
# abundances, the number of aggregation sites, and aggregation sizes, while
# breaking the association between species identity and aggregation membership.
#
# Species and Color are randomized jointly because each species has a fixed
# color phenotype.
#
# IMPORTANT
# ---------
# This script contains only the conditional color-assortativity analyses (P1-P6).
# Overall species-assortment analyses are implemented separately in
# 01_species_assortment_null_model.R.
#
# Figure construction is intentionally kept in a separate figure script.
# ==============================================================================

rm(list = ls())

library(readxl)
library(dplyr)
library(purrr)
library(progress)
library(tibble)
library(tidyr)
library(openxlsx2)


# ==============================================================================
# 1. IMPORT AND CLEAN DATA
# ==============================================================================

df <- read_excel("dados.xlsx")

data_clean <- df %>%
  mutate(
    Site = trimws(as.character(Site)),
    Species = trimws(as.character(Species)),
    Color = trimws(toupper(as.character(Color)))
  )

unexpected_colors <- setdiff(unique(na.omit(data_clean$Color)), c("CRY", "APO"))
if (length(unexpected_colors) > 0) {
  stop(paste("Unexpected Color values:", paste(unexpected_colors, collapse = ", ")))
}

# ==============================================================================
# 2. CLASSIFY SPECIES COMPOSITION
# ==============================================================================

classify_species_composition <- function(data) {
  data %>%
    group_by(Site, Day) %>%
    mutate(
      Species_composition = case_when(
        n() == 1 ~ 0,
        n_distinct(Species) == 1 ~ 1,
        n_distinct(Species) > 1 ~ 2
      )
    ) %>%
    ungroup()
}

data_obs <- classify_species_composition(data_clean)

safe_prop <- function(num, den) {
  if (den > 0) num / den else NA_real_
}

# ==============================================================================
# 3. SIX COLOR-ASSORTATIVITY PREDICTIONS
#
# Aggregation level
# P1 = mixed-color MSA / all MSA
# P2 = mixed-color MSA / MSA containing CRY
# P3 = mixed-color MSA / MSA containing APO
#
# Individual-record level
# P4 = records in mixed-color MSA / all records in MSA
# P5 = CRY records in mixed-color MSA / all CRY records in MSA
# P6 = APO records in mixed-color MSA / all APO records in MSA
# ==============================================================================

calc_color_assortativity <- function(data) {

  agg_data <- data %>%
    group_by(Site, Day) %>%
    summarise(
      Species_composition = first(Species_composition),
      has_CRY = any(Color == "CRY"),
      has_APO = any(Color == "APO"),
      n_records = n(),
      n_CRY = sum(Color == "CRY"),
      n_APO = sum(Color == "APO"),
      .groups = "drop"
    ) %>%
    mutate(
      is_MSA = Species_composition == 2,
      is_mixed_color = has_CRY & has_APO
    )

  msa <- agg_data %>% filter(is_MSA)

  # P1
  p1_num <- sum(msa$is_mixed_color)
  p1_den <- nrow(msa)

  # P2
  msa_CRY <- msa %>% filter(has_CRY)
  p2_num <- sum(msa_CRY$has_APO)
  p2_den <- nrow(msa_CRY)

  # P3
  msa_APO <- msa %>% filter(has_APO)
  p3_num <- sum(msa_APO$has_CRY)
  p3_den <- nrow(msa_APO)

  # P4
  p4_num <- sum(msa$n_records[msa$is_mixed_color])
  p4_den <- sum(msa$n_records)

  # P5
  p5_num <- sum(msa$n_CRY[msa$is_mixed_color])
  p5_den <- sum(msa$n_CRY)

  # P6
  p6_num <- sum(msa$n_APO[msa$is_mixed_color])
  p6_den <- sum(msa$n_APO)

  tibble(
    prediction = paste0("P", 1:6),
    scale = c(rep("Aggregation level", 3), rep("Individual-record level", 3)),
    metric = c(
      "Mixed-color MSA among all MSA",
      "Mixed-color MSA among MSA containing CRY",
      "Mixed-color MSA among MSA containing APO",
      "Records in mixed-color MSA among all MSA records",
      "CRY records in mixed-color MSA among CRY MSA records",
      "APO records in mixed-color MSA among APO MSA records"
    ),
    numerator = c(p1_num, p2_num, p3_num, p4_num, p5_num, p6_num),
    denominator = c(p1_den, p2_den, p3_den, p4_den, p5_den, p6_den),
    value = c(
      safe_prop(p1_num, p1_den),
      safe_prop(p2_num, p2_den),
      safe_prop(p3_num, p3_den),
      safe_prop(p4_num, p4_den),
      safe_prop(p5_num, p5_den),
      safe_prop(p6_num, p6_den)
    )
  )
}

# ==============================================================================
# 4. OBSERVED VALUES
# ==============================================================================

observed_predictions <- calc_color_assortativity(data_obs)
print(observed_predictions)

# ==============================================================================
# 5. TEMPORAL NULL MODEL
#
# Species identity and color are shuffled jointly within Day. Site x Day
# positions and aggregation sizes remain fixed.
# ==============================================================================

generate_null_model_temporal <- function(data) {
  data %>%
    group_by(Day) %>%
    mutate(
      shuffle = sample.int(n()),
      Species = Species[shuffle],
      Color = Color[shuffle]
    ) %>%
    select(-shuffle) %>%
    ungroup() %>%
    classify_species_composition()
}

# ==============================================================================
# 6. RUN 1,000 RANDOMIZATIONS
# ==============================================================================

seed_model <- 123
n_sims <- 1000
set.seed(seed_model)

pb_temporal <- progress_bar$new(
  format = "Temporal [:bar] :percent eta: :eta | :current/:total",
  total = n_sims
)

null_predictions <- map_dfr(seq_len(n_sims), function(i) {
  pb_temporal$tick()

  generate_null_model_temporal(data_clean) %>%
    calc_color_assortativity() %>%
    mutate(run = i)
})

saveRDS(null_predictions, "long_simulations_COLOR_ASSORTATIVITY_P1_P6.rds")

# ==============================================================================
# 7. NULL DISTRIBUTIONS, SES, AND TWO-TAILED MONTE CARLO P-VALUES
# ==============================================================================

calc_results <- function(null_data, observed_data) {

  null_summary <- null_data %>%
    group_by(prediction, scale, metric) %>%
    summarise(
      mean_null = mean(value, na.rm = TRUE),
      sd_null = sd(value, na.rm = TRUE),
      n_valid_simulations = sum(!is.na(value)),
      .groups = "drop"
    )

  p_values <- null_data %>%
    group_by(prediction, scale, metric) %>%
    group_modify(~ {
      null_vec <- .x$value
      null_vec <- null_vec[!is.na(null_vec)]

      obs_val <- observed_data %>%
        filter(prediction == .y$prediction) %>%
        pull(value)

      n_null <- length(null_vec)

      p_ge_obs <- (sum(null_vec >= obs_val) + 1) / (n_null + 1)
      p_le_obs <- (sum(null_vec <= obs_val) + 1) / (n_null + 1)

      tibble(
        observed = obs_val,
        p_ge_obs = p_ge_obs,
        p_le_obs = p_le_obs,
        p_two_sided = min(2 * min(p_ge_obs, p_le_obs), 1)
      )
    }) %>%
    ungroup()

  null_summary %>%
    left_join(p_values, by = c("prediction", "scale", "metric")) %>%
    mutate(SES = (observed - mean_null) / sd_null) %>%
    select(
      prediction, scale, metric,
      observed, mean_null, sd_null, SES,
      p_ge_obs, p_le_obs, p_two_sided,
      n_valid_simulations
    )
}

results_predictions <- calc_results(null_predictions, observed_predictions)
print(results_predictions)

# ==============================================================================
# 8. EXPORT RESULTS
# ==============================================================================

write.csv(
  results_predictions,
  "RESULTS_COLOR_ASSORTATIVITY_P1_P6.csv",
  row.names = FALSE
)

write_xlsx(
  results_predictions,
  "RESULTS_COLOR_ASSORTATIVITY_P1_P6.xlsx"
)

write.csv(
  observed_predictions,
  "OBSERVED_COLOR_ASSORTATIVITY_P1_P6.csv",
  row.names = FALSE
)



sessionInfo()
