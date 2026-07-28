rm(list = ls())


library(readxl)
library(dplyr)
library(purrr)
library(progress)
library(tibble)
library(openxlsx2)
###############################################
# 1. READ DATA
###############################################

df <- read_excel("dados.xlsx")

###############################################
# 2. DATA CLEANING
###############################################

data_clean <- df %>%
  mutate(
    Site = trimws(as.character(Site)),
    Day = Day,
    Species = trimws(as.character(Species)),
    Color = trimws(toupper(as.character(Color)))
  )


###############################################
# 3. RECALCULATE Species_composition
###############################################

data_obs <- data_clean %>%
  group_by(Site, Day) %>%
  mutate(
    Species_composition = case_when(
      n() == 1 ~ 0,
      n_distinct(Species) == 1 ~ 1,
      n_distinct(Species) > 1 ~ 2
    )
  ) %>%
  ungroup()

###############################################
# 4. METRIC CALCULATION FUNCTION
# ALL METRICS ARE CALCULATED AT THE AGGREGATION LEVEL
###############################################

calc_metrics <- function(data) {
  
  ###############################################
  # AGGREGATION TABLE
  # One row represents one aggregation (Site + Day)
  ###############################################
  
  agg_data <- data %>%
    group_by(Site, Day) %>%
    summarise(
      Species_composition = first(Species_composition),
      has_APO = any(Color == "APO"),
      has_CRY = any(Color == "CRY"),
      only_APO = all(Color == "APO"),
      only_CRY = all(Color == "CRY"),
      .groups = "drop"
    )
  
  ###############################################
  # SINGLE-COLOR AGGREGATIONS
  # (MSA or SSA composed exclusively of APO or CRY individuals)
  ###############################################
  
  counts_by_color <- bind_rows(
    agg_data %>%
      filter(only_APO, Species_composition != 0) %>%
      count(Species_composition, name = "n") %>%
      mutate(Color = "APO"),
    
    agg_data %>%
      filter(only_CRY, Species_composition != 0) %>%
      count(Species_composition, name = "n") %>%
      mutate(Color = "CRY")
  ) %>%
    mutate(metric_type = "count_by_color") %>%
    select(metric_type, Color, Species_composition, n)
  
  ###############################################
  # TOTAL NUMBER OF AGGREGATIONS BY Species_composition
  ###############################################
  
  counts_all_colors <- agg_data %>%
    count(Species_composition, name = "n") %>%
    mutate(
      metric_type = "count_all_colors",
      Color = NA_character_
    ) %>%
    select(metric_type, Color, Species_composition, n)
  
  ###############################################
  # TOTAL NUMBER OF AGGREGATIONS CONTAINING EACH COLOR
  ###############################################
  
  counts_color_total <- bind_rows(
    agg_data %>%
      filter(has_APO) %>%
      summarise(n = n()) %>%
      mutate(Color = "APO"),
    
    agg_data %>%
      filter(has_CRY) %>%
      summarise(n = n()) %>%
      mutate(Color = "CRY")
  ) %>%
    mutate(
      metric_type = "count_color_total",
      Species_composition = NA_real_
    ) %>%
    select(metric_type, Color, Species_composition, n)
  
  ###############################################
  # SINGLE-COLOR AGGREGATIONS
  # Exclude solitary individuals (Species_composition != 0)
  ###############################################
  
  mono_metrics <- tibble(
    metric_type = "aggregations_single_color",
    Color = c("APO", "CRY"),
    Species_composition = NA_real_,
    n = c(
      sum(agg_data$Species_composition != 0 & agg_data$only_APO),
      sum(agg_data$Species_composition != 0 & agg_data$only_CRY)
    )
  )
  
  ###############################################
  # MSA AGGREGATIONS CONTAINING BOTH APO AND CRY
  ###############################################
  
  msa_both_colors <- agg_data %>%
    filter(
      Species_composition == 2,
      has_APO,
      has_CRY
    ) %>%
    summarise(n = n()) %>%
    mutate(
      metric_type = "aggregations_MSA_with_both_APO_and_CRY",
      Color = NA_character_,
      Species_composition = 2
    ) %>%
    select(metric_type, Color, Species_composition, n)
  
  bind_rows(
    counts_by_color,
    counts_all_colors,
    counts_color_total,
    mono_metrics,
    msa_both_colors
  )
}

###############################################
# 5. TEMPORAL NULL MODEL
###############################################

generate_null_model_temporal <- function(data) {
  
  data %>%
    group_by(Day) %>%
    mutate(
      shuffle = sample.int(dplyr::n()),
      Species = Species[shuffle],
      Color = Color[shuffle]
    ) %>%
    select(-shuffle) %>%
    ungroup() %>%
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

###############################################
# 6. OBSERVED METRICS
###############################################

obs_metrics <- calc_metrics(data_obs)

###############################################
# 7. SIMULATION PARAMETERS
###############################################

seed_model <- 123
n_sims <- 1000

set.seed(seed_model)

###############################################
# 8. TEMPORAL SIMULATIONS
###############################################

pb_temporal <- progress_bar$new(
  format = "Temporal [:bar] :percent eta: :eta | :current/:total",
  total = n_sims
)

null_temporal <- map_dfr(1:n_sims, function(i) {
  
  pb_temporal$tick()
  
  null_data <- generate_null_model_temporal(data_clean)
  
  calc_metrics(null_data) %>%
    mutate(run = i)
})

saveRDS(null_temporal, "long_simulations_AGGREGATION_TYPES.rds")
###############################################
# 9. RESULT SUMMARY FUNCTION
###############################################

calc_results <- function(null_data, obs_data, model_name) {
  
  summary_null <- null_data %>%
    group_by(metric_type, Color, Species_composition) %>%
    summarise(
      mean_null = mean(n),
      sd_null = sd(n),
      .groups = "drop"
    )
  
  pvals <- null_data %>%
    group_by(metric_type, Color, Species_composition) %>%
    group_modify(~{
      
      null_vec <- .x$n
      
      obs_val <- obs_data %>%
        filter(
          metric_type == .y$metric_type,
          ((is.na(Color) & is.na(.y$Color)) | Color == .y$Color),
          ((is.na(Species_composition) & is.na(.y$Species_composition)) |
             Species_composition == .y$Species_composition)
        ) %>%
        pull(n)
      
      n_null <- sum(!is.na(null_vec))
      
      p_ge_obs <- (sum(null_vec >= obs_val) + 1) / (n_null + 1)
      p_le_obs <- (sum(null_vec <= obs_val) + 1) / (n_null + 1)
      p_two <- min(2 * min(p_ge_obs, p_le_obs), 1)
      
      tibble(
        obs = obs_val,
        p_ge_obs = p_ge_obs,
        p_le_obs = p_le_obs,
        p_two_sided = p_two
      )
    }) %>%
    ungroup()
  
  summary_null %>%
    left_join(
      pvals,
      by = c("metric_type", "Color", "Species_composition")
    ) %>%
    mutate(
      SES = (obs - mean_null) / sd_null,
      model = model_name
    ) %>%
    select(
      model,
      metric_type,
      Color,
      Species_composition,
      obs,
      mean_null,
      sd_null,
      SES,
      p_ge_obs,
      p_le_obs,
      p_two_sided
    )
}

###############################################
# 10. CALCULATE FINAL RESULTS
###############################################

results_final <- calc_results(null_temporal, obs_metrics, "temporal")

###############################################
# 11. EXPORT RESULTS
###############################################

output_file <- "RESULTADOS_AGGREGATION_TYPES_TEMPORAL_103.csv"

write.csv(results_final, output_file, row.names = FALSE)


output_file <- "RESULTADOS_AGGREGATION_TYPES_TEMPORAL_103.xlsx"
write_xlsx(results_final, output_file )

print(results_final)
