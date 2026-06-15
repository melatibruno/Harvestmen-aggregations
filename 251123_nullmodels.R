setwd("C:/Users/bruno/OneDrive/Documentos/trabalhos_extra/Opilhões")

###############################################
# DISPLAY FULL TIBBLES WITHOUT TRUNCATION
###############################################
options(tibble.print_max = Inf, tibble.width = Inf)

#############################################
# 0. LOAD PACKAGES
###############################################
library(readxl)
library(dplyr)
library(purrr)
library(tibble)

###############################################
# 1. LER PLANILHA
###############################################

df <- read_excel("dados.xlsx")

###############################################
# 2. DATA CLEANING
###############################################

data_clean <- df %>%
  mutate(
    Site    = trimws(Site),
    Day     = Day,                     # keeps the original format (number, date, etc.)
    Species = trimws(Species),
    Color   = trimws(toupper(Color)),
    Species_composition = trimws(toupper(Species_composition))
  )

###############################################
# 3. DEFINE SPECIES_COMPOSITION FOR OBSERVED DATA
#     (same logic used in the null models)
###############################################

data_obs <- data_clean %>%
  group_by(Site, Day) %>%
  mutate(
    Species_composition = case_when(
      n() == 1 ~ "SOLITARY",
      n_distinct(Species) == 1 ~ "SSA",
      n_distinct(Species) > 1  ~ "MSA"
    )
  ) %>%
  ungroup()

###############################################
# 4. TEMPORAL NULL MODEL
#    - fixed Day
#    - shuffles Species and Color within each Day
###############################################

generate_null_model_temporal <- function(data, run_id = NULL) {
  data %>%
    group_by(Day) %>% 
    mutate(
      shuffle_idx = sample(n(), replace = FALSE),
      Species     = Species[shuffle_idx],
      Color       = Color[shuffle_idx]
    ) %>%
    select(-shuffle_idx) %>% 
    ungroup() %>%
    
    group_by(Site, Day) %>%
    mutate(
      Species_composition = case_when(
        n() == 1 ~ "SOLITARY",
        n_distinct(Species) == 1 ~ "SSA",
        n_distinct(Species) > 1  ~ "MSA"
      )
    ) %>%
    ungroup() %>%
    mutate(run = run_id)
}


###############################################
# 5. AGGREGATION-LEVEL METRICS (Site + Day)
###############################################

calc_aggregation_metrics <- function(data) {
  
  aggs <- data %>%
    group_by(Site, Day) %>%
    summarise(
      Species_composition = first(Species_composition),
      only_CRY = all(Color == "CRY"),
      only_APO = all(Color == "APO"),
      n_colors = n_distinct(Color),
      .groups = "drop"
    )
  
  tibble(
    n_MSA          = sum(aggs$Species_composition == "MSA"),
    n_SSA          = sum(aggs$Species_composition == "SSA"),
    
    n_MSA_only_CRY = sum(aggs$Species_composition == "MSA" & aggs$only_CRY),
    n_MSA_only_APO = sum(aggs$Species_composition == "MSA" & aggs$only_APO),
    
    n_SSA_only_CRY = sum(aggs$Species_composition == "SSA" & aggs$only_CRY),
    n_SSA_only_APO = sum(aggs$Species_composition == "SSA" & aggs$only_APO),
    
    n_MSA_two_colors = sum(aggs$Species_composition == "MSA" & aggs$n_colors >= 2),
    
    # agregações MSA+SSA só CRY / só APO
    n_MSA_SSA_only_CRY = sum(aggs$Species_composition %in% c("MSA","SSA") & aggs$only_CRY),
    n_MSA_SSA_only_APO = sum(aggs$Species_composition %in% c("MSA","SSA") & aggs$only_APO)
  )
}

###############################################
# 6. INDIVIDUAL-LEVEL METRICS
###############################################

calc_individual_metrics <- function(data) {
  
  tmp <- data %>%
    group_by(Site, Day) %>%
    mutate(
      Species_composition_group = first(Species_composition),
      only_CRY_group  = all(Color == "CRY"),
      only_APO_group  = all(Color == "APO"),
      n_colors_group  = n_distinct(Color)
    ) %>%
    ungroup() %>%
    mutate(
      msa_two_colors_group = Species_composition_group == "MSA" &
                             n_colors_group >= 2
    )
  
  tibble(
    indiv_MSA          = sum(tmp$Species_composition_group == "MSA"),
    indiv_SSA          = sum(tmp$Species_composition_group == "SSA"),
    indiv_SOLITARY     = sum(tmp$Species_composition_group == "SOLITARY"),
    
    indiv_MSA_CRY      = sum(tmp$Species_composition_group == "MSA" &
                             tmp$Color == "CRY"),
    indiv_MSA_APO      = sum(tmp$Species_composition_group == "MSA" &
                             tmp$Color == "APO"),
    
    indiv_SSA_CRY      = sum(tmp$Species_composition_group == "SSA" &
                             tmp$Color == "CRY"),
    indiv_SSA_APO      = sum(tmp$Species_composition_group == "SSA" &
                             tmp$Color == "APO"),
    
    indiv_SOLITARY_CRY = sum(tmp$Species_composition_group == "SOLITARY" &
                             tmp$Color == "CRY"),
    indiv_SOLITARY_APO = sum(tmp$Species_composition_group == "SOLITARY" &
                             tmp$Color == "APO"),
    
    indiv_MSA_two_colors = sum(tmp$msa_two_colors_group),
    
    # individuals in CRY-only or APO-only MSA/SSA aggregations
    indiv_MSA_SSA_only_CRY = sum(tmp$Species_composition_group %in% c("MSA","SSA") &
                                 tmp$only_CRY_group),
    indiv_MSA_SSA_only_APO = sum(tmp$Species_composition_group %in% c("MSA","SSA") &
                                 tmp$only_APO_group)
  )
}

###############################################
# 7. OBSERVED METRICS
###############################################

# 7.1 Aggregation level
obs_aggs  <- calc_aggregation_metrics(data_obs)

# 7.2 Indvidual level
obs_indiv <- calc_individual_metrics(data_obs)

###############################################
# 9. RUN TEMPORAL NULL MODEL SIMULATIONS
###############################################

set.seed(123)
n_sims <- 1000

###############################################
# 9.1 Temporal null model
###############################################

null_temporal <- map_dfr(1:n_sims, function(i) {
  message("Simulação temporal: ", i)
  null_data <- generate_null_model_temporal(data_clean, run_id = i)
  
  aggs_val  <- calc_aggregation_metrics(null_data)
  indiv_val <- calc_individual_metrics(null_data)
  
  tibble(run = i) %>%
    bind_cols(indiv_val) %>%
    bind_cols(aggs_val),
})

head(null_temporal)


###############################################
# 10. AUXILIARY FUNCTION FOR P-VALUES
###############################################

compute_pvals <- function(obs_val, null_vec, metric_name) {
  mean_null <- mean(null_vec)
  sd_null   <- sd(null_vec)
  
  p_lower   <- mean(null_vec >= obs_val)  # hypothesis: observed < expected
  p_upper   <- mean(null_vec <= obs_val)  # hypothesis: observed > expected
  p_two     <- 2 * min(p_lower, p_upper)
  p_two     <- min(p_two, 1)
  
  tibble(
    metric      = metric_name,
    obs         = obs_val,
    mean_null   = mean_null,
    sd_null     = sd_null,
    p_lower     = p_lower,
    p_upper     = p_upper,
    p_two_sided = p_two
  )
}

###############################################
# 11. METRIC LISTS
###############################################

metrics_indiv <- c(
  "indiv_MSA",
  "indiv_SSA",
  "indiv_SOLITARY",
  
  "indiv_MSA_CRY",
  "indiv_MSA_APO",
  "indiv_SSA_CRY",
  "indiv_SSA_APO",
  "indiv_SOLITARY_CRY",
  "indiv_SOLITARY_APO",
  
  "indiv_MSA_two_colors",
  "indiv_MSA_SSA_only_CRY",
  "indiv_MSA_SSA_only_APO"
)

metrics_aggs <- c(
  "n_MSA",
  "n_SSA",
  "n_MSA_only_CRY",
  "n_MSA_only_APO",
  "n_SSA_only_CRY",
  "n_SSA_only_APO",
  "n_MSA_two_colors",
  "n_MSA_SSA_only_CRY",
  "n_MSA_SSA_only_APO"
)



###############################################
# 12. P-VALUES – TEMPORAL MODEL
###############################################

pvals_indiv_temporal <- map_dfr(metrics_indiv, function(m) {
  compute_pvals(
    obs_val     = obs_indiv[[m]][1],
    null_vec    = null_temporal[[m]],
    metric_name = paste0(m, " (temporal)")
  )
})

pvals_aggs_temporal <- map_dfr(metrics_aggs, function(m) {
  compute_pvals(
    obs_val     = obs_aggs[[m]][1],
    null_vec    = null_temporal[[m]],
    metric_name = paste0(m, " (temporal)")
  )
})


pvals_temporal <- bind_rows(
  pvals_indiv_temporal,
  pvals_aggs_temporal
)


###############################################
# 13. FINAL FORMATTING: 3 FIXED DECIMAL PLACES
###############################################

format_3_decimals <- function(df) {
  df %>% 
    mutate(across(
      .cols = where(is.numeric),
      .fns  = ~ sprintf("%.3f", .)
    ))
}

# Apply formatting to the temporal model results
pvals_temporal_fmt   <- format_3_decimals(pvals_temporal)

# View formatted results
pvals_temporal_fmt


library(writexl)
null_temporal = as_tibble(null_temporal)

write.csv(null_temporal, "C:/Users/bruno/OneDrive/Documentos/trabalhos_extra/Opilhões/null_temporal.csv")
write.csv(pvals_temporal_fmt, "C:/Users/bruno/OneDrive/Documentos/trabalhos_extra/Opilhões/pvals_temporal_fmt.csv")     

p = read.csv("C:/Users/bruno/OneDrive/Documentos/trabalhos_extra/Opilhões/pvals_temporal_fmt.csv")
p
