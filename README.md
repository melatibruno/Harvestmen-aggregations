# Harvestmen aggregation assortment analyses

This repository contains the R scripts used to evaluate species and color
assortment in *Prionostemma* aggregations using permutation-based null models.

The analyses are divided into two hierarchical components:

1. **Overall assortment by species identity**
2. **Color assortativity conditional on multispecies aggregation (MSA) formation**

Both analyses use the same general randomization procedure. Within each sampling
day, individual records are randomly reassigned among the observed aggregation
positions. Species identity and color phenotype are shuffled jointly because
each species has a fixed color phenotype.

This procedure preserves daily species abundances, the number of aggregation
sites, and aggregation sizes, while breaking the association between species
identity and aggregation membership.

Each analysis uses 1,000 randomizations with a fixed random seed (`seed = 123`).


## Scripts

### `01_species_assortment_null_model.R`

Tests overall assortment by species identity at two hierarchical levels.

**Aggregation level**

- Number of multispecies aggregations (MSA)
- Number of single-species aggregations (SSA)

**Individual-record level**

- Cryptic records occurring in MSA
- Cryptic records occurring in SSA
- Conspicuous records occurring in MSA
- Conspicuous records occurring in SSA

For each response variable, the observed value is compared with its null
distribution using a standardized effect size (SES) and a two-tailed
Monte Carlo P-value.

This analysis provides the results represented in Figure 2.


### `02_color_assortativity_P1_P6.R`

Tests color assortativity conditional on MSA formation using six response
statistics.

**Aggregation level**

- **P1:** proportion of all MSA classified as mixed-color
- **P2:** among MSA containing cryptic species, proportion also containing
  conspicuous species
- **P3:** among MSA containing conspicuous species, proportion also containing
  cryptic species

**Individual-record level**

- **P4:** proportion of all individual records within MSA occurring in
  mixed-color MSA
- **P5:** among cryptic individual records within MSA, proportion occurring in
  mixed-color MSA
- **P6:** among conspicuous individual records within MSA, proportion occurring
  in mixed-color MSA

For each statistic, the observed value is compared with the corresponding null
distribution using SES and a two-tailed Monte Carlo P-value.

These analyses provide the results represented in Figure 3.


## Null-model inference

For each response variable or prediction, the standardized effect size is
calculated as:

SES = (observed - mean(null)) / SD(null)

Two-tailed Monte Carlo P-values are calculated from the upper and lower tails
of the null distribution using the +1 correction:

p_upper = (number of null values >= observed + 1) / (N + 1)

p_lower = (number of null values <= observed + 1) / (N + 1)

p_two-sided = min(2 × min(p_upper, p_lower), 1)


## Input data

Both scripts use:

`dados.xlsx`

The analyses require the following variables:

- `Site` — aggregation site
- `Day` — sampling day
- `Species` — species identity
- `Color` — color phenotype (`CRY` or `APO`)


## Aggregation classification

For each `Site × Day` combination:

- **Solitary:** one individual record
- **SSA:** more than one individual, all belonging to the same species
- **MSA:** more than one individual and more than one species


## Reproducibility

The analysis scripts save the complete null distributions as `.rds` files.
These files can be used to regenerate figures without rerunning the 1,000
randomizations.

Figure construction is kept in a separate R script and does not alter the
statistical analyses.
