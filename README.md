# Harvestmen-aggregations
Data and analysis scripts supporting the manuscript "Cryptic and conspicuous species segregate in mixed-species aggregations of a Neotropical arachnid". 

The raw data  used in this study are available in data.xlsx. The data are composed by the colluns:

Site <- Sampling site;

Site2 <- Not used;

Day <- Sampling day;

Species <- Species identity;

Color <- Species phenotype (CRY = cryptic, APO = conspicuous/aposematic);

nickname_species <- Species dorsal coloration pattern used for visual classification;

Group_size <- Number of recorded individuals in the aggregation;

Solo_pair_group <- Categorical size of aggregation (Solo = one recorded individual, pair = two recorded individuals, group = three or more recorded individuals);

Species_composition <- Aggregation type (MSA = mixed species aggregation, SSA = same species aggregation).

R scripts used in this study are available in membership_null_model and aggregation_null_model. Scripts used to perform all null model analyses, calculate aggregation and individual-level metrics, estimate transition matrices, and compute p-values for temporal and site-fixed null models for membership and aggregations tests.
