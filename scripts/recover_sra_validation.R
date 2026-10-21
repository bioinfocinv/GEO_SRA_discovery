#!/usr/bin/env Rscript
# Usage:
# Rscript scripts/recover_sra_validation.R alzheimer > outputs/alzheimer/recover_sra.log 2>&1 &

source("R/discover_geo_sra.R")
disease_id <- commandArgs(trailingOnly = TRUE)
result <- run_sra_validation_from_previous("config/discovery.yml", if (length(disease_id)) disease_id[[1]] else NULL)
print(result |> select(gse, srp, title, n_case, n_control, n_srr, n_paired_srr, decision, exclusion_reasons, score), n = Inf)
