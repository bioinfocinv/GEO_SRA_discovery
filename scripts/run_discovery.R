#!/usr/bin/env Rscript
# Usage:
# Rscript scripts/run_discovery.R
# Rscript scripts/run_discovery.R parkinson
# Rscript scripts/run_discovery.R > outputs/log.txt 2>&1 &

source("R/discover_geo_sra.R")
disease_id <- commandArgs(trailingOnly = TRUE)
result <- run_discovery("config/discovery.yml", if (length(disease_id)) disease_id[[1]] else NULL)
print(result |> select(gse, srp, title, n_case, n_control, n_srr, n_paired_srr, decision, exclusion_reasons, score), n = Inf)
