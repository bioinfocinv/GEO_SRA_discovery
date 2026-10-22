#!/usr/bin/env Rscript
source("R/discover_geo_sra.R")
args <- commandArgs(trailingOnly = TRUE)
disease_id <- if (length(args) && args[[1]] != "all") args[[1]] else NULL
all_studies <- "all" %in% args
result <- audit_geo_supplements("config/discovery.yml", disease_id = disease_id, all_studies = all_studies)
print(result |> select(gse, title, has_count_like_file, has_ncbi_named_count_file, has_metadata_like_file, count_like_filenames), n = Inf)
