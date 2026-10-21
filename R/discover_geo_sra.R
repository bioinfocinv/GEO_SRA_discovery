suppressPackageStartupMessages({
  #library(dplyr);
  library(GEOquery);
  library(purrr);
  library(readr);
  library(rentrez);
  #library(stringr);
  #library(tibble)
  library(tidyverse);
  library(yaml)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
clean_chr <- function(x) str_squish(str_to_lower(paste(unlist(x), collapse = " | ")))
first_nonempty <- function(x) { x <- x[!is.na(x) & nzchar(x)]; if (length(x)) x[[1]] else NA_character_ }
has_any <- function(x, patterns) any(str_detect(x, regex(patterns, ignore_case = TRUE)), na.rm = TRUE)
empty_runs <- function() tibble(gse = character(), sra_query_accession = character(), srr = character(), srx = character(), srs = character(), srp = character(), bioproject = character(), biosample = character(), sample_name_sra = character(), gsm_sra = character(), library_layout = character(), library_strategy = character(), assay_type = character(), platform = character(), instrument = character(), release_date = character(), load_date = character(), center_name = character(), spots = character(), bases = character(), run_total_bases = character(), sra_found = logical())
empty_gsms <- function() tibble(gse = character(), gsm = character(), gsm_title = character(), source_name = character(), characteristics = character(), molecule = character(), organism = character(), gsm_text = character(), gse_relation = character(), geo_download_error = character())

# NCBI may occasionally return a transient non-E-utilities response that rentrez
# cannot parse. Retry it, but never let one failed GSE terminate a long discovery run.
ncbi_retry <- function(request, attempts = 3L, pause_seconds = 1) {
  errors <- character()
  for (attempt in seq_len(attempts)) {
    response <- tryCatch(request(), error = function(e) e)
    if (!inherits(response, "error")) return(list(ok = TRUE, value = response, attempts = attempt))
    errors <- c(errors, conditionMessage(response))
    if (attempt < attempts) Sys.sleep(pause_seconds * attempt)
  }
  list(ok = FALSE, value = NULL, attempts = attempts, error = paste(unique(errors), collapse = " | "))
}

ncbi_init <- function(ncbi) {
  if (identical(ncbi$email, "replace-with-your-email@example.org"))
    stop("Set ncbi.email in config/discovery.yml before querying NCBI.")
  # rentrez 1.2.4 does not export global email/tool setters. They are supplied
  # explicitly to every request below, which also works with older releases.
  key <- Sys.getenv(ncbi$api_key_env %||% "NCBI_API_KEY", unset = "")
  if (nzchar(key)) rentrez::set_entrez_key(key)
}

# A disease profile supplies only condition-specific search/classification terms;
# shared sequencing and QC rules remain in config$screen.
resolve_config <- function(cfg) {
  id <- cfg$active_disease %||% "alzheimer"
  profile <- cfg$diseases[[id]]
  if (is.null(profile)) stop("active_disease must name an entry in diseases: ", id)
  required <- c("label", "query_term", "case_patterns")
  missing <- required[!required %in% names(profile)]
  if (length(missing)) stop("Disease profile '", id, "' lacks: ", paste(missing, collapse = ", "))
  cfg$query$disease <- profile$label
  cfg$query$geo_term <- paste(profile$query_term,
    sprintf('"%s"[Organism]', cfg$query$organism),
    '"expression profiling by high throughput sequencing"[DataSet Type]', sep = " AND ")
  cfg$screen$case_patterns <- profile$case_patterns
  cfg$disease_id <- id
  cfg$output_dir <- file.path(cfg$output_dir, id)
  cfg
}

# ESearch returns GEO DataSets (GDS) UIDs; ESummary maps them to GSE accessions.
find_geo_series <- function(query, batch_size, ncbi) {
  hit <- entrez_search(db = "gds", term = query$geo_term, retmax = query$retmax, use_history = TRUE, email = ncbi$email, tool = ncbi$tool)
  if (!length(hit$ids)) return(tibble(geo_uid = character(), gse = character(), title = character(), summary = character(), gds_type = character(), n_samples_geo = integer(), pub_date = character()))
  ids <- hit$ids
  chunks <- split(ids, ceiling(seq_along(ids) / batch_size))
  map_dfr(chunks, function(chunk) {
    sm <- entrez_summary(db = "gds", id = chunk, email = ncbi$email, tool = ncbi$tool)
    map_dfr(sm, function(x) tibble(
      geo_uid = as.character(x$uid), gse = x$accession %||% NA_character_,
      title = x$title %||% NA_character_, summary = x$summary %||% NA_character_,
      gds_type = x$gdstype %||% NA_character_, n_samples_geo = suppressWarnings(as.integer(x$n_samples %||% NA)),
      pub_date = x$PDAT %||% NA_character_
    ))
  }) |> filter(str_detect(gse, "^GSE")) |> distinct(gse, .keep_all = TRUE)
}

gsm_metadata <- function(gse) {
  obj <- tryCatch(getGEO(gse, GSEMatrix = FALSE, getGPL = FALSE), error = function(e) e)
  if (inherits(obj, "error")) return(tibble(gse = gse, gsm = NA_character_, gsm_title = NA_character_, source_name = NA_character_, characteristics = NA_character_, molecule = NA_character_, organism = NA_character_, gsm_text = NA_character_, gse_relation = NA_character_, geo_download_error = conditionMessage(obj)))
  gsms <- GSMList(obj)
  relation <- clean_chr(Meta(obj)$relation)
  imap_dfr(gsms, function(gsm_obj, gsm_id) {
    m <- Meta(gsm_obj)
    tibble(gse = gse, gsm = gsm_id,
      gsm_title = clean_chr(m$title), source_name = clean_chr(m$source_name_ch1),
      characteristics = clean_chr(m$characteristics_ch1),
      molecule = clean_chr(m$molecule_ch1), organism = clean_chr(m$organism_ch1),
      gsm_text = clean_chr(c(m$title, m$source_name_ch1, m$characteristics_ch1)),
      gse_relation = relation,
      geo_download_error = NA_character_)
  })
}

classify_gsms <- function(gsm, screen) {
  case_hit <- map_lgl(gsm$gsm_text, has_any, patterns = screen$case_patterns)
  ctrl_hit <- map_lgl(gsm$gsm_text, has_any, patterns = screen$control_patterns)
  excluded <- map_lgl(gsm$gsm_text, has_any, patterns = screen$exclude_patterns)
  review <- map_lgl(gsm$gsm_text, has_any, patterns = screen$review_patterns)
  phenotype <- case_when(case_hit & !ctrl_hit ~ "case", ctrl_hit & !case_hit ~ "control",
    case_hit & ctrl_hit ~ "ambiguous", TRUE ~ "unclassified")
  evidence <- map2_chr(gsm$gsm_text, phenotype, function(txt, label) {
    pats <- switch(label, case = screen$case_patterns, control = screen$control_patterns, character())
    str_c(pats[map_lgl(pats, ~ str_detect(txt, regex(.x, ignore_case = TRUE)))], collapse = ";")
  })
  gsm |> mutate(
    phenotype = phenotype,
    modality = if_else(excluded, "excluded_keyword", "bulk_candidate"),
    review_flag = review,
    phenotype_evidence = na_if(evidence, "")
  )
}

annotate_tissue <- function(gsm, screen) {
  patterns <- screen$tissue_patterns %||% list()
  if (!length(patterns)) return(gsm |> mutate(tissue_primary = "unclassified", tissue_labels = "unclassified", tissue_evidence = NA_character_))
  hits <- map(gsm$gsm_text, function(txt) {
    names(patterns)[map_lgl(patterns, ~ has_any(txt, .x))]
  })
  tissue_evidence <- map2_chr(gsm$gsm_text, hits, function(txt, labels) {
    if (!length(labels)) return(NA_character_)
    matched_patterns <- unlist(patterns[labels], use.names = FALSE)
    str_c(matched_patterns[map_lgl(matched_patterns, ~ str_detect(txt, regex(.x, ignore_case = TRUE)))], collapse = ";")
  })
  gsm |> mutate(
    tissue_primary = map_chr(hits, ~ if (length(.x)) .x[[1]] else "unclassified"),
    tissue_labels = map_chr(hits, ~ if (length(.x)) str_c(.x, collapse = ";") else "unclassified"),
    tissue_evidence = tissue_evidence
  )
}

# Extract stable SRA/BioProject identifiers from the GEO Series relation field.
# Prefer SRP, which is an SRA study accession, then use BioProject as a fallback.
sra_targets <- function(gsms) {
  gsms |> group_by(gse) |> summarise(gse_relation = first_nonempty(gse_relation), .groups = "drop") |> mutate(
    srp_relation = str_extract(gse_relation, "\\bSRP[0-9]+\\b"),
    bioproject_relation = str_extract(gse_relation, "\\bPRJ(?:NA|EB|DB)[0-9]+\\b"),
    sra_query_accession = coalesce(srp_relation, bioproject_relation),
    sra_link_source = case_when(!is.na(srp_relation) ~ "GEO_SRA_relation", !is.na(bioproject_relation) ~ "GEO_BioProject_relation", TRUE ~ "no_GEO_SRA_or_BioProject_relation")
  )
}

# SRA efetch format=runinfo is the programmatic equivalent of the RunInfo table.
sra_runinfo <- function(gse, sra_query_accession, ncbi) {
  if (is.na(sra_query_accession) || !nzchar(sra_query_accession)) return(tibble(gse = gse, sra_query_accession = NA_character_, sra_found = FALSE, sra_status = "no_geo_sra_relation", sra_attempts = 0L, sra_error = NA_character_))
  q <- paste0(sra_query_accession, "[All Fields]")
  request_config <- httr::timeout(ncbi$sra_timeout_seconds %||% 60)
  attempts <- ncbi$sra_retry_attempts %||% 3L
  retmax <- ncbi$sra_search_retmax %||% 10000L
  searched <- ncbi_retry(function() entrez_search(db = "sra", term = q, retmax = retmax, use_history = FALSE, email = ncbi$email, tool = ncbi$tool, config = request_config), attempts = attempts)
  if (!searched$ok) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "search_failed", sra_attempts = searched$attempts, sra_error = searched$error))
  h <- searched$value
  if (!length(h$ids)) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "no_sra_runs", sra_attempts = searched$attempts, sra_error = NA_character_))
  if (h$count > retmax) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "sra_result_limit_exceeded", sra_attempts = searched$attempts, sra_error = paste0("SRA search found ", h$count, " runs; configured sra_search_retmax is ", retmax)))
  fetched <- ncbi_retry(function() entrez_fetch(db = "sra", id = h$ids, rettype = "runinfo", retmode = "text", email = ncbi$email, tool = ncbi$tool, config = request_config), attempts = attempts)
  if (!fetched$ok) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "runinfo_failed", sra_attempts = searched$attempts + fetched$attempts, sra_error = fetched$error))
  x <- tryCatch(suppressMessages(read_csv(I(fetched$value), show_col_types = FALSE, name_repair = "minimal")), error = function(e) e)
  if (inherits(x, "error")) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "runinfo_parse_failed", sra_attempts = searched$attempts + fetched$attempts, sra_error = conditionMessage(x)))
  if (!nrow(x)) return(tibble(gse = gse, sra_query_accession = sra_query_accession, sra_found = FALSE, sra_status = "no_sra_runs", sra_attempts = searched$attempts + fetched$attempts, sra_error = NA_character_))
  x |> mutate(gse = gse, sra_query_accession = sra_query_accession, sra_found = TRUE, sra_status = "ok", sra_attempts = searched$attempts + fetched$attempts, sra_error = NA_character_, across(everything(), as.character))
}

normalise_runs <- function(runinfo) {
  if (!nrow(runinfo) || !"Run" %in% names(runinfo)) return(empty_runs())
  pick <- function(nm) if (nm %in% names(runinfo)) runinfo[[nm]] else NA_character_
  geo_col <- first_nonempty(intersect(c("GEO_Accession", "GEO Accession", "GEO_Accession (exp)"), names(runinfo)))
  tibble(gse = runinfo$gse, sra_query_accession = pick("sra_query_accession"), srr = pick("Run"), srx = pick("Experiment"), srs = pick("Sample"),
    srp = pick("SRAStudy"), bioproject = pick("BioProject"), biosample = pick("BioSample"),
    sample_name_sra = pick("SampleName"), gsm_sra = if (!is.na(geo_col)) runinfo[[geo_col]] else NA_character_, library_layout = str_to_upper(pick("LibraryLayout")),
    library_strategy = pick("LibraryStrategy"), assay_type = pick("Assay Type"),
    platform = pick("Platform"), instrument = pick("Model"), release_date = pick("ReleaseDate"), load_date = pick("LoadDate"), center_name = pick("CenterName"), spots = pick("spots"), bases = pick("bases"),
    run_total_bases = pick("size_MB"), sra_found = TRUE) |> filter(!is.na(srr), nzchar(srr))
}

match_runs_to_gsm <- function(gsm, runs) {
  if (!nrow(runs)) return(gsm |> mutate(srr = NA_character_, srx = NA_character_, srs = NA_character_, srp = NA_character_, library_layout = NA_character_))
  # Exact GEO GSM identifiers are preferred. Sample-name matching is only a fallback,
  # because GEO and SRA labels often differ. Unmatched runs remain visible in sra_runs.csv.
  with_geo <- gsm |> left_join(runs |> select(-gse) |> filter(!is.na(gsm_sra), nzchar(gsm_sra)), by = c("gsm" = "gsm_sra"), suffix = c("", "_sra"))
  unmatched <- with_geo |> filter(is.na(srr)) |> select(names(gsm))
  by_name <- unmatched |> left_join(runs |> select(-gse), by = c("gsm" = "sample_name_sra"), suffix = c("", "_sra"))
  bind_rows(with_geo |> filter(!is.na(srr)), by_name) |> mutate(match_method = if_else(!is.na(srr), if_else(!is.na(gsm_sra), "GEO_accession", "sample_name"), "unmatched"))
}

screen_study <- function(studies, gsms, runs, cfg) {
  gsm_counts <- gsms |> group_by(gse) |> summarise(
    n_gsm = n(), n_case = sum(phenotype == "case"), n_control = sum(phenotype == "control"),
    n_ambiguous = sum(phenotype == "ambiguous"), n_unclassified = sum(phenotype == "unclassified"),
    n_excluded_keyword = sum(modality != "bulk_candidate"), n_review_flag = sum(review_flag), .groups = "drop")
  tissue_counts <- gsms |> group_by(gse) |> summarise(
    tissue_primary_labels = str_c(sort(unique(tissue_primary)), collapse = ";"),
    n_tissue_primary_labels = n_distinct(tissue_primary), .groups = "drop")
  run_counts <- runs |> group_by(gse) |> summarise(
    n_srr = n(), n_paired_srr = sum(library_layout == "PAIRED", na.rm = TRUE),
    n_single_srr = sum(library_layout == "SINGLE", na.rm = TRUE), n_unknown_layout = sum(is.na(library_layout) | !library_layout %in% c("PAIRED", "SINGLE")),
    srp = first_nonempty(srp), .groups = "drop")
  studies |> left_join(gsm_counts, by = "gse") |> left_join(tissue_counts, by = "gse") |> left_join(run_counts, by = "gse") |> mutate(across(starts_with("n_"), ~coalesce(.x, 0L)),
    pass_groups = n_case >= cfg$screen$min_case_samples & n_control >= cfg$screen$min_control_samples,
    pass_bulk = n_excluded_keyword == 0,
    pass_phenotype = if (isTRUE(cfg$screen$require_complete_phenotype)) n_ambiguous + n_unclassified == 0 else TRUE,
    pass_layout = if (isTRUE(cfg$screen$require_paired_end)) n_paired_srr > 0 & (!isTRUE(cfg$screen$require_all_runs_paired) | (n_single_srr == 0 & n_unknown_layout == 0)) else TRUE,
    decision = case_when(pass_groups & pass_bulk & pass_phenotype & pass_layout ~ "eligible", TRUE ~ "review_or_exclude"),
    exclusion_reasons = pmap_chr(list(pass_groups, pass_bulk, pass_phenotype, pass_layout), ~ str_c(c(if (!..1) "too_few_case_or_control", if (!..2) "single_cell_or_other_exclusion_keyword", if (!..3) "incomplete_or_ambiguous_phenotype", if (!..4) "not_confirmed_paired_end"), collapse = ";")),
    score = 100 + 4 * pmin(n_case, n_control) + 2 * n_paired_srr - 15 * n_review_flag - 20 * (n_ambiguous + n_unclassified) - 30 * n_excluded_keyword) |> arrange(desc(decision == "eligible"), desc(score), desc(n_case + n_control))
}

run_discovery <- function(config_path = "config/discovery.yml", disease_id = NULL) {
  cfg <- read_yaml(config_path)
  if (!is.null(disease_id)) cfg$active_disease <- disease_id
  cfg <- resolve_config(cfg); ncbi_init(cfg$ncbi); dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  studies <- find_geo_series(cfg$query, cfg$query$batch_size, cfg$ncbi)
  gsms <- (if (nrow(studies)) map_dfr(studies$gse, gsm_metadata) else empty_gsms()) |> classify_gsms(cfg$screen) |> annotate_tissue(cfg$screen)
  pre_sra_table <- screen_study(studies, gsms, empty_runs(), cfg)
  targets <- sra_targets(gsms)
  if (isTRUE(cfg$ncbi$sra_pre_screen_only)) targets <- targets |> inner_join(pre_sra_table |> filter(pass_groups, pass_bulk, pass_phenotype) |> select(gse), by = "gse")
  raw_runs <- if (nrow(targets)) pmap_dfr(targets |> select(gse, sra_query_accession), function(gse, sra_query_accession) {
    index <- match(gse, targets$gse)
    message(sprintf("SRA RunInfo: %s via %s (%d of %d GSEs)", gse, sra_query_accession %||% "no GEO SRA relation", index, nrow(targets)))
    sra_runinfo(gse, sra_query_accession, ncbi = cfg$ncbi)
  }) else tibble(); runs <- normalise_runs(raw_runs)
  study_table <- screen_study(studies, gsms, runs, cfg)
  run_table <- match_runs_to_gsm(gsms, runs)
  write_csv(study_table, file.path(cfg$output_dir, "study_screening.csv"))
  write_csv(gsms, file.path(cfg$output_dir, "gsm_phenotype_audit.csv"))
  write_csv(runs, file.path(cfg$output_dir, "sra_runs.csv"))
  sra_failures <- if ("sra_error" %in% names(raw_runs)) raw_runs |> filter(!is.na(sra_error)) |> select(gse, sra_status, sra_attempts, sra_error) |> distinct() else tibble(gse = character(), sra_status = character(), sra_attempts = character(), sra_error = character())
  write_csv(sra_failures, file.path(cfg$output_dir, "sra_query_failures.csv"))
  write_csv(run_table, file.path(cfg$output_dir, "fetch_sra_manifest.csv"))
  write_yaml(cfg, file.path(cfg$output_dir, "discovery_config_used.yml"))
  #saveRDS(list(studies = study_table, gsms = gsms, runs = runs, config = cfg, session = sessionInfo()), file.path(cfg$output_dir, "discovery_bundle.rds"))
  invisible(study_table)
}

# Recovery path: reuse a completed GEO/GSM screen and re-run SRA linkage only
# for the GSEs that passed the phenotype and bulk pre-screen.
run_sra_validation_from_previous <- function(config_path = "config/discovery.yml", disease_id = NULL) {
  cfg <- read_yaml(config_path)
  if (!is.null(disease_id)) cfg$active_disease <- disease_id
  cfg <- resolve_config(cfg); ncbi_init(cfg$ncbi)
  required <- c("study_screening.csv", "gsm_phenotype_audit.csv")
  paths <- file.path(cfg$output_dir, required)
  if (!all(file.exists(paths))) stop("Run the GEO discovery stage first; missing: ", paste(required[!file.exists(paths)], collapse = ", "))
  previous <- read_csv(file.path(cfg$output_dir, "study_screening.csv"), show_col_types = FALSE)
  gsms <- read_csv(file.path(cfg$output_dir, "gsm_phenotype_audit.csv"), show_col_types = FALSE)
  if (!"tissue_primary" %in% names(gsms)) gsms <- annotate_tissue(gsms, cfg$screen)
  base_studies <- previous |> select(any_of(c("geo_uid", "gse", "title", "summary", "gds_type", "n_samples_geo", "pub_date")))
  candidates <- previous |> filter(pass_groups, pass_bulk, pass_phenotype) |> pull(gse)
  message(sprintf("SRA recovery: %d pre-screened GSEs", length(candidates)))
  relation_rows <- map_dfr(candidates, gsm_metadata)
  targets <- sra_targets(relation_rows) |> filter(gse %in% candidates)
  raw_runs <- if (nrow(targets)) pmap_dfr(targets |> select(gse, sra_query_accession), function(gse, sra_query_accession) {
    index <- match(gse, targets$gse)
    message(sprintf("SRA RunInfo: %s via %s (%d of %d GSEs)", gse, sra_query_accession %||% "no GEO SRA relation", index, nrow(targets)))
    sra_runinfo(gse, sra_query_accession, ncbi = cfg$ncbi)
  }) else tibble()
  runs <- normalise_runs(raw_runs)
  study_table <- screen_study(base_studies, gsms, runs, cfg)
  run_table <- match_runs_to_gsm(gsms, runs)
  sra_failures <- if ("sra_error" %in% names(raw_runs)) raw_runs |> filter(!is.na(sra_error)) |> select(gse, sra_query_accession, sra_status, sra_attempts, sra_error) |> distinct() else tibble(gse = character(), sra_query_accession = character(), sra_status = character(), sra_attempts = character(), sra_error = character())
  write_csv(study_table, file.path(cfg$output_dir, "study_screening.csv"))
  write_csv(gsms, file.path(cfg$output_dir, "gsm_phenotype_audit.csv"))
  write_csv(runs, file.path(cfg$output_dir, "sra_runs.csv"))
  write_csv(run_table, file.path(cfg$output_dir, "fetch_sra_manifest.csv"))
  write_csv(sra_failures, file.path(cfg$output_dir, "sra_query_failures.csv"))
  write_yaml(cfg, file.path(cfg$output_dir, "discovery_config_used.yml"))
  #saveRDS(list(studies = study_table, gsms = gsms, runs = runs, config = cfg, session = sessionInfo()), file.path(cfg$output_dir, "discovery_bundle.rds"))
  invisible(study_table)
}
