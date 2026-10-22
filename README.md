# GEO_SRA_discovery
An R pipeline for systematic screening of GEO/SRA studies associated with neurodegenerative diseases.

This workflow discovers candidate human Alzheimer's disease GEO Series, links each GSE to SRA runs, audits GSM-level phenotype evidence, and writes a ranked, reviewable manifest before fetching (`fetch_SRA`) or `edgeR::DGEList` construction.

## Run

1. In `config/discovery.yml`, set `active_disease` to `alzheimer`, `huntington`, or `parkinson`; then adjust the corresponding disease profile, shared control vocabulary, and minimum group sizes. Add another key beneath `diseases:` to support an additional condition.

2. From the project root run:

   ```bash
   Rscript scripts/run_discovery.R
   ```

   Or run a supplied profile without editing the file:

   ```bash
   Rscript scripts/run_discovery.R alzheimer
   ```

Each disease writes to its own directory (for example `outputs/alzheimer/`), so results are not overwritten. Outputs are deterministic for a given NCBI response and configuration. Preserve `discovery_config_used.yml` with run parameters.

### Recovering a completed GEO screen

If an earlier run produced `study_screening.csv` and `gsm_phenotype_audit.csv` but has empty SRA results, do not repeat the full GEO scan. After updating the workflow, run:

```bash
Rscript scripts/recover_sra_validation.R alzheimer
```

This reuses the GSM screening result and retrieves GEO relations/SRA RunInfo only for studies passing the non-SRA gates. It overwrites the SRA-derived tables and recalculates `study_screening.csv`.

### Auditing GEO supplementary counts and metadata

The main discovery and recovery workflows audit every discovered GSE by default. To repeat the check independently—without downloading data—run:

```bash
Rscript scripts/audit_geo_supplementary.R alzheimer
```

It writes `geo_supplementary_files.csv` and adds availability columns to `study_screening.csv`: `has_count_like_file`, `has_ncbi_named_count_file`, `has_metadata_like_file`, and `count_like_filenames`. `has_ncbi_named_count_file` is strict: it recognizes GEO's NCBI-generated pattern `*_raw_counts_<assembly>_NCBI.tsv.gz` (for example, `GSE105041_raw_counts_GRCh38.p13_NCBI.tsv.gz`). Set `supplementary_audit.pre_screen_only: true` to inspect only pre-screened candidates on future runs; that is faster but leaves availability fields as `NA` for all other GSEs. Add `all` to override that setting in an independent audit, for example `Rscript scripts/audit_geo_supplementary.R alzheimer all`.

## What it queries

1. NCBI E-utilities searches the **GEO DataSets** (`gds`) database for the configured disease, organism, and high-throughput expression profiling term. It uses `ESearch` + `ESummary` to obtain candidate GSEs. The configured NCBI email and tool name are sent explicitly with every request (rather than relying on `rentrez` global setters).
2. `GEOquery` downloads each GSE SOFT record and extracts each GSM's title, source name, and `characteristics_ch1` fields.
3. NCBI SRA is searched for the GSE accession; `EFetch(rettype = "runinfo")` retrieves the same run-level family of fields exposed by the SRA Run Selector/RunInfo export. This produces GSE ↔ SRP ↔ SRX ↔ SRS ↔ SRR identifiers, library layout, strategy, and platform.

NCBI documents E-utilities as its programmatic interface and recommends identifying the tool and email (and using an API key for higher request rates). The SRA documentation describes RunInfo/Run Selector metadata exports. GEOquery supports parsing the more complete SOFT GSE records used here. [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), [SRA metadata export](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/), [GEOquery documentation](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html).

## Study screening

`study_screening.csv` has one row per GSE and is ranked with the following explicit gates:

| Gate | Pass condition |
|---|---|
| RNA-seq discovery | GEO data type is high-throughput expression profiling (candidate filter, not proof) |
| Bulk | no configured single-cell/single-nucleus/exclusion keyword in GSM fields |
| Disease/control | at least configured `min_case_samples` and `min_control_samples` from GSM text |
| Phenotype completeness | no ambiguous or unclassified GSMs when `require_complete_phenotype: true` |
| Paired-end | SRA `LibraryLayout == PAIRED`; optionally every run must be paired with no unknown layouts |

`decision = eligible` only if all gates pass. The score ranks eligible studies by balanced group size and number of paired runs, while penalizing review/ambiguity/exclusion flags; it is a triage order, not scientific evidence.

## Output files

| File | Granularity | Key fields/use |
|---|---|---|
| `study_screening.csv` | GSE | `gse`, `srp`, title, inferred tissue labels, group/run counts, gate booleans, `decision`, `exclusion_reasons`, score |
| `gsm_phenotype_audit.csv` | GSM | raw phenotype text, inferred phenotype, tissue labels/evidence, matching-pattern evidence, modality and review flags |
| `sra_runs.csv` | SRR | GSE/SRP/SRX/SRS/SRR IDs, layout, assay, platform, instrument, SRA release/load dates, sequencing center |
| `sra_query_failures.csv` | GSE | SRA search/fetch/parse failures after automatic retries; rerun or inspect these GSEs manually |
| `geo_supplementary_files.csv` | GSE × file | public GEO supplementary-file index and transparent count/metadata filename flags |
| `fetch_sra_manifest.csv` | GSM joined to SRR | reviewable input to `fetch_SRA`; do not download rows until the crosswalk is verified |

## Validation

Before analysis, review each eligible GSE against its GEO record/publication and record the decision:

- Confirm tissue/brain region, species, bulk library preparation, and that no single-cell, sorted-cell, pooled, cell-line, or organoid samples enter the contrast.
- Verify `case`/`control` from the study design, not just keyword matches; distinguish AD pathology, clinical AD, MCI, and other dementias. Explicitly mark exclusions and unresolved labels.
