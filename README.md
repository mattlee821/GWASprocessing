# GWASprocessing

GWASprocessing is a manifest-driven Nextflow workflow for staging GWAS summary statistics, standardising them to GWAS-VCF, and writing persistent outputs that can be resumed across repeated manifest updates.

## Layout

- `src/` contains user-facing helper scripts. Run the workflow with `bash src/standardise.sh` or submit it with `sbatch src/standardise.sh`.
- `pipeline/` contains the Nextflow workflow, pipeline Python utilities, profiles, and tests.
- `manifests/` contains persistent input manifests and optional study YAML overrides.
- `references/` is the root location for reference data generated or downloaded by scripts in `src/reference_data/`.
- `GWAS/<STUDY>/raw/<SOURCE_ID>/` receives staged raw GWAS files.
- `GWAS/<STUDY>/standardise/<SOURCE_ID>/` receives standardised GWAS-VCF outputs, metadata, QC JSON, and QC TIFF plots.
- `temp/` contains Nextflow runtime metadata, `.nextflow/`, `.nextflow.log*`, work directories, and run logs.
- `deprecated/` contains legacy one-off scripts and study directories. It is intentionally not tracked.

## Quick Start

Edit the root `.env` file to set paths/tokens as needed. It is local-only and ignored by git.

```bash
$EDITOR .env
```

Run the default manifest:

```bash
bash src/standardise.sh
```

On a local machine, the wrapper runs only one eligible manifest row per invocation. Re-run the same command to process the next incomplete row. Completed rows are recorded in `GWAS/.standardise_state.tsv` and skipped automatically.

On SLURM, use the `slurm` profile to fan out all eligible rows. Each manifest row is submitted as one row-level Nextflow job that performs staging, standardisation, optional QC plotting, and returns one state record.

Useful options:

```bash
bash src/standardise.sh --manifest manifests/main.tsv
bash src/standardise.sh --manifest manifests/locke_2015_25673413.tsv
bash src/standardise.sh --study-yaml manifests/ferkingstad_2021_34857953.yaml
bash src/standardise.sh --profile slurm
bash src/standardise.sh --row-limit 0
bash src/standardise.sh --only-study locke_2015_25673413
bash src/standardise.sh --force
bash src/standardise.sh --qcplot false
```

`--row-limit 1` is the default for `--profile local`; `--row-limit 0` means no limit and is the default for `--profile slurm`.

The required manifest columns are:

```text
GWAS_location	phenotype	ancestry	author	year	PMID
```

Optional manifest columns can be added only when automatic detection is not enough:

```text
sex	source_type	input_build	population	delimiter	format	study_yaml
```

`ancestry` and `population` have different purposes. `ancestry` describes the GWAS sample ancestry reported by the source or paper, and should normally be populated for every row, for example `EUR`, `African American`, `East Asian`, `multi-ancestry`, or `all`. `population` is only a technical override for the reference-panel lookup used during rsID/SNPID mapping, coordinate filling, and EAF filling. Leave `population` blank unless the ancestry label does not map cleanly to a 1000 Genomes group or you deliberately want a specific reference lookup such as `AFR`, `AMR`, `EAS`, `EUR`, `SAS`, or `ALL`.

The pipeline derives `STUDYID` as `<author>_<year>_<PMID>`, records completed rows in `GWAS/.standardise_state.tsv`, and skips rows whose GWAS-VCF, index, metadata JSON, QC JSON, and required QC plot already exist.

## Outputs

Each manifest row is staged and standardised under:

```text
GWAS/<STUDYID>/raw/<SOURCE_ID>/<ORIGINAL_FILE>
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.gwas.vcf.gz
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.gwas.vcf.gz.tbi
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.metadata.json
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.qc.json
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.tiff
```

`SOURCE_ID` starts with the source name and then a source-specific identifier, for example `GWAScatalog_GCST002783` or `openGWAS_ieu-a-2`.

`--qcplot true` is the default. It creates a compressed TIFF with raw QQ/manhattan plots on the top row and standardised QQ/manhattan plots on the bottom row.

## Source Handling

`GWAS_location` is used to infer the staging method:

- GWAS Catalog study URLs/accessions use the GWAS Catalog downloader and prioritise harmonised FTP summary statistics.
- OpenGWAS dataset URLs/IDs use `ieugwaspy.query.gwasinfo()` and `ieugwaspy.query.gwasinfo_files()`.
- Direct URLs and local files use generic staging, automatic delimiter detection, global column aliases, and reference lookup for missing coordinates where possible.
- Study-specific YAML files in `manifests/` provide overrides for joins, exclusions, allele rules, unusual column names, and P-value derivation. Locke is intentionally handled without a YAML config.

## Conversion

To convert to a normal dataframe:

```r
vcf <- gwasvcf::readVcf(file)
# 1. Convert the VCF object to a tibble
raw_stats <- gwasvcf::vcf_to_tibble(vcf)
# 2. Extract and format into a standard data frame
# Note: P is calculated by reversing the -log10 transformation
df <- data.frame(
  CHR  = as.character(raw_stats$seqnames),
  POS  = as.numeric(raw_stats$start),
  SNP  = as.character(raw_stats$ID),
  REF  = as.character(raw_stats$REF),
  ALT  = as.character(raw_stats$ALT),
  EAF  = as.numeric(raw_stats$AF),
  BETA = as.numeric(raw_stats$ES),
  SE   = as.numeric(raw_stats$SE),
  P    = 10^-(as.numeric(raw_stats$LP)),
  N    = as.numeric(raw_stats$SS),
  stringsAsFactors = FALSE
)
```

## References

Reference data scripts live in `src/reference_data/` and write to `references/` by default. The standardisation workflow expects 1000 Genomes lookup tables under:

```text
references/1000genomes/phase3/lookup/
```

The lookup tables are used for rsID mapping and EAF filling when those values are missing from a GWAS.
