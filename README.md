# GWASprocessing

GWASprocessing is a manifest-driven Nextflow workflow that stages GWAS summary statistics and standardises each one to a GRCh38 GWAS-VCF (specification 1.0.0, as the MRC IEU OpenGWAS tools write and read it). It installs with its tools and reference data on an HPC cluster or a Mac. It brings together in 1 install the steps OpenGWAS spreads over several tools (input checks, build detection and liftover, harmonisation against the reference genome, dbSNP annotation, QC flags), with strand handling from the GWAS Catalog harmoniser and column names from MungeSumstats. Later runs reuse a completed row until its manifest fields, study YAML, reference data or the output contract change.

Principles: GRCh38 only; the output reads like an OpenGWAS file (deliberate differences are listed under [Output](#output)); standardisation does not filter for analysis (variants are removed only when they cannot be placed, aligned or read, or are exact duplicates); every removed variant is counted once against 1 reason; decisions are taken from the data where possible, with metadata as a fallback.

## Layout

- `src/` contains the user-facing scripts: `src/source.sh` installs tools and builds the reference data, `src/standardise.sh` runs the workflow (`bash src/standardise.sh`, or `sbatch src/standardise.sh` from the repository root), and `src/reference_data/` holds the individual reference-data steps.
- `pipeline/` contains the Nextflow workflow (`main.nf`), its profiles (`conf/`), the Python steps (`bin/`) and the test suite (`tests/`).
- `manifests/` contains persistent input manifests and optional study YAML overrides.
- `.env` lists every setting with `/WHERE/...` placeholder paths. It is tracked and stays fixed.
- `server.env` holds your real settings. It is not tracked.

The data are large and stay outside the repository, in the folders the settings point to. Below, `references/`, `GWAS/` and `temp/` stand for those folders:

- `references/` (`GWAS_REFERENCE_ROOT`) receives the built reference data; raw reference downloads go to `GWAS_REFERENCE_RAW_ROOT` (default `references/raw`).
- `GWAS/<STUDY>/raw/<SOURCE_ID>/` (`GWAS_WORK_ROOT`) receives staged raw GWAS files.
- `GWAS/<STUDY>/standardise/<SOURCE_ID>/` receives the standardised GWAS-VCF, its tabix index, metadata and QC JSON, and the HTML report.
- `temp/` (`GWAS_TEMP_ROOT`) contains Nextflow runtime metadata, `.nextflow/`, `.nextflow.log*`, work directories and run logs.

## Quick Start

### 1: settings, tools and `references/`

Create `server.env` in the repository root and set the paths in it, using the names in [.env](.env):

```bash
GWAS_WORK_ROOT=/WHERE/THE/GWAS/FILES/GO
GWAS_REFERENCE_ROOT=/WHERE/THE/REFERENCE/DATA/IS
GWAS_REFERENCE_RAW_ROOT=/WHERE/THE/RAW/REFERENCE/DOWNLOADS/GO
GWAS_TEMP_ROOT=/WHERE/THE/TEMP/FILES/GO
GWAS_LOG_ROOT=/WHERE/THE/TEMP/FILES/GO/logs
GWAS_NEXTFLOW_WORK=/WHERE/THE/TEMP/FILES/GO/nextflow_work
GWAS_NEXTFLOW_LOG=/WHERE/THE/TEMP/FILES/GO/.nextflow.log
GWAS_WORKERS=3
OPENGWAS_JWT=your-token
```

The scripts read `server.env` first and then `.env`, so `server.env` needs only the lines you change, and a run stops with a message while a path is still a placeholder. Variables already set in the environment, and command-line options, take precedence over both files. Do not edit `.env`: it stays fixed so that `git pull` always works, and new settings get a default in the scripts instead.

Keep `GWAS_REFERENCE_ROOT`, `GWAS_WORK_ROOT`, `GWAS_TEMP_ROOT` and `GWAS_NEXTFLOW_WORK` on local or scratch storage: the built references are read at random positions on every run, which is about 100 times slower from a USB drive, and in a cloud-synced folder (OneDrive, Dropbox, iCloud) every output is uploaded; `src/standardise.sh` warns when it detects this. `GWAS_REFERENCE_RAW_ROOT` can be on a slow or external drive: the raw downloads are read once, start to finish, while the references are built.

OpenGWAS sources require authentication. Set `OPENGWAS_JWT` in `server.env` to a valid token (generate or renew one at [https://api.opengwas.io/profile/](https://api.opengwas.io/profile/)). If the token is missing or rejected, the run stops with exit code 78. Manifests without OpenGWAS sources do not need a token.

Run `bash src/source.sh` to install the tools (plink2, Nextflow, bgzip/tabix, the Python packages, and bcftools built with the BCFtools/liftover plugin) and build the reference data end to end. This is the recommended way to set up the project on a new machine, including a cluster. It is idempotent: every step skips work whose output already exists, so it is safe to re-run, including after an interrupted run. Tool locations are written to `server.env` (`PLINK2`, `NEXTFLOW`, `GWAS_VENV`, `GWAS_BCFTOOLS`, `BCFTOOLS_PLUGINS`). If bcftools with the liftover plugin cannot be built, setup continues with a warning and GRCh37 inputs are lifted with pyliftover (positions only).

```bash
bash src/source.sh                  # install tools, then build all reference data
bash src/source.sh --tools-only     # only install/verify tools and Python packages
bash src/source.sh --reference-only # only build/download reference data (tools must already exist)
bash src/source.sh --r-packages     # only the R packages for the R audit test (gwasvcf)
bash src/source.sh --force          # rebuild lookup tables and the rsID coordinate index
```

The reference build downloads about 40 GB (dbSNP 30 GB, 1000 Genomes about 8 GB, 2 FASTAs 2 GB) and builds about 30 GB of lookups. The built references are:

| Path under `references/` | Contents and use | Source |
|---|---|---|
| `genome/GRCh38.fa` (+ `.fai`) | GRCh38 primary assembly: REF check, alignment, indel normalisation, liftover target | Ensembl release 113 |
| `genome/GRCh37.fa` (+ `.fai`) | GRCh37 primary assembly: source genome for BCFtools/liftover only | Ensembl GRCh37 release 113 |
| `genome/hg19ToHg38.over.chain.gz` | UCSC liftover chain | UCSC |
| `dbsnp/GRCh38/` | rsID-sorted positions (for inputs with rsIDs and no positions) and position-sorted rsIDs with allele hashes (for the ID column), memory-mapped | NCBI dbSNP build 157 (`GCF_000001405.40`) |
| `1000genomes/GRCh38/` | ALT-allele frequencies of SNVs per super-population (AFR, AMR, EAS, EUR, SAS, ALL), keyed on GRCh38 position, for the 2,504 unrelated phase 3 samples: palindromic strand, AF fill, frequency check | 1000 Genomes 30x GRCh38 (Byrska-Bishop et al. 2022), plink2 resources |
| `build_detection/` | The same common SNPs placed in GRCh37 and GRCh38, for build detection | 1000 Genomes phase 3 EUR lookup placed in GRCh38 by dbSNP |
| `1000genomes/phase3/` | GRCh37 lookup and rsID index used when the GRCh38 references above are not installed | 1000 Genomes phase 3 |

`python3 src/reference_data/validate_references.py --reference-root /WHERE/THE/REFERENCE/DATA/IS` reports anything missing or in an outdated format.

### 2: manifest

Use [manifests/test.tsv](manifests/test.tsv) as the starting point for how to structure your manifest file:

```bash
bash src/standardise.sh --manifest manifests/test.tsv
```

The required manifest columns are:

```text
GWAS_location	phenotype	ancestry	author	year	PMID
```

Optional manifest columns can be added only when automatic detection is not enough:

```text
sex	source_type	input_build	population	delimiter	format	study_yaml
```

- `ancestry` must resolve to a 1000 Genomes super-population code: `ALL`, `AFR`, `AMR`, `EAS`, `EUR` or `SAS`. Unambiguous labels such as `European`, `East Asian`, `African American`, `Hispanic`, `Mixed` or `multi-ancestry` are accepted; labels such as `Asian`, `American` or `unknown` are rejected so that the reference population is chosen explicitly.
- `population` overrides the reference population used for frequencies. Leave it blank unless the ancestry does not map cleanly.
- `input_build` (`hg19`/`GRCh37` or `hg38`/`GRCh38`) is a fallback: the build is detected from the data, and a declared build is used only when detection is not possible (see [Genome build](#genome-build)).
- `format` (`table` or `vcf`) overrides the automatic choice of reader.

Rows are validated before anything is downloaded. An unrecognised ancestry, an unsupported source, an OpenGWAS location without a dataset ID, a GWAS Catalog location without an accession, or 2 rows that would write to the same output directory stop the run with exit code 78.

The pipeline derives `STUDYID` as `<author>_<year>_<PMID>` and records completed rows in `GWAS/.standardise_state.tsv`.

### 3: run

On a local machine, the wrapper runs 1 eligible manifest row per invocation. Re-run the same command to process the next incomplete row. On SLURM, submit the wrapper from the repository root with `sbatch src/standardise.sh`: the `slurm` profile is selected automatically and all eligible rows fan out, 1 Nextflow job per row. The head job's `#SBATCH --time` (48 h) must stay above the per-row limit (12 h in `pipeline/conf/slurm.config`).

```bash
bash src/standardise.sh --manifest manifests/main.tsv
bash src/standardise.sh --study-yaml manifests/ferkingstad_2021_34857953.yaml
bash src/standardise.sh --profile slurm
bash src/standardise.sh --row-limit 0
bash src/standardise.sh --only-study locke_2015_25673413
bash src/standardise.sh --force
bash src/standardise.sh --report false
```

`--row-limit 1` is the default for `--profile local`; `--row-limit 0` means no limit and is the default for `--profile slurm`. `--report false` skips the per-GWAS HTML report, for batches with thousands of GWAS; a report can be made later for any 1 GWAS (see [Report](#report)).

Reruns. A completed row is skipped only while its outputs exist and were produced with the same settings: the manifest fields, the resolved study YAML, the reference files (size and modification time) and the output-contract version (`OUTPUT_SCHEMA_VERSION` in `pipeline/bin/gwas_pipeline_common.py`). Changing any of these regenerates the row without `--force`. The fingerprint is stored as `settings_hash` in the state file, the QC JSON and the VCF header.

Resources. Tables and GWAS-VCFs are read in pieces of 1 million rows (`GWAS_TABLE_CHUNK`), so memory does not grow with file size. `--workers N` (Nextflow passes the task's `cpus`; under the local profile it comes from `GWAS_WORKERS`) standardises N pieces at once. The references are memory-mapped, so each worker holds only the pages it reads. The output does not depend on the piece size or the number of workers: the genome build is decided once from the first 200,000 rows, and sorting, deduplication in input order, lambda GC and the frequency check run once over the merged pieces. Large intermediate files go to `GWAS_TMPDIR` (default: the system temporary directory, node-local scratch under SLURM).

Failures. Each row runs independently. A row whose input is unusable exits with code 65 and is recorded as failed; other errors are retried once. Configuration errors (exit code 78, for example a build that can be neither detected nor declared) stop the whole run. At the end `src/standardise.sh` prints a summary, lists the rows that failed, and exits with status 1 if any row failed.

## Sources

`GWAS_location` decides how the file is staged:

- GWAS Catalog study URLs and accessions: the raw author-submitted file is downloaded by default, and the URL of the Catalog's harmonised file is recorded. If the raw file cannot be standardised (a data or configuration error, for example unmapped columns or an undetectable build), the harmonised file is downloaded and used instead, and the fallback is recorded in the step log and warnings. Study YAML `gwas_catalog_file: harmonised` chooses the harmonised file from the start. Downloads are verified by MD5 checksum or gzip integrity, and `genome_assembly` is read from the `-meta.yaml` sidecar when there is one. The raw file is preferred because the harmonised one is a processed copy that adds no information we cannot derive, and its processing removes every variant without an rsID, orients all palindromic SNPs by 1 file-wide strand assumption, converts −log10 P to P (which underflows), and does not remove duplicates.
- OpenGWAS dataset URLs and IDs use the authenticated OpenGWAS API with the `OPENGWAS_JWT` token. The GWAS-VCF is read into the same steps as a table (see below), so it is re-harmonised to GRCh38 like any other input.
- Direct URLs are downloaded with retries; `s3://` and `gs://` locations use the `aws` and `gsutil` tools. Local files are copied into the raw directory and copied again when the source changes.
- Study YAML files in `manifests/` provide overrides: `columns`, `joins`, `exclusions`, `transforms`, `extras`, `input_build`, `palindromic_tolerance`, `study_type`, `ncase`, `ncontrol`, `sample_size`, `sample` (which sample of a multi-sample VCF), `gwas_catalog_file`.

## Steps

Every input, table or GWAS-VCF, goes through the same steps in this order. The method and its source:

| Step | Method | Source |
|---|---|---|
| 1. Read | Tables: delimiter detected, columns matched by alias. GWAS-VCFs: records read into the same columns, 1 row per ALT allele | Ours |
| 2. Parse and derive | BETA = log OR; SE from the CI, (log upper − log lower) / (2 × 1.96) for OR intervals or (upper − lower) / (2 × 1.96) for beta intervals; SE = abs(BETA / qnorm(P/2)); LP = −log10 P, from LOG10P directly, or from BETA/SE when P is 0 or missing | GwasDataImport and MungeSumstats formulas; P = 0 recovery is ours |
| 3. Genome build | Detected from the data; declared build only when detection fails | Ours; principle from GwasDataImport and MungeSumstats |
| 4. GRCh38 coordinates | rsID-only inputs: GRCh38 position from dbSNP, no liftover. GRCh37 inputs: BCFtools/liftover | dbSNP lookup as in MungeSumstats and the GWAS Catalog; Genovese et al. 2024 |
| 5. Required fields | Remove rows without a GRCh38 position, 2 different A/C/G/T alleles, ES, SE > 0 or LP; the first cause is recorded | Ours |
| 6. 1000 Genomes annotation | GRCh38 frequency of the matched population; AF filled where the study has none | Ours |
| 7. Palindromic SNPs | Strand inferred per variant from study vs 1000 Genomes frequency (tolerance 0.08); ambiguous SNPs kept and marked | TwoSampleMR action 2 |
| 8. GRCh38 alignment | REF must match the GRCh38 FASTA: keep, swap (negate ES, AF to 1 − AF), complement, or remove; indels checked and left-aligned | gwas2vcf swap, GWAS Catalog strand complementing, bcftools-style normalisation |
| 9. Duplicates | Same GRCh38 position and allele pair: first in input order kept | Ours |
| 10. rsIDs | dbSNP GRCh38 rsID matched on position and alleles | OpenGWAS dbSNP annotation; GWAS Catalog matching rule |
| 11. Write and QC | GWAS-VCF, bgzip, tabix; QC metrics and OpenGWAS QC flags; HTML report | gwas2vcf conventions; opengwas-reports thresholds |

## Column handling

Tables are read with the delimiter from the manifest or YAML, or else detected from the header: tab, comma, or runs of whitespace. Empty cells stay in their own column, `.` reads as missing, and leading `##` comment lines are skipped.

Columns are matched by name against built-in aliases (for example `EA`, `effect_allele`, `ALLELE1`, `A1`; `NCASE` for case counts; `INFO`, `Rsq` for imputation accuracy), and a YAML `columns:` block can name them explicitly. Allele-dependent columns (alleles, effect size, EAF, odds ratio and its confidence interval) are always taken from 1 family: the original columns, or the GWAS Catalog harmonised `hm_*` columns when the original family lacks alleles or an effect size. Mixing families could pair alleles with an effect size oriented to different alleles. Compound identifiers (`rs123:A:G`, `1:12345:A:G`) keep the original text in `INFO/SNPID`; when a table has no chromosome or position columns, `chr:pos[:a1:a2]` identifiers supply them. YAML `extras:` columns are written as `INFO` string fields.

Mislabelled effect alleles are caught by the frequency check: the study EAF of non-palindromic SNPs is correlated with the 1000 Genomes frequency of the same allele, and a negative correlation over at least 1,000 SNPs stops the row with exit code 65 (for example SAIGE reports `Allele2` as the effect allele). A weak correlation is reported in `warnings`.

## Genome build

The build is detected from the data in every case. A fixed-seed sample of up to 20,000 rows with positions, from the first 200,000 rows of the file, is matched on chromosome, position and alleles (either order, either strand) against the same common SNPs placed in GRCh37 and in GRCh38. A build is detected when it has at least 20 matches (as GwasDataImport requires) and at least 10 times the other build's matches. Only when detection fails is a declared build used: the manifest or YAML `input_build`, then source metadata (OpenGWAS `build`, GWAS Catalog `genome_assembly`), then a VCF header's `##contig assembly=` or `##reference`. A declared build that disagrees with the detected one is reported in `warnings` and not used. With neither, the row stops with exit code 78. Inputs with rsIDs and no positions need no build: their GRCh38 positions come from dbSNP.

## Coordinates and liftover

- rsIDs without positions: the GRCh38 chromosome and position come straight from dbSNP build 157, so no liftover is involved. rsIDs at more than 1 position in dbSNP are not used.
- GRCh37 positions are lifted with BCFtools/liftover (Genovese et al. 2024, [BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies](https://doi.org/10.1093/bioinformatics/btae038)), which drops fewer variants than position-only tools, handles indels, and swaps alleles where the reference allele changed between builds. Alleles are first oriented to the GRCh37 FASTA (REF must match the source genome), then lifted with both FASTAs and the UCSC chain; the effect allele is followed through allele swaps and minus-strand blocks. Variants that cannot be lifted are counted by the tool's reason, and so are variants whose alleles do not match GRCh37 and those where GRCh38 has a third allele. Without the plugin or the GRCh37 FASTA, pyliftover lifts positions only.

## Reference alignment

When `references/genome/GRCh38.fa` is present, `REF` is the GRCh38 reference sequence and `ALT` the other allele, with the effect size and frequency flipped to match:

- If the effect allele is the reference base, the alleles are swapped, `ES` and `EZ` are negated and `AF` becomes `1 − AF`.
- A non-palindromic SNP reported on the reverse strand is complemented.
- A SNV whose alleles match neither strand is removed. An indel whose alleles do not match the reference sequence is removed; verified indels are left-aligned and trimmed (as `bcftools norm`) before duplicates are removed, so the same indel written 2 ways is recognised as a duplicate.

## Palindromic SNPs

This follows the [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) `harmonise` allele-frequency approach (its default `action = 2`), with the second frequency from the matched 1000 Genomes population. For each A/T or C/G SNP the study-reported EAF is compared with the 1000 Genomes frequency of the same allele:

- Same side of 0.5: same strand, kept (`HM=palindrome_forward`).
- Opposite sides: the study is on the other strand, so the alleles are complemented (`HM=palindrome_complemented`), with ES and EAF unchanged.
- Ambiguous, kept as reported and marked `HM=palindrome_ambiguous`: either frequency within the tolerance of 0.5 (default `palindromic_tolerance` 0.08, the window 0.42–0.58), no allele-matched 1000 Genomes frequency, or no study-reported EAF (a frequency filled from 1000 Genomes carries no strand information). Downstream tools can exclude these, and MR tools re-resolve palindromes against their own data.

The step log reports the counts per decision.

## Output

Each manifest row is staged and standardised under:

```text
GWAS/<STUDYID>/raw/<SOURCE_ID>/<ORIGINAL_FILE>
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.gwas.vcf.gz
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.gwas.vcf.gz.tbi
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.metadata.json
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.qc.json
GWAS/<STUDYID>/standardise/<SOURCE_ID>/standardised.report.html
```

`SOURCE_ID` starts with the source name and then a source-specific identifier, for example `GWAScatalog_GCST002783` or `openGWAS_ieu-a-2`.

The VCF follows GWAS-VCF specification 1.0.0 as gwas2vcf writes it, so R gwasvcf, pygwasvcf, gwasglue2 and other tools built for OpenGWAS files read it.

| Element | Content |
|---|---|
| `#CHROM` | `1`–`22`, `X`, `Y`, `MT`, no `chr` prefix, in that order |
| `POS` | GRCh38, 1-based; indels left-aligned |
| `ID` | dbSNP rsID matched on position and alleles, else `.` |
| `REF` | non-effect allele, matching the GRCh38 FASTA |
| `ALT` | effect allele; 1 row per allele, no multi-allelic rows |
| `QUAL`, `FILTER` | `.`, `PASS` |
| `INFO/SNPID` | the study's own variant identifier |
| `INFO/HM` | harmonisation action: `forward`, `swapped`, `complemented`, `complemented_swapped`, `unverified` (no reference sequence at the position), or for palindromic SNPs `palindrome_forward`, `palindrome_complemented`, `palindrome_ambiguous`, each with `_swapped` when the alleles were then swapped |
| `FORMAT/ES`, `SE` | effect size relative to ALT (log OR for binary traits), its standard error |
| `FORMAT/LP` | −log10 P; recover `P = 10^(−LP)` |
| `FORMAT/AF` | ALT frequency: study-reported, else the 1000 Genomes frequency of the matched population |
| `FORMAT/AFS` | source of AF: `study`, or `1000G_<POP>` (for example `1000G_EUR`) |
| `FORMAT/SS`, `NC`, `EZ`, `SI` | per-variant sample size, cases, Z-score, imputation accuracy; written only when the data have them |
| `FORMAT/ID` | the rsID, as in the ID column (`.` without one) |
| Sample column | the output ID, for example `GWAScatalog_GCST009004` |

The header declares `##reference` (the FASTA), `##contig` lines with `assembly=GRCh38`, `##settings_hash`, every INFO and FORMAT field, and a `##SAMPLE` line with a `##META` definition for each key:

- `TotalVariants`, `VariantsNotRead`, `VariantsNotHarmonised`, `HarmonisedVariants`, `SwitchedAlleles`, which add up: `TotalVariants = VariantsNotRead + VariantsNotHarmonised + HarmonisedVariants`. `VariantsNotRead` counts rows removed for missing coordinates or required fields; `VariantsNotHarmonised` counts alignment removals and duplicates.
- `AFFilled`: variants whose AF came from 1000 Genomes.
- `StudyType` (`Continuous` or `CaseControl`), `TotalCases` and `TotalControls` as gwas2vcf writes them: case-control studies give both; continuous traits put the total sample size in `TotalControls`. Values come from the study YAML or the source metadata.
- `StudyID`, `Phenotype`, `Ancestry`, `Sex`, `PMID`.

Deliberate differences from an OpenGWAS file: GRCh38 coordinates; `FORMAT/AF` filled from 1000 Genomes with `FORMAT/AFS` marking the source (OpenGWAS leaves it missing); no 1000 Genomes frequencies in INFO; `INFO/SNPID` and `INFO/HM` added; trait metadata and `AFFilled` in `##SAMPLE`; 25 primary contigs only; `SS` typed Float as the specification says. Reader behaviour to know: R gwasvcf replaces per-variant `SS` with `TotalCases + TotalControls` (or `TotalControls` alone) when the header has them.

The QC JSON records what happened to the row: `steps` (the step log), `incomplete_by_cause`, `liftover`, `build_detection`, `palindromic_resolution`, `reference_alignment`, `duplicates_removed`, `rsid_dbsnp`, `af_concordance`, `lambda_gc`, `qc_metrics`, `qc_flags`, `sample_header`, `warnings`, `timings_s`, `settings_hash` and `provenance`.

## Report

The step log is always written to the QC JSON and printed to the task log. With `--report true` (the default) each GWAS also gets `standardised.report.html`: 1 self-contained file with no external files or links, holding the summary, OpenGWAS QC flags and metrics, the `##SAMPLE` counts, the step table with every removal reason, change and note, raw and standardised QQ plots, warnings and timings. Each output ID gets its own report, so a study with several GWAS (for example Pulit, 9) has 1 per GWAS. To report 1 GWAS after a run with `--report false`:

```bash
python3 pipeline/bin/report.py --dir GWAS/<STUDYID>/standardise/<SOURCE_ID>
```

QC flags follow the OpenGWAS report thresholds and never stop a row: absolute frequency correlation with 1000 Genomes below 0.7, lambda GC above 1.2, maximum sample size below 10,000, absolute mean effect above 0.5 or 0.1, more than 1,000 variants at P < 5e-8. The LDSC and clumping flags are not computed yet.

## Conversion

To convert to a data frame in R:

```r
vcf <- gwasvcf::readVcf(file)
raw_stats <- gwasvcf::vcf_to_tibble(vcf)
df <- data.frame(
  CHR  = as.character(raw_stats$seqnames),
  POS  = as.numeric(raw_stats$start),
  SNP  = as.character(raw_stats$rsid),
  REF  = as.character(raw_stats$REF),
  ALT  = as.character(raw_stats$ALT),
  EAF  = as.numeric(raw_stats$AF),
  AFS  = as.character(raw_stats$AFS),
  BETA = as.numeric(raw_stats$ES),
  SE   = as.numeric(raw_stats$SE),
  P    = 10^-(as.numeric(raw_stats$LP)),
  stringsAsFactors = FALSE
)
```

## Testing

Run the test suite from the repository root:

```bash
python3 -m unittest discover -s pipeline/tests -v
```

The suite reads its small synthetic example data from `pipeline/tests/examples/` and builds small GRCh38 reference fixtures itself. It needs the packages in `requirements.txt` and bgzip/tabix on `PATH`. The writer tests also use `bcftools`, the BCFtools/liftover test needs bcftools with the plugin (`GWAS_BCFTOOLS`, `BCFTOOLS_PLUGINS`), the integration test needs `nextflow` and Java, and the R audit test needs R with `gwasvcf` (`bash src/source.sh --r-packages`); tests whose tools are missing are skipped. Tests never use the network: downloads run against local HTTP servers that mimic the GWAS Catalog FTP and OpenGWAS API, every non-local request is routed to a dead proxy, and shell scripts run with stubbed `wget`, `curl` and package managers.

| Test file | Component |
|---|---|
| `test_common.py` | delimiter detection, table reading, identifiers, ancestry and build labels, atomic state and JSON writes |
| `test_stats.py` | −log10 P from P, LOG10P and z-scores; genomic inflation |
| `test_coordinates.py` | pyliftover outcomes, build detection, variant-ID parsing, rsID coordinate index |
| `test_alleles.py` | annotation, EAF fill, palindromic resolution, de-duplication, frequency concordance, FASTA access, alignment and indel normalisation |
| `test_vcf.py` | GWAS-VCF 1.0.0 header, `##SAMPLE`/`##META`, body formatting, chromosome order, bgzip/tabix validated with `bcftools` |
| `test_standardise_table.py`, `test_standardise_passthrough.py` | end-to-end standardisation of tables and GWAS-VCFs |
| `test_standardise_pieces.py` | pieces and workers give the same output; the step log and `##SAMPLE` counts add up |
| `test_grch38_references.py` | dbSNP positions and rsIDs, 1000 Genomes GRCh38 frequencies, build detection, BCFtools/liftover |
| `test_prepare_manifest.py`, `test_stage_source.py`, `test_update_state.py` | manifest validation and skipping, staging routes, state updates |
| `test_downloads.py` | `download_url`, GWAS Catalog (raw first, harmonised fallback) and OpenGWAS clients |
| `test_report.py` | −log10 P extraction, QQ thinning, the HTML report and its regeneration |
| `test_reference_data.py` | lookup and coordinate-index builders, reference validation |
| `test_shell_scripts.py` | `src/standardise.sh`, `src/source.sh` and the reference-data shell scripts |
| `test_audit.py` | Python and R (`gwasvcf`) readers on pipeline output |
| `test_nextflow_integration.py` | the real workflow via `src/standardise.sh`: failures, skipping and settings-driven reruns |
