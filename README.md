# GWASprocessing

GWASprocessing is a manifest-driven Nextflow workflow that stages GWAS summary statistics, standardises each one to an analysis-ready GWAS-VCF, and keeps persistent outputs. Later runs reuse a completed row until its manifest fields, study YAML, analysis settings, reference data or the output contract change.

## Layout

- `src/` contains the user-facing scripts: `src/source.sh` installs tools and builds the reference data, `src/standardise.sh` runs the workflow (`bash src/standardise.sh`, or `sbatch src/standardise.sh` from the repository root), and `src/reference_data/` holds the individual reference-data steps.
- `pipeline/` contains the Nextflow workflow (`main.nf`), its profiles (`conf/`), the Python steps (`bin/`) and the test suite (`tests/`).
- `manifests/` contains persistent input manifests and optional study YAML overrides.
- `references/` is the default root for reference data written by `src/reference_data/`.
- `GWAS/<STUDY>/raw/<SOURCE_ID>/` receives staged raw GWAS files.
- `GWAS/<STUDY>/standardise/<SOURCE_ID>/` receives the standardised GWAS-VCF, its tabix index, metadata and QC JSON, and the QC plot.
- `temp/` contains Nextflow runtime metadata, `.nextflow/`, `.nextflow.log*`, work directories and run logs.
- `deprecated/` contains legacy one-off scripts and study directories. It is intentionally not tracked.

## Quick Start

### 1: `.env`, tools and `references/`

Copy [.env.example](.env.example) to `.env` and set paths and tokens. Values in `.env` are defaults: variables already set in the environment, and command-line options, take precedence.

Keep `GWAS_WORK_ROOT`, `GWAS_REFERENCE_ROOT`, `GWAS_TEMP_ROOT` and `GWAS_NEXTFLOW_WORK` on local or scratch storage. In a cloud-synced folder (OneDrive, Dropbox, iCloud) every output is uploaded and cloud-only reference files are downloaded again on each run; `src/standardise.sh` prints a warning when it detects this.

OpenGWAS sources require authentication. Set `OPENGWAS_JWT` in `.env` to a valid token (generate or renew one at [https://api.opengwas.io/profile/](https://api.opengwas.io/profile/)). If the token is missing or rejected (expired or invalid), the run stops with exit code 78 and a message asking for a valid token, rather than silently skipping the row. Manifests without OpenGWAS sources do not need a token.

Run `bash src/source.sh` to install the tools (plink2, Nextflow, bgzip/tabix and the Python packages) and build the reference data end to end. This is the recommended way to set up the project on a new machine, including the server it will eventually run on. It is idempotent: every step skips work whose output already exists, so it is safe to re-run at any time, including after an interrupted run. The locations of tools it installs are written to `.env` (`PLINK2`, `NEXTFLOW`, `GWAS_VENV`) and picked up by `src/standardise.sh`. If the system Python refuses package installs (PEP 668), the packages go into a project virtual environment, `.venv`, which `src/standardise.sh` puts first on `PATH`.

```bash
bash src/source.sh                  # install tools, then build all reference data
bash src/source.sh --tools-only     # only install/verify tools and Python packages
bash src/source.sh --reference-only # only build/download reference data (tools must already exist)
bash src/source.sh --force          # rebuild lookup tables and the rsID coordinate index
```

The reference build downloads the 1000 Genomes phase 3 dataset (about 4 GB compressed, and about 50 GB of scratch space while `plink2` decompresses and filters it) and the GRCh38 reference FASTA (about 1 GB compressed, 3.1 GB uncompressed), so check free disk space first. Each step in `src/reference_data/` is idempotent and can be run by hand. The outputs are:

| Path under `references/` | Contents and use |
|---|---|
| `1000genomes/phase3/lookup/<POP>.tsv.gz` | Common biallelic SNPs on chromosomes 1–22 and X per super-population (and `ALL`): `CHR`, `SNP`, `POS19` (GRCh37), `A1` (ALT), `A2` (REF), `AF1` (frequency of `A1`). Used to detect the genome build, fill rsIDs and missing EAF, check effect-allele frequencies, and resolve palindromic SNPs. |
| `1000genomes/phase3/coordinates/` | rsID index of every 1000 Genomes variant with its chromosome and GRCh37 position, for inputs that report rsIDs without positions. |
| `genome/hg19ToHg38.over.chain.gz`, `genome/hg38ToHg19.over.chain.gz` | UCSC liftover chains for the dual-build coordinates. |
| `genome/GRCh38.fa` (+ `.fai`) | GRCh38 primary assembly, used to make `REF` the GRCh38 reference base. |

`python3 src/reference_data/validate_references.py` reports anything missing or in an outdated format. Lookup tables that already carry `A1`, `A2` and `AF1` keep working; rebuilding them with `--force` adds chromosome X. The 1000 Genomes download and processing is based on:
- [https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3](https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3)
- [https://dougspeed.com/reference-panel/](https://dougspeed.com/reference-panel/)

Note on genome build: 1000 Genomes phase 3 is only available on GRCh37 (there is no phase 3 remap to GRCh38; the only native-GRCh38 1000 Genomes release is the newer, differently QC'd 3,202-sample high-coverage callset, which would be a separate undertaking to adopt). This does not limit GWASprocessing: the lookup tables provide real GRCh37 positions, and the pipeline's liftover stage (see [Coordinate build](#coordinate-build)) derives the GRCh38 position from them, so every standardised output is dual-build regardless of the reference's native build.

### 2: manifest

Use [manifests/test.tsv](manifests/test.tsv) as the starting point for how to structure your manifest file. You can run the pipeline using the test manifest with:

```bash
bash src/standardise.sh --manifest manifests/test.tsv
```

The test manifest includes 3 GWAS. These are the same GWAS from 3 different sources. Their standardised outputs differ slightly because the data are filtered and processed slightly differently by OpenGWAS and the GWAS Catalog.

The required manifest columns are:

```text
GWAS_location	phenotype	ancestry	author	year	PMID
```

Optional manifest columns can be added only when automatic detection is not enough:

```text
sex	source_type	input_build	population	delimiter	format	study_yaml
```

- `ancestry` describes the GWAS sample ancestry reported by the source or paper and must be populated for every row with a 1000 Genomes super-population code: `ALL`, `AFR`, `AMR`, `EAS`, `EUR` or `SAS`. Unambiguous labels such as `European`, `East Asian`, `South Asian`, `African American`, `Hispanic`, `Mixed` or `multi-ancestry` are accepted. Labels such as `Asian`, `American`, `Middle East` or `unknown` are rejected so that the reference population is chosen explicitly.
- `population` is only a technical override for the reference-panel lookup used for rsID mapping, EAF filling and strand checks. Leave it blank unless the ancestry label does not map cleanly to a 1000 Genomes group or you deliberately want a specific reference lookup.
- `input_build` (`hg19`/`GRCh37` or `hg38`/`GRCh38`) is optional. When it is blank the build is detected from the data (see [Coordinate build](#coordinate-build)).
- `format` (`table` or `vcf`) overrides the automatic choice between the table and GWAS-VCF paths.

Rows are validated before anything is downloaded. An unrecognised ancestry, an unsupported source (Synapse is not supported yet), an OpenGWAS location without a dataset ID, a GWAS Catalog location without an accession, or 2 rows that would write to the same output directory stop the run with exit code 78 and a list of the offending manifest rows.

The pipeline derives `STUDYID` as `<author>_<year>_<PMID>` and records completed rows in `GWAS/.standardise_state.tsv`.

### 3: run

On a local machine, the wrapper runs 1 eligible manifest row per invocation. Re-run the same command to process the next incomplete row.

On SLURM, submit the wrapper from the repository root with `sbatch src/standardise.sh`. Inside an sbatch job the `slurm` profile is selected automatically (an explicit `--profile` still wins) and all eligible rows fan out, 1 Nextflow job per row. The head job's `#SBATCH --time` (48 h) must stay above the per-row limit (12 h in `pipeline/conf/slurm.config`); adjust both, and `params.slurm_queue` in `pipeline/nextflow.config`, for your cluster.

Useful options:

```bash
bash src/standardise.sh --manifest manifests/main.tsv
bash src/standardise.sh --manifest manifests/test.tsv
bash src/standardise.sh --study-yaml manifests/ferkingstad_2021_34857953.yaml
bash src/standardise.sh --profile slurm
bash src/standardise.sh --row-limit 0
bash src/standardise.sh --only-study locke_2015_25673413
bash src/standardise.sh --force
bash src/standardise.sh --qcplot false
bash src/standardise.sh --analysis-qc true
```

`--row-limit 1` is the default for `--profile local`; `--row-limit 0` means no limit and is the default for `--profile slurm`.

`--analysis-qc` (default `false`) applies standard GWAS analysis QC before writing (see [Analysis QC](#analysis-qc)).

Reruns. A completed row is skipped only while its outputs exist and were produced with the same settings: the manifest fields, the resolved study YAML, `--analysis-qc`, the reference files (size and modification time) and the output-contract version (`OUTPUT_SCHEMA_VERSION` in `pipeline/bin/gwas_pipeline_common.py`). Changing any of these regenerates the row on the next run without `--force`; `--force` regenerates rows regardless. The fingerprint is stored as `settings_hash` in the state file and the QC JSON.

Resources. A row needs about 1.6 GB for the reference lookup plus about 0.45 GB per million input variants (measured on a 2.56-million-variant GWAS Catalog file, which took about 2.5 minutes), so the 16 GB per-row limit in `pipeline/conf/slurm.config` covers tables of roughly 30 million variants. Pre-built GWAS-VCFs are processed in chunks of 1 million records (`GWAS_PASSTHROUGH_CHUNK`), so their memory does not grow with file size. Large intermediate files are written to `GWAS_TMPDIR` (default: the system temporary directory, which is node-local scratch under SLURM), never next to the outputs.

Failures. Each row runs independently. A row whose input is unusable exits with code 65 and is recorded as failed; other errors are retried once. Configuration errors (exit code 78, for example a missing chain file or an undeterminable genome build) stop the whole run. At the end `src/standardise.sh` prints a summary, lists the rows that failed, and exits with status 1 if any row failed.

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

`--qcplot true` is the default. It creates a compressed TIFF with 2 QQ plots side by side, raw (left) and standardised (right), each titled with its genomic inflation factor. Plots work in −log10 P space, so P-values beyond double precision are shown; every point with −log10 P ≥ 2 is drawn and the dense region below is thinned to plotting resolution. A plot is regenerated whenever the standardised VCF is newer.

The QC JSON records what happened to the row. The fields most worth checking are:

| Field | Meaning |
|---|---|
| `warnings` | Anything that needs attention, for example retaining fewer than 50% of the input rows. |
| `build_detection` | Declared and detected build, matches per build, and whether the build was declared or detected. |
| `liftover` | Per direction: variants lifted (`ok`), `unmapped`, `multiple`, `other_chrom`, and `minus_strand`. |
| `coordinate_source`, `reference_fills` | Where positions came from, and how many rows received `CHR`, `POS19`, an rsID or EAF from the reference. |
| `missing_required_by_field`, `rows_removed_missing_required` | Rows dropped because a required field was missing or invalid. |
| `derivations`, `p_value_source` | Rows whose BETA came from an odds ratio, whose SE came from a confidence interval or P, and whose −log10 P came from LOG10P, P or the z-score. |
| `af_concordance` | Correlation of study EAF with the reference frequency of the same allele. |
| `qc_removed`, `palindromic_resolution` | Analysis-QC drop counts (when enabled). |
| `reference_alignment` | Alleles swapped, strand-flipped or dropped by alignment to the GRCh38 FASTA. |
| `lambda_gc` | Genomic inflation factor of the output. |
| `settings_hash`, `provenance` | The rerun fingerprint, and the output-contract version, package versions and git commit. |

## Source Handling

`GWAS_location` is used to infer the staging method:

- GWAS Catalog study URLs and accessions use the GWAS Catalog downloader, which prioritises harmonised FTP summary statistics, verifies MD5 checksums or gzip integrity, and records the genome assembly from the harmonised `-meta.yaml` sidecar when the Catalog provides one.
- OpenGWAS dataset URLs and IDs use the authenticated OpenGWAS API (`POST /api/gwasinfo` and `POST /api/gwasinfo/files`) with the `OPENGWAS_JWT` token. Dataset IDs keep their case and suffix (for example `ieu-b-5129`, `ukb-d-30780_irnt`, `ebi-a-GCST90000123`, `finn-b-I9_HEARTFAIL`). The downloaded GWAS-VCF is passed through: its `FORMAT`/`INFO` fields are preserved, trait headers are re-injected, coordinates are lifted to GRCh38 (see [Coordinate build](#coordinate-build)), and `FORMAT/SS` is added from the OpenGWAS sample size to records that lack it (OpenGWAS headers often declare SS while their records omit it).
- Direct URLs are downloaded with retries; `s3://` and `gs://` locations use the `aws` and `gsutil` command-line tools. Local files are copied into the raw directory and copied again when the source file changes.
- Study-specific YAML files in `manifests/` provide overrides for joins, exclusions, column names, transforms, extra columns and analysis-QC settings. Locke is intentionally handled without a YAML config.

## Column handling

Tables are read with the delimiter from the manifest or YAML, or else detected from the header line: tab when it contains a tab, comma when it contains a comma, and runs of whitespace otherwise. Empty cells stay in their own column, `.` reads as missing, and leading `##` comment lines are skipped.

Columns are matched by name against built-in aliases (for example `EA`, `effect_allele`, `ALLELE1`, `A1`), and a YAML `columns:` block can name them explicitly. Allele-dependent columns (alleles, effect size, EAF, odds ratio and its confidence interval) are always taken from 1 family: the original columns, or the GWAS Catalog harmonised `hm_*` columns when the original family lacks alleles or an effect size. Mixing families could pair alleles with an effect size oriented to different alleles.

- An odds ratio becomes `BETA = log(OR)`, and its 95% confidence interval gives `SE = (log(upper) − log(lower)) / (2 × 1.96)`. A confidence interval for a beta gives `SE = (upper − lower) / (2 × 1.96)`. With `BETA` and `P` but no SE, `SE = |BETA| / z`, where z is the normal quantile of P/2.
- −log10 P is carried natively: a LOG10P column is used as it is (and negated if it holds log10 P), P converts directly, and a P that underflowed to 0 is recovered from the z-score `|BETA/SE|`. Strong cis-pQTLs therefore keep their −log10 P instead of becoming missing.
- Compound identifiers (`rs123:A:G`, `1:12345:A:G`) keep their rsID in `ID` and the original text in `INFO/SNPID`. When a table has no chromosome or position columns, `chr:pos[:a1:a2]` identifiers supply them.
- YAML `extras:` columns are written as `INFO` string fields. Recognised `transforms:` are `strip_hash_header`, `uppercase_alleles` and `derive_oa_from_ref_alt`; unknown YAML keys and transforms are reported in `warnings`.

Rows need a primary chromosome (1–22, X, Y, MT), a GRCh38 position, both alleles, `BETA`, `SE > 0` and a valid −log10 P. Rows without them are dropped and counted per field. A row with no usable variants at all exits with code 65 and names the missing fields.

## Coordinate build

Every output VCF is dual-build: `#CHROM`/`POS` are GRCh38 and `INFO/POS19` is GRCh37, and both are genuinely those builds. `coordinate_build` in `standardised.metadata.json` and `standardised.qc.json` is therefore always `GRCh38`.

- Build. Columns named for a build (`POS19`, `POS38`, `POS_hg19`, ...) are used as such. A generic position column needs the input build, which is either declared (manifest or YAML `input_build`, or source metadata such as the OpenGWAS `build` or the GWAS Catalog `genome_assembly`) or detected. Detection samples up to 20,000 rows and counts matches of chromosome, position and alleles against the GRCh37 lookup, directly (GRCh37) and after lifting GRCh38 to GRCh37 (GRCh38). A build is detected when its matches reach 5% of the sample and 10 times the other build's matches. A declared build contradicted by at least 10 matches stops the row with exit code 78, and so does a build that can be neither declared nor detected; there is no default build.
- Rows without positions get them from their rsID: the coordinate index covers every 1000 Genomes variant; without it, the population lookup covers only common SNPs and a warning is recorded.
- Liftover uses `pyliftover` with the UCSC chains. A variant is lifted only when it maps to exactly 1 place on the same chromosome. A variant without a GRCh38 position is dropped and counted; a missing GRCh37 position is written as `POS19=.`. When a segment maps to the minus strand, SNV alleles are complemented so they stay on the forward strand, and other alleles are dropped.
- Passthrough sources (OpenGWAS or pre-built GWAS-VCF) are lifted the same way. Chromosome names such as `chr1` or `23` are normalised. The source's `##contig` and `##reference` lines are replaced by GRCh38 declarations.
- The `hg19ToHg38` chain is required whenever a source has GRCh37 coordinates; the run stops with exit code 78 if it is missing.

## Reference alignment

When `references/genome/GRCh38.fa` is present, alleles are oriented so that `REF` is the GRCh38 reference base and `ALT` is the other allele, with the effect size and frequency flipped to match. This follows the GWAS-VCF convention and makes outputs from different sources directly comparable on `(CHROM, POS, REF, ALT)`.

- If the effect allele is the reference base, the alleles are swapped, `ES` is negated and `AF` becomes `1 − AF` (`EZ` and `INFO/AF` too in passthrough files).
- A non-palindromic SNP reported on the reverse strand is complemented.
- An SNV whose alleles cannot match the reference base is dropped. Indels are checked against the reference sequence and kept when they do not match. Alleles that are not A, C, G or T are kept unaligned.

Without the FASTA, `REF` is the study's other allele, `reference_alignment.performed` is `false`, and a warning is recorded.

## Effect allele frequency

`FORMAT/AF` is the study's effect allele frequency when the input reports one (`INFO/EAF_SOURCE=input`). Otherwise it is filled from the population lookup as the 1000 Genomes frequency of the effect allele (`EAF_SOURCE=reference`), matched on position and alleles in GRCh37 orientation; a non-palindromic SNP may match on the reverse strand. Reference-filled frequencies are never used to resolve strand.

As a check on allele labelling, the study EAF of non-palindromic SNPs is correlated with the reference frequency of the same allele. A negative correlation over at least 1,000 SNPs stops the row with exit code 65, because it almost always means the effect and other allele columns are swapped (for example SAIGE reports `Allele2` as the effect allele). A weak correlation is reported in `warnings`.

## Analysis QC

`--analysis-qc` (default `false`) enables a standard GWAS analysis-QC stage applied before the VCF is written. It is a general, opt-in capability; with it off, no rows are filtered for analysis reasons. When on, these filters run in order and their drop counts are recorded in `standardised.qc.json` under `qc_removed`:

1. Biallelic single-base SNPs only: keep rows where `EA` and `OA` are each 1 base in {A,C,G,T} and `EA != OA` (drops indels and multiallelic records).
2. Palindromic (strand-ambiguous) SNPs: `A/T`, `T/A`, `C/G`, `G/C`. Rather than dropping them all, the pipeline resolves strand from allele frequencies against the 1000 Genomes reference (see below) and keeps the resolvable ones; only the genuinely ambiguous are dropped.
3. MAF filter: `MAF = min(EAF, 1−EAF)`; require `EAF ∈ [0,1]` and `MAF >= 0.01` (drops missing or invalid EAF).
4. Positional de-duplication: key on `(CHR, POS38)` and the unordered allele pair, keeping the first record. The same variant coded `A/G` and `G/A` therefore counts as a duplicate; rsIDs are never used.

Each filter (and the MAF threshold) can be overridden per study via an `analysis_qc:` block in the study YAML, for example `analysis_qc: {maf_threshold: 0.005, remove_palindromic: false}`.

### Palindromic strand resolution

This follows the [TwoSampleMR](https://github.com/MRCIEU/TwoSampleMR/blob/master/R/harmonise.R) `harmonise` allele-frequency approach (equivalent to its default `action = 2`), except that the second frequency comes from the ancestry-matched 1000 Genomes reference instead of a second dataset. For each palindromic SNP the pipeline compares the study-reported effect allele frequency (`EAF`) with the reference frequency of that allele, read in GRCh37 orientation:

- Genuinely ambiguous, dropped: either frequency is within the tolerance of 0.5, i.e. `0.5 − tolerance < AF < 0.5 + tolerance`. With the default `palindromic_tolerance` of 0.08 (matching TwoSampleMR) this is the window `0.42 < AF < 0.58`. In that band the 2 strand hypotheses are indistinguishable.
- Otherwise strand is inferred from which side of 0.5 the 2 frequencies fall on: the same side means the same strand (kept unchanged); opposite sides mean the study is on the other strand, so the alleles are complemented (A↔T, C↔G) and kept, with effect size and `EAF` unchanged.
- Also dropped: SNPs with no allele-matched reference (`dropped_no_reference`) and SNPs whose EAF was not reported by the study (`dropped_no_study_eaf`); a frequency filled from the reference carries no strand information.

Per-decision counts are recorded in `standardised.qc.json` under `palindromic_resolution` (`kept_same_strand`, `flipped`, `dropped_ambiguous`, `dropped_no_reference`, `dropped_no_study_eaf`, `dropped_resolution_disabled`). Table-derived and passthrough (OpenGWAS) outputs are resolved the same way. Tune via the study YAML, for example `analysis_qc: {palindromic_tolerance: 0.10, resolve_palindromic_strand: false}`; `resolve_palindromic_strand: false` drops all palindromic SNPs (TwoSampleMR's `action = 3`).

Caveat: frequency-based strand inference approximates the GWAS population with the 1000 Genomes panel, so it is least reliable for ancestries poorly represented there and for SNPs with MAF near 0.5, hence the conservative ambiguity window.

## VCF header metadata

Every output declares its assembly with `##reference=GRCh38` and `##contig` lines carrying GRCh38 lengths. Trait-level metadata travels inside every VCF as simple `##key=value` header lines (the contract downstream consumers parse):

| Header | Meaning |
|---|---|
| `##study_id=` | `<author>_<year>_<PMID>` |
| `##phenotype=` | trait/phenotype (never blank; falls back to the study id) |
| `##ancestry=` | super-population code: `ALL`, `AFR`, `AMR`, `EAS`, `EUR`, `SAS` |
| `##sex=` | reported sex, or `NA` when unknown (never blank) |
| `##PMID=` | source PMID |

Header values have tabs and line breaks replaced by spaces, and the sample column of table-derived outputs is the phenotype with whitespace replaced by `_`. `ancestry` is mandatory: a manifest row whose `ancestry` does not resolve to a super-population code is rejected before the run starts (exit 78) rather than producing an empty `##ancestry=`. Table-derived and passthrough outputs carry the same lines.

## GWAS-VCF field contract

The single-sample VCF body uses these semantics:

| VCF element | Meaning |
|---|---|
| `#CHROM` | chromosome (`1`–`22`, `X`, `Y`, `MT`), in that order |
| `POS` | GRCh38 position (integer) |
| `ID` | rsID, or `.` |
| `REF` | other allele (OA); the GRCh38 reference base when the FASTA is installed |
| `ALT` | effect allele (EA); `ES` and `AF` refer to it |
| `INFO POS19=` | GRCh37 position (`.` if unavailable) |
| `INFO SNPID=` | source variant identifier (`.` if none) |
| `INFO EAF_SOURCE=` | `input` or `reference` origin of `FORMAT/AF` (`.` if no AF) |
| `FORMAT ES` | effect size (beta) |
| `FORMAT SE` | standard error |
| `FORMAT LP` | −log10(P); recover `P = 10^(−LP)` |
| `FORMAT AF` | effect-allele frequency (EAF) |
| `FORMAT SS` | sample size (N) |

With `--analysis-qc true`, `ALT` is never multiallelic. Passthrough outputs keep their source `FORMAT` fields (for example `EZ`, `NC`, `ID`) and source `INFO` fields alongside these, and gain `SS` from the OpenGWAS sample size where a record lacks it.

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

## Testing

Run the test suite from the repository root:

```bash
python3 -m unittest discover -s pipeline/tests -v
```

The suite needs the packages in `requirements.txt` and bgzip/tabix on `PATH`. The writer tests also use `bcftools`, the integration test needs `nextflow` and Java, and the R audit test needs R with the `gwasvcf` package; tests whose tools are missing are skipped. Tests never use the network: downloads run against local HTTP servers that mimic the GWAS Catalog FTP and OpenGWAS API, every non-local request is routed to a dead proxy, and shell scripts run with stubbed `wget`, `curl` and package managers.

| Test file | Component |
|---|---|
| `test_common.py` | delimiter detection, table reading, identifiers, ancestry and build labels, atomic state and JSON writes |
| `test_stats.py` | −log10 P from P, LOG10P and z-scores; genomic inflation (values checked against R) |
| `test_coordinates.py` | liftover outcomes, build detection, variant-ID parsing, rsID coordinate index |
| `test_alleles.py` | lookup loading, positional annotation, EAF fill, palindromic resolution, de-duplication, frequency concordance, FASTA access and alignment |
| `test_vcf.py` | header contract, body formatting, chromosome order, bgzip/tabix output validated with `bcftools` |
| `test_standardise_table.py`, `test_standardise_passthrough.py` | end-to-end standardisation of tables and pre-built GWAS-VCFs |
| `test_prepare_manifest.py`, `test_stage_source.py`, `test_update_state.py` | manifest validation and skipping, staging routes, state updates |
| `test_downloads.py` | `download_url`, GWAS Catalog and OpenGWAS clients |
| `test_qc_plot.py` | −log10 P extraction, QQ thinning, inflation factors, plot regeneration |
| `test_reference_data.py` | lookup and coordinate-index builders, reference validation |
| `test_shell_scripts.py` | `src/standardise.sh`, `src/source.sh` and the reference-data shell scripts |
| `test_audit.py` | Python and R (`gwasvcf`) readers on pipeline output |
| `test_nextflow_integration.py` | the real workflow via `src/standardise.sh`: failures, skipping and settings-driven reruns |
