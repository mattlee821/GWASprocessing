# GWASprocessing

GWASprocessing turns GWAS summary statistics into a standardised GRCh38 GWAS-VCF (specification 1.0.0). It takes a manifest of studies from the GWAS Catalog, OpenGWAS, a URL, or a local file, detects the genome build, places every variant on GRCh38, aligns alleles to the reference genome, resolves strand, annotates dbSNP rsIDs, and writes 1 VCF per GWAS. Standardisation does not QC the data; variants are removed only when they cannot be placed, aligned, read, or are exact duplicates. A log file accounts for all actions. It runs under Nextflow on a Mac or a SLURM cluster, and 1 setup script installs its tools and builds its reference data.

## Installation

You need bash, curl, and git (on a Mac, install the Command Line Tools with `xcode-select --install`), and about 75 GB of free disk for the reference data. Step 4 installs everything else.

#### 1. Clone the repository

```bash
git clone https://github.com/mattlee821/GWASprocessing.git
cd GWASprocessing
```

#### 2. Choose local tools or a container

Both routes use the same repository, setup script, manifests and commands; they differ only in where the tools that process the data come from. Choose 1 and use its command in step 4.

| Route | Use it when | What step 4 installs |
|---|---|---|
| Local tools (default) | on a Mac, or on a cluster without Apptainer | Java, Nextflow, Python, bcftools with the liftover plugin, plink2, bgzip and tabix, in 1 private environment under `.tools/` (or `GWAS_TOOLS_ROOT`); no administrator rights needed |
| Container | on a cluster with Apptainer (`module load apptainer` first if your cluster uses modules), or on a machine with Docker | the published image `ghcr.io/mattlee821/gwasprocessing`, inside which every pipeline step runs, plus a small local environment with Nextflow and Java to launch them |

#### 3. Choose where the data go (optional)

The reference data (about 41 GB of downloads and 30 GB of lookups), the GWAS outputs and the temporary files live outside the repository. By default they go in folders inside the repository, or in `~/GWASprocessing-data` when the repository is in a cloud-synced folder. To put them somewhere else, start the setup command in step 4 with `GWAS_DATA_HOME=/path`, or edit the paths in `.env` after setup:

| Setting in `.env` | Holds |
|---|---|
| `GWAS_WORK_ROOT` | the GWAS outputs and `manifest.html` |
| `GWAS_REFERENCE_ROOT` | the built reference lookups, read at random positions on every run |
| `GWAS_TEMP_ROOT` (with `GWAS_LOG_ROOT`, `GWAS_NEXTFLOW_WORK`, `GWAS_NEXTFLOW_LOG`) | logs and Nextflow work files |
| `GWAS_REFERENCE_RAW_ROOT` | the raw reference downloads, read once while the references are built, so a slow or external drive will do |

Keep the first 3 on local or scratch storage, and avoid cloud-synced folders.

#### 4. Run setup

```bash
bash src/source.sh                        # local tools
bash src/source.sh --container apptainer  # container, with Apptainer
bash src/source.sh --container docker     # container, with Docker
```

Setup creates `.env` from [.env.example](.env.example), installs the tools of the chosen route (a few minutes) and writes their locations to `.env`, then downloads and builds the reference data (several hours). It is safe to rerun: every step skips work that is already done. To leave the container route later, empty `GWAS_CONTAINER` in `.env`. Parts of setup can be run alone:

```bash
bash src/source.sh --env-only        # only create or complete .env
bash src/source.sh --tools-only      # only install or update the tools
bash src/source.sh --reference-only  # only download and build the reference data
bash src/source.sh --force           # rebuild the reference lookups (and pull the image again)
bash src/source.sh --r-packages      # gwasvcf for the optional R reader test (needs R installed)
```

#### 5. Add your OpenGWAS token and review the settings

Add your OpenGWAS token to `.env`; it is needed only for OpenGWAS sources (generate one at [api.opengwas.io/profile](https://api.opengwas.io/profile/)):

```bash
OPENGWAS_JWT=your-token
```

`GWAS_WORKERS` sets how many pieces of a GWAS are standardised at once on 1 machine, about 2 GB of memory each – essentially how many parallel workers you want to use per GWAS. Setup fills in everything else, including the tool locations (`GWAS_ENV`, `NEXTFLOW`, `PLINK2`, `GWAS_BCFTOOLS`, `BCFTOOLS_PLUGINS`, and `GWAS_CONTAINER` and `GWAS_CONTAINER_IMAGE` for the container route); a value you set yourself is kept. `.env` is not tracked by git.

#### 6. Run the tests

About 4 minutes; they must all pass. There is 1 test for each of the 12 [pipeline steps](#pipeline-steps) (`pipeline/tests/test_step01_stage.py` to `test_step12_write_report.py`) and 1 each for the manifest, the run, the manifest page and setup. Each test checks every rule this README states for its component, and a failure names the rule. On the container route the tests run inside the image, so they test the image's tools.

```bash
bash src/test.sh
```

To run 1 test: `bash src/test.sh test_step08_palindromes`.

#### 7. Run the benchmark

About 50 minutes the first time; it must end with `Result: PASS`. It standardises the 4 rows of `manifests/test.tsv` and compares them with the values stored in `src/benchmark/expected/` (see [GWAS comparison](#gwas-comparison)).

```bash
bash src/benchmark.sh
```

## Usage

#### 1. Write a manifest

A manifest is a tab-separated file with 1 row per GWAS. Copy [manifests/test.tsv](manifests/test.tsv), the benchmark, which fills in every column:

```text
GWAS_location	phenotype	ancestry	sex	author	year	PMID	source_type	input_build	population	delimiter	format	study_yaml
https://www.ebi.ac.uk/gwas/studies/GCST002783	BMI	ALL	combined	Locke	2015	25673413	gwas_catalog	hg19	EUR	tab	table	
https://opengwas.io/datasets/ebi-a-GCST002783	BMI	ALL	combined	Locke	2015	25673413	opengwas	hg19	EUR		vcf	
https://www.ebi.ac.uk/gwas/studies/GCST006464	endometrial cancer	EUR	female	OMara	2018	30093612	gwas_catalog	hg19	EUR	tab	table	manifests/omara_2018_30093612.yaml
https://opengwas.io/datasets/ebi-a-GCST006464	endometrial cancer	EUR	female	OMara	2018	30093612	opengwas	hg19	EUR		vcf	manifests/omara_2018_30093612.yaml
```

| Column | Required | Values | Meaning |
|---|---|---|---|
| `GWAS_location` | yes | a GWAS Catalog study URL or accession (`GCST002783`); an OpenGWAS dataset URL or ID (`ebi-a-GCST002783`); a URL (`https://`, `ftp://`); a local file path | where the summary statistics are |
| `phenotype` | yes | text | the trait, written to the VCF header |
| `ancestry` | yes | `AFR`, `AMR`, `EAS`, `EUR`, `SAS` or `ALL` (`European`, `East Asian` and similar labels are accepted; `ALL` for a multi-ancestry GWAS) | the ancestry of the GWAS sample, written to the VCF header; it also chooses the 1000 Genomes population unless `population` is set |
| `sex` | no | `combined`, `male` or `female` | the sex of the GWAS sample, written to the VCF header (`NA` when empty) |
| `author`, `year`, `PMID` | yes | first author's surname, publication year, PubMed ID | written to the VCF header, and together name the study folder (`locke_2015_25673413`) |
| `source_type` | no | `gwas_catalog`, `opengwas`, `direct_url`, `local` | how the file is fetched. It is worked out from `GWAS_location`; set it only when that goes wrong, for example a local file whose name looks like an OpenGWAS ID |
| `input_build` | no | `hg19` (GRCh37) or `hg38` (GRCh38); labels such as `GRCh37` are accepted, and any other value fails that row | the genome build of the file's positions. The build is detected from the data; this value is used only when detection fails, and a value that contradicts the data is reported and ignored |
| `population` | no | `AFR`, `AMR`, `EAS`, `EUR`, `SAS` or `ALL`; any other value fails that row | the 1000 Genomes population whose allele frequencies are compared with the study's, to resolve the strand of palindromic SNPs, fill missing frequencies and check the effect allele. It defaults to `ancestry`; set it when the file's frequency column comes from another population. Locke is an all-ancestry meta-analysis (`ancestry` `ALL`) whose frequencies come from HapMap CEU, so its `population` is `EUR` |
| `delimiter` | no | `tab`, `comma`, `space` (any run of spaces or tabs), or the character itself | the column separator of a table. It is detected from the header when empty; leave it empty for a VCF |
| `format` | no | `table` or `vcf` | whether the file is a table of summary statistics or a GWAS-VCF (OpenGWAS files are GWAS-VCFs). It is detected from the file when empty |
| `study_yaml` | no | a path, relative to the repository or to `manifests/` | a study YAML with settings for this GWAS (step 2); empty when the file needs none |
| `output_name` | no | letters, digits, `_` and `-`, starting with a letter or digit | the name of the output folder and of every output file (`standardise/<output_name>/<output_name>.gwas.vcf.gz`, `.qc.json`, `.report.html`); when empty the folder is named from the source and the files are `standardised.*`. Each row of a manifest needs its own name; a name with other characters fails that row, and 2 rows with the same name stop the run |

#### 2. Add a study YAML when a file needs help

A study YAML holds settings that a manifest row cannot: column names the pipeline does not recognise, case and control counts, joins with other files, rows to exclude, extra columns to keep. A join or exclude file too large to read for every piece (the deCODE annotation file has 27 million rows) is marked `index: true`: the first GWAS that needs it builds an indexed copy beside it (`<file>.index/`, needs pyarrow, which `src/source.sh` installs), and every GWAS after that reads only the rows it needs. `columns: [a, b]` on a join limits it to those columns of the joined file, which makes an indexed join faster still. The index is rebuilt when the file changes. A study that writes a marker instead of the other allele at multiallelic sites (deCODE writes `!`: the effect allele was tested against every other allele combined) is handled with `multiallelic: marker: "!"`. Where the variant ID (`SNPID`, written `chrom:pos:allele:allele`) names the effect allele and 1 other allele, that allele becomes the other allele and the record carries `INFO/MA`, which says that `ES` and `AF` are for the effect allele against all other alleles combined and not against `REF` alone; `bcftools view -e 'INFO/MA=1'` drops these records, and `-i 'INFO/MA=1'` selects them. A marked row whose ID gives no other allele is removed with the invalid alleles. Every key is optional, so a YAML holds only what its study needs. [manifests/omara_2018_30093612.yaml](manifests/omara_2018_30093612.yaml) is the example: it sets the case and control counts of the endometrial cancer GWAS, which the file does not carry, and lists every other key as a commented line with what it does.

#### 3. Run it

```bash
bash src/standardise.sh --manifest manifests/test.tsv --row-limit 0   # every row not yet complete
bash src/standardise.sh --manifest manifests/test.tsv                 # 1 row at a time (the local default)
bash src/standardise.sh --only-study locke_2015_25673413              # 1 study
bash src/standardise.sh --force                                       # rerun rows that are complete
bash src/standardise.sh --report false                                # no HTML report per GWAS
sbatch src/standardise.sh --manifest manifests/test.tsv               # SLURM: every row, 1 job each
```

#### 4. Check what has been processed

Open `<GWAS_WORK_ROOT>/manifest.html`. It lists every GWAS this install has processed or tried to, with its manifest fields, status, variants in and out, sample size, QC flags, warnings and links to its folder, VCF and report, and is updated at the end of every run. Each GWAS has its own folder, `<GWAS_WORK_ROOT>/<STUDYID>/standardise/<SOURCE_ID>/` (for example `locke_2015_25673413/standardise/GWAScatalog_GCST002783/`):

```text
standardised.gwas.vcf.gz       the GWAS-VCF, with its tabix index (.tbi)
standardised.qc.json           step log, QC metrics, flags and warnings
standardised.metadata.json     study metadata
standardised.report.html       self-contained report: steps, flags, QQ and Manhattan plots
```

When the manifest row has an `output_name`, the folder is `<GWAS_WORK_ROOT>/<STUDYID>/standardise/<output_name>/` and the files are `<output_name>.gwas.vcf.gz`, `<output_name>.qc.json`, `<output_name>.metadata.json` and `<output_name>.report.html`.

#### 5. Rerun when something changes

A completed row is skipped while its outputs exist and its settings are unchanged (manifest row, study YAML, reference files, pipeline output version); a change reruns it. A problem with 1 row fails only that row, and the run lists failed rows, with the reason, at the end: a manifest value that fails a check (ancestry, sex, build, population, format, source, a missing study YAML) sets the row aside before it starts, and an unusable file or a build that cannot be detected fails it during processing. Only a problem with the setup or the manifest as a whole stops the run: a missing or expired OpenGWAS token, missing tools or references, a missing manifest column, or 2 rows writing to the same folder. To make a report for 1 GWAS processed with `--report false`: `python3 pipeline/bin/report.py --dir <GWAS folder>`.

#### 6. Read a VCF in R

```r
vcf <- gwasvcf::readVcf(file)
stats <- gwasvcf::vcf_to_tibble(vcf)
df <- data.frame(CHR = as.character(stats$seqnames), POS = stats$start, SNP = stats$rsid, REF = stats$REF,
                 ALT = stats$ALT, EAF = stats$AF, AFS = stats$AFS, BETA = stats$ES, SE = stats$SE, P = 10^-stats$LP)
```

## Pipeline steps

| Step | What happens | Method |
|---|---|---|
| 1. Stage | Download or copy the file; GWAS Catalog raw file first, harmonised file as fallback; checksum or gzip check | ours |
| 2. Read | Detect the delimiter and map columns by alias, or read GWAS-VCF records into the same columns | ours; aliases extended from MungeSumstats |
| 3. Parse and derive | BETA from OR; SE from a confidence interval or from BETA and P; −log10 P (LP) from P, LOG10P or BETA/SE | GwasDataImport and MungeSumstats formulas |
| 4. Genome build | Detect GRCh37 or GRCh38 from the data; a declared build only as fallback | ours, following GwasDataImport |
| 5. GRCh38 coordinates | dbSNP positions for files with rsIDs only (including retired rsIDs); BCFtools/liftover for GRCh37 positions | dbSNP build 157; BCFtools/liftover |
| 6. Required fields | Remove rows without a GRCh38 position, 2 valid alleles, an effect, a positive SE or a P value | ours |
| 7. 1000 Genomes annotation | Frequencies of the matched population; AF filled where the study has none; frequency check | 1000 Genomes 30x GRCh38 |
| 8. Palindromic SNPs | Strand of each A/T and C/G SNP inferred from study against 1000 Genomes frequency | TwoSampleMR action 2 |
| 9. Reference alignment | REF must match the GRCh38 genome: keep, swap, complement or remove; indels left-aligned | gwas2vcf, GWAS Catalog harmoniser, bcftools norm |
| 10. Duplicates | Remove repeated position and allele pairs | ours |
| 11. rsIDs | dbSNP rsID matched on position and both alleles | dbSNP build 157 |
| 12. Write and report | GWAS-VCF, index, step log, QC flags, HTML report, manifest page | gwas2vcf conventions; OpenGWAS QC thresholds |


### 1. Stage

GWAS Catalog sources are downloaded as the raw author-submitted file, and the URL of the Catalog's harmonised file is recorded. When the raw file cannot be standardised (columns that cannot be mapped, a build that cannot be detected), the harmonised file is used instead and the fallback is recorded; `gwas_catalog_file: harmonised` in a study YAML chooses it from the start. The raw file is preferred because the harmonised file is a processed copy: it removes every variant without an rsID, orients palindromic SNPs by 1 file-wide strand assumption, converts −log10 P to P (which underflows) and keeps duplicates. OpenGWAS sources use the authenticated API, and their GWAS-VCF goes through the same steps as any table, so it is re-harmonised to GRCh38. URLs are downloaded with retries, and a local file is copied into the row's raw folder (and copied again when it changes), except a file that already sits in the study's raw folder (`<GWAS_WORK_ROOT>/<STUDYID>/raw/`), which is read where it is, so the study keeps 1 raw file and 1 standardised output. GWAS Catalog downloads are checked against the Catalog's md5 checksum, or by gzip integrity where there is none; URL and OpenGWAS downloads are checked by gzip integrity when compressed. A local file with a checksum file beside it (`X.txt.gz.md5sum`, `X.txt.md5sum` or the same with `.md5`) is checked against it before it is copied; the checksum may be of the file or, for a `.gz`, of the text inside it, which is what deCODE's md5sum files hold. A local file that does not match its checksum, or a local `.gz` that fails the gzip check, fails that row. A file that fails, including one left damaged by an earlier run, is downloaded once more, and a second failure stops the row.

### 2. Read

The delimiter comes from the manifest or YAML, or is detected from the header (tab, comma or whitespace); `.` reads as missing and leading `##` lines are skipped. Columns are matched against built-in aliases (for example `EA`, `effect_allele`, `ALLELE1` and `A1`; `NCASE`; `INFO` and `Rsq`), and a YAML `columns:` block names any the aliases miss. Allele-dependent columns are taken from 1 family only, the original columns or the GWAS Catalog `hm_*` columns, so alleles are never paired with an effect oriented to other alleles. Compound identifiers (`rs123:A:G`, `1:12345:A:G`) keep their text in `INFO/SNPID` and supply positions when a file has none. A multi-allelic GWAS-VCF record becomes 1 row per ALT allele.

### 3. Parse and derive

BETA is log OR for binary traits. SE comes from a confidence interval, (log upper − log lower) / (2 × 1.96) for odds ratios or (upper − lower) / (2 × 1.96) for betas, or else from abs(BETA / qnorm(P / 2)). LP (−log10 P) is read from a LOG10P column or computed from P, and recovered from BETA/SE when P is 0 or missing, so the most significant variants keep their precision.

### 4. Genome build

The build is detected from the data in every case. A fixed-seed sample of up to 20,000 rows from the first 200,000 is matched on chromosome, position and alleles (either order, either strand) against 2 million common 1000 Genomes SNPs placed in both GRCh37 and GRCh38. A build is called when it has at least 20 matches and at least 10 times the other build's matches. Only when detection fails is a declared build used: the manifest or YAML `input_build`, then the source metadata, then the VCF header. A declared build that disagrees with the detected one is reported and not used; with neither, the row stops. Files with rsIDs and no positions need no build.

### 5. GRCh38 coordinates

Files with rsIDs and no positions take GRCh38 positions straight from dbSNP build 157, so no liftover is involved; rsIDs at more than 1 position are not used. An rsID that dbSNP has retired by merging it into another takes the position of the current rsID, from a map of 21.6 million retired rsIDs: Locke's lead SNP rs12016871, for example, became rs9581854 in dbSNP build 131. `INFO/SNPID` keeps the study's rsID.

GRCh37 positions are lifted with BCFtools/liftover ([BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies](https://doi.org/10.1093/bioinformatics/btae038)), which loses fewer variants than position-only tools, lifts indels, and swaps alleles where the reference base changed between builds. Alleles are first oriented to the GRCh37 genome, then lifted with both genomes and the UCSC chain, and the effect allele is followed through swaps and through chain blocks on the minus strand, where the alleles are complemented. Rows that cannot be lifted are counted by reason: alleles matching neither strand of GRCh37, the tool's reject reason, a move to another chromosome, or a third allele at the site in GRCh38.

### 6. Required fields

Rows without a GRCh38 position, 2 different alleles made of A, C, G and T, an effect, a positive SE or a valid P value are removed, each counted against the first cause that applies. A missing P is first recovered from BETA and SE (step 3), so only a P outside 0 to 1, or a missing P with no SE, removes a row. Nothing else is filtered: frequency, imputation quality, chromosome and biallelic filters are analysis choices left to the user.

### 7. 1000 Genomes annotation

Each variant is matched on GRCh38 position and alleles (either strand) to the frequencies of the manifest population, from the 2,504 unrelated phase 3 samples of the 1000 Genomes 30x GRCh38 release. Where the study reports no frequency, `FORMAT/AF` is filled from 1000 Genomes and `FORMAT/AFS` names the source. The frequency check correlates the study frequencies of non-palindromic SNPs with 1000 Genomes: a negative correlation over at least 1,000 SNPs means the effect-allele column is mislabelled (as when a SAIGE file reports `Allele2` as the effect allele) and stops the row; a correlation below 0.5 over at least 100 SNPs is a warning. If the frequencies of the chosen population are missing from the references, those of ALL are used and a warning says so.

### 8. Palindromic SNPs

A/T and C/G SNPs follow the TwoSampleMR `harmonise_data` frequency approach (its default, action 2), with 1000 Genomes as the second frequency. When the study frequency and the 1000 Genomes frequency of the same allele are on the same side of 0.5, the SNP is kept; on opposite sides, the study is on the other strand and the alleles are complemented, with ES and EAF unchanged. When either frequency is within 0.08 of 0.5 (`palindromic_tolerance`), or either is missing (a filled frequency carries no strand information), the SNP is kept as reported and marked ambiguous. MR tools resolve palindromes again against their own data and can exclude these.

### 9. Reference alignment

`REF` is the GRCh38 reference base and `ALT` the other allele. When the effect allele is the reference base, the alleles are swapped, `ES` and `EZ` negated and `AF` becomes 1 − AF. A non-palindromic SNP reported on the reverse strand is complemented. A SNV matching neither strand is removed. Indels are checked against the reference, left-aligned and trimmed, so the same indel written 2 ways is recognised in step 10. Where the reference has no sequence at the position, or an `N`, the variant is kept as reported and marked `unverified`.

### 10. Duplicates

Records with the same GRCh38 position, `REF` and `ALT` are duplicates; as alignment has already put the alleles in `REF`/`ALT` order, the same variant reported either way round matches, while an insertion and a deletion at the same position (G>GCA, GCA>G) stay separate. Variants that could not be aligned are matched on the unordered allele pair. The first in input order is kept, except that a record under a current rsID is kept over one under a retired rsID.

### 11. rsIDs

The ID column is the dbSNP build 157 rsID matched on GRCh38 position and both alleles, or `.` when there is none; variants without an rsID are kept.

### 12. Write and report

Records are sorted, compressed with bgzip and indexed with tabix. The step log in `standardised.qc.json`, also printed to the run log, gives every step's rows in, rows removed by reason and changes made, and the counts reconcile from input rows to written rows. The QC flags follow the OpenGWAS report thresholds and never stop a row: frequency correlation with 1000 Genomes below 0.7 (a negative correlation included: OpenGWAS uses the absolute value, but after alignment a negative correlation means mislabelled alleles), lambda GC above 1.2, maximum sample size below 10,000, mean effect above 0.5 or 0.1 in absolute value, and more than 1,000 variants at P < 5e-8. The HTML report holds the step table, flags, header counts, raw and standardised QQ plots, raw and standardised Manhattan plots (raw on the left, standardised on the right, with the genome-wide significance line at P = 5e-8), warnings and timings in 1 self-contained file. The raw Manhattan panel uses the raw file's own chromosome and position, so it shows the file's build, and says so when the raw file has no chromosome and position (rsID-only files). `manifest.html` in the work root summarises every GWAS of the install, and keeps failed rows listed until they complete.

## Output columns

The VCF follows GWAS-VCF 1.0.0 as gwas2vcf writes it, so `gwasvcf`,` pygwasvcf`, `gwasglue2` and other tools built for OpenGWAS files read it.

| Column | Content |
|---|---|
| `#CHROM` | `1`–`22`, `X`, `Y`, `MT`, in that order, without a `chr` prefix |
| `POS` | GRCh38 position, 1-based; indels left-aligned |
| `ID` | dbSNP rsID matched on position and alleles, else `.` |
| `REF`, `ALT` | non-effect allele (the GRCh38 reference base) and effect allele; 1 row per allele |
| `QUAL`, `FILTER` | `.` and `PASS` |
| `INFO/SNPID` | the study's own variant identifier |
| `INFO/HM` | what was done to the study's alleles (below) |
| `FORMAT/ES`, `SE` | effect size relative to ALT (log OR for binary traits) and its standard error |
| `FORMAT/LP` | −log10 P |
| `FORMAT/AF` | ALT frequency: the study's (as given), else 1000 Genomes (6 significant digits) |
| `FORMAT/AFS` | source of AF: `study` or `1000G_<population>`, for example `1000G_EUR` |
| `FORMAT/SS`, `NC`, `EZ`, `SI` | sample size, cases, Z-score and imputation accuracy per variant, when the file has them |
| `FORMAT/ID` | the rsID, as in the ID column |
| Sample column | the output ID, for example `GWAScatalog_GCST002783` |

`INFO/HM` (harmonisation) records, for each variant, what the pipeline did to the alleles the study reported. Harmonisation can exchange the effect and other allele, which reverses the sign of the effect, or move the alleles to the other DNA strand; without a record, a variant whose effect sign differs from the paper's cannot be told apart from an error. HM makes every change traceable, and lets you exclude, for example, the palindromic SNPs whose strand could not be inferred. The name follows the GWAS Catalog harmoniser, whose files carry the same per-variant record as `hm_code`; OpenGWAS files have none. `FORMAT/AFS` does the same for the allele frequency.

| Value | Meaning |
|---|---|
| `forward` | as reported |
| `swapped` | the effect allele was the reference base: alleles exchanged, ES negated, AF to 1 − AF |
| `complemented` | alleles now on the other strand from the study's, because the study reported the reverse strand or liftover crossed a minus-strand block (2 complements cancel) |
| `complemented_swapped` | both |
| `unverified` | no reference sequence at the position, or an `N` in it |
| `palindrome_forward`, `palindrome_complemented` | A/T or C/G SNP whose strand was inferred from frequency |
| `palindrome_ambiguous` | A/T or C/G SNP kept as reported, strand not inferable |
| `palindrome_ambiguous_complemented` | as above, but liftover complemented the alleles |

Each value takes `_swapped` when the alleles were then swapped. The header declares the reference genome, `##contig` lines with `assembly=GRCh38`, the settings fingerprint, every INFO and FORMAT field, and a `##SAMPLE` line with a `##META` definition for each key: `TotalVariants`, `VariantsNotRead`, `VariantsNotHarmonised` and `HarmonisedVariants` (which add up), `SwitchedAlleles`, `AFFilled`, `StudyType`, `TotalCases` and `TotalControls` (a continuous trait puts its total sample size in `TotalControls`, as gwas2vcf does), `StudyID`, `Phenotype`, `Ancestry`, `Sex` and `PMID`.

Compared with an OpenGWAS file, the output is on GRCh38, fills `FORMAT/AF` from 1000 Genomes and marks it in `FORMAT/AFS`, adds `INFO/SNPID` and `INFO/HM`, carries the trait metadata in `##SAMPLE`, declares the 25 primary contigs only, and types `SS` as Float, as the specification says. R gwasvcf replaces a per-variant `SS` with `TotalCases + TotalControls` when the header has them.

## Reference data

| Under `GWAS_REFERENCE_ROOT` | Use | Source |
|---|---|---|
| `genome/GRCh38.fa`, `GRCh37.fa` | alignment; liftover source and target | Ensembl release 113 primary assemblies |
| `genome/hg19ToHg38.over.chain.gz`, `hg38ToHg19.over.chain.gz` | liftover of GRCh37 inputs; placing the build-detection sample in GRCh37 | UCSC |
| `dbsnp/GRCh38/` | positions for rsIDs, rsIDs for positions and alleles, 21.6 million retired rsIDs mapped to current ones | NCBI dbSNP build 157 and its merge history |
| `1000genomes/GRCh38/` | ALT frequencies for AFR, AMR, EAS, EUR, SAS and ALL | 1000 Genomes 30x on GRCh38, 2,504 unrelated phase 3 samples |
| `build_detection/` | 2 million common SNPs placed in both builds | 1000 Genomes GRCh38 EUR, fixed seed, lifted with BCFtools/liftover |

## Comparison with other pipelines and benchmark

### Step comparison

We read the code of MungeSumstats (Al-Murphy/MungeSumstats, November 2026), gwas2vcf and GwasDataImport (MRCIEU) and the GWAS Catalog sumstats-harmoniser (v1.1.12).

| Step | OpenGWAS (GwasDataImport, gwas2vcf) | MungeSumstats | GWAS Catalog harmoniser | GWASprocessing |
|---|---|---|---|---|
| Input | column positions given by the uploader | tables, VCF and OpenGWAS IDs; 410-entry header dictionary | GWAS-SSF columns only | tables and GWAS-VCF; alias dictionary and study YAML |
| Build | HapMap3 rsID match (GwasDataImport); taken as declared (gwas2vcf) | 10,000 rsIDs against dbSNP 155, unseeded | from metadata or file name | paired GRCh37/GRCh38 sample, fixed seed; declared build as fallback |
| Output build | GRCh37 | GRCh37 or GRCh38 | GRCh38 | GRCh38 |
| Coordinates | UCSC chain before upload | Ensembl chain, alleles not rechecked | Ensembl positions by rsID, then position-only liftover (X, Y and MT fail) | dbSNP positions for rsID-only files; BCFtools/liftover following the alleles |
| Strand | swap only; reverse-strand SNPs lost | swap only | both strands; palindromes by 1 file-wide strand assumption | swap and complement; palindromes per variant by frequency |
| rsIDs | dbSNP 153 by position | dbSNP 155 by position; variants without an rsID removed | Ensembl (dbSNP 151) by alleles; study rsIDs overwritten | dbSNP 157 by position and alleles; retired rsIDs mapped; variants without an rsID kept |
| Filters | P ≤ 0 or above 1, non-ACGT alleles, REF mismatch | INFO below 0.9, rsIDs with 3 or more alleles, X, Y and MT, missing SE or N | every row not harmonised | only rows that cannot be placed or read, and duplicates |
| P = 0 | row removed (gwas2vcf) | kept as 0 | −log10 P converted to P, which underflows | LP recovered from BETA/SE |
| Duplicates | kept | by rsID, then by position | kept | by position and allele pair |
| Output | GWAS-VCF | table (GWAS-VCF optional) | GWAS-SSF table | GWAS-VCF |
| Accounting | report and flags for a human reviewer | console messages | counts per harmonisation code | step log reconciling every row, flags, report, manifest page |

### GWAS comparison

- 2 public GWAS with all 5 pipelines
- compared the outputs with the raw study files
- ran the same Mendelian randomisation (MR) on each pipeline's output
- both GWAS have a raw file, an OpenGWAS GWAS-VCF, and a GWAS Catalog harmonised file
  - Locke et al. 2015 BMI ([Genetic studies of body mass index yield new insights for obesity biology]) 
    - NOTE: rsIDs without positions
  - O'Mara et al. 2018 endometrial cancer ([Identification of nine new susceptibility loci for endometrial cancer]) 
    - NOTE: gives GRCh37 positions, 1,014,582 indels, 12,906 cases and 108,979 controls. 
- `gwas2vcf` (1.4.3) was run on the Locke harmonised file against GRCh38 and on the endometrial cancer raw file against GRCh37
- `MungeSumstats` ran with its defaults and dbSNP 155
- The sampled GWASprocessing values are stored in `src/benchmark/expected/`, and `bash src/benchmark.sh` checks a new install against them.

Variants were matched by rsID: the lead SNPs plus 10,000 rsIDs drawn with seed 821. Locke's 8 leads are known BMI loci (FTO, MC4R, TMEM18, GNPDA2, SEC16B, TFAP2B, BDNF, NEGR1); the 16 endometrial cancer leads are the raw file's genome-wide significant SNPs clumped against 1000 Genomes EUR (r² 0.001, 10,000 kb). Each record was aligned to the raw effect allele and compared on alleles, effect sign, effect, and SE (within 1e-4), −log10 P (within 0.01) and frequency (within 0.02).

| Pipeline | Locke: present of 10,008 | Locke: sign agrees | Locke leads | Endometrial: present of 10,016 | Endometrial: sign agrees | Endometrial leads |
|---|---|---|---|---|---|---|
| GWASprocessing | 10,008 | 100% | 8 of 8 | 10,016 | 99.97% | 16 of 16 |
| OpenGWAS | 9,921 | 100% | 8 of 8 | 9,973 | 100% | 15 of 16 |
| gwas2vcf | 9,930 | 100% | 8 of 8 | 10,016 | 100% | 16 of 16 |
| GWAS Catalog harmoniser | 9,985 | 100% | 8 of 8 | 9,984 | 100% | 16 of 16 |
| MungeSumstats | 4,618 | 100% | 7 of 8 | 4,304 | 100% | 4 of 16 |

- The 3 GWASprocessing sign differences for endometrial cancer are A/T or C/G SNPs on minus-strand chain blocks, where liftover complements the alleles; for these SNPs a complement looks like a swap to a comparison by allele letter, and `INFO/HM` records it. 
- The OpenGWAS files have 27,303 (Locke) and 68,464 (endometrial cancer) fewer rows than the raw files, and the endometrial file lacks the top lead, rs11263761 at HNF1B (P = 3.2e-20). 
- `gwas2vcf` cannot read 24,229 Locke harmonised rows that have no values, and the harmonised Locke file has no effect for 55 sampled variants. 
- `MungeSumstats` removes every rsID that dbSNP 155 lists with 3 or more alleles (1,290,677 Locke rows, including the TMEM18 lead; 3,317,703 endometrial cancer rows, including 11 leads) and every variant with INFO below 0.9 (1,966,820 endometrial cancer rows, including 1 lead). 
- No pipeline disagreed with the raw study on the effect, SE, or P of a variant it kept

For the MR, Locke BMI was the exposure and endometrial cancer the outcome, each processed by the same pipeline, so differences between the estimates come only from processing. Instruments were exposure variants at P < 5e-8 clumped against 1000 Genomes EUR (r² 0.001, 10,000 kb) and harmonised with `TwoSampleMR` (0.7.13) under action 2 (palindromes resolved by frequency) and action 3 (palindromes removed).

| Pipeline | Variants at P < 5e-8 | Instruments | IVW, action 2: b (se), SNPs | IVW, action 3: b (se), SNPs | MR-Egger, action 2 |
|---|---|---|---|---|---|
| GWASprocessing | 2,042 | 78 | 0.612 (0.078), 78 | 0.615 (0.079), 75 | 0.579 |
| OpenGWAS | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| gwas2vcf | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| GWAS Catalog harmoniser | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| MungeSumstats | 915 | 60 | 0.634 (0.085), 59 | 0.639 (0.088), 56 | 0.383 |

- `GWASprocessing`, OpenGWAS, `gwas2vcf`, and the GWAS Catalog select the same 78 instruments, 3 of them palindromic. Only the GWASprocessing exposure carries the study's allele frequency, so only there can `TwoSampleMR` resolve those 3 under action 2; with palindromes removed, the 4 pipelines agree exactly. 
- `MungeSumstats` keeps 45% of the genome-wide significant exposure variants and shares 38 of its 60 instruments with the others.

The scripts used for the comparison are in `src/benchmark/`: `compare.py` (variant comparison), `to_mr_table.py` and `mr.R` (MR), `run_mungesumstats.R` and the gwas2vcf settings (`gwas2vcf_*.json`).
