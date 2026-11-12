# GWASprocessing

## 1. Description

GWASprocessing turns GWAS summary statistics into a standardised GRCh38 GWAS-VCF (specification 1.0.0, as the MRC IEU OpenGWAS tools write and read it). It takes a manifest of studies from the GWAS Catalog, OpenGWAS, a URL or a local file, detects the genome build, places every variant on GRCh38, aligns alleles to the reference genome, resolves strand, annotates dbSNP rsIDs and writes 1 VCF per GWAS with a step log that accounts for every removed row. It runs under Nextflow on a Mac or a SLURM cluster, and 1 setup script installs its tools and builds its reference data.

The pipeline brings together in 1 install the steps OpenGWAS spreads over GwasDataImport, gwas2vcf, dbSNP annotation and its QC reports, with per-variant strand handling following TwoSampleMR and the GWAS Catalog harmoniser. It follows 5 principles:

- GRCh38 only.
- The output reads like an OpenGWAS file; the few deliberate differences are listed in [Output](#output).
- Standardisation does not filter for analysis: variants are removed only when they cannot be placed, aligned or read, or are exact duplicates.
- Every removed variant is counted once, against 1 reason.
- Decisions are taken from the data where possible, with metadata as a fallback.

## 2. Installation

You need bash, curl, git, Python 3.9 or later, Java 17 or later (for Nextflow), and a C compiler with make and zlib (to build bcftools with the BCFtools/liftover plugin). bgzip and tabix are installed with conda or mamba when they are missing. R is needed only for the optional R reader test.

1. Clone the repository and move into it:

```bash
git clone https://github.com/mattlee821/GWASprocessing.git
cd GWASprocessing
```

2. Choose where the data go. The references (about 30 GB built, plus about 41 GB of raw downloads) and the GWAS outputs live outside the repository. By default they go in folders under the repository, or under `~/GWASprocessing-data` when the repository is in a cloud-synced folder. To put them elsewhere, set `GWAS_DATA_HOME` before step 3, or edit the paths in `.env` after it. Keep the built references, outputs and Nextflow work folder on local or scratch storage: the references are read at random positions on every run, which is about 14 times slower from a USB drive, and a cloud-synced folder uploads every file. The raw downloads (`GWAS_REFERENCE_RAW_ROOT`) can be on a slow or external drive, because they are read once.

3. Run the setup script. It creates `.env` from [.env.example](.env.example) and fills in the data paths, installs plink2, Nextflow, bgzip/tabix, a Python environment with the packages in `requirements.txt`, and bcftools built with the liftover plugin, writing each tool's location to `.env`. It then downloads and builds the reference data. Every step skips work that is already done, so it is safe to re-run after an interruption.

```bash
bash src/source.sh
```

Setup stops with an error if bcftools with the liftover plugin cannot be built; set `GWAS_BCFTOOLS` and `BCFTOOLS_PLUGINS` in `.env` to use an existing build instead. Parts of setup can be run alone:

```bash
bash src/source.sh --env-only       # only create or complete .env
bash src/source.sh --tools-only     # only install or check the tools and Python packages
bash src/source.sh --reference-only # only download and build the reference data
bash src/source.sh --force          # rebuild the built references
bash src/source.sh --r-packages     # R packages for the optional R reader test (gwasvcf)
```

4. Add your OpenGWAS token to `.env`; it is the only value setup cannot fill. Generate or renew one at [api.opengwas.io/profile](https://api.opengwas.io/profile/). Manifests without OpenGWAS sources do not need it.

```bash
OPENGWAS_JWT=your-token
```

Review `GWAS_WORKERS` in `.env` too: the number of 1-million-row pieces standardised at once under the local profile, each using about 2 GB of memory.

5. Check the reference data:

```bash
python3 src/reference_data/validate_references.py --reference-root <GWAS_REFERENCE_ROOT from .env>
```

6. Run the test suite (about 3 minutes):

```bash
python3 -m unittest discover -s pipeline/tests
```

7. Check the install against the benchmark. This standardises Locke et al. 2015 BMI and O'Mara et al. 2018 endometrial cancer and compares 10,008 and 10,016 sampled variants with the values stored in `src/benchmark/expected/`; each check must report 0 differences.

```bash
bash src/standardise.sh --manifest src/benchmark/benchmark.tsv --row-limit 0
python3 src/benchmark/check_install.py --expected src/benchmark/expected/locke_2015.tsv \
    --output <GWAS_WORK_ROOT>/locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.gwas.vcf.gz
python3 src/benchmark/check_install.py --expected src/benchmark/expected/omara_2018.tsv \
    --output <GWAS_WORK_ROOT>/omara_2018_30093612/standardise/GWAScatalog_GCST006464/standardised.gwas.vcf.gz
```

`.env` is not tracked; `.env.example` is the tracked template with `/WHERE/...` placeholders and is never read by the scripts. Variables set in the environment and command-line options take precedence over `.env`, and a run stops while a path is still a placeholder.

## 3. Usage

1. Write a manifest: a tab-separated file with 1 row per GWAS. [manifests/test.tsv](manifests/test.tsv) and [src/template_manifest.tsv](src/template_manifest.tsv) are examples.

```text
GWAS_location	phenotype	ancestry	author	year	PMID
https://www.ebi.ac.uk/gwas/studies/GCST002783	BMI	EUR	Locke	2015	25673413
https://opengwas.io/datasets/ebi-a-GCST006464	endometrial cancer	EUR	OMara	2018	30093612
```

`GWAS_location` is a GWAS Catalog study URL or accession, an OpenGWAS dataset URL or ID, a direct URL (`https://`, `s3://`, `gs://`) or a local path. `ancestry` must resolve to a 1000 Genomes super-population (`AFR`, `AMR`, `EAS`, `EUR`, `SAS` or `ALL`; labels such as `European` or `East Asian` are accepted, ambiguous ones such as `Asian` are rejected). Choose the population of the file's allele-frequency column, because palindromic SNPs are resolved against it: Locke's frequencies come from HapMap CEU, so Locke is run as `EUR`. Optional columns, needed only when automatic detection is not enough: `sex`, `source_type`, `input_build`, `population` (overrides the reference population), `delimiter`, `format` (`table` or `vcf`) and `study_yaml`.

2. Add a study YAML when the file needs help, and name it in the `study_yaml` column. It can name columns, add joins and exclusion files, set `input_build`, `study_type`, `ncase`, `ncontrol`, `sample_size`, `palindromic_tolerance`, `extras` (source columns written to INFO) or `gwas_catalog_file: harmonised`. [src/benchmark/endometrial_omara_2018.yaml](src/benchmark/endometrial_omara_2018.yaml) is a minimal case-control example; the YAMLs in `manifests/` are larger ones.

3. Run the manifest:

```bash
bash src/standardise.sh --manifest manifests/main.tsv                 # 1 incomplete row (local default)
bash src/standardise.sh --manifest manifests/main.tsv --row-limit 0   # every incomplete row
bash src/standardise.sh --only-study locke_2015_25673413              # 1 study
bash src/standardise.sh --force                                       # rerun rows that are complete
bash src/standardise.sh --report false                                # no HTML report per GWAS
sbatch src/standardise.sh --manifest manifests/main.tsv               # SLURM: every row, 1 job each
```

Rows are checked before anything is downloaded; an unrecognised ancestry or source, or 2 rows writing to the same folder, stop the run with exit code 78. Under SLURM the head job's `#SBATCH --time` (48 h) must stay above the per-row limit (12 h in `pipeline/conf/slurm.config`).

4. Find the outputs. Each row is written to `<GWAS_WORK_ROOT>/<STUDYID>/standardise/<SOURCE_ID>/`, where `STUDYID` is `<author>_<year>_<PMID>` and `SOURCE_ID` names the source, for example `GWAScatalog_GCST002783` or `openGWAS_ebi-a-GCST006464`. The raw download is kept in `<STUDYID>/raw/<SOURCE_ID>/`.

```text
standardised.gwas.vcf.gz       the GWAS-VCF
standardised.gwas.vcf.gz.tbi   tabix index
standardised.metadata.json     study metadata
standardised.qc.json           step log, QC metrics, flags and warnings
standardised.report.html       self-contained HTML report (with --report true)
```

The step log is also printed to the task log. To make a report later for 1 GWAS:

```bash
python3 pipeline/bin/report.py --dir <GWAS_WORK_ROOT>/<STUDYID>/standardise/<SOURCE_ID>
```

5. Rerun when something changes. A completed row is skipped while its outputs exist and its settings fingerprint matches: the manifest fields, the study YAML, the reference files and the output version. Changing any of them reruns the row without `--force`. A row whose input is unusable exits with code 65 and is recorded as failed while the other rows continue; at the end the wrapper lists failed rows and exits with status 1.

6. Read the VCF in R:

```r
vcf <- gwasvcf::readVcf(file)
stats <- gwasvcf::vcf_to_tibble(vcf)
df <- data.frame(CHR = as.character(stats$seqnames), POS = stats$start, SNP = stats$rsid, REF = stats$REF,
                 ALT = stats$ALT, EAF = stats$AF, AFS = stats$AFS, BETA = stats$ES, SE = stats$SE, P = 10^-stats$LP)
```

## 4. Pipeline steps

Every input, table or GWAS-VCF, goes through the same steps in this order:

1. Stage: download or copy the file; GWAS Catalog raw file first, harmonised file as fallback; checksum or gzip check.
2. Read: detect the delimiter and map columns by alias, or read GWAS-VCF records into the same columns.
3. Parse and derive: BETA from OR, SE from a confidence interval or from BETA and P, −log10 P (LP) from P, LOG10P or BETA/SE.
4. Detect the genome build from the data, with a declared build as fallback.
5. Place every variant on GRCh38: dbSNP positions for files with rsIDs only (including rsIDs retired by dbSNP merges), BCFtools/liftover for GRCh37 positions.
6. Remove rows that lack a GRCh38 position, 2 valid alleles, an effect, a positive SE or a P value, recording the first cause.
7. Annotate with 1000 Genomes GRCh38 frequencies of the matched population; fill AF where the study has none; check study against reference frequencies.
8. Resolve the strand of palindromic (A/T, C/G) SNPs per variant by frequency.
9. Align to the GRCh38 reference: REF must match the FASTA; swap, complement or remove; left-align indels.
10. Remove duplicates by GRCh38 position and allele pair.
11. Annotate dbSNP GRCh38 rsIDs matched on position and alleles.
12. Write the GWAS-VCF, index it, and write the step log, QC metrics, OpenGWAS QC flags and the HTML report.

## 5. Pipeline in detail and comparison with other pipelines

### Sources

- GWAS Catalog: the raw author-submitted file is downloaded by default, and the URL of the harmonised file is recorded. If the raw file cannot be standardised (unmapped columns, an undetectable build), the harmonised file is used instead and the fallback is recorded. `gwas_catalog_file: harmonised` in a study YAML chooses it from the start. The raw file is preferred because the harmonised file is a processed copy: its processing removes every variant without an rsID, orients palindromic SNPs by 1 file-wide strand assumption, converts −log10 P to P (which underflows) and does not remove duplicates.
- OpenGWAS: the authenticated API, with `OPENGWAS_JWT`. The GWAS-VCF is read into the same steps as a table, so it is re-harmonised to GRCh38 like any other input.
- Direct URLs are downloaded with retries; `s3://` and `gs://` use `aws` and `gsutil`; local files are copied and copied again when they change.

### Columns and statistics

The delimiter comes from the manifest or YAML, or is detected from the header (tab, comma or whitespace). Columns are matched against built-in aliases (for example `EA`, `effect_allele`, `ALLELE1`, `A1`; `NCASE`; `INFO`, `Rsq`), and a YAML `columns:` block can name them. Allele-dependent columns are taken from 1 family only, the original columns or the GWAS Catalog `hm_*` columns, so alleles are never paired with an effect oriented to other alleles. Compound identifiers (`rs123:A:G`, `1:12345:A:G`) keep the original text in `INFO/SNPID`, and supply positions when a table has none.

BETA is log OR for binary traits. SE comes from a confidence interval, (log upper − log lower) / (2 × 1.96) for odds ratios or (upper − lower) / (2 × 1.96) for betas, or from abs(BETA / qnorm(P / 2)). LP is read from LOG10P directly or computed as −log10 P, and recovered from BETA/SE when P is 0 or missing.

Mislabelled effect alleles are caught by the frequency check: the study EAF of non-palindromic SNPs is correlated with the 1000 Genomes frequency of the same allele, and a negative correlation over at least 1,000 SNPs stops the row with exit code 65 (as when a SAIGE file reports `Allele2` as the effect allele). A weak correlation is reported as a warning.

### Genome build

The build is detected from the data in every case. A fixed-seed sample of up to 20,000 rows from the first 200,000 is matched on chromosome, position and alleles (either order, either strand) against 2 million common 1000 Genomes SNPs placed in both GRCh37 and GRCh38. A build is called when it has at least 20 matches and at least 10 times the other build's matches. Only when detection fails is a declared build used: the manifest or YAML `input_build`, then source metadata, then the VCF header. A declared build that disagrees with the detected one is reported and not used; with neither, the row stops with exit code 78. Files with rsIDs and no positions need no build.

### Coordinates and liftover

Files with rsIDs and no positions take GRCh38 positions straight from dbSNP build 157, so no liftover is involved. rsIDs at more than 1 position are not used. An rsID that dbSNP has retired by merging it into another takes the position of the current rsID, from a map of 21.6 million retired rsIDs built from dbSNP's merge history: Locke's lead SNP rs12016871, for example, was merged into rs9581854 in dbSNP build 131. `INFO/SNPID` keeps the study's rsID and the ID column carries the current one. If the file also reports the variant under its current rsID, that record is kept.

GRCh37 positions are lifted with BCFtools/liftover ([BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies](https://doi.org/10.1093/bioinformatics/btae038)), which drops fewer variants than position-only tools, lifts indels, and swaps alleles where the reference allele changed between builds. Alleles are first oriented to the GRCh37 FASTA, then lifted with both FASTAs and the UCSC chain, and the effect allele is followed through swaps and minus-strand chain blocks. Rows that cannot be lifted are counted by reason: alleles matching neither strand of GRCh37, the tool's reject reason, another chromosome, or a third allele at the site in GRCh38.

### Required fields

Rows without a GRCh38 position, 2 different alleles made of A, C, G and T, an effect, a positive SE or a P value are removed, each counted against the first cause that applies. Nothing else is filtered: frequency, imputation quality, chromosome and biallelic filters are analysis choices left to the user.

### 1000 Genomes annotation and allele frequency

Each variant is matched by GRCh38 position and alleles (either strand) to the 1000 Genomes frequencies of the manifest population, from the 2,504 unrelated phase 3 samples in the 30x GRCh38 release. Where the study reports no frequency, `FORMAT/AF` is filled from 1000 Genomes (6 significant digits, the precision it is stored at) and `FORMAT/AFS` says so (`1000G_EUR`); study frequencies are written as given, with `AFS=study`.

### Palindromic SNPs

A/T and C/G SNPs follow the TwoSampleMR `harmonise_data` allele-frequency approach (its default, action 2), with the second frequency from 1000 Genomes. The study EAF is compared with the 1000 Genomes frequency of the same allele:

- Same side of 0.5: same strand, kept (`HM=palindrome_forward`).
- Opposite sides: the study is on the other strand, so the alleles are complemented (`HM=palindrome_complemented`) with ES and EAF unchanged.
- Either frequency within 0.08 of 0.5 (`palindromic_tolerance`), no matching 1000 Genomes frequency, or no study frequency (a filled frequency carries no strand information): kept as reported and marked `HM=palindrome_ambiguous`. MR tools re-resolve palindromes against their own data and can exclude these.

### Reference alignment

`REF` is the GRCh38 reference base and `ALT` the other allele. If the effect allele is the reference base, the alleles are swapped, `ES` and `EZ` negated and `AF` becomes 1 − AF. A non-palindromic SNP on the reverse strand is complemented. A SNV matching neither strand is removed; indels are checked against the reference, left-aligned and trimmed, so the same indel written 2 ways is recognised as a duplicate.

### Duplicates and rsIDs

Records with the same GRCh38 position and allele pair are duplicates; the first in input order is kept, except that a record under a current rsID is kept over one under a retired rsID. The ID column is the dbSNP build 157 rsID matched on GRCh38 position and both alleles, or `.` when there is none; variants without an rsID are kept.

### Output

The VCF follows GWAS-VCF 1.0.0 as gwas2vcf writes it, so R gwasvcf, pygwasvcf, gwasglue2 and other tools built for OpenGWAS files read it.

| Element | Content |
|---|---|
| `#CHROM` | `1`–`22`, `X`, `Y`, `MT`, no `chr` prefix, in that order |
| `POS` | GRCh38, 1-based; indels left-aligned |
| `ID` | dbSNP rsID matched on position and alleles, else `.` |
| `REF`, `ALT` | non-effect allele (the GRCh38 reference base) and effect allele; 1 row per allele |
| `QUAL`, `FILTER` | `.`, `PASS` |
| `INFO/SNPID` | the study's own variant identifier |
| `INFO/HM` | what was done to the study's alleles (below) |
| `FORMAT/ES`, `SE` | effect size relative to ALT (log OR for binary traits) and its standard error |
| `FORMAT/LP` | −log10 P |
| `FORMAT/AF`, `AFS` | ALT frequency and its source: `study`, or `1000G_<POP>` |
| `FORMAT/SS`, `NC`, `EZ`, `SI` | per-variant sample size, cases, Z-score and imputation accuracy, written only when the data have them |
| `FORMAT/ID` | the rsID, as in the ID column |
| Sample column | the output ID, for example `GWAScatalog_GCST002783` |

`INFO/HM` is relative to the alleles as the study reported them: `forward`; `swapped` (effect allele was the reference base); `complemented` (alleles now on the other strand from the study's, because the study reported the reverse strand or liftover crossed a minus-strand chain block; 2 complements cancel); `complemented_swapped`; `unverified` (no reference sequence at the position); and for palindromic SNPs `palindrome_forward`, `palindrome_complemented`, `palindrome_ambiguous` and `palindrome_ambiguous_complemented` (strand not inferable from frequency, but liftover complemented the alleles). Each takes `_swapped` when the alleles were then swapped.

The header declares the FASTA (`##reference`), `##contig` lines with `assembly=GRCh38`, the settings fingerprint, every INFO and FORMAT field, and a `##SAMPLE` line with a `##META` definition for each key: `TotalVariants`, `VariantsNotRead`, `VariantsNotHarmonised`, `HarmonisedVariants` (these add up), `SwitchedAlleles`, `AFFilled`, `StudyType`, `TotalCases` and `TotalControls` (as gwas2vcf writes them: a continuous trait puts its total sample size in `TotalControls`), `StudyID`, `Phenotype`, `Ancestry`, `Sex` and `PMID`.

Deliberate differences from an OpenGWAS file: GRCh38; `FORMAT/AF` filled from 1000 Genomes and marked in `FORMAT/AFS`; no 1000 Genomes frequencies in INFO; `INFO/SNPID` and `INFO/HM`; trait metadata and `AFFilled` in `##SAMPLE`; 25 primary contigs only; `SS` typed Float, as the specification says. R gwasvcf replaces per-variant `SS` with `TotalCases + TotalControls` (or `TotalControls`) when the header has them.

### Step log, QC and report

`standardised.qc.json` records every step with its rows in, rows removed by reason and changes made without removal, and the counts reconcile from input rows to written rows. It also holds the liftover and alignment counts, build detection, palindromic decisions, the frequency check, lambda GC, QC metrics, warnings, timings and the settings fingerprint. QC flags follow the OpenGWAS report thresholds and never stop a row: frequency correlation with 1000 Genomes below 0.7, lambda GC above 1.2, maximum sample size below 10,000, mean effect above 0.5 or 0.1 in absolute value, and more than 1,000 variants at P < 5e-8. The OpenGWAS LDSC and clumping flags are not computed, because the pipeline runs neither. The HTML report is 1 self-contained file per GWAS with the step table, flags, `##SAMPLE` counts, raw and standardised QQ plots, warnings and timings.

### Reference data

| Under `GWAS_REFERENCE_ROOT` | Use | Source |
|---|---|---|
| `genome/GRCh38.fa`, `GRCh37.fa` (+ `.fai`) | alignment and liftover target; liftover source | Ensembl release 113 primary assemblies |
| `genome/hg19ToHg38.over.chain.gz`, `hg38ToHg19.over.chain.gz` | liftover of GRCh37 inputs; placing the build-detection sample in GRCh37 | UCSC |
| `dbsnp/GRCh38/` | positions for rsIDs, rsIDs for positions and alleles, and 21.6 million retired rsIDs mapped to current ones; memory-mapped | NCBI dbSNP build 157 and its merge history |
| `1000genomes/GRCh38/` | ALT frequencies per super-population (AFR, AMR, EAS, EUR, SAS, ALL) of the 2,504 unrelated phase 3 samples | 1000 Genomes 30x on GRCh38, via plink2 resources |
| `build_detection/` | 2 million common biallelic SNPs (frequency 0.05–0.95) placed in both builds | 1000 Genomes GRCh38 EUR, fixed seed 821, lifted with BCFtools/liftover |

All are required except the hg38ToHg19 chain (needed only to rebuild the detection sample) and the merge map (without it, retired rsIDs are removed with a warning).

### Comparison with other pipelines

The comparison follows the code of MungeSumstats (Al-Murphy/MungeSumstats, November 2026), gwas2vcf and GwasDataImport (MRCIEU), and the GWAS Catalog sumstats-harmoniser (v1.1.12).

| Step | OpenGWAS (GwasDataImport, gwas2vcf) | MungeSumstats | GWAS Catalog harmoniser | GWASprocessing |
|---|---|---|---|---|
| Input | column positions supplied by the uploader | tables, VCF, OpenGWAS IDs; 410-entry header dictionary | GWAS-SSF columns only | tables and GWAS-VCF; alias dictionary and study YAML |
| Build | GwasDataImport matches HapMap3 rsIDs; gwas2vcf takes it as declared | 10,000 rsIDs against dbSNP 155, unseeded | from metadata or file name, no inference | paired GRCh37/GRCh38 sample, fixed seed; declared build as fallback |
| Target build | GRCh37 | GRCh37 or GRCh38 | GRCh38 | GRCh38 |
| Coordinates | UCSC chain before upload | rtracklayer with the Ensembl chain, alleles not rechecked | Ensembl positions by rsID, then pyliftover by position only (X, Y, MT fail) | dbSNP positions for rsID-only files; BCFtools/liftover with allele tracking |
| Strand | swap only; reverse-strand SNPs lost | swap only | both strands; palindromes by 1 file-wide strand assumption | swap and complement; palindromes per variant by frequency |
| rsIDs | dbSNP 153 by position | dbSNP 155 by position; variants without an rsID removed | Ensembl (dbSNP 151) by alleles; input rsIDs overwritten | dbSNP 157 by position and alleles; retired rsIDs mapped; variants without an rsID kept |
| Filters | P ≤ 0 or > 1, non-ACGT alleles, REF mismatch | INFO < 0.9, rsIDs with 3 or more alleles, X/Y/MT, missing SE or N | every unharmonised row | only unplaceable or unreadable rows and duplicates |
| P = 0 | row removed (gwas2vcf) | kept as 0 | −log10 P converted to P, underflows | LP recovered from BETA/SE |
| Duplicates | not removed | by rsID, then by position | not removed | by position and allele pair |
| Output | GWAS-VCF | table (GWAS-VCF optional) | GWAS-SSF table | GWAS-VCF |
| Accounting | HTML report and flags for a human reviewer | console messages | counts per harmonisation code | step log reconciling every row, flags, HTML report |

### Benchmark

2 public GWAS, each with a raw file, an OpenGWAS GWAS-VCF and a GWAS Catalog harmonised file, were processed by all 5 pipelines on 12 November 2026: Locke et al. 2015 BMI ([Genetic studies of body mass index yield new insights for obesity biology]; rsIDs without positions, continuous) and O'Mara et al. 2018 endometrial cancer ([Identification of nine new susceptibility loci for endometrial cancer]; GRCh37 positions, 1,014,582 indels, 12,906 cases and 108,979 controls). gwas2vcf 1.4.3 was run on the Locke harmonised file against GRCh38 and on the endometrial cancer raw file against GRCh37; MungeSumstats ran with its defaults and dbSNP 155, per chromosome for endometrial cancer because the whole file needs more than 32 GB of R memory. Edge cases neither GWAS has (P = 0, the OpenGWAS LP = 999 sentinel, the GWAS Catalog fallback) are covered by the test suite.

Variants were matched by rsID: the lead SNPs plus 10,000 rsIDs drawn from the GWASprocessing output with seed 821. Locke's 8 leads are known BMI loci (FTO, MC4R, TMEM18, GNPDA2, SEC16B, TFAP2B, BDNF, NEGR1); the 16 endometrial cancer leads are the raw file's genome-wide significant SNPs clumped against 1000 Genomes EUR (r² 0.001, 10,000 kb). Each record was aligned to the raw effect allele and compared on presence, alleles, effect sign, effect and SE (within 1e-4), −log10 P (within 0.01) and frequency (within 0.02).

| Pipeline | Locke: present of 10,008 | Locke: sign agrees | Locke leads | Endometrial: present of 10,016 | Endometrial: sign agrees | Endometrial leads |
|---|---|---|---|---|---|---|
| GWASprocessing | 10,008 | 100% | 8 of 8 | 10,016 | 99.97% | 16 of 16 |
| OpenGWAS | 9,921 | 100% | 8 of 8 | 9,973 | 100% | 15 of 16 |
| gwas2vcf | 9,930 | 100% | 8 of 8 | 10,016 | 100% | 16 of 16 |
| GWAS Catalog harmoniser | 9,985 | 100% | 8 of 8 | 9,984 | 100% | 16 of 16 |
| MungeSumstats | 4,618 | 100% | 7 of 8 | 4,304 | 100% | 4 of 16 |

Every difference is explained. The 3 GWASprocessing sign differences for endometrial cancer are palindromic SNPs on minus-strand chain blocks, where liftover complements the alleles; for an A/T or C/G SNP that looks like a swap to a comparison by allele letter, and `INFO/HM` records it. The OpenGWAS files have 27,303 (Locke) and 68,464 (endometrial cancer) fewer rows than the raw files, and the endometrial file lacks the top lead, rs11263761 at HNF1B (P = 3.2e-20). gwas2vcf cannot read 24,229 Locke harmonised rows that have no values, and the GWAS Catalog harmonised Locke file has no effect for 55 sampled variants. MungeSumstats removes every rsID that dbSNP 155 lists with 3 or more alleles (1,290,677 Locke rows, including the TMEM18 lead; 3,317,703 endometrial cancer rows, including 11 leads) and every variant with INFO below 0.9 (1,966,820 endometrial cancer rows, including 1 lead). No pipeline disagreed with the raw study on the effect, SE or P of any variant it kept.

The same 2-sample MR was then run on each pipeline's own files, Locke BMI on endometrial cancer, so that differences come only from processing; the estimates compare pipelines and are not a causal analysis. Instruments were exposure variants at P < 5e-8 clumped against 1000 Genomes EUR (r² 0.001, 10,000 kb), harmonised with TwoSampleMR 0.7.13 under action 2 (palindromes resolved by frequency) and action 3 (palindromes removed).

| Pipeline | Variants at P < 5e-8 | Instruments | IVW, action 2: b (se), SNPs | IVW, action 3: b (se), SNPs | MR-Egger, action 2 |
|---|---|---|---|---|---|
| GWASprocessing | 2,042 | 78 | 0.612 (0.078), 78 | 0.615 (0.079), 75 | 0.579 |
| OpenGWAS | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| gwas2vcf | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| GWAS Catalog harmoniser | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| MungeSumstats | 915 | 60 | 0.634 (0.085), 59 | 0.639 (0.088), 56 | 0.383 |

GWASprocessing, OpenGWAS, gwas2vcf and the GWAS Catalog select the same 78 instruments, 3 of them palindromic. Only the GWASprocessing exposure carries the study's allele frequency, so only there can TwoSampleMR resolve those 3 under action 2; with palindromes removed the 4 pipelines agree exactly. MungeSumstats keeps 45% of the genome-wide significant exposure variants and shares 38 of its 60 instruments with the others.

To rerun the benchmark, process `src/benchmark/benchmark.tsv`, run the other pipelines as above (`gwas2vcf_*.json`, `run_mungesumstats.R`), then compare and run the MR:

```bash
python3 src/benchmark/compare.py --raw <raw file> --raw-columns <rsID,EA,OA,EAF,BETA,SE,P column names> \
    --ours <ours>.gwas.vcf.gz --opengwas <ebi-a>.vcf.gz --gwas2vcf <gwas2vcf>.vcf.gz --catalog <harmonised>.h.tsv.gz \
    --mungesumstats <mungesumstats>.tsv.gz --leads <comma-separated rsIDs> --out-dir <folder> --expected src/benchmark/expected/<name>.tsv
python3 src/benchmark/to_mr_table.py --format vcf|catalog|mungesumstats --input <file> --output <table>.tsv.gz
Rscript src/benchmark/mr.R <pairs.tsv> <1000 Genomes EUR bfile> <plink 1.9> <results folder>
```

`pairs.tsv` has the columns pipeline, outcome, exposure_table and outcome_table.

### Tests

`python3 -m unittest discover -s pipeline/tests` runs 208 tests on small synthetic data in `pipeline/tests/examples/` and GRCh38 reference fixtures built by the tests. They cover column mapping and statistics, build detection, dbSNP positions and merged rsIDs, BCFtools/liftover and the strand recorded in `INFO/HM`, palindromic resolution, alignment and indel normalisation, the VCF header and body, reconciliation of the step log across pieces and workers, the downloaders (against local servers that mimic the GWAS Catalog and OpenGWAS), the setup and run scripts, readers in Python and R, and the Nextflow workflow end to end. Tests never use the network, and tests whose tools are missing (bcftools with the liftover plugin, Nextflow, R with gwasvcf) are skipped.
