# Cross-pipeline benchmark

This is a check of the pipeline, not part of it. The same public GWAS is standardised by GWASprocessing and by 4 other pipelines, and a sample of variants is compared with the raw study values. The sampled GWASprocessing values are stored in `expected/`, so a new install (a new HPC, a new version of a dependency) only has to run GWASprocessing on the benchmark GWAS and compare with them.

## Benchmark GWAS

| Case | Why | Routes it tests |
|---|---|---|
| Yengo et al. 2018 BMI (primary) | The raw file is available directly (cnsgenomics.com), from the GWAS Catalog (GCST006900, the same file) and processed by OpenGWAS (ieu-b-40), so every source can be compared with the raw values | GRCh37 positions: build detection from the data, BCFtools/liftover; direct URL, GWAS Catalog and OpenGWAS staging; GWAS-VCF input |
| Locke et al. 2015 BMI (second) | The GWAS Catalog has harmonised it, so the GWAS Catalog harmoniser can be compared | rsIDs without positions: GRCh38 positions from dbSNP |

## Check a new install

1. Standardise Yengo 2018 BMI from its raw file, for example with the manifest row
   `https://cnsgenomics.com/data/yengo_et_al_2018_hmg/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz	BMI	EUR	Yengo	2018	30124842`,
   and Locke 2015 BMI from the GWAS Catalog (`GCST002783	BMI	EUR	Locke	2015	25673413`).
2. Compare each output with the stored values:

```bash
python3 src/benchmark/check_install.py --expected src/benchmark/expected/yengo_2018.tsv --output <Yengo standardised.gwas.vcf.gz>
python3 src/benchmark/check_install.py --expected src/benchmark/expected/locke_2015.tsv --output <Locke standardised.gwas.vcf.gz>
```

Each of the 10,008 expected variants must be present with the same GRCh38 position, alleles, ES, SE, LP, AF and AFS. The script exits 1 and lists the differences otherwise. A change in a reference (a new dbSNP build) can legitimately change positions or rsIDs; rerun the full benchmark below and update `expected/` when that is intended.

## Yengo 2018: 3 sources, 1 result (11 November 2026)

The cnsgenomics and GWAS Catalog downloads are byte-identical (MD5 ccca6c58bc45a53aec16619497da3f00). Standardised from the raw file (either source) and from the OpenGWAS ieu-b-40 GWAS-VCF, all 3 outputs hold the same 2,336,125 variants (of 2,336,269 raw rows) with identical positions, alleles, rsIDs, ES, SE, AF, AFS and SS; LP agrees within the rounding of the OpenGWAS file. The build was detected from the data (GRCh37) and lifted with BCFtools/liftover in each case; 5,408 variants had their alleles swapped because the reference allele differs between GRCh37 and GRCh38.

Variants compared with the raw study values (8 known BMI lead SNPs plus 10,000 rsIDs, seed 821; `expected/yengo_2018_summary.tsv`):

| Pipeline | Found by rsID (of 10,008) | Sign agrees | Effect, SE, LP equal | Leads present |
|---|---|---|---|---|
| GWASprocessing | 10,008 | 99.98% | 99.98%, 100%, 100% | 8 of 8 |
| OpenGWAS ieu-b-40 | 9,432 | 100% | 100% | 8 of 8 |
| gwas2vcf on GRCh37 | 10,007 | 100% | 100% | 7 of 8 |
| MungeSumstats (GRCh37 to GRCh38) | 4,830 | 100% | 100% | 7 of 8 |

Explained differences:
- GWASprocessing, 2 sign disagreements: both palindromic SNPs (rs10910842, rs7480344) in GRCh37-to-GRCh38 chain blocks on the minus strand, where liftover complements the alleles. For an A/T or C/G SNP that is indistinguishable from a swap, so the comparison, which matches palindromes by exact allele order, reports a sign change; the association itself is unchanged (the effect is reported for the complementary allele on the GRCh38 strand).
- OpenGWAS, 576 not found: all are present, but under a different rsID. OpenGWAS annotated rsIDs from dbSNP by position only, so for 5.8% of variants it took another dbSNP record at the same position (for example rs1335002413 instead of rs10503747). GWASprocessing matches dbSNP on position and alleles and kept the study rsID for all 10,008.
- OpenGWAS also stores LP = 999 for P = 0 (FTO variants); GWASprocessing recovers LP from ES/SE in that case (LP 458.9 for rs1558902), whichever source it reads.
- gwas2vcf drops rows with P = 0, including the FTO lead rs1558902.
- MungeSumstats removes every rsID that dbSNP 155 lists with 3 or more alleles (1,208,942 of 2,336,269 rows), including the TMEM18 lead rs13021737.

## Rerun the full benchmark

Locke 2015 inputs, all public (Yengo: the raw file, ieu-b-40, gwas2vcf on GRCh37 with `gwas2vcf_yengo_grch37.json`, and MungeSumstats with `GRCh37 GRCh38`):

| Pipeline | Output compared | How it was made |
|---|---|---|
| GWASprocessing | GCST002783 raw file standardised | this pipeline |
| OpenGWAS | `ieu-a-2.vcf.gz` (GRCh37) | downloaded from OpenGWAS (gwas2vcf plus dbSNP annotation) |
| gwas2vcf on GRCh38 | the GWAS Catalog harmonised file converted with gwas2vcf 1.4.3 against the GRCh38 FASTA | `gwas2vcf_locke_grch38.json` |
| GWAS Catalog harmoniser | `25673413-GCST002783-EFO_0004340.h.tsv.gz` | downloaded from the GWAS Catalog FTP |
| MungeSumstats | the raw file with `ref_genome = "GRCh38"`, dbSNP 155, defaults otherwise | `run_mungesumstats.R` |

```bash
python3 src/benchmark/compare.py --raw All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz \
    --raw-columns SNP,A1,A2,Freq1.Hapmap,b,se,p --ours standardised.gwas.vcf.gz --opengwas ieu-a-2.vcf.gz \
    --gwas2vcf locke_gwas2vcf_grch38.vcf.gz --catalog 25673413-GCST002783-EFO_0004340.h.tsv.gz \
    --mungesumstats locke_mss_grch38.tsv.gz \
    --leads rs1558902,rs6567160,rs13021737,rs10938397,rs543874,rs2207139,rs11030104,rs3101336 \
    --out-dir compare --expected src/benchmark/expected/locke_2015.tsv
```

Variants are matched by rsID. The sample is 8 known BMI lead SNPs (FTO, MC4R, TMEM18, GNPDA2, SEC16B, TFAP2B, BDNF, NEGR1) plus 10,000 rsIDs drawn from the GWASprocessing output with seed 821. Each pipeline's record is aligned to the raw effect allele (palindromic SNPs by exact allele order only) and compared on presence, allele pair, effect sign, effect (within 1e-4), SE (within 1e-4), −log10 P (within 0.01) and allele frequency (within 0.02).

## Locke 2015 result, 11 November 2026 (`expected/locke_2015_summary.tsv`)

| Pipeline | Present (of 10,008) | Sign agrees | Effect equal | SE, LP equal | Leads present |
|---|---|---|---|---|---|
| GWASprocessing | 10,008 | 100% | 100% | 100% | 8 of 8 |
| OpenGWAS | 10,008 | 100% | 100% | 100% | 8 of 8 |
| gwas2vcf on GRCh38 | 9,930 | 100% | 100% | 100% | 8 of 8 |
| GWAS Catalog harmoniser | 9,985 | 100% | 99.45% | 100% | 8 of 8 |
| MungeSumstats | 4,618 | 100% | 100% | 100% | 7 of 8 |

Every difference is explained: gwas2vcf could not read 24,229 harmonised rows that have no values; the GWAS Catalog harmoniser leaves the effect empty for 55 of the sampled variants; MungeSumstats removes every rsID that dbSNP 155 lists with 3 or more alleles (1,290,677 of 2,555,510 rows, including TMEM18 rs13021737). No pipeline disagreed with the raw study on the effect allele, sign, SE or P of any variant it kept.
