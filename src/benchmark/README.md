# Cross-pipeline benchmark

This is a check of the pipeline, not part of it. The same public GWAS is standardised by GWASprocessing and by 4 other pipelines, and a sample of variants is compared with the raw study values. The sampled GWASprocessing values are stored in `expected/`, so a new install (a new HPC, a new version of a dependency) only has to run GWASprocessing on the benchmark GWAS and compare with them.

## Check a new install

1. Standardise Locke et al. 2015 BMI (GWAS Catalog GCST002783, raw file) with the pipeline, for example with a manifest row `GCST002783	BMI	EUR	Locke	2015	25673413`.
2. Compare the output with the stored values:

```bash
python3 src/benchmark/check_install.py --expected src/benchmark/expected/locke_2015.tsv \
    --output GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.gwas.vcf.gz
```

Every one of the 10,008 expected variants must be present with the same GRCh38 position, alleles, ES, SE, LP, AF and AFS. The script exits 1 and lists the differences otherwise. A change in a reference (a new dbSNP build) can legitimately change positions or rsIDs; rerun the full benchmark below and update `expected/` when that is intended.

## Rerun the full benchmark

Inputs, all public:

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

## Result, 11 November 2026 (`expected/locke_2015_summary.tsv`)

| Pipeline | Present (of 10,008) | Sign agrees | Effect equal | SE, LP equal | Leads present |
|---|---|---|---|---|---|
| GWASprocessing | 10,008 | 100% | 100% | 100% | 8 of 8 |
| OpenGWAS | 10,008 | 100% | 100% | 100% | 8 of 8 |
| gwas2vcf on GRCh38 | 9,930 | 100% | 100% | 100% | 8 of 8 |
| GWAS Catalog harmoniser | 9,985 | 100% | 99.45% | 100% | 8 of 8 |
| MungeSumstats | 4,618 | 100% | 100% | 100% | 7 of 8 |

Every difference is explained: gwas2vcf could not read 24,229 harmonised rows that have no values; the GWAS Catalog harmoniser leaves the effect empty for 55 of the sampled variants; MungeSumstats removes every rsID that dbSNP 155 lists with 3 or more alleles (1,290,677 of 2,555,510 rows, including TMEM18 rs13021737). No pipeline disagreed with the raw study on the effect allele, sign, SE or P of any variant it kept.
