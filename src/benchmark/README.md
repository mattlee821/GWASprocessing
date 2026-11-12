# Cross-pipeline benchmark

This is a check of the pipeline, kept outside it. 2 public GWAS are standardised by GWASprocessing and by 4 other pipelines; a sample of variants is compared with the raw study values, and the same Mendelian randomisation is run on each pipeline's output. The sampled GWASprocessing values are stored in `expected/`, so a new install (a new HPC, a new version of a dependency) only has to run GWASprocessing on the 2 GWAS and compare with them.

## Benchmark GWAS

Both GWAS have a raw author file, an OpenGWAS GWAS-VCF and a GWAS Catalog harmonised file, so every pipeline can be compared on the same study.

| GWAS | Raw file | OpenGWAS | GWAS Catalog harmonised | What it tests in GWASprocessing |
|---|---|---|---|---|
| Locke et al. 2015 BMI, [Genetic studies of body mass index yield new insights for obesity biology] | GCST002783 (GIANT, rsIDs without positions, HapMap CEU frequencies) | ebi-a-GCST002783 | 25673413-GCST002783-EFO_0004340.h.tsv.gz | GRCh38 positions from dbSNP, rsIDs retired by dbSNP merges, a continuous trait; OpenGWAS VCF input with GRCh37 positions |
| O'Mara et al. 2018 endometrial cancer, [Identification of nine new susceptibility loci for endometrial cancer] | GCST006464 (GRCh37 positions, imputation INFO, 1,014,582 indels) | ebi-a-GCST006464 | 30093612-GCST006464-EFO_1001512.h.tsv.gz | build detection, BCFtools/liftover (including minus-strand blocks), indels, a case-control trait (12,906 cases, 108,979 controls), INFO as FORMAT/SI |

Edge cases that neither GWAS has (P = 0, the OpenGWAS LP = 999 sentinel, falling back to the GWAS Catalog harmonised file) are covered by the unit tests.

## Check a new install

1. Standardise both GWAS with the benchmark manifest (it names the endometrial cancer study YAML):

```bash
bash src/standardise.sh --manifest src/benchmark/benchmark.tsv --row-limit 0
```

2. Compare the 2 GWAS Catalog outputs with the stored values:

```bash
python3 src/benchmark/check_install.py --expected src/benchmark/expected/locke_2015.tsv \
    --output <work root>/locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.gwas.vcf.gz
python3 src/benchmark/check_install.py --expected src/benchmark/expected/omara_2018.tsv \
    --output <work root>/omara_2018_30093612/standardise/GWAScatalog_GCST006464/standardised.gwas.vcf.gz
```

Every expected variant (10,008 for Locke, 10,016 for endometrial cancer) must be present with the same GRCh38 position, alleles, ES, SE, LP, AF and AFS; the script exits 1 and lists the differences otherwise. A new reference (a new dbSNP build) can legitimately change positions or rsIDs: rerun the full benchmark below and update `expected/` when that is intended.

Locke is run with ancestry EUR because its frequency column comes from HapMap CEU. With ALL, 19 of the sampled palindromic SNPs were resolved to the other strand, because the study frequencies were being compared with 1000 Genomes frequencies of a different population. The ancestry in the manifest should match the population of the file's frequency column.

## Pipelines compared

| Pipeline | Locke output | Endometrial cancer output |
|---|---|---|
| GWASprocessing | GCST002783 raw file, ancestry EUR | GCST006464 raw file, ancestry EUR |
| OpenGWAS | ebi-a-GCST002783.vcf.gz (GRCh37) | ebi-a-GCST006464.vcf.gz (GRCh37) |
| gwas2vcf 1.4.3 | the GWAS Catalog harmonised file against the GRCh38 FASTA (`gwas2vcf_locke_grch38.json`) | the raw file against the GRCh37 FASTA (`gwas2vcf_endometrial_grch37.json`) |
| GWAS Catalog harmoniser | the published harmonised file (GRCh38) | the published harmonised file (GRCh38) |
| MungeSumstats | raw file, `ref_genome = "GRCh38"`, dbSNP 155 (`run_mungesumstats.R`) | raw file, GRCh37 converted to GRCh38, dbSNP 155, run 1 chromosome at a time |

MungeSumstats needs more than 32 GB of R vector memory for the 9.5 million endometrial cancer rows, so on a laptop it is run per chromosome and the outputs are joined; its checks are per variant apart from duplicate removal, so the result matches a whole-file run except for any rsID reported on 2 chromosomes.

## Variant comparison

```bash
python3 src/benchmark/compare.py --raw All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz \
    --raw-columns SNP,A1,A2,Freq1.Hapmap,b,se,p --ours <ours>.gwas.vcf.gz --opengwas ebi-a-GCST002783.vcf.gz \
    --gwas2vcf locke_gwas2vcf_grch38.vcf.gz --catalog 25673413-GCST002783-EFO_0004340.h.tsv.gz \
    --mungesumstats locke_mss_grch38.tsv.gz \
    --leads rs1558902,rs6567160,rs13021737,rs10938397,rs543874,rs2207139,rs11030104,rs3101336 \
    --out-dir compare/locke --expected src/benchmark/expected/locke_2015.tsv
python3 src/benchmark/compare.py --raw GCST006464_GRCh37.tsv.gz \
    --raw-columns variant_id,effect_allele,other_allele,effect_allele_frequency,beta,standard_error,p_value \
    --ours <ours>.gwas.vcf.gz --opengwas ebi-a-GCST006464.vcf.gz --gwas2vcf endometrial_gwas2vcf_grch37.vcf.gz \
    --catalog 30093612-GCST006464-EFO_1001512.h.tsv.gz --mungesumstats endometrial_mss_grch38.tsv.gz \
    --leads rs11263761,rs7981863,rs1740828,rs17601876,rs4733613,rs3184504,rs2747716,rs9668337,rs4767257,rs882380,rs937213,rs1679014,rs10835920,rs148261157,rs113998067,rs139584729 \
    --out-dir compare/endometrial --expected src/benchmark/expected/omara_2018.tsv
```

Variants are matched by rsID, because the outputs are in different builds. The sample is the lead SNPs plus 10,000 rsIDs drawn from the GWASprocessing output with seed 821. Locke's 8 leads are known BMI loci (FTO, MC4R, TMEM18, GNPDA2, SEC16B, TFAP2B, BDNF, NEGR1). The 16 endometrial cancer leads are the raw file's genome-wide significant SNPs clumped against the 1000 Genomes EUR reference (r² 0.001, 10,000 kb). Each record is aligned to the raw effect allele (palindromic SNPs by exact allele order only) and compared on presence, allele pair, effect sign, effect (within 1e-4), SE (within 1e-4), −log10 P (within 0.01) and allele frequency (within 0.02).

### Locke 2015, 12 November 2026 (`expected/locke_2015_summary.tsv`)

| Pipeline | Present (of 10,008) | Sign agrees | Effect equal | SE, LP equal | Leads present |
|---|---|---|---|---|---|
| GWASprocessing | 10,008 | 100% | 100% | 100% | 8 of 8 |
| OpenGWAS | 9,921 | 100% | 100% | 100% | 8 of 8 |
| gwas2vcf | 9,930 | 100% | 100% | 100% | 8 of 8 |
| GWAS Catalog harmoniser | 9,985 | 100% | 99.45% | 100% | 8 of 8 |
| MungeSumstats | 4,618 | 100% | 100% | 100% | 7 of 8 |

### O'Mara 2018, 12 November 2026 (`expected/omara_2018_summary.tsv`)

| Pipeline | Present (of 10,016) | Sign agrees | Effect equal | SE, LP equal | Leads present |
|---|---|---|---|---|---|
| GWASprocessing | 10,016 | 99.97% | 99.97% | 100% | 16 of 16 |
| OpenGWAS | 9,973 | 100% | 100% | 100% | 15 of 16 |
| gwas2vcf | 10,016 | 100% | 100% | 100% | 16 of 16 |
| GWAS Catalog harmoniser | 9,984 | 100% | 100% | 100% | 16 of 16 |
| MungeSumstats | 4,304 | 100% | 100% | 100% | 4 of 16 |

Every difference is explained:
- GWASprocessing, endometrial cancer, 3 sign disagreements: all palindromic SNPs on GRCh37-to-GRCh38 chain blocks on the minus strand (rs8039713, rs12146369, rs4495303). Liftover complements their alleles, which for an A/T or C/G SNP looks like a swap to a comparison by allele letter; the association itself is unchanged. `INFO/HM` records the complement (`palindrome_complemented`, `palindrome_ambiguous_complemented`).
- OpenGWAS, Locke: 87 sampled rsIDs are not in ebi-a-GCST002783, which has 27,303 fewer rows than the raw file. Endometrial cancer: the OpenGWAS file has 68,464 fewer rows than the raw file, and the top lead, rs11263761 at HNF1B (P = 3.2e-20), is not in it.
- gwas2vcf, Locke: 24,229 rows of the harmonised file have no values and cannot be read.
- GWAS Catalog harmoniser, Locke: the effect is empty for 55 of the sampled variants.
- MungeSumstats removes every rsID that dbSNP 155 lists with 3 or more alleles (Locke 1,290,677 of 2,555,510 rows, including the TMEM18 lead rs13021737; endometrial cancer 3,317,703 rows, including 11 of the 16 leads), and by default every variant with INFO below 0.9 (endometrial cancer 1,966,820 rows, including the lead rs148261157, INFO 0.893).

No pipeline disagreed with the raw study on the effect, SE or P of any variant it kept and could align by allele.

## Mendelian randomisation

The same 2-sample MR is run on each pipeline's own exposure and outcome, Locke BMI on endometrial cancer, so that differences come only from processing. The estimates compare pipelines; they are not a causal analysis.

```bash
python3 src/benchmark/to_mr_table.py --format vcf|catalog|mungesumstats --input FILE --output TABLE.tsv.gz
Rscript src/benchmark/mr.R pairs.tsv <1000 Genomes EUR bfile> <plink 1.9> results/
```

`pairs.tsv` has columns pipeline, outcome, exposure_table and outcome_table. Instruments are exposure variants with P < 5e-8, clumped locally with plink 1.9 against the 1000 Genomes EUR reference (r² 0.001, 10,000 kb, the TwoSampleMR defaults). The outcome is harmonised with TwoSampleMR 0.7.13 `harmonise_data`, action 2 (palindromes resolved by allele frequency, the default) and action 3 (all palindromes removed); estimates by IVW, MR-Egger, weighted median and weighted mode.

### Result, 12 November 2026

| Pipeline | Exposure variants at P < 5e-8 | Instruments | IVW, action 2: b (se), SNPs | IVW, action 3: b (se), SNPs | MR-Egger, action 2 |
|---|---|---|---|---|---|
| GWASprocessing | 2,042 | 78 | 0.612 (0.078), 78 | 0.615 (0.079), 75 | 0.579 |
| OpenGWAS | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| gwas2vcf | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| GWAS Catalog harmoniser | 2,042 | 78 | 0.615 (0.079), 75 | 0.615 (0.079), 75 | 0.591 |
| MungeSumstats | 915 | 60 | 0.634 (0.085), 59 | 0.639 (0.088), 56 | 0.383 |

GWASprocessing, OpenGWAS, gwas2vcf and the GWAS Catalog harmoniser select the same 78 instruments. 3 are palindromic (rs14810, rs17001654, rs9304665). Only the GWASprocessing exposure carries the study's allele frequency (the other 3 have none for Locke), so only there can TwoSampleMR resolve them under action 2; with all palindromes removed (action 3) the 4 pipelines give the same estimate. Before the dbSNP merge map was added, GWASprocessing removed the Locke lead rs12016871 (merged into rs9581854 in dbSNP build 131) and clumping chose a neighbour instead; it now keeps it. MungeSumstats keeps 45% of the genome-wide significant exposure variants, shares 38 of its 60 instruments with the other pipelines, and its MR-Egger estimate moves furthest.
