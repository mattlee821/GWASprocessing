#!/usr/bin/env Rscript
# Run MungeSumstats on a raw GWAS file with its defaults, for the benchmark.
# Usage: Rscript src/benchmark/run_mungesumstats.R <raw file> <output .tsv.gz> [threads] [ref_genome] [convert_ref_genome]
#   rsID-only input (Locke): ref_genome GRCh38 (default), no conversion.
#   GRCh37 positions (Yengo): ref_genome GRCh37, convert_ref_genome GRCh38 (MungeSumstats' default Ensembl chain).
# Needs MungeSumstats and the SNPlocs.Hsapiens.dbSNP155 and BSgenome packages for the builds used.
args <- commandArgs(trailingOnly = TRUE)
threads <- if (length(args) >= 3) as.integer(args[3]) else 1L
ref <- if (length(args) >= 4) args[4] else "GRCh38"
convert <- if (length(args) >= 5) args[5] else NULL
suppressPackageStartupMessages(library(MungeSumstats))
out <- format_sumstats(path = args[1], ref_genome = ref, convert_ref_genome = convert, dbSNP = 155, save_path = args[2],
                       nThread = threads, force_new = TRUE, return_data = FALSE, log_mungesumstats_msgs = FALSE)
cat("written:", out, "\n")
