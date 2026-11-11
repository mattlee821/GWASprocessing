#!/usr/bin/env Rscript
# Run MungeSumstats on a raw GWAS file with GRCh38 references and its defaults otherwise, for the benchmark.
# Usage: Rscript src/benchmark/run_mungesumstats.R <raw file> <output .tsv.gz> [threads]
# Needs MungeSumstats, SNPlocs.Hsapiens.dbSNP155.GRCh38 and BSgenome.Hsapiens.NCBI.GRCh38.
args <- commandArgs(trailingOnly = TRUE)
threads <- if (length(args) >= 3) as.integer(args[3]) else 1L
suppressPackageStartupMessages(library(MungeSumstats))
out <- format_sumstats(path = args[1], ref_genome = "GRCh38", dbSNP = 155, save_path = args[2], nThread = threads,
                       force_new = TRUE, return_data = FALSE, log_mungesumstats_msgs = FALSE)
cat("written:", out, "\n")
