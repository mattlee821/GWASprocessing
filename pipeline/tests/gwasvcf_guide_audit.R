#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
}))

default_vcf <- file.path(
  "GWAS", "locke_2015_25673413", "standardise", "GWAScatalog_GCST002783",
  "standardised.gwas.vcf.gz"
)
default_output <- file.path(
  "GWAS", "locke_2015_25673413", "standardise", "GWAScatalog_GCST002783",
  "standardised.flat.r.tsv"
)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- if (length(args) >= 1) args[[1]] else default_vcf
output_file <- if (length(args) >= 2) args[[2]] else default_output

if (!file.exists(vcf_file)) {
  stop("VCF not found: ", vcf_file, call. = FALSE)
}

vcf <- readVcf(vcf_file)
if (!inherits(vcf, "CollapsedVCF")) {
  stop("readVcf did not return a CollapsedVCF object", call. = FALSE)
}

required_format <- c("ES", "SE", "LP", "SS")
missing_format <- setdiff(required_format, names(geno(vcf)))
if (length(missing_format) > 0) {
  stop("Missing FORMAT fields: ", paste(missing_format, collapse = ", "), call. = FALSE)
}

gr <- vcf_to_granges(vcf)
rr <- rowRanges(vcf)
variant_ids <- names(rr)
if (is.null(variant_ids)) {
  variant_ids <- rep(".", length(rr))
}

first_sample <- 1L
extract_geno <- function(name) {
  if (!name %in% names(geno(vcf))) {
    return(rep(NA, nrow(vcf)))
  }
  value <- geno(vcf)[[name]]
  if (length(dim(value)) >= 2L) {
    return(as.vector(value[, first_sample]))
  }
  as.vector(value)
}

extract_integer <- function(name) {
  value <- suppressWarnings(as.numeric(extract_geno(name)))
  rounded <- floor(value + 0.5)
  bad <- !is.na(value) & abs(value - rounded) > 1e-8
  if (any(bad)) {
    stop("Sample size is not an integer after VCF standardisation", call. = FALSE)
  }
  as.integer(rounded)
}

lp <- suppressWarnings(as.numeric(extract_geno("LP")))
flat <- data.frame(
  CHR = as.character(seqnames(gr)),
  POS = as.integer(start(gr)),
  SNP = as.character(variant_ids),
  REF = as.character(mcols(rr)$REF),
  ALT = as.character(unlist(mcols(rr)$ALT, use.names = FALSE)),
  EAF = suppressWarnings(as.numeric(extract_geno("AF"))),
  BETA = suppressWarnings(as.numeric(extract_geno("ES"))),
  SE = suppressWarnings(as.numeric(extract_geno("SE"))),
  P = 10 ^ (-lp),
  N = extract_integer("SS"),
  stringsAsFactors = FALSE
)

required_columns <- c("CHR", "POS", "SNP", "REF", "ALT", "EAF", "BETA", "SE", "P", "N")
missing_columns <- setdiff(required_columns, names(flat))
if (length(missing_columns) > 0) {
  stop("Flat output is missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (nrow(flat) != nrow(vcf)) {
  stop("Flat output row count does not match VCF row count", call. = FALSE)
}
if (any(is.na(flat$CHR) | is.na(flat$POS) | is.na(flat$REF) | is.na(flat$ALT))) {
  stop("Flat output has missing required variant coordinates or alleles", call. = FALSE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(flat, output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("VCF rows:", nrow(vcf), "\n")
cat("Flat rows:", nrow(flat), "\n")
cat("Samples:", paste(samples(header(vcf)), collapse = ","), "\n")
cat("Output:", output_file, "\n")
