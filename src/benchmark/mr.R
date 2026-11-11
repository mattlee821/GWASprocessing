#!/usr/bin/env Rscript
# Two-sample MR run identically on each pipeline's processed exposure and outcome, to compare the pipelines downstream.
#
# Usage: Rscript src/benchmark/mr.R <pairs.tsv> <bfile prefix (1kg.v3 EUR)> <plink binary> <output folder>
#   pairs.tsv columns: pipeline, outcome, exposure_table, outcome_table (tables from to_mr_table.py)
#
# For every pair: instruments are exposure variants with P < 5e-8, clumped locally with plink 1.9 against the same
# LD reference (r2 0.001, 10,000 kb, the TwoSampleMR defaults) so differences come only from processing; the outcome
# is harmonised with TwoSampleMR::harmonise_data, action 2 (palindromes resolved by allele frequency, the default)
# and action 3 (all palindromes dropped, sensitivity); MR by IVW, MR-Egger, weighted median and weighted mode, with
# Cochran's Q and the MR-Egger intercept. Writes instruments.tsv, harmonised.tsv, results.tsv, heterogeneity.tsv,
# pleiotropy.tsv. The estimates compare pipelines; they are not a causal analysis.
suppressPackageStartupMessages({ library(TwoSampleMR); library(ieugwasr); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
pairs <- fread(args[1]); bfile <- args[2]; plink <- args[3]; out <- args[4]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

read_side <- function(path, type, snps = NULL) {
  d <- fread(path, colClasses = list(character = c("SNP", "effect_allele", "other_allele", "chr")))
  if (!is.null(snps)) d <- d[SNP %in% snps]
  format_data(as.data.frame(d), type = type, snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
              effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval",
              chr_col = "chr", pos_col = "pos")
}

instruments_cache <- list()
all_inst <- list(); all_harm <- list(); all_res <- list(); all_het <- list(); all_pleio <- list()
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  key <- paste(p$pipeline, p$exposure_table)
  if (is.null(instruments_cache[[key]])) {
    exposure <- read_side(p$exposure_table, "exposure")
    significant <- exposure[exposure$pval.exposure < 5e-8, ]
    clumped <- ld_clump_local(data.frame(rsid = significant$SNP, pval = significant$pval.exposure), clump_kb = 10000,
                              clump_r2 = 0.001, clump_p = 1, bfile = bfile, plink_bin = plink)
    instruments_cache[[key]] <- significant[significant$SNP %in% clumped$rsid, ]
    message(p$pipeline, ": ", nrow(significant), " exposure variants at P < 5e-8, ", nrow(instruments_cache[[key]]),
            " instruments after clumping")
  }
  inst <- instruments_cache[[key]]
  all_inst[[length(all_inst) + 1]] <- cbind(pipeline = p$pipeline, inst[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                                         "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")])
  outcome <- read_side(p$outcome_table, "outcome", inst$SNP)
  for (action in c(2, 3)) {
    h <- harmonise_data(inst, outcome, action = action)
    tag <- data.frame(pipeline = p$pipeline, outcome_trait = p$outcome, action = action)
    all_harm[[length(all_harm) + 1]] <- cbind(tag[rep(1, nrow(h)), ], h[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                              "effect_allele.outcome", "other_allele.outcome", "beta.exposure", "beta.outcome",
                                              "eaf.exposure", "eaf.outcome", "palindromic", "ambiguous", "mr_keep")])
    res <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
    all_res[[length(all_res) + 1]] <- cbind(tag[rep(1, nrow(res)), ], res[, c("method", "nsnp", "b", "se", "pval")],
                                            instruments = nrow(inst), in_outcome = nrow(outcome))
    het <- mr_heterogeneity(h)
    if (nrow(het)) all_het[[length(all_het) + 1]] <- cbind(tag[rep(1, nrow(het)), ], het[, c("method", "Q", "Q_df", "Q_pval")])
    pl <- mr_pleiotropy_test(h)
    if (nrow(pl)) all_pleio[[length(all_pleio) + 1]] <- cbind(tag, pl[, c("egger_intercept", "se", "pval")])
  }
}
fwrite(rbindlist(all_inst), file.path(out, "instruments.tsv"), sep = "\t")
fwrite(rbindlist(all_harm), file.path(out, "harmonised.tsv"), sep = "\t")
fwrite(rbindlist(all_res), file.path(out, "results.tsv"), sep = "\t")
fwrite(rbindlist(all_het), file.path(out, "heterogeneity.tsv"), sep = "\t")
fwrite(rbindlist(all_pleio), file.path(out, "pleiotropy.tsv"), sep = "\t")
res <- rbindlist(all_res)
print(res[method == "Inverse variance weighted", .(pipeline, outcome_trait, action, nsnp, b = round(b, 4), se = round(se, 4), pval = signif(pval, 3))])
