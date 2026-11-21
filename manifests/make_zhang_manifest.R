#!/usr/bin/env Rscript
# Build the GWASprocessing manifest for the Zhang et al. 2022 (PMID 35501419) ARIC SomaScan cis-pGWAS.
#
# Usage:
#   Rscript manifests/make_zhang_manifest.R --files DIR_OR_LIST --out manifests/zhang_EUR.tsv
#       [--ancestry EUR] [--study-id ID] [--location-dir DIR] [--usefuldata PATH] [--st01 seqid.txt]
#       [--yaml manifests/zhang_2022_35501419.yaml] [--limit N]
#
#   --files         a folder of the per-protein files (SeqId_<a>_<b>.PHENO1.glm.linear), or a text file listing their names
#   --ancestry      manifest ancestry: EUR (or EA / european, default) or AFR (or AA / african)
#   --study-id      output study folder name (default zhang_2022_35501419_<POP>, e.g. zhang_2022_35501419_EUR or zhang_2022_35501419_AFR)
#   --location-dir  folder written into GWAS_location (default: --files when it is a folder)
#   --usefuldata    path of the usefuldata package (default: the installed package; install it with remotes::install_github("mattlee821/usefuldata"))
#   --st01          the study's seqid.txt (seqid_in_sample, uniprot_id, entrezgenesymbol, ...), to report where its UniProt differs
#   --limit         the first N files in SeqId order (a subset to test)
#
# Each file is identified by its SeqId (2 numeric parts, 10000_28), which is how deCODE names the files. The name of the
# standardised output is UNIPROT_SEQID: "_" separates the 2 values and "-" separates the parts of a value, so the
# UniProt accession(s) from usefuldata::somalogic_seqid_uniprot_current (the most recent SomaLogic annotation, retired
# accessions replaced by current ones; several accessions joined with "-") and the SeqId written with "-" (10000-28),
# for example P43320_10000-28 and P01011-P07288_4153-11. A file whose SeqId has no current accession is left out of the manifest and listed in <out>.not_written.tsv.
# Also written: <out>.names.tsv (every file with its name and the mapping behind it) and <out>.provenance.txt.

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NA_character_) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1L]
}
files_arg <- option("files")
out <- option("out")
if (is.na(files_arg) || is.na(out)) stop("Usage: --files DIR_OR_LIST --out MANIFEST.tsv (see the header of this script)", call. = FALSE)
location_dir <- option("location-dir", if (dir.exists(files_arg)) normalizePath(files_arg) else NA_character_)
if (is.na(location_dir)) stop("--location-dir is needed when --files is a list", call. = FALSE)
yaml <- option("yaml", "manifests/zhang_2022_35501419.yaml")
ancestry_arg <- option("ancestry", "EUR")
pop_code <- switch(tolower(ancestry_arg),
                   "eur" = "EUR", "ea" = "EUR", "european" = "EUR",
                   "afr" = "AFR", "aa" = "AFR", "african" = "AFR", "african-american" = "AFR",
                   stop("--ancestry for Zhang 2022 must be EUR (EA) or AFR (AA), got: ", ancestry_arg, call. = FALSE))
study_id <- option("study-id", paste0("zhang_2022_35501419_", pop_code))
st01_path <- option("st01")
limit <- as.integer(option("limit", NA_character_))

# the mapping ====
package_path <- option("usefuldata", NA_character_)
mapping <- NULL
if (is.na(package_path) && requireNamespace("usefuldata", quietly = TRUE)) {
  mapping <- usefuldata::somalogic_seqid_uniprot_current
  package_version <- as.character(utils::packageVersion("usefuldata"))
} else {
  if (is.na(package_path)) stop("The R package usefuldata is needed: remotes::install_github(\"mattlee821/usefuldata\"), or give --usefuldata PATH to a checkout of it", call. = FALSE)
  rda <- file.path(package_path, "data", "somalogic_seqid_uniprot_current.rda")
  if (!file.exists(rda)) stop("usefuldata mapping not found: ", rda, call. = FALSE)
  environment_ <- new.env()
  load(rda, envir = environment_)
  mapping <- environment_$somalogic_seqid_uniprot_current
  description <- read.dcf(file.path(package_path, "DESCRIPTION"), fields = "Version")
  commit <- tryCatch(system2("git", c("-C", package_path, "log", "-1", "--format=%h"), stdout = TRUE, stderr = FALSE), error = function(e) "")
  package_version <- paste0(description[1, 1], if (length(commit) && nzchar(commit[1])) paste0(" (commit ", commit[1], ")") else "")
}
mapping$key <- gsub("-", "_", mapping$SeqId)

# the files ====
file_names <- if (dir.exists(files_arg)) list.files(files_arg, pattern = "\\.glm\\.linear$") else
  basename(trimws(readLines(files_arg)[nzchar(trimws(readLines(files_arg)))]))
file_names <- file_names[grepl("\\.glm\\.linear$", file_names)]
parsed <- regmatches(file_names, regexec("^SeqId_([0-9]+_[0-9]+)\\.PHENO1\\.glm\\.linear$", file_names))
has_seqid <- lengths(parsed) == 2L
not_written <- data.frame(file = file_names[!has_seqid], reason = rep("no SeqId at the start of the file name", sum(!has_seqid)))
files <- data.frame(file = file_names[has_seqid], key = vapply(parsed[has_seqid], `[`, character(1), 2L), stringsAsFactors = FALSE)
files <- files[order(as.integer(sub("_.*", "", files$key)), as.integer(sub(".*_", "", files$key))), , drop = FALSE]
if (!is.na(limit)) files <- head(files, limit)
if (anyDuplicated(files$key)) stop("Two files share a SeqId: ", paste(files$key[duplicated(files$key)], collapse = ", "), call. = FALSE)

names_table <- merge(files, mapping[c("key", "SeqId", "Target", "TargetFullName", "ExampleADAT", "UNIPROT_SOURCE", "UNIPROT_CURRENT",
                                      "UNIPROT_STATUS", "UNIPROT_CHANGES", "UNIPROT_RELEASE", "UNIPROT_CHECKED")],
                     by = "key", all.x = TRUE, sort = FALSE)
names_table <- names_table[match(files$key, names_table$key), , drop = FALSE]
names_table$output_name <- ifelse(is.na(names_table$UNIPROT_CURRENT), NA_character_, paste0(names_table$UNIPROT_CURRENT, "_", names_table$SeqId))
reason <- ifelse(is.na(names_table$SeqId), "SeqId is not in the usefuldata SomaScan mapping",
                 ifelse(is.na(names_table$UNIPROT_CURRENT), "UniProt has no current accession for this SeqId", NA_character_))
not_written <- rbind(not_written, data.frame(file = names_table$file[!is.na(reason)], reason = reason[!is.na(reason)]))
ok <- is.na(reason)
if (anyDuplicated(names_table$output_name[ok])) stop("Two files get the same output name: ", paste(unique(names_table$output_name[ok][duplicated(names_table$output_name[ok])]), collapse = ", "), call. = FALSE)

# the manifest ====
manifest <- data.frame(
  GWAS_location = file.path(location_dir, names_table$file[ok]),
  phenotype = paste0(ifelse(is.na(names_table$TargetFullName[ok]), names_table$Target[ok], names_table$TargetFullName[ok]),
                     " plasma protein level, cis-pGWAS (SomaScan aptamer ", names_table$SeqId[ok], ")"),
  ancestry = pop_code, sex = "combined", author = "Zhang", year = "2022", PMID = "35501419",
  source_type = "local", input_build = "hg38", population = pop_code, delimiter = "tab", format = "table",
  study_yaml = yaml, output_name = names_table$output_name[ok], study_id = study_id,
  stringsAsFactors = FALSE
)
write.table(manifest, out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(names_table, paste0(out, ".names.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(not_written, paste0(out, ".not_written.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# the cross-check against the paper's table, and the provenance ====
st01_note <- "not compared (no --st01)"
if (!is.na(st01_path)) {
  st01 <- read.delim(st01_path, stringsAsFactors = FALSE, quote = "")
  st01$key <- sub("^SeqId_", "", st01$seqid_in_sample)
  st01$UniProt <- st01$uniprot_id
  paper <- st01$UniProt[match(names_table$key, st01$key)]
  differs <- !is.na(paper) & !is.na(names_table$UNIPROT_CURRENT) & paper != names_table$UNIPROT_CURRENT
  names_table$seqid.txt_UniProt <- paper
  write.table(names_table, paste0(out, ".names.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  st01_note <- sprintf("%d of %d files have the same UniProt accession(s) as the paper's seqid.txt; %d differ%s", sum(!differs & !is.na(paper)),
                       nrow(names_table), sum(differs),
                       if (any(differs)) paste0(" (", paste(head(names_table$key[differs], 20), collapse = ", "), ")") else "")
}
writeLines(c(
  paste("manifest:", out),
  paste("written on:", Sys.Date()),
  paste("files:", nrow(files), "| in the manifest:", nrow(manifest), "| not written:", nrow(not_written)),
  paste("usefuldata:", package_version),
  paste("UniProt release:", paste(unique(stats::na.omit(names_table$UNIPROT_RELEASE)), collapse = ","),
        "| checked:", paste(unique(stats::na.omit(names_table$UNIPROT_CHECKED)), collapse = ",")),
  paste("names whose UniProt was updated from the SomaLogic annotation:", sum(names_table$UNIPROT_STATUS == "updated", na.rm = TRUE)),
  paste("seqid.txt:", st01_note)
), paste0(out, ".provenance.txt"))
cat(readLines(paste0(out, ".provenance.txt")), sep = "\n")
