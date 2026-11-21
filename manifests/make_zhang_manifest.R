#!/usr/bin/env Rscript
# Build the GWASprocessing manifest for the Zhang et al. 2022 (PMID 35501419) ARIC SomaScan cis-pGWAS.
#
# Usage:
#   Rscript manifests/make_zhang_manifest.R [--ancestry EUR|AFR] [--archive-url URL | --files DIR_OR_LIST]
#       --out manifests/zhang_EUR.tsv [--study-id ID] [--location-dir DIR] [--usefuldata PATH]
#       [--st01 seqid.txt] [--yaml manifests/zhang_2022_35501419.yaml] [--limit N]
#
#   --ancestry      manifest ancestry: EUR (or EA / european, default) or AFR (or AA / african)
#   --archive-url   remote .zip archive URL holding the per-protein .glm.linear files (default when --files is omitted:
#                   https://jh-pwas.s3.amazonaws.com/results/EA.zip for EUR, or .../AA.zip for AFR; each row gets
#                   GWAS_location = <archive-url>#<member> with source_type = direct_url so no pre-download is needed)
#   --files         optional local folder of SeqId_<a>_<b>.PHENO1.glm.linear files, or a text file listing names/URLs
#   --study-id      output study folder name (default zhang_2022_35501419_<POP>, e.g. zhang_2022_35501419_EUR or zhang_2022_35501419_AFR)
#   --location-dir  folder written into GWAS_location for local plain filename lists (default: --files when it is a folder)
#   --usefuldata    path of the usefuldata package (default: the installed package; install it with remotes::install_github("mattlee821/usefuldata"))
#   --st01          the study's seqid.txt (seqid_in_sample, uniprot_id, entrezgenesymbol, ...), or URL; fetched automatically
#                   from https://jh-pwas.s3.amazonaws.com/results/seqid.txt when using the default archive URL
#   --limit         the first N files in SeqId order (a subset to test)
#
# Each file is identified by its SeqId (2 numeric parts, 10000_28). The name of the standardised output is UNIPROT_SEQID:
# "_" separates the 2 values and "-" separates the parts of a value, using usefuldata::somalogic_seqid_uniprot_current
# (for example P43320_10000-28 and P01011-P07288_4153-11). A file whose SeqId has no current accession is left out of the
# manifest and listed in <out>.not_written.tsv. Also written: <out>.names.tsv and <out>.provenance.txt.

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NA_character_) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1L]
}
out <- option("out")
if (is.na(out)) stop("Usage: [--ancestry EUR|AFR] [--archive-url URL | --files DIR_OR_LIST] --out MANIFEST.tsv (see script header)", call. = FALSE)
yaml <- option("yaml", "manifests/zhang_2022_35501419.yaml")
ancestry_arg <- option("ancestry", "EUR")
pop_code <- switch(tolower(ancestry_arg),
                   "eur" = "EUR", "ea" = "EUR", "european" = "EUR",
                   "afr" = "AFR", "aa" = "AFR", "african" = "AFR", "african-american" = "AFR",
                   stop("--ancestry for Zhang 2022 must be EUR (EA) or AFR (AA), got: ", ancestry_arg, call. = FALSE))
study_id <- option("study-id", paste0("zhang_2022_35501419_", pop_code))
files_arg <- option("files")
default_archive <- if (pop_code == "EUR") "https://jh-pwas.s3.amazonaws.com/results/EA.zip" else "https://jh-pwas.s3.amazonaws.com/results/AA.zip"
archive_url <- option("archive-url", if (is.na(files_arg)) default_archive else NA_character_)
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
urls <- NULL
if (!is.na(archive_url)) {
  if (file.exists(".env")) readRenviron(".env")
  env_dir <- Sys.getenv("GWAS_ENV", "")
  py <- if (nzchar(env_dir) && file.exists(file.path(env_dir, "bin", "python3"))) file.path(env_dir, "bin", "python3") else
    if (file.exists(".tools/env/bin/python3")) ".tools/env/bin/python3" else Sys.which("python3")
  tmp_members <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp_members), add = TRUE)
  status <- system2(py, c("pipeline/bin/download_archive.py", "--list-members", archive_url, "--out", tmp_members))
  if (status != 0L) stop("Failed to list archive members from ", archive_url, call. = FALSE)
  members <- trimws(readLines(tmp_members))
  members <- members[nzchar(members)]
  urls <- paste0(archive_url, "#", members)
  file_names <- basename(members)
  if (is.na(st01_path) && grepl("^https?://jh-pwas\\.s3\\.amazonaws\\.com/results/", archive_url)) {
    tmp_seqid <- tempfile(fileext = ".txt")
    on.exit(unlink(tmp_seqid), add = TRUE)
    ok_dl <- tryCatch({ utils::download.file("https://jh-pwas.s3.amazonaws.com/results/seqid.txt", tmp_seqid, quiet = TRUE); TRUE },
                      error = function(e) FALSE)
    if (ok_dl && file.exists(tmp_seqid) && file.info(tmp_seqid)$size > 0) st01_path <- tmp_seqid
  }
} else if (dir.exists(files_arg)) {
  file_names <- list.files(files_arg, pattern = "\\.glm\\.linear$")
} else {
  raw_lines <- trimws(readLines(files_arg))
  raw_lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]
  if (length(raw_lines) > 0L && all(grepl("^https?://", raw_lines, ignore.case = TRUE))) {
    urls <- raw_lines
    file_names <- basename(sub("^.*#", "", raw_lines))
  } else {
    file_names <- basename(raw_lines)
  }
}

is_glm <- grepl("\\.glm\\.linear$", file_names)
file_names <- file_names[is_glm]
if (!is.null(urls)) urls <- urls[is_glm]
location_dir <- option("location-dir", if (!is.na(files_arg) && dir.exists(files_arg)) normalizePath(files_arg) else NA_character_)
if (is.null(urls) && is.na(location_dir)) stop("--location-dir is needed when --files is a plain filename list", call. = FALSE)

parsed <- regmatches(file_names, regexec("^SeqId_([0-9]+_[0-9]+)\\.PHENO1\\.glm\\.linear$", file_names))
has_seqid <- lengths(parsed) == 2L
not_written <- data.frame(file = file_names[!has_seqid], reason = rep("no SeqId at the start of the file name", sum(!has_seqid)))
files <- data.frame(file = file_names[has_seqid], key = vapply(parsed[has_seqid], `[`, character(1), 2L), stringsAsFactors = FALSE)
if (!is.null(urls)) files$url <- urls[has_seqid]
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
  GWAS_location = if ("url" %in% names(names_table)) names_table$url[ok] else file.path(location_dir, names_table$file[ok]),
  phenotype = paste0(ifelse(is.na(names_table$TargetFullName[ok]), names_table$Target[ok], names_table$TargetFullName[ok]),
                     " plasma protein level, cis-pGWAS (SomaScan aptamer ", names_table$SeqId[ok], ")"),
  ancestry = pop_code, sex = "combined", author = "Zhang", year = "2022", PMID = "35501419",
  source_type = if ("url" %in% names(names_table)) "direct_url" else "local",
  input_build = "hg38", population = pop_code, delimiter = "tab", format = "table",
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
