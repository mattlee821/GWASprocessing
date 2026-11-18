#!/usr/bin/env Rscript
# Build the GWASprocessing manifest for the Sun et al. 2023 (PMID 37794186) UKB-PPP Olink pGWAS.
#
# Usage:
#   Rscript manifests/make_sun_manifest.R (--synapse-folder SYN_ID | --files DIR_OR_LIST) --out manifests/sun_combined.tsv --map-dir DIR
#       [--location-dir DIR] [--yaml manifests/sun_2023_37794186.yaml] [--ancestry combined] [--limit N]
#
#   --synapse-folder  Synapse folder ID holding the per-protein .tar files (uses SYNAPSE_TOKEN from .env or environment)
#   --files           a folder of the per-protein tar files (<Assay>_<UniProt>_<OlinkID>_v1_<Panel>.tar), a 2-column TSV of
#                     synapse_id and file name (from download_synapse.py --list-children), or a text file listing file names
#   --map-dir         folder with olink_protein_map_1.5k_v1.tsv and olink_protein_map_3k_v1.tsv (Synapse syn51396703, Protein annotation)
#   --location-dir    folder written into GWAS_location for local file lists (default: --files when it is a folder)
#   --ancestry        manifest ancestry (default combined = the multi-ancestry GWAS; european, african, ... for the others)
#   --limit           the first N files in OlinkID order
#
# The name of the standardised output is UNIPROT_OLINKID, for example P16442_OID30675. "_" separates the 2 values and "-"
# joins the parts of a value: the UniProt accessions of a complex are joined with "-" (Q14213-Q8NEV9_OID21389), and an
# isoform keeps its own "-" suffix (O43521-2_OID20513): an accession part made only of digits after a "-" is an isoform
# number, never a second accession. The OlinkID is unique to each file (6 proteins, such as TNF and IL6, are on all 4
# panels with 4 OlinkIDs), so names never repeat. Files whose OlinkID is not in the protein maps are left out and listed
# in <out>.not_written.tsv. Also written: <out>.names.tsv and <out>.provenance.txt.

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default = NA_character_) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[i + 1L]
}
files_arg <- option("files"); synapse_folder <- option("synapse-folder"); out <- option("out"); map_dir <- option("map-dir")
if ((is.na(files_arg) && is.na(synapse_folder)) || is.na(out) || is.na(map_dir)) {
  stop("Usage: (--synapse-folder SYN_ID | --files DIR_OR_LIST) --out MANIFEST.tsv --map-dir DIR (see the header of this script)", call. = FALSE)
}
yaml <- option("yaml", "manifests/sun_2023_37794186.yaml")
ancestry <- option("ancestry", "combined")
limit <- as.integer(option("limit", NA_character_))

read_map <- function(name) read.delim(file.path(map_dir, name), stringsAsFactors = FALSE, quote = "", check.names = FALSE)
map <- rbind(read_map("olink_protein_map_1.5k_v1.tsv")[c("UKBPPP_ProteinID", "olink_target_fullname", "OlinkID", "UniProt", "Assay", "Panel")],
             read_map("olink_protein_map_3k_v1.tsv")[c("UKBPPP_ProteinID", "olink_target_fullname", "OlinkID", "UniProt", "Assay", "Panel")])
map <- map[!duplicated(map$UKBPPP_ProteinID), , drop = FALSE]
if (anyDuplicated(map$OlinkID)) stop("OlinkID is not unique in the protein maps", call. = FALSE)
map$UNIPROT <- gsub("_", "-", map$UniProt)   # a complex "Q14213_Q8NEV9" -> "Q14213-Q8NEV9"; an isoform "O43521-2" is kept

if (!is.na(synapse_folder)) {
  if (file.exists(".env") && !nzchar(Sys.getenv("SYNAPSE_TOKEN")) && !nzchar(Sys.getenv("SYNAPSE_AUTH_TOKEN"))) readRenviron(".env")
  env_dir <- Sys.getenv("GWAS_ENV", "")
  py <- if (nzchar(env_dir) && file.exists(file.path(env_dir, "bin", "python3"))) file.path(env_dir, "bin", "python3") else
    if (file.exists(".tools/env/bin/python3")) ".tools/env/bin/python3" else Sys.which("python3")
  tmp_children <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_children), add = TRUE)
  status <- system2(py, c("pipeline/bin/download_synapse.py", "--list-children", synapse_folder, "--out", tmp_children))
  if (status != 0L) stop("Failed to list Synapse folder ", synapse_folder, call. = FALSE)
  files_arg <- tmp_children
}

synapse_ids <- NULL
if (dir.exists(files_arg)) {
  file_names <- list.files(files_arg, pattern = "\\.tar$")
} else {
  raw_lines <- trimws(readLines(files_arg))
  raw_lines <- raw_lines[nzchar(raw_lines)]
  parts <- strsplit(raw_lines, "[\t ]+")
  if (length(parts) > 0L && all(lengths(parts) >= 2L)) {
    c1 <- vapply(parts, `[`, character(1), 1L)
    c2 <- vapply(parts, `[`, character(1), 2L)
    if (tolower(c1[1L]) %in% c("synapse_id", "id", "synid")) { c1 <- c1[-1L]; c2 <- c2[-1L] }
    if (length(c1) > 0L && all(grepl("^syn[0-9]+(\\.[0-9]+)?$", c1, ignore.case = TRUE))) {
      synapse_ids <- tolower(c1)
      file_names <- basename(c2)
    }
  }
  if (is.null(synapse_ids)) file_names <- basename(raw_lines)
}
is_tar <- grepl("\\.tar$", file_names)
file_names <- file_names[is_tar]
if (!is.null(synapse_ids)) synapse_ids <- synapse_ids[is_tar]
location_dir <- option("location-dir", if (!is.na(files_arg) && dir.exists(files_arg)) normalizePath(files_arg) else NA_character_)
if (is.null(synapse_ids) && is.na(location_dir)) stop("--location-dir is needed when --files is a plain filename list", call. = FALSE)

parsed <- regmatches(file_names, regexec("_(OID[0-9]+)_v1_", file_names))
has_id <- lengths(parsed) == 2L
not_written <- data.frame(file = file_names[!has_id], reason = rep("no OlinkID in the file name", sum(!has_id)))
files <- data.frame(file = file_names[has_id], OlinkID = vapply(parsed[has_id], `[`, character(1), 2L), stringsAsFactors = FALSE)
if (!is.null(synapse_ids)) files$synapse_id <- synapse_ids[has_id]
files <- files[order(files$OlinkID), , drop = FALSE]
if (!is.na(limit)) files <- head(files, limit)
if (anyDuplicated(files$OlinkID)) stop("Two files share an OlinkID", call. = FALSE)

names_table <- merge(files, map, by = "OlinkID", all.x = TRUE, sort = FALSE)
names_table <- names_table[match(files$OlinkID, names_table$OlinkID), , drop = FALSE]
names_table$output_name <- ifelse(is.na(names_table$UNIPROT), NA_character_, paste0(names_table$UNIPROT, "_", names_table$OlinkID))
# the file name carries the same assay, UniProt and panel as the map: report any file that disagrees
expected <- paste0(names_table$Assay, "_", gsub("-", "_", names_table$UniProt), "_", names_table$OlinkID, "_v1_", names_table$Panel, ".tar")
disagree <- !is.na(names_table$Assay) & expected != names_table$file
reason <- ifelse(is.na(names_table$UNIPROT), "OlinkID is not in the protein maps", NA_character_)
not_written <- rbind(not_written, data.frame(file = names_table$file[!is.na(reason)], reason = reason[!is.na(reason)]))
ok <- is.na(reason)
if (anyDuplicated(names_table$output_name[ok])) stop("Two files get the same output name", call. = FALSE)

manifest <- data.frame(
  GWAS_location = if ("synapse_id" %in% names(names_table)) names_table$synapse_id[ok] else file.path(location_dir, names_table$file[ok]),
  phenotype = paste0(names_table$olink_target_fullname[ok], " plasma protein level (Olink Explore, ", sub("_", " ", names_table$Panel[ok]), " panel, ", names_table$OlinkID[ok], ")"),
  ancestry = ancestry, sex = "combined", author = "Sun", year = "2023", PMID = "37794186",
  source_type = if ("synapse_id" %in% names(names_table)) "synapse" else "local",
  input_build = "hg38", population = toupper(if (ancestry == "combined") "ALL" else substr(ancestry, 1, 3)), delimiter = "space", format = "table",
  study_yaml = yaml, output_name = names_table$output_name[ok], stringsAsFactors = FALSE
)
write.table(manifest, out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(names_table, paste0(out, ".names.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(not_written, paste0(out, ".not_written.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  paste("manifest:", out), paste("written on:", Sys.Date()),
  paste("files:", nrow(files), "| in the manifest:", nrow(manifest), "| not written:", nrow(not_written)),
  paste("file names that disagree with the protein map:", sum(disagree)),
  paste("names with a complex (2 or more accessions):", sum(grepl("^[A-Z0-9]+-[A-Z][A-Z0-9]+_", manifest$output_name))),
  paste("names with an isoform suffix:", sum(grepl("^[A-Z0-9]{6,10}-[0-9]+(-|_)", manifest$output_name)))
), paste0(out, ".provenance.txt"))
cat(readLines(paste0(out, ".provenance.txt")), sep = "\n")
