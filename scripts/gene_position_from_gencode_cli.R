#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript gene_position_from_gencode_cli.R --project_name <name> [options]\n\n",
      "Options:\n",
      "  --project_name <name>   Overlap project name, for example: TWB1_LAA\n",
      "  --base_dir <dir>        Base directory. Default: /mnt/data/ai_agent/gene_analysis\n",
      "  --gene_file <path>      Input gene CSV. Default: overlap/<project_name>/overlap_gene_name_pairwise_union.csv\n",
      "  --version <v19|v26>     GENCODE version. v19 for GRCh37/hg19, v26 for GRCh38/hg38. Default: v26\n\n",
      "Examples:\n",
      "  Rscript gene_position_from_gencode_cli.R --project_name TWB1_LAA\n",
      "  Rscript gene_position_from_gencode_cli.R --project_name TWB1_LAA --version v19\n",
      "  Rscript gene_position_from_gencode_cli.R --project_name TWB1_LAA --gene_file /path/to/gene_list.csv --base_dir /mnt/data/ai_agent/gene_analysis\n"
    )
  )
}

parse_args <- function(args) {
  parsed <- list(
    project_name = NA_character_,
    base_dir = "/mnt/data/ai_agent/gene_analysis",
    gene_file = NA_character_,
    version = "v26"
  )

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (!startsWith(arg, "--")) {
      stop("All arguments must be named. Unknown argument: ", arg, call. = FALSE)
    }

    if (i == length(args)) {
      stop("Missing value after ", arg, call. = FALSE)
    }

    value <- args[[i + 1]]

    if (arg == "--project_name") {
      parsed$project_name <- value
    } else if (arg == "--base_dir") {
      parsed$base_dir <- value
    } else if (arg == "--gene_file") {
      parsed$gene_file <- value
    } else if (arg == "--version") {
      parsed$version <- tolower(value)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }

    i <- i + 2
  }

  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rtracklayer))

project_name <- args$project_name
base_dir <- args$base_dir
version <- args$version

if (is.na(project_name) || project_name == "") {
  stop("project_name is required.", call. = FALSE)
}

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

if (!version %in% c("v19", "v26")) {
  stop("version must be either 'v19' or 'v26'.", call. = FALSE)
}

gene_file <- args$gene_file
if (is.na(gene_file) || gene_file == "") {
  gene_file <- file.path(
    base_dir,
    "overlap",
    project_name,
    "overlap_gene_name_pairwise_union.csv"
  )
}

gtf_file <- if (version == "v19") {
  file.path(base_dir, "tools", "GENCODE", "gencode.v19.annotation.gtf")
} else {
  file.path(base_dir, "tools", "GENCODE", "gencode.v26.annotation.gtf")
}

if (!file.exists(gene_file)) {
  stop("gene_file not found: ", gene_file, call. = FALSE)
}

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file, call. = FALSE)
}

gene_dir <- dirname(gene_file)
gene_basename <- basename(gene_file)
gene_stem <- sub("\\.csv$", "", gene_basename, ignore.case = TRUE)
out_file <- file.path(gene_dir, paste0(gene_stem, "_with_position.csv"))

message("Project name: ", project_name)
message("Base dir: ", base_dir)
message("GENCODE version: ", version)
message("Gene file: ", gene_file)
message("GTF file: ", gtf_file)
message("Output file: ", out_file)

gene <- fread(gene_file)
gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

required_cols <- c("ensembl_gene_id_clean", "gene_name")
missing_cols <- setdiff(required_cols, colnames(gene))
if (length(missing_cols) > 0) {
  stop(
    "gene_file is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

gene_anno <- gtf_df %>%
  filter(type == "gene") %>%
  transmute(
    ensembl_gene_id = gene_id,
    ensembl_gene_id_clean = sub("\\..*$", "", gene_id),
    gene_name = gene_name,
    chr = as.character(seqnames),
    start = start,
    end = end,
    gene_type = gene_type
  )

result <- gene %>%
  left_join(
    gene_anno,
    by = c("ensembl_gene_id_clean" = "ensembl_gene_id_clean", "gene_name" = "gene_name")
  ) %>%
  distinct()

result_simple <- result[, c("ensembl_gene_id_clean", "gene_name", "chr", "start", "end", "gene_type")]

fwrite(result_simple, out_file)

message("Matched genes: ", sum(!is.na(result_simple$start)))
message("Unmatched genes: ", sum(is.na(result_simple$start)))
message("Done.")
