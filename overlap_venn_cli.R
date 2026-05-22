#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript overlap_venn_cli.R --project_name_overlap <name> --base_dir <dir> [options]\n\n",
      "Options:\n",
      "  --project_name_v7 <name|None>      Project name for v7 TWAS tables\n",
      "  --project_name_v8 <name|None>      Project name for v8 TWAS tables\n",
      "  --project_name_cojo <name|None>    Project name for COJO table\n",
      "  --project_name_overlap <name>      Output folder name under overlap/\n",
      "  --base_dir <dir>                   Base directory\n",
      "  --ver <v7|v8|all>                 TWAS version selection. Default: all\n",
      "  --twas_method <fusion|spredixcan|all>  TWAS method selection. Default: all\n",
      "  --fdr_cutoff <value>              Default: 0.15\n",
      "  --cojo_p_cutoff <value>           Default: 1e-6\n\n",
      "Example:\n",
      "  Rscript overlap_venn_cli.R \\\n",
      "    --project_name_v7 TWB1_LAA \\\n",
      "    --project_name_v8 TWB1_LAA_hg38 \\\n",
      "    --project_name_cojo TWB1_LAA_hg38 \\\n",
      "    --project_name_overlap TWB1_LAA_overlap \\\n",
      "    --base_dir /mnt/data/ai_agent/gene_analysis \\\n",
      "    --ver all \\\n",
      "    --twas_method all \\\n",
      "    --fdr_cutoff 0.15 \\\n",
      "    --cojo_p_cutoff 1e-6\n"
    )
  )
}

parse_args <- function(args) {
  defaults <- list(
    project_name_v7 = "None",
    project_name_v8 = "None",
    project_name_cojo = "None",
    project_name_overlap = NA_character_,
    base_dir = "/mnt/data/ai_agent/gene_analysis",
    ver = "all",
    twas_method = "all",
    FDR_cutoff = 0.15,
    cojo_p_cutoff = 1e-6
  )

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  parsed <- defaults
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

    if (arg == "--project_name_v7") {
      parsed$project_name_v7 <- value
    } else if (arg == "--project_name_v8") {
      parsed$project_name_v8 <- value
    } else if (arg == "--project_name_cojo") {
      parsed$project_name_cojo <- value
    } else if (arg == "--project_name_overlap") {
      parsed$project_name_overlap <- value
    } else if (arg == "--base_dir") {
      parsed$base_dir <- value
    } else if (arg == "--ver") {
      parsed$ver <- tolower(value)
    } else if (arg == "--twas_method") {
      parsed$twas_method <- tolower(value)
    } else if (arg == "--fdr_cutoff") {
      parsed$FDR_cutoff <- as.numeric(value)
    } else if (arg == "--cojo_p_cutoff") {
      parsed$cojo_p_cutoff <- as.numeric(value)
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
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggvenn))

project_name_v7 <- args$project_name_v7
project_name_v8 <- args$project_name_v8
project_name_cojo <- args$project_name_cojo
project_name_overlap <- args$project_name_overlap
base_dir <- args$base_dir
ver <- args$ver
twas_method <- args$twas_method
FDR_cutoff <- as.numeric(args$FDR_cutoff)
cojo_p_cutoff <- as.numeric(args$cojo_p_cutoff)

if (is.na(project_name_overlap) || project_name_overlap == "") {
  stop("project_name_overlap is required.", call. = FALSE)
}

if (!ver %in% c("v7", "v8", "all")) {
  stop("ver must be v7, v8, or all.", call. = FALSE)
}

if (!twas_method %in% c("fusion", "spredixcan", "all")) {
  stop("twas_method must be fusion, spredixcan, or all.", call. = FALSE)
}

if (is.na(FDR_cutoff)) {
  stop("FDR_cutoff must be numeric.", call. = FALSE)
}

if (is.na(cojo_p_cutoff)) {
  stop("cojo_p_cutoff must be numeric.", call. = FALSE)
}

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

message("project_name_v7: ", project_name_v7)
message("project_name_v8: ", project_name_v8)
message("project_name_cojo: ", project_name_cojo)
message("project_name_overlap: ", project_name_overlap)
message("base_dir: ", base_dir)
message("version selection: ", ver)
message("twas method: ", twas_method)
message("FDR cutoff: ", FDR_cutoff)
message("COJO p cutoff: ", cojo_p_cutoff)

overlap_output_folder <- file.path(base_dir, "overlap", project_name_overlap)
dir.create(overlap_output_folder, recursive = TRUE, showWarnings = FALSE)

strip_version <- function(x) sub("\\..*", "", x)

convert_p_candidates <- function(p) {
  p_sci <- format(p, scientific = TRUE)
  p_compact <- gsub("e-0*", "e", p_sci)
  p_preserve <- gsub("e-", "e", p_sci)
  unique(c(paste0("p", p_compact), paste0("p", p_preserve)))
}

read_if_exists <- function(path, label) {
  if (!file.exists(path)) {
    warning("File not found for ", label, ": ", path)
    return(NULL)
  }
  fread(path)
}

get_fusion_tables <- function() {
  out <- list(v7 = NULL, v8 = NULL)

  if (twas_method %in% c("fusion", "all")) {
    if (ver %in% c("v7", "all") && project_name_v7 != "None") {
      path_v7 <- file.path(base_dir, "FUSION/result_v7", project_name_v7, "table", "FUSION_v7_TWAS.csv")
      out$v7 <- read_if_exists(path_v7, "FUSION v7")
    }
    if (ver %in% c("v8", "all") && project_name_v8 != "None") {
      path_v8 <- file.path(base_dir, "FUSION/result_v8", project_name_v8, "table", "FUSION_v8_TWAS.csv")
      out$v8 <- read_if_exists(path_v8, "FUSION v8")
    }
  }

  out
}

get_spredixcan_tables <- function() {
  out <- list(v7 = NULL, v8 = NULL)

  if (twas_method %in% c("spredixcan", "all")) {
    if (ver %in% c("v7", "all") && project_name_v7 != "None") {
      path_v7 <- file.path(base_dir, "S_PrediXcan/result_v7", project_name_v7, "table", "SPrediXcan_v7_TWAS.csv")
      out$v7 <- read_if_exists(path_v7, "S-PrediXcan v7")
    }
    if (ver %in% c("v8", "all") && project_name_v8 != "None") {
      path_v8 <- file.path(base_dir, "S_PrediXcan/result_v8", project_name_v8, "table", "SPrediXcan_v8_TWAS.csv")
      out$v8 <- read_if_exists(path_v8, "S-PrediXcan v8")
    }
  }

  out
}

get_cojo_table <- function() {
  if (project_name_cojo == "None") {
    return(NULL)
  }
  cojo_table_dir <- file.path(base_dir, "COJO", project_name_cojo, "table")
  candidate_files <- file.path(
    cojo_table_dir,
    paste0(convert_p_candidates(cojo_p_cutoff), "_annotation_genes.csv")
  )

  for (path in candidate_files) {
    if (file.exists(path)) {
      return(fread(path))
    }
  }

  warning(
    "File not found for COJO. Tried: ",
    paste(candidate_files, collapse = ", ")
  )
  NULL
}

fusion_tables <- get_fusion_tables()
spredixcan_tables <- get_spredixcan_tables()
cojo_table <- get_cojo_table()

available_methods <- c()
if (!is.null(fusion_tables$v7) || !is.null(fusion_tables$v8)) available_methods <- c(available_methods, "FUSION")
if (!is.null(spredixcan_tables$v7) || !is.null(spredixcan_tables$v8)) available_methods <- c(available_methods, "SPrediXcan")
if (!is.null(cojo_table)) available_methods <- c(available_methods, "GWAS")

if (length(available_methods) < 2) {
  stop("At least two methods with valid input tables are required to calculate overlap.", call. = FALSE)
}

get_fusion_set <- function(threshold) {
  sets <- c()
  if (!is.null(fusion_tables$v7)) {
    sets <- c(sets, strip_version(unique(fusion_tables$v7[FDR_tissue <= threshold, ensembl_gene_id])))
  }
  if (!is.null(fusion_tables$v8)) {
    sets <- c(sets, strip_version(unique(fusion_tables$v8[FDR_tissue <= threshold, ensembl_gene_id])))
  }
  unique(sets)
}

get_spredixcan_set <- function(threshold) {
  sets <- c()
  if (!is.null(spredixcan_tables$v7)) {
    sets <- c(sets, strip_version(unique(spredixcan_tables$v7[FDR_tissue <= threshold, ensembl_gene_id])))
  }
  if (!is.null(spredixcan_tables$v8)) {
    sets <- c(sets, strip_version(unique(spredixcan_tables$v8[FDR_tissue <= threshold, ensembl_gene_id])))
  }
  unique(sets)
}

get_cojo_set <- function() {
  if (is.null(cojo_table)) return(character(0))
  unique(strip_version(cojo_table$ensembl_gene_id_version))
}

build_method_sets <- function(threshold) {
  sets <- list()
  if ("FUSION" %in% available_methods) sets$FUSION <- get_fusion_set(threshold)
  if ("SPrediXcan" %in% available_methods) sets$SPrediXcan <- get_spredixcan_set(threshold)
  if ("GWAS" %in% available_methods) sets$GWAS <- get_cojo_set()
  sets
}

calc_overlap_summary <- function(method_sets, threshold) {
  summary_row <- data.table(FDR_tissue = threshold)

  for (nm in names(method_sets)) {
    summary_row[[paste0(nm, "_n")]] <- length(method_sets[[nm]])
  }

  if (length(method_sets) >= 2) {
    pair_names <- combn(names(method_sets), 2, simplify = FALSE)
    for (pair in pair_names) {
      overlap_n <- length(intersect(method_sets[[pair[1]]], method_sets[[pair[2]]]))
      summary_row[[paste(pair, collapse = "_")]] <- overlap_n
    }
  }

  all_overlap <- Reduce(intersect, method_sets)
  summary_row[["ALL_METHODS"]] <- length(all_overlap)

  list(summary = summary_row, all_overlap = all_overlap)
}

thresholds <- c(0.10, 0.15, 0.20, 0.25, FDR_cutoff)
thresholds <- sort(unique(thresholds))

all_results <- list()
all_gene_lists <- list()

for (t in thresholds) {
  method_sets <- build_method_sets(t)
  overlap_result <- calc_overlap_summary(method_sets, t)
  all_results[[as.character(t)]] <- overlap_result$summary
  all_gene_lists[[as.character(t)]] <- method_sets
}

final_overlap_table <- rbindlist(all_results, fill = TRUE)
write.csv(
  final_overlap_table,
  file.path(overlap_output_folder, "overlap_summary.csv"),
  row.names = FALSE
)

target_key <- as.character(FDR_cutoff)
target_sets <- all_gene_lists[[target_key]]
target_intersection <- Reduce(intersect, target_sets)

gene_map_list <- list()

if (!is.null(fusion_tables$v7)) {
  gene_map_list[["fusion_v7"]] <- data.table(
    ensembl_gene_id_clean = strip_version(fusion_tables$v7$ensembl_gene_id),
    gene_name = fusion_tables$v7$gene_name
  )
}
if (!is.null(fusion_tables$v8)) {
  gene_map_list[["fusion_v8"]] <- data.table(
    ensembl_gene_id_clean = strip_version(fusion_tables$v8$ensembl_gene_id),
    gene_name = fusion_tables$v8$gene_name
  )
}
if (!is.null(spredixcan_tables$v7)) {
  gene_map_list[["spredixcan_v7"]] <- data.table(
    ensembl_gene_id_clean = strip_version(spredixcan_tables$v7$ensembl_gene_id),
    gene_name = spredixcan_tables$v7$gene_name
  )
}
if (!is.null(spredixcan_tables$v8)) {
  gene_map_list[["spredixcan_v8"]] <- data.table(
    ensembl_gene_id_clean = strip_version(spredixcan_tables$v8$ensembl_gene_id),
    gene_name = spredixcan_tables$v8$gene_name
  )
}
if (!is.null(cojo_table)) {
  gene_map_list[["cojo"]] <- data.table(
    ensembl_gene_id_clean = strip_version(cojo_table$ensembl_gene_id_version),
    gene_name = cojo_table$hgnc_symbol
  )
}

all_gene_map <- unique(rbindlist(gene_map_list, fill = TRUE))
all_gene_map <- all_gene_map[gene_name != "" & !is.na(gene_name)]

overlap_gene_table <- all_gene_map[ensembl_gene_id_clean %in% target_intersection]
write.csv(
  overlap_gene_table,
  file.path(overlap_output_folder, "overlap_gene_name.csv"),
  row.names = FALSE
)

twas_gene_map_list <- list()

if (!is.null(fusion_tables$v7)) {
  twas_gene_map_list[["fusion_v7"]] <- data.table(
    ensembl_gene_id_clean = strip_version(fusion_tables$v7$ensembl_gene_id),
    gene_name = fusion_tables$v7$gene_name
  )
}
if (!is.null(fusion_tables$v8)) {
  twas_gene_map_list[["fusion_v8"]] <- data.table(
    ensembl_gene_id_clean = strip_version(fusion_tables$v8$ensembl_gene_id),
    gene_name = fusion_tables$v8$gene_name
  )
}
if (!is.null(spredixcan_tables$v7)) {
  twas_gene_map_list[["spredixcan_v7"]] <- data.table(
    ensembl_gene_id_clean = strip_version(spredixcan_tables$v7$ensembl_gene_id),
    gene_name = spredixcan_tables$v7$gene_name
  )
}
if (!is.null(spredixcan_tables$v8)) {
  twas_gene_map_list[["spredixcan_v8"]] <- data.table(
    ensembl_gene_id_clean = strip_version(spredixcan_tables$v8$ensembl_gene_id),
    gene_name = spredixcan_tables$v8$gene_name
  )
}

if (length(twas_gene_map_list) > 0) {
  twas_detected_gene_table <- unique(rbindlist(twas_gene_map_list, fill = TRUE))
  twas_detected_gene_table <- twas_detected_gene_table[
    gene_name != "" & !is.na(gene_name)
  ]
  twas_detected_gene_table <- unique(
    twas_detected_gene_table[order(ensembl_gene_id_clean, gene_name)]
  )

  write.csv(
    twas_detected_gene_table,
    file.path(overlap_output_folder, "twas_detected_gene_name.csv"),
    row.names = FALSE
  )
}

write_gene_set_csv <- function(gene_ids, output_path) {
  gene_table <- all_gene_map[ensembl_gene_id_clean %in% gene_ids]
  gene_table <- unique(gene_table[order(ensembl_gene_id_clean, gene_name)])
  write.csv(gene_table, output_path, row.names = FALSE)
}

if (length(target_sets) >= 2) {
  pair_names <- combn(names(target_sets), 2, simplify = FALSE)
  pairwise_overlap_union <- character(0)

  for (pair in pair_names) {
    pair_label <- paste(pair, collapse = "_")
    pair_intersection <- intersect(target_sets[[pair[1]]], target_sets[[pair[2]]])
    pair_union <- union(target_sets[[pair[1]]], target_sets[[pair[2]]])
    pairwise_overlap_union <- union(pairwise_overlap_union, pair_intersection)

    write_gene_set_csv(
      pair_intersection,
      file.path(
        overlap_output_folder,
        paste0("overlap_gene_name_", pair_label, ".csv")
      )
    )

    write_gene_set_csv(
      pair_union,
      file.path(
        overlap_output_folder,
        paste0("union_gene_name_", pair_label, ".csv")
      )
    )
  }

  write_gene_set_csv(
    pairwise_overlap_union,
    file.path(
      overlap_output_folder,
      "overlap_gene_name_pairwise_union.csv"
    )
  )
}

venn_input <- target_sets

venn_file <- file.path(
  overlap_output_folder,
  paste0(
    "venn_",
    paste(names(venn_input), collapse = "_"),
    "_FDR",
    gsub("\\.", "", format(FDR_cutoff, nsmall = 2)),
    ".png"
  )
)

png(
  filename = venn_file,
  width = 3000,
  height = 3000,
  res = 600
)

print(
  ggvenn(
    venn_input,
    fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")[seq_along(venn_input)],
    stroke_size = 0.5,
    set_name_size = 4
  )
)

dev.off()

message("Available methods: ", paste(available_methods, collapse = ", "))
message("Overlap summary: ", file.path(overlap_output_folder, "overlap_summary.csv"))
message("Overlap gene names: ", file.path(overlap_output_folder, "overlap_gene_name.csv"))
message("Venn plot: ", venn_file)
