#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript pathway_enrichment_cli.R <project_name_overlap> [base_dir] [fdr_cutoff]\n",
      "  Rscript pathway_enrichment_cli.R --project_name_overlap <name> [--base_dir <dir>] [--fdr_cutoff <value>] [--gobp_fdr_cutoff <value>] [--kegg_fdr_cutoff <value>] [--gene_file <name|path>] [--universe_file <name|path>]\n\n",
      "Arguments:\n",
      "  project_name_overlap  Overlap project name under overlap/\n",
      "  base_dir              Base directory. Default: /mnt/data/ai_agent/gene_analysis\n",
      "  fdr_cutoff            Default cutoff used for output-file naming. Default: 0.15\n",
      "  gobp_fdr_cutoff       GOBP-specific cutoff for output-file naming. Default: same as fdr_cutoff\n",
      "  kegg_fdr_cutoff       KEGG-specific cutoff for output-file naming. Default: same as fdr_cutoff\n",
      "  gene_file             Candidate gene file. Default: overlap_gene_name_pairwise_union.csv\n",
      "                        - Bare filename (no '/') is resolved under overlap/<project_name_overlap>/\n",
      "                        - A path containing '/' is used as-is (absolute or relative)\n",
      "  universe_file         Universe gene file. Default: twas_detected_gene_name.csv\n",
      "                        - Same resolution rules as gene_file\n\n",
      "Examples:\n",
      "  Rscript pathway_enrichment_cli.R TWB1_LAA /mnt/data/ai_agent/gene_analysis 0.2\n",
      "  Rscript pathway_enrichment_cli.R --project_name_overlap TWB1_LAA --base_dir /mnt/data/ai_agent/gene_analysis --fdr_cutoff 0.2\n",
      "  # Use a different candidate set inside overlap/<project>/:\n",
      "  Rscript pathway_enrichment_cli.R --project_name_overlap TWB1_LAA --base_dir /your/base --gene_file coloc_mr_pass_gene.csv\n",
      "  # Specify both files by filename:\n",
      "  Rscript pathway_enrichment_cli.R --project_name_overlap TWB1_LAA --base_dir /your/base --gene_file coloc_mr_pass_gene.csv --universe_file twas_detected_gene_name.csv\n",
      "  # Different FDR tags for GOBP vs KEGG outputs:\n",
      "  Rscript pathway_enrichment_cli.R --project_name_overlap TWB1_LAA --base_dir /your/base --gobp_fdr_cutoff 0.05 --kegg_fdr_cutoff 0.2\n"
    )
  )
}

parse_args <- function(args) {
  default_base_dir <- "/mnt/data/ai_agent/gene_analysis"

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  if (any(grepl("^--", args))) {
    parsed <- list(
      project_name_overlap = NA_character_,
      base_dir = default_base_dir,
      fdr_cutoff = 0.15,
      gobp_fdr_cutoff = NA_real_,
      kegg_fdr_cutoff = NA_real_,
      gene_file = NA_character_,
      universe_file = NA_character_
    )
    i <- 1

    while (i <= length(args)) {
      arg <- args[[i]]

      if (arg %in% c("--project_name_overlap", "--project-name-overlap")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$project_name_overlap <- args[[i + 1]]
      } else if (arg %in% c("--base_dir", "--base-dir")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$base_dir <- args[[i + 1]]
      } else if (arg %in% c("--fdr_cutoff", "--fdr-cutoff")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$fdr_cutoff <- as.numeric(args[[i + 1]])
      } else if (arg %in% c("--gobp_fdr_cutoff", "--gobp-fdr-cutoff")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$gobp_fdr_cutoff <- as.numeric(args[[i + 1]])
      } else if (arg %in% c("--kegg_fdr_cutoff", "--kegg-fdr-cutoff")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$kegg_fdr_cutoff <- as.numeric(args[[i + 1]])
      } else if (arg %in% c("--gene_file", "--gene-file")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$gene_file <- args[[i + 1]]
      } else if (arg %in% c("--universe_file", "--universe-file")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$universe_file <- args[[i + 1]]
      } else {
        stop("Unknown argument: ", arg, call. = FALSE)
      }

      i <- i + 2
    }

    if (is.na(parsed$project_name_overlap) || parsed$project_name_overlap == "") {
      stop("project_name_overlap is required.", call. = FALSE)
    }

    return(parsed)
  }

  list(
    project_name_overlap = args[[1]],
    base_dir = ifelse(length(args) >= 2, args[[2]], default_base_dir),
    fdr_cutoff = ifelse(length(args) >= 3, as.numeric(args[[3]]), 0.15),
    gobp_fdr_cutoff = NA_real_,
    kegg_fdr_cutoff = NA_real_,
    gene_file = NA_character_,
    universe_file = NA_character_
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(KEGGREST))

project_name_overlap <- args$project_name_overlap
base_dir <- args$base_dir
fdr_cutoff <- as.numeric(args$fdr_cutoff)

if (is.na(fdr_cutoff)) {
  stop("fdr_cutoff must be numeric.", call. = FALSE)
}

# Per-analysis FDR cutoffs: if not specified by the user, fall back to fdr_cutoff.
gobp_fdr_cutoff <- as.numeric(args$gobp_fdr_cutoff)
if (is.na(gobp_fdr_cutoff)) gobp_fdr_cutoff <- fdr_cutoff
kegg_fdr_cutoff <- as.numeric(args$kegg_fdr_cutoff)
if (is.na(kegg_fdr_cutoff)) kegg_fdr_cutoff <- fdr_cutoff

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

overlap_input_folder <- file.path(base_dir, "overlap", project_name_overlap)
pathway_root_folder <- file.path(base_dir, "pathway")
gobp_output_folder <- file.path(pathway_root_folder, "GOBP", project_name_overlap)
kegg_output_folder <- file.path(pathway_root_folder, "KEGG", project_name_overlap)
log_folder <- file.path(pathway_root_folder, "log", project_name_overlap)

dir.create(gobp_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(kegg_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

# Resolve a user-supplied --gene_file / --universe_file argument.
#   - If the user omitted it, fall back to the default filename under overlap/<project>/.
#   - If the user passed a bare filename (no '/'), resolve it under overlap/<project>/.
#   - If the user passed a path containing '/' (absolute or relative), use it as-is.
resolve_overlap_file <- function(user_input, default_filename) {
  if (is.na(user_input) || user_input == "") {
    return(file.path(overlap_input_folder, default_filename))
  }
  if (!grepl("/", user_input, fixed = TRUE)) {
    return(file.path(overlap_input_folder, user_input))
  }
  user_input
}

gene_file     <- resolve_overlap_file(args$gene_file,     "overlap_gene_name_pairwise_union.csv")
universe_file <- resolve_overlap_file(args$universe_file, "twas_detected_gene_name.csv")

fdr_tag      <- gsub("\\.", "", format(fdr_cutoff,      nsmall = 2))
gobp_fdr_tag <- gsub("\\.", "", format(gobp_fdr_cutoff, nsmall = 2))
kegg_fdr_tag <- gsub("\\.", "", format(kegg_fdr_cutoff, nsmall = 2))

log_file <- file.path(
  log_folder,
  paste0("process_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

write_log <- function(message) {
  cat(message, file = log_file, append = TRUE)
}

select_quiet <- function(keys, keytype, columns) {
  suppressWarnings(
    suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys = keys,
        keytype = keytype,
        columns = columns
      )
    )
  )
}

cat(
  paste0(
    "========== Pipeline Start ==========\n",
    "Time: ", Sys.time(), "\n",
    "Project name overlap: ", project_name_overlap, "\n",
    "Base dir: ", base_dir, "\n",
    "Gene file: ", gene_file, "\n",
    "Universe file: ", universe_file, "\n",
    "FDR cutoff (default tag): ", fdr_cutoff, "\n",
    "GOBP FDR cutoff: ", gobp_fdr_cutoff, "\n",
    "KEGG FDR cutoff: ", kegg_fdr_cutoff, "\n\n"
  ),
  file = log_file
)

message("Project name overlap: ", project_name_overlap)
message("Base dir: ", base_dir)
message("Gene file: ", gene_file)
message("Universe file: ", universe_file)
message("GOBP output folder: ", gobp_output_folder)
message("KEGG output folder: ", kegg_output_folder)
message("Log file: ", log_file)

if (!file.exists(gene_file)) {
  stop("gene_file not found: ", gene_file, call. = FALSE)
}

if (!file.exists(universe_file)) {
  stop("universe_file not found: ", universe_file, call. = FALSE)
}

gene <- fread(gene_file)
universe_gene <- fread(universe_file)

if (!"gene_name" %in% colnames(gene)) {
  stop("gene_file must contain column: gene_name", call. = FALSE)
}

if (!"gene_name" %in% colnames(universe_gene)) {
  stop("universe_file must contain column: gene_name", call. = FALSE)
}

candidate_symbols <- unique(gene$gene_name[!is.na(gene$gene_name) & gene$gene_name != ""])
universe_symbols <- unique(universe_gene$gene_name[!is.na(universe_gene$gene_name) & universe_gene$gene_name != ""])

if (length(candidate_symbols) == 0) {
  stop("No valid gene_name found in gene_file.", call. = FALSE)
}

if (length(universe_symbols) == 0) {
  stop("No valid gene_name found in universe_file.", call. = FALSE)
}

gene_map <- suppressWarnings(
  bitr(
    candidate_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

universe_map <- suppressWarnings(
  bitr(
    universe_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

if (nrow(gene_map) == 0) {
  stop("No candidate genes could be mapped by bitr().", call. = FALSE)
}

if (nrow(universe_map) == 0) {
  stop("No universe genes could be mapped by bitr().", call. = FALSE)
}

candidate_mapped_n <- length(unique(gene_map$SYMBOL))
universe_mapped_n <- length(unique(universe_map$SYMBOL))
candidate_failed_n <- length(setdiff(candidate_symbols, unique(gene_map$SYMBOL)))
universe_failed_n <- length(setdiff(universe_symbols, unique(universe_map$SYMBOL)))

mapping_text <- paste0(
  "Candidate genes input: ", length(candidate_symbols), "\n",
  "Candidate genes mapped: ", candidate_mapped_n, "\n",
  "Candidate genes failed: ", candidate_failed_n, "\n",
  "Universe genes input: ", length(universe_symbols), "\n",
  "Universe genes mapped: ", universe_mapped_n, "\n",
  "Universe genes failed: ", universe_failed_n, "\n"
)

message(mapping_text)
write_log(mapping_text)

ego_bp <- enrichGO(
  gene = unique(gene_map$SYMBOL),
  universe = unique(universe_map$SYMBOL),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

ego_union_df <- as.data.frame(ego_bp) %>%
  filter(
    Count > 1,
    p.adjust <= gobp_fdr_cutoff
  ) %>%
  arrange(p.adjust)

ego_union_df$pathway_gene_num <- rep(0, nrow(ego_union_df))
if (nrow(ego_union_df) > 0) {
  for (path_index in seq_len(nrow(ego_union_df))) {
    go_id <- ego_union_df$ID[path_index]
    offspring <- c(go_id, unlist(GOBPOFFSPRING[[go_id]]))
    go2gene <- suppressWarnings(select_quiet(
      keys = offspring,
      keytype = "GOALL",
      columns = c("ENTREZID")
    ))
    ego_union_df$pathway_gene_num[path_index] <- length(unique(go2gene$ENTREZID))
  }
}

if (nrow(ego_union_df) > 0) {
  plot_df_go <- ego_union_df %>%
    arrange(desc(RichFactor)) %>%
    head(min(20, nrow(.))) %>%
    mutate(geneID = gsub("/", ", ", geneID))

  p_go <- ggplot(
    plot_df_go,
    aes(x = RichFactor, y = reorder(Description, RichFactor))
  ) +
    geom_point(aes(color = p.adjust), size = 4) +
    geom_text(aes(label = geneID), hjust = -0.1, size = 3) +
    scale_color_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
    labs(
      x = "Rich Factor",
      y = "GO Biological Process",
      title = paste0("GO Biological Process Enrichment of ", project_name_overlap, " Candidate Genes")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    ) +
    coord_cartesian(xlim = c(0, max(plot_df_go$RichFactor) * 1.35))

  ggsave(
    file.path(gobp_output_folder, paste0(project_name_overlap, "_GOBP_ORA_FDR", gobp_fdr_tag, ".png")),
    plot = p_go,
    width = 15,
    height = 8,
    dpi = 600
  )
}

write.csv(
  ego_union_df %>% arrange(desc(RichFactor)),
  file.path(gobp_output_folder, paste0(project_name_overlap, "_GOBP_ORA_FDR", gobp_fdr_tag, "_summary.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_bp) %>% arrange(desc(RichFactor)),
  file.path(gobp_output_folder, paste0(project_name_overlap, "_GOBP_ORA_full.csv")),
  row.names = FALSE
)

ekegg <- enrichKEGG(
  gene = unique(gene_map$ENTREZID),
  universe = unique(universe_map$ENTREZID),
  organism = "hsa",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1
)

ekegg <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

ekegg_union_df <- as.data.frame(ekegg) %>%
  filter(
    Count > 1,
    p.adjust <= kegg_fdr_cutoff
  ) %>%
  arrange(desc(RichFactor)) %>%
  mutate(geneID = gsub("/", ", ", geneID))

ekegg_union_df$pathway_gene_num <- rep(0, nrow(ekegg_union_df))
if (nrow(ekegg_union_df) > 0) {
  for (path_index in seq_len(nrow(ekegg_union_df))) {
    pid <- ekegg_union_df$ID[path_index]
    genes <- KEGGREST::keggLink("hsa", pid)
    ekegg_union_df$pathway_gene_num[path_index] <- length(unique(genes))
  }
}

if (nrow(ekegg_union_df) > 0) {
  plot_df_kegg <- ekegg_union_df %>%
    arrange(desc(RichFactor)) %>%
    head(min(20, nrow(.)))

  p_kegg <- ggplot(
    plot_df_kegg,
    aes(x = RichFactor, y = reorder(Description, RichFactor))
  ) +
    geom_point(aes(color = p.adjust), size = 4) +
    geom_text(aes(label = geneID), hjust = -0.1, size = 3) +
    scale_color_gradient(low = "#d73027", high = "#4575b4", name = "p.adjust") +
    labs(
      x = "Rich Factor",
      y = "KEGG Pathway",
      title = paste0("KEGG Pathway Enrichment of ", project_name_overlap, " Candidate Genes")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10)
    ) +
    coord_cartesian(xlim = c(0, max(plot_df_kegg$RichFactor) * 1.35))

  ggsave(
    file.path(kegg_output_folder, paste0(project_name_overlap, "_KEGG_ORA_FDR", kegg_fdr_tag, ".png")),
    plot = p_kegg,
    width = 15,
    height = 8,
    dpi = 600
  )
}

write.csv(
  ekegg_union_df %>% arrange(desc(RichFactor)),
  file.path(kegg_output_folder, paste0(project_name_overlap, "_KEGG_ORA_FDR", kegg_fdr_tag, "_summary.csv")),
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg) %>% arrange(desc(RichFactor)),
  file.path(kegg_output_folder, paste0(project_name_overlap, "_KEGG_ORA_full.csv")),
  row.names = FALSE
)

summary_text <- paste0(
  "\n========== Summary ==========\n",
  "Candidate gene count: ", length(candidate_symbols), "\n",
  "Universe gene count: ", length(universe_symbols), "\n",
  "Mapped candidate genes: ", candidate_mapped_n, "\n",
  "Mapped universe genes: ", universe_mapped_n, "\n",
  "GO BP rows: ", nrow(ego_union_df), "\n",
  "KEGG rows: ", nrow(ekegg_union_df), "\n",
  "GOBP output folder: ", gobp_output_folder, "\n",
  "KEGG output folder: ", kegg_output_folder, "\n",
  "Pipeline finished at: ", Sys.time(), "\n",
  "=============================\n"
)

cat(summary_text)
write_log(summary_text)
