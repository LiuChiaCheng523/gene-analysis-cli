#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript spredixcan_manhattan_heatmap_cli.R <project_name> [base_dir] [spredixcan_ver] [fdr_cutoff]\n",
      "  Rscript spredixcan_manhattan_heatmap_cli.R --project_name <project_name> [--base_dir <base_dir>] [--spredixcan_ver <v7|v8>] [--fdr_cutoff <value>]\n",
      "    [--gencode_v26_gtf <path>] [--gencode_v19_gtf <path>]\n\n",
      "Arguments:\n",
      "  project_name      Project name, for example: TWB1_LAA_test1\n",
      "  base_dir          Base directory. Default: /mnt/data/ai_agent/gene_analysis\n",
      "  spredixcan_ver    S-PrediXcan version. Default: v8\n",
      "  fdr_cutoff        FDR cutoff. Default: 0.15\n",
      "  gencode_v26_gtf   Default: <base_dir>/tools/GENCODE/gencode.v26.annotation.gtf\n",
      "  gencode_v19_gtf   Default: <base_dir>/tools/GENCODE/gencode.v19.annotation.gtf\n\n",
      "Example:\n",
      "  Rscript spredixcan_manhattan_heatmap_cli.R TWB1_LAA_test1 /mnt/data/ai_agent/gene_analysis v8 0.15\n"
    )
  )
}

parse_args <- function(args) {
  defaults <- list(
    project_name = NA_character_,
    base_dir = "/mnt/data/ai_agent/gene_analysis",
    spredixcan_ver = "v8",
    fdr_cutoff = 0.15,
    gencode_v26_gtf = NA_character_,
    gencode_v19_gtf = NA_character_
  )

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  if (any(grepl("^--", args))) {
    parsed <- defaults
    i <- 1

    while (i <= length(args)) {
      arg <- args[[i]]

      if (arg %in% c("--project_name", "--project-name", "-p")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$project_name <- args[[i + 1]]
        i <- i + 2
      } else if (arg %in% c("--base_dir", "--base-dir", "-b")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$base_dir <- args[[i + 1]]
        i <- i + 2
      } else if (arg %in% c("--spredixcan_ver", "--spredixcan-ver", "-v")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$spredixcan_ver <- args[[i + 1]]
        i <- i + 2
      } else if (arg %in% c("--fdr_cutoff", "--fdr-cutoff", "-f")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$fdr_cutoff <- as.numeric(args[[i + 1]])
        i <- i + 2
      } else if (arg == "--gencode_v26_gtf") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$gencode_v26_gtf <- args[[i + 1]]
        i <- i + 2
      } else if (arg == "--gencode_v19_gtf") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$gencode_v19_gtf <- args[[i + 1]]
        i <- i + 2
      } else {
        stop("Unknown argument: ", arg, call. = FALSE)
      }
    }

    if (is.na(parsed$project_name) || parsed$project_name == "") {
      stop("project_name is required.", call. = FALSE)
    }

    return(parsed)
  }

  parsed <- defaults
  parsed$project_name <- args[[1]]
  if (length(args) >= 2) parsed$base_dir <- args[[2]]
  if (length(args) >= 3) parsed$spredixcan_ver <- args[[3]]
  if (length(args) >= 4) parsed$fdr_cutoff <- as.numeric(args[[4]])
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))

project_name <- args$project_name
base_dir <- args$base_dir
spredixcan_ver <- args$spredixcan_ver
FDR_cutoff <- as.numeric(args$fdr_cutoff)
gencode_v26_gtf <- args$gencode_v26_gtf
gencode_v19_gtf <- args$gencode_v19_gtf

# If not explicitly provided, default GENCODE GTFs to base_dir/tools/GENCODE/
if (is.na(gencode_v26_gtf) || gencode_v26_gtf == "") {
  gencode_v26_gtf <- file.path(base_dir, "tools/GENCODE/gencode.v26.annotation.gtf")
}
if (is.na(gencode_v19_gtf) || gencode_v19_gtf == "") {
  gencode_v19_gtf <- file.path(base_dir, "tools/GENCODE/gencode.v19.annotation.gtf")
}

if (!spredixcan_ver %in% c("v7", "v8")) {
  stop("spredixcan_ver must be either 'v7' or 'v8'.", call. = FALSE)
}

if (is.na(FDR_cutoff)) {
  stop("fdr_cutoff must be numeric.", call. = FALSE)
}

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

spredixcan_result_v8_folder <- file.path(base_dir, "S_PrediXcan/result_v8", project_name)
spredixcan_result_v7_folder <- file.path(base_dir, "S_PrediXcan/result_v7", project_name)

spredixcan_result_v8_table_output_folder <- file.path(spredixcan_result_v8_folder, "table")
spredixcan_result_v7_table_output_folder <- file.path(spredixcan_result_v7_folder, "table")
spredixcan_result_v8_manhanttan_output_folder <- file.path(spredixcan_result_v8_folder, "manhanttan")
spredixcan_result_v7_manhanttan_output_folder <- file.path(spredixcan_result_v7_folder, "manhanttan")
spredixcan_result_v8_heatmap_output_folder <- file.path(spredixcan_result_v8_folder, "heatmap")
spredixcan_result_v7_heatmap_output_folder <- file.path(spredixcan_result_v7_folder, "heatmap")

if (spredixcan_ver == "v8" && !dir.exists(spredixcan_result_v8_folder)) {
  stop("S-PrediXcan v8 result folder not found: ", spredixcan_result_v8_folder, call. = FALSE)
}

if (spredixcan_ver == "v7" && !dir.exists(spredixcan_result_v7_folder)) {
  stop("S-PrediXcan v7 result folder not found: ", spredixcan_result_v7_folder, call. = FALSE)
}

message("Project name: ", project_name)
message("Base dir: ", base_dir)
message("S-PrediXcan version: ", spredixcan_ver)
message("FDR cutoff: ", FDR_cutoff)

if (!file.exists(gencode_v26_gtf)) {
  stop("gencode_v26_gtf not found: ", gencode_v26_gtf, call. = FALSE)
}

if (!file.exists(gencode_v19_gtf)) {
  stop("gencode_v19_gtf not found: ", gencode_v19_gtf, call. = FALSE)
}

gtf <- import(gencode_v26_gtf)
gene_gr <- gtf[gtf$type == "gene"]

gene_dt <- data.table(
  ensembl_gene_id = gene_gr$gene_id,
  gene_name = gene_gr$gene_name,
  chr = as.character(seqnames(gene_gr)),
  start = start(gene_gr),
  end = end(gene_gr),
  gene_type = gene_gr$gene_type
)
gene_dt[, chr := sub("^chr", "", chr)]
gene_dt <- gene_dt[chr %in% c(as.character(1:22))]

gtf_v19 <- import(gencode_v19_gtf)
gene_gr_v19 <- gtf_v19[gtf_v19$type == "gene"]
gene_dt_v19 <- data.table(
  ensembl_gene_id = gene_gr_v19$gene_id,
  gene_name = gene_gr_v19$gene_name,
  chr = as.character(seqnames(gene_gr_v19)),
  start = start(gene_gr_v19),
  end = end(gene_gr_v19),
  gene_type = gene_gr_v19$gene_type
)
gene_dt_v19[, chr := sub("^chr", "", chr)]
gene_dt_v19 <- gene_dt_v19[chr %in% c(as.character(1:22))]

if (spredixcan_ver == "v8") {
  dir.create(spredixcan_result_v8_table_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(spredixcan_result_v8_manhanttan_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(spredixcan_result_v8_heatmap_output_folder, recursive = TRUE, showWarnings = FALSE)

  spredixcan_v8_tissue_file <- file.path(base_dir, "S_PrediXcan/tissues_v8_elastic_net.txt")
  if (!file.exists(spredixcan_v8_tissue_file)) {
    stop("S-PrediXcan v8 tissue list not found: ", spredixcan_v8_tissue_file, call. = FALSE)
  }

  SP_v8_tissue <- fread(spredixcan_v8_tissue_file, header = FALSE)
  SP_v8_tissue_list <- SP_v8_tissue$V1

  SP_v8 <- data.frame()
  for (tissue in SP_v8_tissue_list) {
    twas_file <- paste0(spredixcan_result_v8_folder, "/", tissue, "_SPrediXcan_v8_en.csv")
    if (!file.exists(twas_file)) {
      warning("File not found: ", twas_file)
      next
    }
    twas <- fread(twas_file)
    twas$tissue <- tissue
    SP_v8 <- rbind(SP_v8, twas)
  }

  if (nrow(SP_v8) == 0) {
    stop("No S-PrediXcan v8 result found.", call. = FALSE)
  }

  SP_v8 <- SP_v8[order(SP_v8$pvalue), ]
  setnames(SP_v8, "gene", "ensembl_gene_id")

  SP_v8_annot <- merge(
    SP_v8,
    gene_dt,
    by.x = c("ensembl_gene_id", "gene_name"),
    by.y = c("ensembl_gene_id", "gene_name")
  )

  SP_v8_annot <- SP_v8_annot[order(SP_v8_annot$pvalue), ]
  setnames(SP_v8_annot, "chr", "CHR")
  SP_v8_annot[, BP := start]

  SP_v8_manhattan <- SP_v8_annot %>%
    filter(CHR %in% 1:22) %>%
    mutate(P = pvalue)

  SP_v8_manhattan <- SP_v8_manhattan %>% filter(!is.na(P))
  SP_v8_manhattan <- as.data.table(SP_v8_manhattan)
  SP_v8_manhattan[, FDR := p.adjust(P, method = "BH")]
  SP_v8_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]

  write.csv(
    SP_v8_manhattan,
    paste0(spredixcan_result_v8_table_output_folder, "/SPrediXcan_v8_TWAS.csv"),
    row.names = FALSE
  )

  write.csv(
    SP_v8_manhattan[SP_v8_manhattan$FDR_tissue < FDR_cutoff, ],
    paste0(spredixcan_result_v8_table_output_folder, "/SPrediXcan_v8_TWAS_FDR", FDR_cutoff, ".csv"),
    row.names = FALSE
  )

  SP_v8_manhattan_clean <- SP_v8_manhattan %>%
    dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)

  SP_v8_manhattan_clean <- SP_v8_manhattan_clean %>%
    filter(!is.na(CHR), !is.na(BP))

  SP_v8_plot_df <- SP_v8_manhattan_clean %>%
    mutate(
      CHR = as.numeric(CHR),
      BP = as.numeric(BP)
    ) %>%
    arrange(CHR, BP)

  SP_v8_chr_max <- SP_v8_plot_df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

  SP_v8_plot_df <- SP_v8_plot_df %>%
    left_join(SP_v8_chr_max, by = "CHR") %>%
    mutate(BP_cum = BP + offset)

  SP_v8_chr_centers <- SP_v8_plot_df %>%
    group_by(CHR) %>%
    summarise(center = mean(range(BP_cum)), .groups = "drop")

  SP_v8_gene_fdr <- SP_v8_plot_df[SP_v8_plot_df$FDR_tissue < FDR_cutoff, ]
  SP_v8_label_genes <- unique(SP_v8_gene_fdr$gene_name)
  SP_v8_label_df <- SP_v8_plot_df %>%
    filter(
      gene_name %in% SP_v8_label_genes,
      FDR_tissue < FDR_cutoff
    ) %>%
    group_by(gene_name) %>%
    slice_min(P, n = 1) %>%
    ungroup()

  SP_v8_plot_df <- SP_v8_plot_df %>%
    mutate(
      chr_color = factor(CHR %% 2),
      logP = -log10(pmax(P, .Machine$double.xmin))
    )

  SP_v8_label_df <- SP_v8_label_df %>%
    mutate(logP = -log10(pmax(P, .Machine$double.xmin)))

  SP_v8_y_max <- ceiling(max(SP_v8_plot_df$logP, na.rm = TRUE)) + 1

  SP_v8_p <- ggplot(SP_v8_plot_df, aes(x = BP_cum, y = logP)) +
    geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
    scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
    geom_point(
      data = SP_v8_label_df,
      aes(x = BP_cum, y = logP),
      color = "red",
      size = 2
    ) +
    geom_text_repel(
      data = SP_v8_label_df,
      aes(label = gene_name),
      color = "red",
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    labs(
      title = "S-PrediXcan GTEx_v8",
      x = "Chromosome",
      y = expression(-log[10](P))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    scale_x_continuous(
      breaks = SP_v8_chr_centers$center,
      labels = SP_v8_chr_centers$CHR
    ) +
    scale_y_continuous(
      limits = c(0, SP_v8_y_max),
      breaks = seq(0, SP_v8_y_max, by = 1)
    )

  ggsave(
    paste0(spredixcan_result_v8_manhanttan_output_folder, "/SPrediXcan_GTEx_v8_manhanttan_FDR", FDR_cutoff, ".png"),
    SP_v8_p,
    width = 12,
    height = 6,
    dpi = 600
  )

  SP_v8_fdr <- SP_v8_manhattan %>%
    filter(FDR_tissue < FDR_cutoff)

  heat_df <- SP_v8_fdr %>%
    mutate(sig_value = -log10(FDR_tissue))

  heat_mat <- heat_df %>%
    dplyr::select(gene_name, tissue, sig_value) %>%
    pivot_wider(
      names_from = tissue,
      values_from = sig_value,
      values_fill = 0
    ) %>%
    column_to_rownames("gene_name") %>%
    as.matrix()

  sig_cutoff <- -log10(FDR_cutoff)
  gene_count <- colSums(heat_mat > sig_cutoff)
  tissue_count <- rowSums(heat_mat > sig_cutoff)

  tissue_order <- order(colnames(heat_mat))
  heat_mat <- heat_mat[, tissue_order, drop = FALSE]
  gene_count <- gene_count[tissue_order]

  col_fun <- colorRamp2(
    c(0, median(heat_mat[heat_mat > 0]), max(heat_mat)),
    c("grey", "#6baed6", "#d73027")
  )

  ht <- Heatmap(
    heat_mat,
    name = "Significance",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 9),
    column_names_rot = 60,
    top_annotation = HeatmapAnnotation(
      GeneCount = anno_barplot(
        gene_count,
        gp = grid::gpar(fill = "grey40"),
        height = grid::unit(2, "cm")
      )
    ),
    left_annotation = rowAnnotation(
      TissueCount = anno_barplot(
        tissue_count,
        gp = grid::gpar(fill = "grey40"),
        width = grid::unit(2, "cm")
      )
    )
  )

  png(
    filename = paste0(spredixcan_result_v8_heatmap_output_folder, "/SPrediXcan_GTEx_v8_heatmap_FDR", FDR_cutoff, ".png"),
    width = 4200,
    height = 3000,
    res = 300
  )

  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )

  dev.off()

  message("S-PrediXcan v8 processing complete.")
} else if (spredixcan_ver == "v7") {
  dir.create(spredixcan_result_v7_table_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(spredixcan_result_v7_manhanttan_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(spredixcan_result_v7_heatmap_output_folder, recursive = TRUE, showWarnings = FALSE)

  spredixcan_v7_tissue_file <- file.path(base_dir, "S_PrediXcan/tissues.txt")
  if (!file.exists(spredixcan_v7_tissue_file)) {
    stop("S-PrediXcan v7 tissue list not found: ", spredixcan_v7_tissue_file, call. = FALSE)
  }

  SP_v7_tissue <- fread(spredixcan_v7_tissue_file, header = FALSE)
  SP_v7_tissue_list <- SP_v7_tissue$V1

  SP_v7 <- data.frame()
  for (tissue in SP_v7_tissue_list) {
    twas_file <- paste0(spredixcan_result_v7_folder, "/", tissue, "_SPrediXcan_v7_en.csv")
    if (!file.exists(twas_file)) {
      warning("File not found: ", twas_file)
      next
    }
    twas <- fread(twas_file)
    twas$tissue <- tissue
    SP_v7 <- rbind(SP_v7, twas)
  }

  if (nrow(SP_v7) == 0) {
    stop("No S-PrediXcan v7 result found.", call. = FALSE)
  }

  SP_v7 <- SP_v7[order(SP_v7$pvalue), ]
  setnames(SP_v7, "gene", "ensembl_gene_id")

  SP_v7_annot <- merge(
    SP_v7,
    gene_dt_v19,
    by.x = c("ensembl_gene_id", "gene_name"),
    by.y = c("ensembl_gene_id", "gene_name")
  )

  SP_v7_annot <- SP_v7_annot[order(SP_v7_annot$pvalue), ]
  setnames(SP_v7_annot, "chr", "CHR")
  SP_v7_annot[, BP := start]

  SP_v7_manhattan <- SP_v7_annot %>%
    filter(CHR %in% 1:22) %>%
    mutate(P = pvalue)

  SP_v7_manhattan <- SP_v7_manhattan %>% filter(!is.na(P))
  SP_v7_manhattan <- as.data.table(SP_v7_manhattan)
  SP_v7_manhattan[, FDR := p.adjust(P, method = "BH")]
  SP_v7_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]

  write.csv(
    SP_v7_manhattan,
    paste0(spredixcan_result_v7_table_output_folder, "/SPrediXcan_v7_TWAS.csv"),
    row.names = FALSE
  )

  write.csv(
    SP_v7_manhattan[SP_v7_manhattan$FDR_tissue < FDR_cutoff, ],
    paste0(spredixcan_result_v7_table_output_folder, "/SPrediXcan_v7_TWAS_FDR", FDR_cutoff, ".csv"),
    row.names = FALSE
  )

  SP_v7_manhattan_clean <- SP_v7_manhattan %>%
    dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)

  SP_v7_manhattan_clean <- SP_v7_manhattan_clean %>%
    filter(!is.na(CHR), !is.na(BP))

  SP_v7_plot_df <- SP_v7_manhattan_clean %>%
    mutate(
      CHR = as.numeric(CHR),
      BP = as.numeric(BP)
    ) %>%
    arrange(CHR, BP)

  SP_v7_chr_max <- SP_v7_plot_df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

  SP_v7_plot_df <- SP_v7_plot_df %>%
    left_join(SP_v7_chr_max, by = "CHR") %>%
    mutate(BP_cum = BP + offset)

  SP_v7_chr_centers <- SP_v7_plot_df %>%
    group_by(CHR) %>%
    summarise(center = mean(range(BP_cum)), .groups = "drop")

  SP_v7_gene_fdr <- SP_v7_plot_df[SP_v7_plot_df$FDR_tissue < FDR_cutoff, ]
  SP_v7_label_genes <- unique(SP_v7_gene_fdr$gene_name)
  SP_v7_label_df <- SP_v7_plot_df %>%
    filter(
      gene_name %in% SP_v7_label_genes,
      FDR_tissue < FDR_cutoff
    ) %>%
    group_by(gene_name) %>%
    slice_min(P, n = 1) %>%
    ungroup()

  SP_v7_plot_df <- SP_v7_plot_df %>%
    mutate(
      chr_color = factor(CHR %% 2),
      logP = -log10(pmax(P, .Machine$double.xmin))
    )

  SP_v7_label_df <- SP_v7_label_df %>%
    mutate(logP = -log10(pmax(P, .Machine$double.xmin)))

  SP_v7_y_max <- ceiling(max(SP_v7_plot_df$logP, na.rm = TRUE)) + 1

  SP_v7_p <- ggplot(SP_v7_plot_df, aes(x = BP_cum, y = logP)) +
    geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
    scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
    geom_point(
      data = SP_v7_label_df,
      aes(x = BP_cum, y = logP),
      color = "red",
      size = 2
    ) +
    geom_text_repel(
      data = SP_v7_label_df,
      aes(label = gene_name),
      color = "red",
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    labs(
      title = "S-PrediXcan GTEx_v7",
      x = "Chromosome",
      y = expression(-log[10](P))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    ) +
    scale_x_continuous(
      breaks = SP_v7_chr_centers$center,
      labels = SP_v7_chr_centers$CHR
    ) +
    scale_y_continuous(
      limits = c(0, SP_v7_y_max),
      breaks = seq(0, SP_v7_y_max, by = 1)
    )

  ggsave(
    paste0(spredixcan_result_v7_manhanttan_output_folder, "/SPrediXcan_GTEx_v7_manhanttan_FDR", FDR_cutoff, ".png"),
    SP_v7_p,
    width = 12,
    height = 6,
    dpi = 600
  )

  SP_v7_fdr <- SP_v7_manhattan %>%
    filter(FDR_tissue < FDR_cutoff)

  heat_df <- SP_v7_fdr %>%
    mutate(sig_value = -log10(FDR_tissue))

  heat_mat <- heat_df %>%
    dplyr::select(gene_name, tissue, sig_value) %>%
    pivot_wider(
      names_from = tissue,
      values_from = sig_value,
      values_fill = 0
    ) %>%
    column_to_rownames("gene_name") %>%
    as.matrix()

  sig_cutoff <- -log10(FDR_cutoff)
  gene_count <- colSums(heat_mat > sig_cutoff)
  tissue_count <- rowSums(heat_mat > sig_cutoff)

  tissue_order <- order(colnames(heat_mat))
  heat_mat <- heat_mat[, tissue_order, drop = FALSE]
  gene_count <- gene_count[tissue_order]

  col_fun <- colorRamp2(
    c(0, median(heat_mat[heat_mat > 0]), max(heat_mat)),
    c("grey", "#6baed6", "#d73027")
  )

  ht <- Heatmap(
    heat_mat,
    name = "Significance",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 9),
    column_names_rot = 60,
    top_annotation = HeatmapAnnotation(
      GeneCount = anno_barplot(
        gene_count,
        gp = grid::gpar(fill = "grey40"),
        height = grid::unit(2, "cm")
      )
    ),
    left_annotation = rowAnnotation(
      TissueCount = anno_barplot(
        tissue_count,
        gp = grid::gpar(fill = "grey40"),
        width = grid::unit(2, "cm")
      )
    )
  )

  png(
    filename = paste0(spredixcan_result_v7_heatmap_output_folder, "/SPrediXcan_GTEx_v7_heatmap_FDR", FDR_cutoff, ".png"),
    width = 4200,
    height = 3000,
    res = 300
  )

  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )

  dev.off()

  message("S-PrediXcan v7 processing complete.")
}
