#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript fusion_manhattan_heatmap_cli.R <project_name> [base_dir] [fusion_ver] [fdr_cutoff]\n",
      "  Rscript fusion_manhattan_heatmap_cli.R --project_name <project_name> [--base_dir <base_dir>] [--fusion_ver <v7|v8>] [--fdr_cutoff <value>]\n",
      "    [--gencode_v26_gtf <path>] [--gencode_v19_gtf <path>] [--weight_dir_v7 <path>]\n\n",
      "Arguments:\n",
      "  project_name      Project name, for example: TWB1_LAA_test1\n",
      "  base_dir          Base directory. Default: /mnt/data/ai_agent/gene_analysis\n",
      "  fusion_ver        FUSION version. Default: v8\n",
      "  fdr_cutoff        FDR cutoff. Default: 0.15\n",
      "  gencode_v26_gtf   Default: <base_dir>/tools/GENCODE/gencode.v26.annotation.gtf\n",
      "  gencode_v19_gtf   Default: <base_dir>/tools/GENCODE/gencode.v19.annotation.gtf\n",
      "  weight_dir_v7     Default: <base_dir>/FUSION/WEIGHTS_v7/GTEx.ALL\n\n",
      "Example:\n",
      "  Rscript fusion_manhattan_heatmap_cli.R TWB1_LAA_test1 /mnt/data/ai_agent/gene_analysis v8 0.15\n"
    )
  )
}

parse_args <- function(args) {
  defaults <- list(
    project_name = NA_character_,
    base_dir = "/mnt/data/ai_agent/gene_analysis",
    fusion_ver = "v8",
    fdr_cutoff = 0.15,
    gencode_v26_gtf = NA_character_,
    gencode_v19_gtf = NA_character_,
    weight_dir_v7 = NA_character_
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
      } else if (arg %in% c("--fusion_ver", "--fusion-ver", "-v")) {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$fusion_ver <- args[[i + 1]]
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
      } else if (arg == "--weight_dir_v7") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$weight_dir_v7 <- args[[i + 1]]
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
  if (length(args) >= 3) parsed$fusion_ver <- args[[3]]
  if (length(args) >= 4) parsed$fdr_cutoff <- as.numeric(args[[4]])
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))

project_name <- args$project_name
base_dir <- args$base_dir
fusion_ver <- args$fusion_ver
FDR_cutoff <- as.numeric(args$fdr_cutoff)
gencode_v26_gtf <- args$gencode_v26_gtf
gencode_v19_gtf <- args$gencode_v19_gtf
weight_dir <- args$weight_dir_v7

# If not explicitly provided, default paths relative to base_dir
if (is.na(gencode_v26_gtf) || gencode_v26_gtf == "") {
  gencode_v26_gtf <- file.path(base_dir, "tools/GENCODE/gencode.v26.annotation.gtf")
}
if (is.na(gencode_v19_gtf) || gencode_v19_gtf == "") {
  gencode_v19_gtf <- file.path(base_dir, "tools/GENCODE/gencode.v19.annotation.gtf")
}
if (is.na(weight_dir) || weight_dir == "") {
  weight_dir <- file.path(base_dir, "FUSION/WEIGHTS_v7/GTEx.ALL")
}

if (!fusion_ver %in% c("v7", "v8")) {
  stop("fusion_ver must be either 'v7' or 'v8'.", call. = FALSE)
}

if (is.na(FDR_cutoff)) {
  stop("fdr_cutoff must be numeric.", call. = FALSE)
}

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

fusion_result_v8_folder <- file.path(base_dir, "FUSION/result_v8", project_name)
fusion_result_v7_folder <- file.path(base_dir, "FUSION/result_v7", project_name)

fusion_result_v8_table_output_folder <- file.path(fusion_result_v8_folder, "table")
fusion_result_v7_table_output_folder <- file.path(fusion_result_v7_folder, "table")
fusion_result_v8_manhanttan_output_folder <- file.path(fusion_result_v8_folder, "manhanttan")
fusion_result_v7_manhanttan_output_folder <- file.path(fusion_result_v7_folder, "manhanttan")
fusion_result_v8_heatmap_output_folder <- file.path(fusion_result_v8_folder, "heatmap")
fusion_result_v7_heatmap_output_folder <- file.path(fusion_result_v7_folder, "heatmap")

message("Project name: ", project_name)
message("Base dir: ", base_dir)
message("FUSION version: ", fusion_ver)
message("FDR cutoff: ", FDR_cutoff)

if (!file.exists(gencode_v26_gtf)) {
  stop("gencode_v26_gtf not found: ", gencode_v26_gtf, call. = FALSE)
}

if (!file.exists(gencode_v19_gtf)) {
  stop("gencode_v19_gtf not found: ", gencode_v19_gtf, call. = FALSE)
}

if (!dir.exists(weight_dir)) {
  stop("weight_dir_v7 not found: ", weight_dir, call. = FALSE)
}

if (fusion_ver == "v8" && !dir.exists(fusion_result_v8_folder)) {
  stop("FUSION v8 result folder not found: ", fusion_result_v8_folder, call. = FALSE)
}

if (fusion_ver == "v7" && !dir.exists(fusion_result_v7_folder)) {
  stop("FUSION v7 result folder not found: ", fusion_result_v7_folder, call. = FALSE)
}

# load gencode database -----
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

# load FUSION v7 wgt filename
wgt_files <- list.files(
  weight_dir,
  pattern = "\\.wgt\\.RDat$",
  full.names = TRUE,
  recursive = TRUE
)

gene_map_wgt <- data.table(filename = basename(wgt_files))
gene_map_wgt[, ensembl_gene_id := str_extract(filename, "ENSG[0-9]+\\.[0-9]+")]
gene_map_wgt[, gene_id := sub("\\..*$", "", ensembl_gene_id)]
gene_map_wgt[, gene_name := sub(".*ENSG[0-9]+\\.[0-9]+\\.", "", filename)]
gene_map_wgt[, gene_name := sub("\\.wgt\\.RDat$", "", gene_name)]
gene_map_wgt <- unique(gene_map_wgt[, .(gene_name, ensembl_gene_id, gene_id)])

if (fusion_ver == "v8") {
  
  dir.create(fusion_result_v8_table_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(fusion_result_v8_manhanttan_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(fusion_result_v8_heatmap_output_folder, recursive = TRUE, showWarnings = FALSE)
  
  fusion_v8_tissue_file <- file.path(base_dir, "FUSION/WEIGHTS/tissue_list.txt")
  if (!file.exists(fusion_v8_tissue_file)) {
    stop("FUSION v8 tissue list not found: ", fusion_v8_tissue_file, call. = FALSE)
  }

  FU_v8_tissue <- fread(fusion_v8_tissue_file, header = FALSE)
  FU_v8_tissue_list <- FU_v8_tissue$V1

  FU_v8 <- data.frame()
  for (num in 1:22) {
    for (tissue in FU_v8_tissue_list) {
      twas_file <- paste0(fusion_result_v8_folder, "/chr", num, "/", tissue, "_FUSION_v8.dat")
      if (!file.exists(twas_file)) {
        warning("File not found: ", twas_file)
        next
      }

      twas <- fread(twas_file)
      twas$tissue <- str_remove(twas$PANEL, "^GTExv8.EUR\\.")
      twas <- twas[, -c(1, 2)]
      FU_v8 <- rbind(FU_v8, twas)
    }
  }

  if (nrow(FU_v8) == 0) {
    stop("No FUSION v8 result found.", call. = FALSE)
  }

  setnames(FU_v8, "ID", "ensembl_gene_id")
  FU_v8[, CHR := as.character(CHR)]

  FU_v8_annot <- merge(
    FU_v8,
    gene_dt,
    by.x = c("ensembl_gene_id", "CHR"),
    by.y = c("ensembl_gene_id", "chr")
  )

  FU_v8_annot_dup_symbol <- FU_v8_annot[, uniqueN(ensembl_gene_id), by = gene_name][V1 > 1]
  FU_v8_annot <- FU_v8_annot[order(FU_v8_annot$TWAS.P), ]
  FU_v8_annot[, BP := start]

  FU_v8_manhattan <- FU_v8_annot %>%
    filter(CHR %in% 1:22) %>%
    mutate(P = TWAS.P)

  FU_v8_manhattan <- FU_v8_manhattan %>% filter(!is.na(P))
  FU_v8_manhattan <- as.data.table(FU_v8_manhattan)
  FU_v8_manhattan[, FDR := p.adjust(P, method = "BH")]
  FU_v8_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]

  write.csv(
    FU_v8_manhattan,
    paste0(fusion_result_v8_table_output_folder, "/FUSION_v8_TWAS.csv"),
    row.names = FALSE
  )

  write.csv(
    FU_v8_manhattan[FU_v8_manhattan$FDR_tissue < FDR_cutoff, ],
    paste0(fusion_result_v8_table_output_folder, "/FUSION_v8_TWAS_FDR", FDR_cutoff, ".csv"),
    row.names = FALSE
  )

  FU_v8_manhattan_clean <- FU_v8_manhattan %>%
    dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)

  FU_v8_manhattan_clean <- FU_v8_manhattan_clean %>%
    filter(
      !is.na(CHR),
      !is.na(BP)
    )

  FU_v8_plot_df <- FU_v8_manhattan_clean %>%
    mutate(
      CHR = as.numeric(CHR),
      BP = as.numeric(BP)
    ) %>%
    arrange(CHR, BP)

  FU_v8_chr_max <- FU_v8_plot_df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

  FU_v8_plot_df <- FU_v8_plot_df %>%
    left_join(FU_v8_chr_max, by = "CHR") %>%
    mutate(BP_cum = BP + offset)

  FU_v8_chr_centers <- FU_v8_plot_df %>%
    group_by(CHR) %>%
    summarise(center = mean(range(BP_cum)), .groups = "drop")

  FU_v8_gene_fdr <- FU_v8_plot_df[FU_v8_plot_df$FDR_tissue < FDR_cutoff, ]
  FU_v8_label_genes <- unique(FU_v8_gene_fdr$gene_name)
  FU_v8_label_df <- FU_v8_plot_df %>%
    filter(
      gene_name %in% FU_v8_label_genes,
      FDR_tissue < FDR_cutoff
    ) %>%
    group_by(gene_name) %>%
    slice_min(P, n = 1) %>%
    ungroup()

  FU_v8_plot_df <- FU_v8_plot_df %>%
    mutate(
      chr_color = factor(CHR %% 2),
      logP = -log10(pmax(P, .Machine$double.xmin))
    )

  FU_v8_label_df <- FU_v8_label_df %>%
    mutate(logP = -log10(pmax(P, .Machine$double.xmin)))

  FU_v8_y_max <- ceiling(max(FU_v8_plot_df$logP, na.rm = TRUE)) + 1

  FU_v8_p <- ggplot(FU_v8_plot_df, aes(x = BP_cum, y = logP)) +
    geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
    scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
    geom_point(
      data = FU_v8_label_df,
      aes(x = BP_cum, y = logP),
      color = "red",
      size = 2
    ) +
    geom_text_repel(
      data = FU_v8_label_df,
      aes(label = gene_name),
      color = "red",
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    labs(
      title = "FUSION GTEx_v8",
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
      breaks = FU_v8_chr_centers$center,
      labels = FU_v8_chr_centers$CHR
    ) +
    scale_y_continuous(
      limits = c(0, FU_v8_y_max),
      breaks = seq(0, FU_v8_y_max, by = 1)
    )

  ggsave(
    paste0(fusion_result_v8_manhanttan_output_folder, "/FUSION_GTEx_v8_manhanttan_FDR", FDR_cutoff, ".png"),
    FU_v8_p,
    width = 12,
    height = 6,
    dpi = 600
  )

  FU_v8_fdr <- FU_v8_manhattan %>%
    filter(FDR_tissue < FDR_cutoff)

  heat_df <- FU_v8_fdr %>%
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
    filename = paste0(fusion_result_v8_heatmap_output_folder, "/FUSION_GTEx_v8_heatmap_FDR", FDR_cutoff, ".png"),
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

  message("FUSION v8 processing complete.")
} else if (fusion_ver == "v7") {
  
  dir.create(fusion_result_v7_table_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(fusion_result_v7_manhanttan_output_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(fusion_result_v7_heatmap_output_folder, recursive = TRUE, showWarnings = FALSE)
  
  fusion_v7_tissue_file <- file.path(base_dir, "FUSION/WEIGHTS_v7/GTEx.ALL/tissue_list_v7.txt")
  if (!file.exists(fusion_v7_tissue_file)) {
    stop("FUSION v7 tissue list not found: ", fusion_v7_tissue_file, call. = FALSE)
  }

  FU_v7_tissue <- fread(fusion_v7_tissue_file, header = FALSE)
  FU_v7_tissue_list <- FU_v7_tissue$V1

  FU_v7 <- data.frame()
  for (num in 1:22) {
    for (tissue in FU_v7_tissue_list) {
      twas_file <- paste0(fusion_result_v7_folder, "/chr", num, "/", tissue, "_FUSION_v7.dat")
      if (!file.exists(twas_file)) {
        warning("File not found: ", twas_file)
        next
      }

      twas <- fread(twas_file)
      twas$tissue <- tissue
      twas <- twas[, -c(1, 2)]
      FU_v7 <- rbind(FU_v7, twas)
    }
  }

  if (nrow(FU_v7) == 0) {
    stop("No FUSION v7 result found.", call. = FALSE)
  }

  FU_v7 <- FU_v7[order(FU_v7$TWAS.P), ]

  setnames(FU_v7, "ID", "gene_name")
  FU_v7[, CHR := as.character(CHR)]

  FU_v7_annot <- merge(
    FU_v7,
    gene_map_wgt[, .(gene_name, ensembl_gene_id)],
    by = "gene_name",
    all.x = TRUE
  )

  FU_v7_annot <- merge(
    FU_v7_annot,
    gene_dt_v19,
    by.x = c("ensembl_gene_id", "gene_name", "CHR"),
    by.y = c("ensembl_gene_id", "gene_name", "chr")
  )

  FU_v7_annot <- FU_v7_annot[
    P0 >= start &
      P1 <= end
  ]

  FU_v7_dup_groups <- FU_v7_annot[, .N, by = .(gene_name, CHR, P0, P1, tissue)][N > 1]

  FU_v7_annot_clean <- FU_v7_annot[
    !FU_v7_annot[, .I[.N > 1], by = .(gene_name, CHR, P0, P1, tissue)]$V1
  ]

  FU_v7_annot_dup_symbol <- FU_v7_annot_clean[, uniqueN(ensembl_gene_id), by = gene_name][V1 > 1]
  FU_v7_annot_clean <- FU_v7_annot_clean[order(FU_v7_annot_clean$TWAS.P), ]
  FU_v7_annot_clean[, BP := start]

  FU_v7_manhattan <- FU_v7_annot_clean %>%
    filter(CHR %in% 1:22) %>%
    mutate(P = TWAS.P)

  FU_v7_manhattan <- FU_v7_manhattan %>% filter(!is.na(P))
  FU_v7_manhattan <- as.data.table(FU_v7_manhattan)
  FU_v7_manhattan[, FDR := p.adjust(P, method = "BH")]
  FU_v7_manhattan[, FDR_tissue := p.adjust(P, method = "BH"), by = tissue]

  write.csv(
    FU_v7_manhattan,
    paste0(fusion_result_v7_table_output_folder, "/FUSION_v7_TWAS.csv"),
    row.names = FALSE
  )

  write.csv(
    FU_v7_manhattan[FU_v7_manhattan$FDR_tissue < FDR_cutoff, ],
    paste0(fusion_result_v7_table_output_folder, "/FUSION_v7_TWAS_FDR", FDR_cutoff, ".csv"),
    row.names = FALSE
  )

  FU_v7_manhattan_clean <- FU_v7_manhattan %>%
    dplyr::select(CHR, BP, P, FDR_tissue, gene_name, tissue)

  FU_v7_manhattan_clean <- FU_v7_manhattan_clean %>%
    filter(
      !is.na(CHR),
      !is.na(BP)
    )

  FU_v7_plot_df <- FU_v7_manhattan_clean %>%
    mutate(
      CHR = as.numeric(CHR),
      BP = as.numeric(BP)
    ) %>%
    arrange(CHR, BP)

  FU_v7_chr_max <- FU_v7_plot_df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(offset = lag(cumsum(chr_len), default = 0))

  FU_v7_plot_df <- FU_v7_plot_df %>%
    left_join(FU_v7_chr_max, by = "CHR") %>%
    mutate(BP_cum = BP + offset)

  FU_v7_chr_centers <- FU_v7_plot_df %>%
    group_by(CHR) %>%
    summarise(center = mean(range(BP_cum)), .groups = "drop")

  FU_v7_gene_fdr <- FU_v7_plot_df[FU_v7_plot_df$FDR_tissue < FDR_cutoff, ]
  FU_v7_label_genes <- unique(FU_v7_gene_fdr$gene_name)
  FU_v7_label_df <- FU_v7_plot_df %>%
    filter(
      gene_name %in% FU_v7_label_genes,
      FDR_tissue < FDR_cutoff
    ) %>%
    group_by(gene_name) %>%
    slice_min(P, n = 1) %>%
    ungroup()

  FU_v7_plot_df <- FU_v7_plot_df %>%
    mutate(
      chr_color = factor(CHR %% 2),
      logP = -log10(pmax(P, .Machine$double.xmin))
    )

  FU_v7_label_df <- FU_v7_label_df %>%
    mutate(logP = -log10(pmax(P, .Machine$double.xmin)))

  FU_v7_y_max <- ceiling(max(FU_v7_plot_df$logP, na.rm = TRUE)) + 1

  FU_v7_p <- ggplot(FU_v7_plot_df, aes(x = BP_cum, y = logP)) +
    geom_point(aes(color = chr_color), size = 0.6, alpha = 0.8) +
    scale_color_manual(values = c("#2C7BB6", "#F28E2B"), guide = "none") +
    geom_point(
      data = FU_v7_label_df,
      aes(x = BP_cum, y = logP),
      color = "red",
      size = 2
    ) +
    geom_text_repel(
      data = FU_v7_label_df,
      aes(label = gene_name),
      color = "red",
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      max.overlaps = Inf
    ) +
    labs(
      title = "FUSION GTEx_v7",
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
      breaks = FU_v7_chr_centers$center,
      labels = FU_v7_chr_centers$CHR
    ) +
    scale_y_continuous(
      limits = c(0, FU_v7_y_max),
      breaks = seq(0, FU_v7_y_max, by = 1)
    )

  ggsave(
    paste0(fusion_result_v7_manhanttan_output_folder, "/FUSION_GTEx_v7_manhanttan_FDR", FDR_cutoff, ".png"),
    FU_v7_p,
    width = 12,
    height = 6,
    dpi = 600
  )

  FU_v7_fdr <- FU_v7_manhattan %>%
    filter(FDR_tissue < FDR_cutoff)

  heat_df <- FU_v7_fdr %>%
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
    filename = paste0(fusion_result_v7_heatmap_output_folder, "/FUSION_GTEx_v7_heatmap_FDR", FDR_cutoff, ".png"),
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

  message("FUSION v7 processing complete.")
}
