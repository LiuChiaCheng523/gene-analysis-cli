#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript cojo_annotation_manhattan_cli.R <project_name> [base_dir] [genome_build]\n",
      "  Rscript cojo_annotation_manhattan_cli.R --project_name <project_name> [--base_dir <base_dir>] [--genome_build <grch37|grch38>] [--gtf_file <path>]\n\n",
      "Arguments:\n",
      "  project_name  Project name, for example: TWB1_LAA_test1\n",
      "  base_dir      Base directory. Default: /mnt/data/ai_agent/gene_analysis\n\n",
      "  genome_build  grch37 or grch38. Default: grch37\n",
      "  gtf_file      Optional explicit GTF path. If not given, use tools/ensembl default by genome_build\n\n",
      "Example:\n",
      "  Rscript cojo_annotation_manhattan_cli.R TWB1_LAA_test1 /mnt/data/ai_agent/gene_analysis grch37\n",
      "  Rscript cojo_annotation_manhattan_cli.R --project_name TWB1_LAA_test1 --base_dir /mnt/data/ai_agent/gene_analysis --genome_build grch38\n"
    )
  )
}

parse_args <- function(args) {
  default_base_dir <- "/mnt/data/ai_agent/gene_analysis"
  default_genome_build <- "grch37"

  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  if (any(grepl("^--", args))) {
    parsed <- list(
      project_name = NA_character_,
      base_dir = default_base_dir,
      genome_build = default_genome_build,
      gtf_file = NA_character_
    )
    i <- 1

    while (i <= length(args)) {
      arg <- args[[i]]

      if (arg %in% c("--project_name", "--project-name", "-p")) {
        if (i == length(args)) {
          stop("Missing value after ", arg, call. = FALSE)
        }
        parsed$project_name <- args[[i + 1]]
        i <- i + 2
      } else if (arg %in% c("--base_dir", "--base-dir", "-b")) {
        if (i == length(args)) {
          stop("Missing value after ", arg, call. = FALSE)
        }
        parsed$base_dir <- args[[i + 1]]
        i <- i + 2
      } else if (arg %in% c("--genome_build", "--genome-build", "-g")) {
        if (i == length(args)) {
          stop("Missing value after ", arg, call. = FALSE)
        }
        parsed$genome_build <- tolower(args[[i + 1]])
        i <- i + 2
      } else if (arg == "--gtf_file") {
        if (i == length(args)) {
          stop("Missing value after ", arg, call. = FALSE)
        }
        parsed$gtf_file <- args[[i + 1]]
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

  list(
    project_name = args[[1]],
    base_dir = ifelse(length(args) >= 2, args[[2]], default_base_dir),
    genome_build = ifelse(length(args) >= 3, tolower(args[[3]]), default_genome_build),
    gtf_file = NA_character_
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(rtracklayer))

project_name <- args$project_name
base_dir <- args$base_dir
genome_build <- args$genome_build

if (!(genome_build %in% c("grch37", "grch38"))) {
  stop("genome_build must be grch37 or grch38.", call. = FALSE)
}

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

cojo_result_input_folder <- file.path(base_dir, "COJO", project_name, "result")
glm_logistic_input_folder <- file.path(base_dir, "PLINK/imputed/glm_logistic", project_name)
cojo_table_output_folder <- file.path(base_dir, "COJO", project_name, "table")
cojo_manhattan_output_folder <- file.path(base_dir, "COJO", project_name, "manhattan")
log_folder <- file.path(base_dir, "COJO", project_name, "log")
gtf_file <- if (!is.na(args$gtf_file) && args$gtf_file != "") {
  args$gtf_file
} else if (genome_build == "grch37") {
  file.path(base_dir, "tools/ensembl/Homo_sapiens.GRCh37.87.gtf")
} else {
  file.path(base_dir, "tools/ensembl/Homo_sapiens.GRCh38.115.gtf")
}

if (!dir.exists(cojo_result_input_folder)) {
  stop("COJO result input folder not found: ", cojo_result_input_folder, call. = FALSE)
}

if (!dir.exists(glm_logistic_input_folder)) {
  stop("GWAS input folder not found: ", glm_logistic_input_folder, call. = FALSE)
}

dir.create(cojo_table_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(cojo_manhattan_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_folder,
  paste0("process_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

write_log <- function(message) {
  cat(message, file = log_file, append = TRUE)
}

cat(
  paste0(
    "========== Pipeline Start ==========\n",
    "Time: ", Sys.time(), "\n",
    "Project name: ", project_name, "\n",
    "Base dir: ", base_dir, "\n\n"
  ),
  file = log_file
)

message("Project name: ", project_name)
message("Base dir: ", base_dir)
message("Genome build: ", genome_build)
message("COJO result input folder: ", cojo_result_input_folder)
message("GWAS input folder: ", glm_logistic_input_folder)
message("GTF file: ", gtf_file)
message("Log file: ", log_file)

if (!file.exists(gtf_file)) {
  stop("GTF file not found: ", gtf_file, call. = FALSE)
}

# load local Ensembl GRCh37 GTF -----
gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

gene_anno <- gtf_df %>%
  filter(type == "gene") %>%
  transmute(
    hgnc_symbol = gene_name,
    ensembl_gene_id_version = gene_id,
    chromosome_name = as.character(seqnames),
    start_position = start,
    end_position = end,
    gene_biotype = gene_biotype
  )

# LAA gwas fine-mapping (COJO) -----
annotate_one_snp <- function(chr, bp, snp_id, window, gene_anno) {
  genes <- gene_anno %>%
    filter(
      chromosome_name == as.character(chr),
      start_position <= bp + window,
      end_position >= bp - window
    )

  if (nrow(genes) == 0) return(NULL)

  genes$SNP <- snp_id
  genes$SNP_bp <- bp
  genes$distance <- ifelse(
    bp >= genes$start_position & bp <= genes$end_position,
    0,
    pmin(
      abs(genes$start_position - bp),
      abs(genes$end_position - bp)
    )
  )
  genes$location <- ifelse(
    bp >= genes$start_position & bp <= genes$end_position,
    "intragenic",
    ifelse(
      bp < genes$start_position,
      "upstream",
      "downstream"
    )
  )

  genes
}

chr_list <- 1:22
window <- 500000

all_cojo_snps <- data.frame()
for (chr in chr_list) {
  cojo_file <- paste0(
    cojo_result_input_folder, "/chr", chr, "_p1e5.jma.cojo"
  )

  if (!file.exists(cojo_file)) {
    warning_msg <- paste0("File not found: ", cojo_file, "\n")
    warning(warning_msg)
    write_log(warning_msg)
    next
  }

  cojo_result <- fread(cojo_file)
  cojo_snps <- data.frame(
    SNP = cojo_result$SNP,
    chr = cojo_result$Chr,
    bp = cojo_result$bp
  )
  all_cojo_snps <- rbind(all_cojo_snps, cojo_snps)
}

all_anno_list <- list()
for (chr in chr_list) {
  message("Annotating chr", chr)
  write_log(paste0("[", Sys.time(), "] Annotating chr", chr, "\n"))

  cojo_file <- paste0(
    cojo_result_input_folder, "/chr", chr, "_p1e5.jma.cojo"
  )

  if (!file.exists(cojo_file)) {
    warning_msg <- paste0("File not found: ", cojo_file, "\n")
    warning(warning_msg)
    write_log(warning_msg)
    next
  }

  cojo_result <- fread(cojo_file)

  if (nrow(cojo_result) == 0) next

  cojo_snps <- data.frame(
    SNP = cojo_result$SNP,
    chr = cojo_result$Chr,
    bp = cojo_result$bp
  )

  chr_anno <- lapply(
    seq_len(nrow(cojo_snps)),
    function(i) {
      annotate_one_snp(
        chr = cojo_snps$chr[i],
        bp = cojo_snps$bp[i],
        snp_id = cojo_snps$SNP[i],
        window = window,
        gene_anno = gene_anno
      )
    }
  )

  chr_anno <- do.call(rbind, chr_anno)

  if (!is.null(chr_anno)) {
    all_anno_list[[paste0("chr", chr)]] <- chr_anno
  }
}

anno <- do.call(rbind, all_anno_list)
if (is.null(anno) || nrow(anno) == 0) {
  stop("No annotation result found. Please check COJO result files.", call. = FALSE)
}

write.csv(
  anno,
  paste0(cojo_table_output_folder, "/p1e5_annotation_genes.csv"),
  row.names = FALSE
)

protein_nearest_gene <- anno %>%
  filter(
    gene_biotype == "protein_coding",
    hgnc_symbol != ""
  ) %>%
  group_by(SNP) %>%
  slice_min(distance, n = 1, with_ties = FALSE) %>%
  ungroup()

write.csv(
  protein_nearest_gene,
  paste0(cojo_table_output_folder, "/p1e5_protein_coding_nearest_gene.csv"),
  row.names = FALSE
)

# GWAS Manhattan plot -----
all_gwas <- data.frame()
for (chr_num in 1:22) {
  gwas_path <- file.path(
    glm_logistic_input_folder,
    paste0("chr", chr_num, ".PHENO1.glm.logistic.hybrid")
  )

  if (!file.exists(gwas_path)) {
    warning_msg <- paste0("File not found: ", gwas_path, "\n")
    warning(warning_msg)
    write_log(warning_msg)
    next
  }

  gwas <- fread(gwas_path)
  gwas <- gwas[!is.na(gwas$P), ]
  all_gwas <- rbind(all_gwas, gwas)
}

if (nrow(all_gwas) == 0) {
  stop("No GWAS data found. Please check PLINK glm logistic input files.", call. = FALSE)
}

all_gwas_plot <- as.data.frame(all_gwas)[, c("#CHROM", "POS", "ID", "P"), drop = FALSE]
colnames(all_gwas_plot) <- c("CHR", "BP", "SNP", "P")

gwas_df <- all_gwas_plot %>%
  filter(
    CHR %in% 1:22,
    !is.na(P),
    is.finite(P)
  ) %>%
  mutate(
    CHR = as.numeric(CHR),
    BP = as.numeric(BP)
  ) %>%
  arrange(CHR, BP)

chr_max_gwas <- gwas_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

gwas_df <- gwas_df %>%
  left_join(chr_max_gwas, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

label_df_gwas <- protein_nearest_gene %>%
  distinct(SNP, hgnc_symbol, SNP_bp) %>%
  left_join(
    gwas_df,
    by = c("SNP" = "SNP")
  ) %>%
  filter(!is.na(P))

gwas_df <- gwas_df %>%
  mutate(chr_color = factor(CHR %% 2))

p_gwas <- ggplot(gwas_df, aes(x = BP_cum, y = -log10(P))) +
  geom_point(
    aes(color = chr_color),
    size = 0.4,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("grey60", "black"),
    guide = "none"
  ) +
  geom_hline(
    yintercept = -log10(1e-5),
    color = "blue",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_point(
    data = label_df_gwas,
    aes(x = BP_cum, y = -log10(P)),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = label_df_gwas,
    aes(label = hgnc_symbol),
    color = "red",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  labs(
    title = "GWAS Manhattan Plot",
    x = "Chromosome",
    y = expression(-log[10](P))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

chr_centers_gwas <- gwas_df %>%
  group_by(CHR) %>%
  summarise(center = mean(range(BP_cum)), .groups = "drop")

p_gwas <- p_gwas +
  scale_x_continuous(
    breaks = chr_centers_gwas$center,
    labels = chr_centers_gwas$CHR
  )

p_gwas <- p_gwas +
  scale_y_continuous(
    breaks = seq(0, ceiling(max(-log10(gwas_df$P))), by = 1)
  )

ggsave(
  file.path(cojo_manhattan_output_folder, "cojo_protein_coding_manhanttan.png"),
  p_gwas,
  width = 12,
  height = 6,
  dpi = 600
)

summary_text <- paste0(
  "\n========== Summary ==========\n",
  "Annotation table: ",
  file.path(cojo_table_output_folder, "p1e5_annotation_genes.csv"),
  "\n",
  "Protein-coding nearest gene table: ",
  file.path(cojo_table_output_folder, "p1e5_protein_coding_nearest_gene.csv"),
  "\n",
  "Manhattan plot: ",
  file.path(cojo_manhattan_output_folder, "cojo_protein_coding_manhanttan.png"),
  "\n",
  "Pipeline finished at: ",
  Sys.time(),
  "\n=============================\n"
)

cat(summary_text)
write_log(summary_text)
