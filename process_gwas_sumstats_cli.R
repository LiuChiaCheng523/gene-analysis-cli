#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript process_gwas_sumstats_cli.R <project_name> [base_dir]\n",
      "  Rscript process_gwas_sumstats_cli.R --project_name <project_name> [--base_dir <base_dir>]\n\n",
      "Arguments:\n",
      "  project_name  Project name, for example: TWB1_LAA\n",
      "  base_dir      Base directory. Default: /mnt/data/ai_agent/gene_analysis\n\n",
      "Example:\n",
      "  Rscript process_gwas_sumstats_cli.R TWB1_LAA /mnt/data/ai_agent/gene_analysis\n"
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
    parsed <- list(project_name = NA_character_, base_dir = default_base_dir)
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
    base_dir = ifelse(length(args) >= 2, args[[2]], default_base_dir)
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

project_name <- args$project_name
base_dir <- args$base_dir

glm_logistic_input_folder <- file.path(base_dir, "PLINK/imputed/glm_logistic", project_name)
cojo_output_folder <- file.path(base_dir, "COJO", project_name)
fusion_output_folder <- file.path(base_dir, "FUSION/GWAS", project_name)
spredixcan_output_folder <- file.path(base_dir, "S_PrediXcan/GWAS", project_name)
log_folder <- file.path(fusion_output_folder, "log")

dir.create(cojo_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(fusion_output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(spredixcan_output_folder, recursive = TRUE, showWarnings = FALSE)
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
message("Input folder: ", glm_logistic_input_folder)
message("Log file: ", log_file)

spredixcan_twas_format <- data.frame()

chr_list <- 1:22
chr_with_data <- c()
chr_without_data <- c()

for (chr_num in chr_list) {
  input_file <- file.path(
    glm_logistic_input_folder,
    paste0("chr", chr_num, ".PHENO1.glm.logistic.hybrid")
  )

  write_log(paste0("[", Sys.time(), "] Checking chr", chr_num, "\n"))

  if (!file.exists(input_file)) {
    warning_msg <- paste0(
      "[", Sys.time(), "] WARNING: chr",
      chr_num,
      " file does not exist. Skip.\n"
    )

    warning(warning_msg)
    write_log(warning_msg)

    chr_without_data <- c(chr_without_data, chr_num)
    next
  }

  chr_with_data <- c(chr_with_data, chr_num)

  write_log(paste0("[", Sys.time(), "] Reading chr", chr_num, "\n"))

  gwas <- fread(input_file)
  gwas <- gwas[!is.na(P)]

  cojo_ma <- data.frame(
    SNP = gwas$ID,
    A1 = gwas$A1,
    A2 = ifelse(gwas$A1 == gwas$REF, gwas$ALT, gwas$REF),
    freq = gwas$A1_FREQ,
    b = log(gwas$OR),
    se = gwas$`LOG(OR)_SE`,
    p = gwas$P,
    N = gwas$OBS_CT
  )

  write.table(
    cojo_ma,
    file = file.path(cojo_output_folder, paste0("chr", chr_num, ".cojo.ma")),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )

  fusion_twas_format <- data.frame(
    SNP = gwas$ID,
    A1 = gwas$A1,
    A2 = ifelse(gwas$A1 == gwas$REF, gwas$ALT, gwas$REF),
    Z = gwas$Z_STAT,
    CHR = gwas$`#CHROM`,
    p = gwas$P
  )

  fwrite(
    fusion_twas_format,
    file.path(fusion_output_folder, paste0("chr", chr_num, ".sumstats")),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  spredixcan_twas_format <- rbind(
    spredixcan_twas_format,
    fusion_twas_format
  )

  write_log(paste0("[", Sys.time(), "] chr", chr_num, " completed.\n"))
}

if (nrow(spredixcan_twas_format) > 0) {
  fwrite(
    spredixcan_twas_format,
    file.path(spredixcan_output_folder, "chr1tochr22.sumstats"),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  write_log(paste0("[", Sys.time(), "] SPrediXcan combined sumstats file created.\n"))
} else {
  warning_msg <- paste0(
    "[", Sys.time(),
    "] WARNING: No chromosome data found. SPrediXcan combined file was not created.\n"
  )

  warning(warning_msg)
  write_log(warning_msg)
}

summary_text <- paste0(
  "\n========== Summary ==========\n",
  "Chromosomes with data: ",
  paste(chr_with_data, collapse = ", "),
  "\n",
  "Chromosomes without data: ",
  paste(chr_without_data, collapse = ", "),
  "\n",
  "Pipeline finished at: ",
  Sys.time(),
  "\n=============================\n"
)

cat(summary_text)
write_log(summary_text)
