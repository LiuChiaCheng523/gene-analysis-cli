#!/usr/bin/env Rscript

usage <- function() {
  cat(
    paste0(
      "Usage:\n",
      "  Rscript make_twb2_pheno_cli.R <var_name> <var_type> [case_label] [control_label] [base_dir]\n",
      "  Rscript make_twb2_pheno_cli.R --var_name <var_name> --var_type <categorical|continuous> [--case_label <value>] [--control_label <value>] [--base_dir <dir>]\n\n",
      "Arguments:\n",
      "  var_name        Variable name, for example: VERTIGO_SELF\n",
      "  var_type        categorical or continuous\n",
      "  case_label      Case label for categorical phenotype, for example: 1\n",
      "  control_label   Control label for categorical phenotype, for example: 0\n\n",
      "  base_dir        Data base directory. Default: /mnt/SP-siliconpower/TWB20250806download\n\n",
      "Example:\n",
      "  Rscript make_twb2_pheno_cli.R VERTIGO_SELF categorical 1 0\n",
      "  Rscript make_twb2_pheno_cli.R BMI continuous\n",
      "  Rscript make_twb2_pheno_cli.R BMI continuous NA NA /mnt/SP-siliconpower/TWB20250806download\n"
    )
  )
}

parse_args <- function(args) {
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1, 0), save = "no")
  }

  parsed <- list(
    var_name = NA_character_,
    var_type = NA_character_,
    case_label = NA_character_,
    control_label = NA_character_,
    base_dir = "/mnt/SP-siliconpower/TWB20250806download"
  )

  if (any(grepl("^--", args))) {
    i <- 1
    while (i <= length(args)) {
      arg <- args[[i]]

      if (arg == "--var_name") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$var_name <- args[[i + 1]]
        i <- i + 2
      } else if (arg == "--var_type") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$var_type <- args[[i + 1]]
        i <- i + 2
      } else if (arg == "--case_label") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$case_label <- args[[i + 1]]
        i <- i + 2
      } else if (arg == "--control_label") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$control_label <- args[[i + 1]]
        i <- i + 2
      } else if (arg == "--base_dir") {
        if (i == length(args)) stop("Missing value after ", arg, call. = FALSE)
        parsed$base_dir <- args[[i + 1]]
        i <- i + 2
      } else {
        stop("Unknown argument: ", arg, call. = FALSE)
      }
    }
  } else {
    parsed$var_name <- args[[1]]
    parsed$var_type <- args[[2]]
    if (length(args) >= 3) parsed$case_label <- args[[3]]
    if (length(args) >= 4) parsed$control_label <- args[[4]]
    if (length(args) >= 5) parsed$base_dir <- args[[5]]
  }

  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

var_name <- args$var_name
var_type <- args$var_type
case_label <- args$case_label
control_label <- args$control_label
base_dir <- args$base_dir

if (!dir.exists(base_dir)) {
  stop("base_dir not found: ", base_dir, call. = FALSE)
}

release_cols_path <- file.path(base_dir, "survey/release_list_colnames.txt")
if (!file.exists(release_cols_path)) {
  stop("release_list_colnames.txt not found: ", release_cols_path, call. = FALSE)
}

release_cols <- fread(
  release_cols_path,
  header = FALSE
)$V1

var_type_list <- c("categorical", "continuous")

if (!(var_name %in% release_cols)) {
  message(paste0("Cannot find variable: ", var_name, " in release_list_colnames.txt!"))
  message("Please try another variable name!")
  stop(call. = FALSE)
}

if (!(var_type %in% var_type_list)) {
  message("var_type should be categorical or continuous")
  message("Please try again!")
  stop(call. = FALSE)
}

if (var_type == "categorical" && (is.na(case_label) || is.na(control_label))) {
  message("case_label and control_label are required for categorical phenotype")
  stop(call. = FALSE)
}

release_list_path <- file.path(base_dir, "survey/release_list_survey.csv")
lab_info_path <- file.path(base_dir, "lab_info/lab_info.csv")
twb2_fam_path <- file.path(base_dir, "Imputed.120161.TWB2/imputed_120161/TWB2.hg38.impu.v4.fam")

if (!file.exists(release_list_path)) {
  stop("release_list_survey.csv not found: ", release_list_path, call. = FALSE)
}

if (!file.exists(lab_info_path)) {
  stop("lab_info.csv not found: ", lab_info_path, call. = FALSE)
}

if (!file.exists(twb2_fam_path)) {
  stop("TWB2.hg38.impu.v4.fam not found: ", twb2_fam_path, call. = FALSE)
}

release_list <- fread(release_list_path)
lab_info <- fread(lab_info_path)
TWB2_fam <- fread(twb2_fam_path)

release_list_var <- release_list[, c(
  "Release_No", "FOLLOW", "CRF_NAME_QN", "SURVEY_DATE",
  "ID_BIRTH", "AGE", "SEX", var_name
), with = FALSE]

release_list_var <- release_list_var %>%
  filter(FOLLOW == "Baseline")

lab_info_TWB2 <- lab_info %>%
  filter(TWB2_ID != "") %>%
  dplyr::select("Release_No", "FOLLOW", "PROJECT", "TWB2_ID", "TWB2_Batch")

TWB2_info <- merge(lab_info_TWB2, release_list_var, by = "Release_No", all.x = FALSE)

pheno <- copy(TWB2_fam[, .(FID = V1, IID = V2)])

output_dir <- file.path(base_dir, "Imputed.120161.TWB2/imputed_120161")
pheno_file <- file.path(output_dir, paste0("TWB2_", var_name, "_pheno.txt"))
keep_file <- file.path(output_dir, paste0("TWB2_", var_name, "_keep.txt"))

if (var_type == "categorical") {
  TWB2_var_case <- TWB2_info$TWB2_ID[as.character(TWB2_info[[var_name]]) == case_label]
  TWB2_var_control <- TWB2_info$TWB2_ID[as.character(TWB2_info[[var_name]]) == control_label]

  if (length(TWB2_var_case) == 0) {
    message(paste0("There is no ", var_name, " case sample in TWB2!"))
    stop(call. = FALSE)
  } else if (length(TWB2_var_control) == 0) {
    message(paste0("There is no ", var_name, " control sample in TWB2!"))
    stop(call. = FALSE)
  }

  pheno[, (var_name) := NA_integer_]
  pheno[IID %in% TWB2_var_control, (var_name) := 1L]
  pheno[IID %in% TWB2_var_case, (var_name) := 2L]

  pheno2 <- pheno[!is.na(get(var_name))]
  keep <- pheno2[, .(FID, IID)]

  fwrite(pheno2, pheno_file, sep = "\t", quote = FALSE)
  fwrite(keep, keep_file, sep = "\t", quote = FALSE)

  message("case in fam: ", sum(pheno$IID %in% TWB2_var_case))
  message("control in fam: ", sum(pheno$IID %in% TWB2_var_control))
  message("total: ", nrow(pheno2))
} else if (var_type == "continuous") {
  TWB2_var <- TWB2_info %>%
    dplyr::select(TWB2_ID, all_of(var_name)) %>%
    filter(!is.na(.data[[var_name]]))

  if (nrow(TWB2_var) == 0) {
    message(paste0("There is no non-missing ", var_name, " sample in TWB2!"))
    stop(call. = FALSE)
  }

  pheno_value_map <- data.table(
    IID = TWB2_var$TWB2_ID,
    pheno_value = TWB2_var[[var_name]]
  )

  pheno <- merge(pheno, pheno_value_map, by = "IID", all.x = TRUE)
  setcolorder(pheno, c("FID", "IID", "pheno_value"))
  setnames(pheno, "pheno_value", var_name)

  pheno2 <- pheno[!is.na(get(var_name))]
  keep <- pheno2[, .(FID, IID)]

  fwrite(pheno2, pheno_file, sep = "\t", quote = FALSE)
  fwrite(keep, keep_file, sep = "\t", quote = FALSE)

  message("non-missing in fam: ", nrow(pheno2))
  message("total with value in survey: ", nrow(TWB2_var))
}

message("pheno file: ", pheno_file)
message("keep file: ", keep_file)
