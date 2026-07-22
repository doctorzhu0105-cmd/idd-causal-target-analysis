#!/usr/bin/env Rscript

mrpresso_lib <- Sys.getenv("MRPRESSO_LIB", unset = "")
if (nzchar(mrpresso_lib)) {
  .libPaths(c(mrpresso_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  result <- list(repo_root = ".", nb_distribution = 10000L)
  index <- 1L
  while (index <= length(args)) {
    if (args[[index]] == "--repo-root") {
      result$repo_root <- args[[index + 1L]]
      index <- index + 2L
    } else if (args[[index]] == "--nb-distribution") {
      result$nb_distribution <- as.integer(args[[index + 1L]])
      index <- index + 2L
    } else {
      stop(sprintf("Unknown argument: %s", args[[index]]))
    }
  }
  result
}

extract_scalar <- function(value, name) {
  if (is.null(value)) return(NA_real_)
  if (is.data.frame(value) && name %in% names(value)) return(as.numeric(value[[name]][1]))
  if (is.matrix(value) && name %in% colnames(value)) return(as.numeric(value[1, name]))
  if (is.list(value) && !is.null(value[[name]])) return(as.numeric(value[[name]][1]))
  NA_real_
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- normalizePath(args$repo_root, winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(repo_root, "manuscript", ".r_lib"), .libPaths()))
input_path <- file.path(repo_root, "artifacts", "core", "genetics", "harmonized_instruments.tsv")
output_dir <- file.path(repo_root, "manuscript", "supplement")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  stop("MRPRESSO 1.0 is required. The analysis is pinned to GitHub commit 3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87.")
}

set.seed(20260716)
instruments <- fread(input_path)
forward <- instruments[direction == "expression_to_disease"]
groups <- unique(forward[, .(gene_symbol, ensembl_gene_id, outcome_id)])
setorder(groups, gene_symbol, outcome_id)

rows <- vector("list", nrow(groups))
for (index in seq_len(nrow(groups))) {
  key <- groups[index]
  subset <- forward[
    gene_symbol == key$gene_symbol &
      ensembl_gene_id == key$ensembl_gene_id &
      outcome_id == key$outcome_id
  ]
  n_instruments <- nrow(subset)
  row <- data.table(
    gene_symbol = key$gene_symbol,
    ensembl_gene_id = key$ensembl_gene_id,
    outcome_id = key$outcome_id,
    n_instruments = n_instruments,
    global_rssobs = NA_real_,
    global_pvalue = NA_real_,
    outlier_count = NA_integer_,
    distortion_coefficient = NA_real_,
    distortion_pvalue = NA_real_,
    status = "not_estimable",
    reason = "fewer_than_4_instruments"
  )
  if (n_instruments >= 4L) {
    result <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta_outcome",
        BetaExposure = "beta_exposure",
        SdOutcome = "se_outcome",
        SdExposure = "se_exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = as.data.frame(subset),
        NbDistribution = args$nb_distribution,
        SignifThreshold = 0.05
      ),
      error = function(error) error
    )
    if (inherits(result, "error")) {
      row$reason <- paste0("mr_presso_error:", conditionMessage(result))
    } else {
      tests <- result[["MR-PRESSO results"]]
      global <- tests[["Global Test"]]
      outlier <- tests[["Outlier Test"]]
      distortion <- tests[["Distortion Test"]]
      row$global_rssobs <- extract_scalar(global, "RSSobs")
      row$global_pvalue <- extract_scalar(global, "Pvalue")
      row$outlier_count <- if (is.null(outlier)) 0L else nrow(as.data.frame(outlier))
      row$distortion_coefficient <- extract_scalar(distortion, "Distortion Coefficient")
      row$distortion_pvalue <- extract_scalar(distortion, "Pvalue")
      row$status <- "estimable"
      row$reason <- ""
    }
  }
  rows[[index]] <- row
}

results <- rbindlist(rows, fill = TRUE)
results[, global_fdr := NA_real_]
results[status == "estimable", global_fdr := p.adjust(global_pvalue, method = "BH")]
results[, `:=`(
  seed = 20260716L,
  nb_distribution = args$nb_distribution,
  package_version = as.character(utils::packageVersion("MRPRESSO")),
  package_commit = "3e3c92d7eda6dce0d1d66077373ec0f7ff4f7e87"
)]
fwrite(results, file.path(output_dir, "mr_presso_results.tsv"), sep = "\t", na = "not_estimable")

summary <- results[, .(
  total_combinations = .N,
  estimable_combinations = sum(status == "estimable"),
  bh_significant_global_tests = sum(global_fdr < 0.05, na.rm = TRUE),
  priority_target_estimable_combinations = sum(
    status == "estimable" & gene_symbol %chin% c("AKR1C2", "NT5E", "CCL28", "P2RY13")
  )
)]
summary[, interpretation := "A non-significant global test is reported as no detected evidence, not proof of no horizontal pleiotropy."]
fwrite(summary, file.path(output_dir, "mr_presso_summary.tsv"), sep = "\t", na = "not_estimable")
cat(sprintf("MR-PRESSO complete: %d/%d combinations estimable\n", summary$estimable_combinations, summary$total_combinations))
