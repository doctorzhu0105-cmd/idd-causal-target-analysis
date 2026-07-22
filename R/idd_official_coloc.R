#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function(args) {
  result <- list(repo_root = normalizePath(".", winslash = "/", mustWork = TRUE))
  index <- 1L
  while (index <= length(args)) {
    key <- args[[index]]
    if (key == "--repo-root" && index < length(args)) {
      result$repo_root <- normalizePath(args[[index + 1L]], winslash = "/", mustWork = TRUE)
      index <- index + 2L
    } else {
      stop("Unknown or incomplete argument: ", key)
    }
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- args$repo_root
local_library <- file.path(repo_root, "env", "R-library")
.libPaths(c(local_library, .libPaths()))

if (!requireNamespace("coloc", quietly = TRUE)) {
  stop("The official coloc package is not installed in the registered R library")
}

coloc_version <- as.character(utils::packageVersion("coloc"))
focus_genes <- c("AKR1C2", "AKR1C1", "NT5E", "CCL28", "P2RY13")
outcomes <- c("M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA")
outcome_counts <- data.table(
  outcome_id = outcomes,
  cases = c(53270L, 45696L, 28094L),
  controls = c(347768L, 347768L, 347768L)
)

core_dir <- file.path(repo_root, "artifacts", "core", "genetics")
processed_dir <- file.path(repo_root, "data", "processed", "idd", "genetics")
output_path <- file.path(repo_root, "HANDOFF", "COLOC_COMPARISON.tsv")
qc_path <- file.path(repo_root, "manuscript", "qc", "official_coloc_results.tsv")

old_coloc <- fread(file.path(core_dir, "coloc_results.tsv"), na.strings = c("", "NA"))
mr <- fread(file.path(core_dir, "mr_results.tsv"), na.strings = c("", "NA"))
full_eqtl <- fread(file.path(processed_dir, "candidate_full_eqtl.tsv.gz"), na.strings = c("", "NA"))

setorder(full_eqtl, ensembl_gene_id, snp, pvalue)
full_eqtl <- full_eqtl[, .SD[1L], by = .(ensembl_gene_id, snp)]
setkey(full_eqtl, ensembl_gene_id)

outcome_tables <- setNames(vector("list", length(outcomes)), outcomes)
for (outcome_id in outcomes) {
  path <- file.path(processed_dir, paste0(outcome_id, ".candidate_variants.tsv.gz"))
  table <- fread(path, na.strings = c("", "NA"))
  setorder(table, snp, pvalue)
  outcome_tables[[outcome_id]] <- table[, .SD[1L], by = snp]
}

mr_lookup <- mr[
  direction == "expression_to_disease",
  .(
    ensembl_gene_id,
    gene_symbol,
    outcome_id,
    mr_status = status,
    mr_pvalue = pvalue,
    formal_set_fdr = fdr
  )
]

run_one <- function(index) {
  archived <- old_coloc[index]
  gene_id <- archived$ensembl_gene_id[[1L]]
  outcome_id <- archived$outcome_id[[1L]]
  exposure <- full_eqtl[.(gene_id)]
  outcome <- outcome_tables[[outcome_id]]
  shared <- merge(
    exposure,
    outcome,
    by = "snp",
    suffixes = c("_expression", "_outcome"),
    all = FALSE,
    sort = TRUE
  )
  shared <- shared[
    is.finite(beta_expression) & is.finite(se_expression) & se_expression > 0 &
      is.finite(beta_outcome) & is.finite(se_outcome) & se_outcome > 0
  ]

  base <- data.table(
    ensembl_gene_id = gene_id,
    gene_symbol = archived$gene_symbol[[1L]],
    outcome_id = outcome_id,
    method = "coloc::coloc.abf",
    coloc_version = coloc_version,
    expression_type = "quant",
    expression_sdY = 1,
    outcome_type = "cc",
    window_definition = "eQTLGen variants within +/-1,000,000 bp of the registered gene position",
    requested_p1 = 1e-4,
    requested_p2 = 1e-4,
    requested_p12 = 1e-5,
    expression_effect_prior_sd = 0.15,
    outcome_effect_prior_sd = 0.20,
    shared_variants = nrow(shared),
    archived_shared_variants = archived$shared_variants[[1L]],
    archived_pp0 = archived$pp0[[1L]],
    archived_pp1 = archived$pp1[[1L]],
    archived_pp2 = archived$pp2[[1L]],
    archived_pp3 = archived$pp3[[1L]],
    archived_pp4 = archived$pp4[[1L]],
    official_pp0 = NA_real_,
    official_pp1 = NA_real_,
    official_pp2 = NA_real_,
    official_pp3 = NA_real_,
    official_pp4 = NA_real_,
    effective_p1 = NA_real_,
    effective_p2 = NA_real_,
    effective_p12 = NA_real_,
    pp4_absolute_difference = NA_real_,
    pp4_difference_gt_0_05 = NA,
    status = "not_estimable",
    reason = "insufficient_shared_variants",
    warning = ""
  )
  if (nrow(shared) < 50L) {
    return(base)
  }

  current_outcome_id <- outcome_id
  counts <- outcome_counts[outcome_id == current_outcome_id]
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      {
        capture.output(
          fit_value <- coloc::coloc.abf(
            dataset1 = list(
              beta = shared$beta_expression,
              varbeta = shared$se_expression^2,
              type = "quant",
              sdY = 1,
              snp = shared$snp
            ),
            dataset2 = list(
              beta = shared$beta_outcome,
              varbeta = shared$se_outcome^2,
              type = "cc",
              N = counts$cases[[1L]] + counts$controls[[1L]],
              s = counts$cases[[1L]] / (counts$cases[[1L]] + counts$controls[[1L]]),
              snp = shared$snp
            ),
            p1 = 1e-4,
            p2 = 1e-4,
            p12 = 1e-5
          )
        )
        fit_value
      },
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  if (inherits(fit, "condition")) {
    base$reason <- paste0("official_coloc_error: ", conditionMessage(fit))
    base$warning <- gsub("[[:space:]]+", " ", trimws(paste(unique(warnings), collapse = " | ")))
    return(base)
  }

  summary <- fit$summary
  base[, `:=`(
    official_pp0 = unname(summary[["PP.H0.abf"]]),
    official_pp1 = unname(summary[["PP.H1.abf"]]),
    official_pp2 = unname(summary[["PP.H2.abf"]]),
    official_pp3 = unname(summary[["PP.H3.abf"]]),
    official_pp4 = unname(summary[["PP.H4.abf"]]),
    effective_p1 = unname(fit$priors[["p1"]]),
    effective_p2 = unname(fit$priors[["p2"]]),
    effective_p12 = unname(fit$priors[["p12"]]),
    status = "estimable",
    reason = "",
    warning = gsub("[[:space:]]+", " ", trimws(paste(unique(warnings), collapse = " | ")))
  )]
  base[, pp4_absolute_difference := abs(official_pp4 - archived_pp4)]
  base[, pp4_difference_gt_0_05 := pp4_absolute_difference > 0.05]
  base
}

result <- rbindlist(lapply(seq_len(nrow(old_coloc)), run_one), fill = TRUE)
result <- merge(result, mr_lookup, by = c("ensembl_gene_id", "gene_symbol", "outcome_id"), all.x = TRUE)
result[, focus_gene := gene_symbol %chin% focus_genes]
result[, mr_raw_p_lt_0_05 := is.finite(mr_pvalue) & mr_pvalue < 0.05]
result[, formal_set_fdr_lt_0_05 := is.finite(formal_set_fdr) & formal_set_fdr < 0.05]
setorder(result, -focus_gene, gene_symbol, outcome_id)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(qc_path), recursive = TRUE, showWarnings = FALSE)
fwrite(result, output_path, sep = "\t", na = "not_estimable", quote = FALSE)
fwrite(result, qc_path, sep = "\t", na = "not_estimable", quote = FALSE)

cat("official coloc version:", coloc_version, "\n")
cat("rows:", nrow(result), "\n")
cat("estimable:", result[status == "estimable", .N], "\n")
cat("MR raw P < 0.05 rows:", result[mr_raw_p_lt_0_05 == TRUE, .N], "\n")
cat("PP4 absolute difference > 0.05:", result[pp4_difference_gt_0_05 == TRUE, .N], "\n")
cat("output:", output_path, "\n")
