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

sanitize <- function(values) {
  gsub("[[:space:]]+", " ", trimws(paste(unique(values[nzchar(values)]), collapse = " | ")))
}

is_palindromic <- function(a1, a2) {
  pairs <- paste0(pmin(a1, a2), pmax(a1, a2))
  pairs %chin% c("AT", "CG")
}

complement <- function(allele) {
  replacements <- c(A = "T", T = "A", C = "G", G = "C")
  unname(replacements[allele])
}

align_beta_to_a1 <- function(beta, effect_allele, other_allele, a1, a2) {
  if (!all(c(effect_allele, other_allele, a1, a2) %chin% c("A", "C", "G", "T"))) {
    return(NA_real_)
  }
  if (effect_allele == a1 && other_allele == a2) {
    return(beta)
  }
  if (effect_allele == a2 && other_allele == a1) {
    return(-beta)
  }
  if (complement(effect_allele) == a1 && complement(other_allele) == a2) {
    return(beta)
  }
  if (complement(effect_allele) == a2 && complement(other_allele) == a1) {
    return(-beta)
  }
  NA_real_
}

run_command <- function(command, arguments, stdout_path, stderr_path) {
  status <- system2(
    command,
    args = arguments,
    stdout = stdout_path,
    stderr = stderr_path,
    wait = TRUE
  )
  if (!identical(status, 0L)) {
    stop("Command failed with exit status ", status, ": ", command)
  }
}

run_susie <- function(dataset, label) {
  messages <- character()
  warnings <- character()
  started <- Sys.time()
  fit <- tryCatch(
    withCallingHandlers(
      coloc::runsusie(
        dataset,
        maxit = 1000,
        repeat_until_convergence = FALSE
      ),
      message = function(condition) {
        messages <<- c(messages, conditionMessage(condition))
        invokeRestart("muffleMessage")
      },
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (inherits(fit, "condition")) {
    return(list(
      status = "failed",
      reason = paste0(label, ": ", conditionMessage(fit)),
      warning = sanitize(warnings),
      message = sanitize(messages),
      elapsed_seconds = elapsed,
      fit = NULL
    ))
  }
  list(
    status = "success",
    reason = "",
    warning = sanitize(warnings),
    message = sanitize(messages),
    elapsed_seconds = elapsed,
    fit = fit
  )
}

run_or_load_susie <- function(dataset, label, cache_path) {
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }
  result <- run_susie(dataset, label)
  saveRDS(result, cache_path, compress = "xz")
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- args$repo_root
.libPaths(c(file.path(repo_root, "env", "R-library"), .libPaths()))

if (!requireNamespace("coloc", quietly = TRUE) || !requireNamespace("susieR", quietly = TRUE)) {
  stop("coloc and susieR must be installed in the registered R library")
}

genes <- data.table(
  gene_symbol = c("AKR1C1", "AKR1C2"),
  ensembl_gene_id = c("ENSG00000187134", "ENSG00000151632")
)
outcome_id <- "M13_INTERVERTEB"
outcome_cases <- 53270L
outcome_controls <- 347768L

processed_dir <- file.path(repo_root, "data", "processed", "idd", "genetics")
ld_source_prefix <- file.path(repo_root, "data", "raw", "idd", "genetics", "ld", "1000G_EUR")
plink <- file.path(repo_root, "tools", "idd", "plink", "plink.exe")
work_dir <- file.path(processed_dir, "coloc_susie_akr1c_locus")
qc_dir <- file.path(repo_root, "manuscript", "qc")
summary_path <- file.path(qc_dir, "coloc_susie_summary.tsv")
pair_path <- file.path(qc_dir, "coloc_susie_signal_pairs.tsv")
handoff_path <- file.path(repo_root, "HANDOFF", "COLOC_SUSIE.tsv")
comparison_path <- file.path(repo_root, "HANDOFF", "COLOC_COMPARISON.tsv")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

full_eqtl <- fread(file.path(processed_dir, "candidate_full_eqtl.tsv.gz"), na.strings = c("", "NA"))
full_eqtl <- full_eqtl[ensembl_gene_id %chin% genes$ensembl_gene_id]
setorder(full_eqtl, ensembl_gene_id, snp, pvalue)
full_eqtl <- full_eqtl[, .SD[1L], by = .(ensembl_gene_id, snp)]

make_gene_table <- function(gene_id, label) {
  result <- full_eqtl[ensembl_gene_id == gene_id]
  setnames(
    result,
    c("beta", "se", "n", "effect_allele", "other_allele", "pos_b37", "gene_pos_b37"),
    paste0(c("beta_", "se_", "n_", "effect_allele_", "other_allele_", "pos_b37_", "gene_pos_b37_"), label)
  )
  result[, c("ensembl_gene_id", "gene_symbol") := NULL]
  result
}

akr1c1 <- make_gene_table(genes[gene_symbol == "AKR1C1", ensembl_gene_id], "akr1c1")
akr1c2 <- make_gene_table(genes[gene_symbol == "AKR1C2", ensembl_gene_id], "akr1c2")
outcome <- fread(
  file.path(processed_dir, paste0(outcome_id, ".candidate_variants.tsv.gz")),
  na.strings = c("", "NA")
)
setorder(outcome, snp, pvalue)
outcome <- outcome[, .SD[1L], by = snp]
setnames(
  outcome,
  c("beta", "se", "effect_allele", "other_allele", "pos_b38"),
  c("beta_outcome", "se_outcome", "effect_allele_outcome", "other_allele_outcome", "pos_b38_outcome")
)

bim <- fread(
  paste0(ld_source_prefix, ".bim"),
  header = FALSE,
  col.names = c("chr", "snp", "cm", "pos_b37_ld", "a1", "a2")
)

shared <- Reduce(
  function(left, right) merge(left, right, by = "snp", all = FALSE, sort = FALSE),
  list(akr1c1, akr1c2, outcome, bim)
)
shared_initial <- nrow(shared)
shared <- shared[
  !is_palindromic(a1, a2) &
    is.finite(beta_akr1c1) & is.finite(se_akr1c1) & se_akr1c1 > 0 &
    is.finite(beta_akr1c2) & is.finite(se_akr1c2) & se_akr1c2 > 0 &
    is.finite(beta_outcome) & is.finite(se_outcome) & se_outcome > 0
]
palindromic_or_missing_excluded <- shared_initial - nrow(shared)

shared[, beta_akr1c1_aligned := mapply(
  align_beta_to_a1,
  beta_akr1c1,
  effect_allele_akr1c1,
  other_allele_akr1c1,
  a1,
  a2
)]
shared[, beta_akr1c2_aligned := mapply(
  align_beta_to_a1,
  beta_akr1c2,
  effect_allele_akr1c2,
  other_allele_akr1c2,
  a1,
  a2
)]
shared[, beta_outcome_aligned := mapply(
  align_beta_to_a1,
  beta_outcome,
  effect_allele_outcome,
  other_allele_outcome,
  a1,
  a2
)]
allele_mismatch_excluded <- shared[
  !is.finite(beta_akr1c1_aligned) |
    !is.finite(beta_akr1c2_aligned) |
    !is.finite(beta_outcome_aligned),
  .N
]
shared <- shared[
  is.finite(beta_akr1c1_aligned) &
    is.finite(beta_akr1c2_aligned) &
    is.finite(beta_outcome_aligned)
]

extract_path <- file.path(work_dir, "common_nonpalindromic_snps.txt")
fwrite(shared[, .(snp)], extract_path, col.names = FALSE, quote = FALSE)

subset_prefix <- file.path(work_dir, "1000G_EUR_AKR1C_locus")
run_command(
  plink,
  c(
    "--bfile", ld_source_prefix,
    "--extract", extract_path,
    "--keep-allele-order",
    "--make-bed",
    "--out", subset_prefix
  ),
  paste0(subset_prefix, ".makebed.stdout.log"),
  paste0(subset_prefix, ".makebed.stderr.log")
)

subset_bim <- fread(
  paste0(subset_prefix, ".bim"),
  header = FALSE,
  col.names = c("chr", "snp", "cm", "pos_b37_ld", "a1_subset", "a2_subset")
)
shared <- shared[match(subset_bim$snp, shared$snp)]
if (anyNA(shared$snp) || !identical(shared$snp, subset_bim$snp)) {
  stop("PLINK subset variant order could not be matched to the summary statistics")
}

ld_prefix <- file.path(work_dir, "1000G_EUR_AKR1C_locus_r")
run_command(
  plink,
  c(
    "--bfile", subset_prefix,
    "--r", "square", "bin4",
    "--out", ld_prefix
  ),
  paste0(ld_prefix, ".stdout.log"),
  paste0(ld_prefix, ".stderr.log")
)

ld_candidates <- c(paste0(ld_prefix, ".ld.bin"), paste0(ld_prefix, ".ld.bin4"))
ld_path <- ld_candidates[file.exists(ld_candidates)][1L]
if (is.na(ld_path)) {
  stop("PLINK did not create the expected binary LD matrix")
}

variant_count <- nrow(subset_bim)
expected_bytes <- as.double(variant_count) * as.double(variant_count) * 4
actual_bytes <- file.info(ld_path)$size
if (!identical(as.double(actual_bytes), expected_bytes)) {
  stop("Unexpected binary LD size: expected ", expected_bytes, " bytes, observed ", actual_bytes)
}

connection <- file(ld_path, open = "rb")
ld_values <- readBin(connection, what = numeric(), n = variant_count * variant_count, size = 4L)
close(connection)
ld <- matrix(ld_values, nrow = variant_count, ncol = variant_count, byrow = TRUE)
dimnames(ld) <- list(subset_bim$snp, subset_bim$snp)

bad_ld <- which(rowSums(!is.finite(ld)) > 0L | !is.finite(diag(ld)))
if (length(bad_ld)) {
  keep <- setdiff(seq_len(nrow(ld)), bad_ld)
  ld <- ld[keep, keep, drop = FALSE]
  shared <- shared[keep]
}
ld_nonfinite_excluded <- length(bad_ld)

if (nrow(ld) < 50L) {
  stop("Fewer than 50 variants remained after LD and allele quality control")
}
if (max(abs(ld - t(ld))) > 1e-6) {
  stop("LD matrix is not symmetric within tolerance")
}
diag(ld) <- 1

exposure_n_akr1c1 <- as.integer(round(stats::median(shared$n_akr1c1)))
exposure_n_akr1c2 <- as.integer(round(stats::median(shared$n_akr1c2)))
outcome_n <- outcome_cases + outcome_controls

make_dataset <- function(beta, se, n, type, sd_y = NULL) {
  result <- list(
    beta = beta,
    varbeta = se^2,
    type = type,
    snp = shared$snp,
    position = shared$pos_b37_ld,
    N = n,
    LD = ld
  )
  if (!is.null(sd_y)) {
    result$sdY <- sd_y
  }
  result
}

dataset_outcome <- make_dataset(
  shared$beta_outcome_aligned,
  shared$se_outcome,
  outcome_n,
  "cc"
)
dataset_akr1c1 <- make_dataset(
  shared$beta_akr1c1_aligned,
  shared$se_akr1c1,
  exposure_n_akr1c1,
  "quant",
  1
)
dataset_akr1c2 <- make_dataset(
  shared$beta_akr1c2_aligned,
  shared$se_akr1c2,
  exposure_n_akr1c2,
  "quant",
  1
)

stage_log <- file.path(work_dir, "run_stages.log")
writeLines(paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "starting outcome SuSiE"), stage_log)
outcome_susie <- run_or_load_susie(
  dataset_outcome,
  "FinnGen outcome SuSiE",
  file.path(work_dir, "finngen_outcome_susie.rds")
)
write(
  paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "outcome SuSiE", outcome_susie$status),
  stage_log,
  append = TRUE
)
gene_runs <- list(
  AKR1C1 = run_or_load_susie(
    dataset_akr1c1,
    "AKR1C1 eQTL SuSiE",
    file.path(work_dir, "akr1c1_eqtl_susie.rds")
  ),
  AKR1C2 = run_or_load_susie(
    dataset_akr1c2,
    "AKR1C2 eQTL SuSiE",
    file.path(work_dir, "akr1c2_eqtl_susie.rds")
  )
)
write(
  paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    "gene SuSiE",
    paste(names(gene_runs), vapply(gene_runs, `[[`, character(1L), "status"), collapse = ";")
  ),
  stage_log,
  append = TRUE
)

summary_rows <- list()
pair_rows <- list()
for (current_gene in names(gene_runs)) {
  gene_run <- gene_runs[[current_gene]]
  result_status <- "not_estimable"
  result_reason <- ""
  signal_pairs <- 0L
  max_pp4 <- NA_real_
  pp4_ge_0_8_pairs <- 0L
  expression_signals <- if (is.null(gene_run$fit$sets$cs)) 0L else length(gene_run$fit$sets$cs)
  outcome_signals <- if (is.null(outcome_susie$fit$sets$cs)) 0L else length(outcome_susie$fit$sets$cs)

  if (gene_run$status != "success") {
    result_status <- "failed"
    result_reason <- gene_run$reason
  } else if (outcome_susie$status != "success") {
    result_status <- "failed"
    result_reason <- outcome_susie$reason
  } else if (expression_signals == 0L || outcome_signals == 0L) {
    result_status <- "not_estimable"
    result_reason <- "SuSiE produced no credible set for at least one trait"
  } else {
    coloc_fit <- tryCatch(
      coloc::coloc.susie(
        gene_run$fit,
        outcome_susie$fit,
        p1 = 1e-4,
        p2 = 1e-4,
        p12 = 1e-5
      ),
      error = function(condition) condition
    )
    if (inherits(coloc_fit, "condition")) {
      result_status <- "failed"
      result_reason <- paste0("coloc.susie: ", conditionMessage(coloc_fit))
    } else if (!is.list(coloc_fit) || is.null(coloc_fit$summary) || !nrow(coloc_fit$summary)) {
      result_status <- "not_estimable"
      result_reason <- "coloc.susie returned no signal-pair summary"
    } else {
      pairs <- as.data.table(coloc_fit$summary)
      current_gene_id <- genes[gene_symbol == current_gene, ensembl_gene_id]
      pairs[, `:=`(
        gene_symbol = current_gene,
        ensembl_gene_id = current_gene_id,
        outcome_id = outcome_id
      )]
      pair_rows[[current_gene]] <- pairs
      signal_pairs <- nrow(pairs)
      max_pp4 <- max(pairs[["PP.H4.abf"]], na.rm = TRUE)
      pp4_ge_0_8_pairs <- sum(pairs[["PP.H4.abf"]] >= 0.8, na.rm = TRUE)
      result_status <- "estimable"
    }
  }

  current_gene_id <- genes[gene_symbol == current_gene, ensembl_gene_id]
  summary_rows[[current_gene]] <- data.table(
    gene_symbol = current_gene,
    ensembl_gene_id = current_gene_id,
    outcome_id = outcome_id,
    method = "coloc::runsusie + coloc::coloc.susie",
    coloc_version = as.character(utils::packageVersion("coloc")),
    susieR_version = as.character(utils::packageVersion("susieR")),
    ld_reference = "1000 Genomes Phase 3 EUR, 503 participants, GRCh37",
    locus_definition = "Common AKR1C1/AKR1C2 eQTLGen cis variants in both registered +/-1 Mb windows and FinnGen M13_INTERVERTEB",
    shared_before_allele_qc = shared_initial,
    palindromic_or_missing_excluded = palindromic_or_missing_excluded,
    allele_mismatch_excluded = allele_mismatch_excluded,
    ld_nonfinite_excluded = ld_nonfinite_excluded,
    variants_analyzed = nrow(ld),
    expression_n = if (current_gene == "AKR1C1") exposure_n_akr1c1 else exposure_n_akr1c2,
    expression_n_min = if (current_gene == "AKR1C1") min(shared$n_akr1c1) else min(shared$n_akr1c2),
    expression_n_max = if (current_gene == "AKR1C1") max(shared$n_akr1c1) else max(shared$n_akr1c2),
    outcome_n = outcome_n,
    requested_p1 = 1e-4,
    requested_p2 = 1e-4,
    requested_p12 = 1e-5,
    expression_susie_status = gene_run$status,
    outcome_susie_status = outcome_susie$status,
    expression_credible_sets = expression_signals,
    outcome_credible_sets = outcome_signals,
    signal_pairs = signal_pairs,
    max_pp4 = max_pp4,
    pp4_ge_0_8_pairs = pp4_ge_0_8_pairs,
    status = result_status,
    reason = result_reason,
    expression_warning = gene_run$warning,
    outcome_warning = outcome_susie$warning,
    expression_elapsed_seconds = gene_run$elapsed_seconds,
    outcome_elapsed_seconds = outcome_susie$elapsed_seconds
  )
}

summary_table <- rbindlist(summary_rows, fill = TRUE)
pair_table <- if (length(pair_rows)) {
  rbindlist(pair_rows, fill = TRUE)
} else {
  data.table(
    gene_symbol = character(),
    ensembl_gene_id = character(),
    outcome_id = character(),
    idx1 = integer(),
    idx2 = integer(),
    nsnps = integer(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric()
  )
}
fwrite(summary_table, summary_path, sep = "\t", na = "not_estimable", quote = FALSE)
fwrite(summary_table, handoff_path, sep = "\t", na = "not_estimable", quote = FALSE)
fwrite(pair_table, pair_path, sep = "\t", na = "not_estimable", quote = FALSE)

comparison <- fread(comparison_path, na.strings = "not_estimable")
comparison[, `:=`(
  susie_attempted = FALSE,
  susie_status = "not_attempted",
  susie_variants_analyzed = NA_integer_,
  susie_expression_credible_sets = NA_integer_,
  susie_outcome_credible_sets = NA_integer_,
  susie_signal_pairs = NA_integer_,
  susie_max_pp4 = NA_real_,
  susie_pp4_ge_0_8_pairs = NA_integer_,
  susie_reason = ""
)]
for (current_gene in summary_table$gene_symbol) {
  source <- summary_table[gene_symbol == current_gene]
  comparison[
    gene_symbol == current_gene & outcome_id == "M13_INTERVERTEB",
    `:=`(
      susie_attempted = TRUE,
      susie_status = source$status,
      susie_variants_analyzed = source$variants_analyzed,
      susie_expression_credible_sets = source$expression_credible_sets,
      susie_outcome_credible_sets = source$outcome_credible_sets,
      susie_signal_pairs = source$signal_pairs,
      susie_max_pp4 = source$max_pp4,
      susie_pp4_ge_0_8_pairs = source$pp4_ge_0_8_pairs,
      susie_reason = source$reason
    )
  ]
}
fwrite(comparison, comparison_path, sep = "\t", na = "not_estimable", quote = FALSE)
write(
  paste(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "outputs written"),
  stage_log,
  append = TRUE
)

cat("coloc version:", as.character(utils::packageVersion("coloc")), "\n")
cat("susieR version:", as.character(utils::packageVersion("susieR")), "\n")
cat("variants analyzed:", nrow(ld), "\n")
print(summary_table)
