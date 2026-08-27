#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table); library(jsonlite) })

parse_args <- function(args) {
  out <- list(repo_root = ".", pairs = NULL, output = NULL)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--repo-root", "--pairs", "--output") && i < length(args)) {
      out[[gsub("-", "_", sub("^--", "", key))]] <- args[[i + 1L]]; i <- i + 2L
    } else stop("Unknown or incomplete argument: ", key)
  }
  if (is.null(out$pairs) || is.null(out$output)) stop("--pairs and --output are required")
  out$repo_root <- gsub("\\\\", "/", out$repo_root)
  out$pairs <- gsub("\\\\", "/", out$pairs)
  out$output <- gsub("\\\\", "/", out$output)
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo_root <- args$repo_root
cfg <- fromJSON(file.path(repo_root, "config", "coloc_smr_authoritative_raw.json"), simplifyVector = TRUE)
expected_root <- normalizePath(cfg$execution$required_repo_root, winslash = "/", mustWork = TRUE)
actual_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
if (!identical(actual_root, expected_root)) stop("Runtime root does not match raw authoritative config")
if (as.integer(cfg$coloc$cis_window_bp) != 1000000L) stop("Raw coloc window is not +/-1 Mb")
if (!identical(as.logical(cfg$coloc$outcome_input_policy$depends_on_eqtl_selection), FALSE)) stop("Outcome input depends on eQTL selection")
if (!identical(as.logical(cfg$coloc$outcome_input_policy$eQTL_SNP_whitelist), FALSE)) stop("Outcome input uses an eQTL SNP whitelist")
local_library <- if (dir.exists("D:/idd_r_lib_ascii")) "D:/idd_r_lib_ascii" else file.path(dirname(repo_root), "env", "R-library")
.libPaths(c(local_library, .libPaths()))
if (!requireNamespace("coloc", quietly = TRUE)) { library(coloc, lib.loc = local_library) }
coloc_version <- as.character(utils::packageVersion("coloc"))
if (!identical(coloc_version, as.character(cfg$coloc$package_version_required))) stop("Unexpected coloc version: ", coloc_version)

pairs <- fread(args$pairs, na.strings = c("", "NA"))
eqtl_path <- file.path(repo_root, cfg$inputs$eqtl_input)
eqtl <- fread(eqtl_path, na.strings = c("", "NA"))
setorderv(eqtl, c("ensembl_gene_id", "snp", "pvalue"))
eqtl <- eqtl[, .SD[1L], by = .(ensembl_gene_id, snp)]
counts <- data.table(outcome_id = c("M13_INTERVERTEB", "M13_LOWBACKPAIN", "M13_SCIATICA"), cases = c(53270L,45696L,28094L), controls = c(347768L,347768L,347768L))
# Registered p12 is 1e-5.  This fixed grid was set before evaluating PP4 and
# stays within p12 <= min(p1,p2) for p1=p2=1e-4.
p12_grid <- c(1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4)

run_one <- function(pair, p12) {
  gene <- pair$gene_symbol; gene_id <- pair$ensembl_gene_id; outcome_id <- pair$outcome_id
  outcome_rel <- cfg$outputs$outcome_region_template
  outcome_rel <- gsub("\\{gene\\}", gene, outcome_rel)
  outcome_rel <- gsub("\\{outcome\\}", outcome_id, outcome_rel)
  outcome_path <- file.path(repo_root, outcome_rel)
  oy <- fread(outcome_path, na.strings = c("", "NA"))
  setorder(oy, snp, fg_se); oy <- oy[, .SD[1L], by = snp]
  shared <- merge(eqtl[ensembl_gene_id == gene_id], oy, by = "snp", all = FALSE, sort = TRUE)
  shared <- shared[is.finite(beta) & is.finite(se) & se > 0 & is.finite(fg_beta) & is.finite(fg_se) & fg_se > 0]
  base <- data.table(gene_symbol = gene, ensembl_gene_id = gene_id, outcome_id = outcome_id, p12 = p12, shared_variants = nrow(shared), pp0=NA_real_,pp1=NA_real_,pp2=NA_real_,pp3=NA_real_,pp4=NA_real_, status="not_estimable", reason="insufficient_shared_variants", p12_registered=1e-5, p4_ge_0_80=NA, eqtl_input_path=eqtl_path, outcome_input_path=outcome_path, coloc_version=coloc_version)
  if (nrow(shared) < 50L) return(base)
  current_outcome <- outcome_id
  n <- counts[outcome_id == current_outcome, cases + controls][[1L]]; s <- counts[outcome_id == current_outcome, cases][[1L]] / n
  fit <- tryCatch(coloc::coloc.abf(dataset1=list(beta=shared$beta,varbeta=shared$se^2,type="quant",sdY=1,snp=shared$snp), dataset2=list(beta=shared$fg_beta,varbeta=shared$fg_se^2,type="cc",N=n,s=s,snp=shared$snp), p1=cfg$coloc$p1,p2=cfg$coloc$p2,p12=p12), error=function(e)e)
  if (inherits(fit,"condition")) { base$reason <- conditionMessage(fit); return(base) }
  sm <- fit$summary
  base[, `:=`(pp0=unname(sm[["PP.H0.abf"]]),pp1=unname(sm[["PP.H1.abf"]]),pp2=unname(sm[["PP.H2.abf"]]),pp3=unname(sm[["PP.H3.abf"]]),pp4=unname(sm[["PP.H4.abf"]]),status="estimable",reason="",p4_ge_0_80=unname(sm[["PP.H4.abf"]])>=0.80)]
  base
}

res <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) rbindlist(lapply(p12_grid, function(p12) run_one(pairs[i], p12)))), fill=TRUE)
dir.create(dirname(args$output), recursive=TRUE, showWarnings=FALSE)
fwrite(res, args$output, sep="\t", na="not_estimable", quote=FALSE)
cat("rows:", nrow(res), "estimable:", res[status=="estimable",.N], "\n")
