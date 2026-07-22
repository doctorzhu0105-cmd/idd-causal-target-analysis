#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
if (length(args) == 2L && identical(args[[1L]], "--repo-root")) repo_root <- args[[2L]]

source(file.path(repo_root, "R", "idd_utils.R"))
local_library <- normalizePath(
  file.path(repo_root, "env", "R-library"),
  winslash = "/",
  mustWork = TRUE
)
.libPaths(c(local_library, .libPaths()))
idd_require(c("Biobase", "GEOquery", "limma", "WGCNA", "jsonlite", "AnnotationDbi", "org.Hs.eg.db"))
set.seed(IDD_SEED)
options(stringsAsFactors = FALSE)
WGCNA::disableWGCNAThreads()
cor <- WGCNA::cor

config <- jsonlite::fromJSON(
  file.path(repo_root, "config", "idd_analysis.json"),
  simplifyVector = TRUE
)
output_dir <- file.path(repo_root, "manuscript", "supplement")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primary <- idd_read_geo(repo_root, "GSE70362")
metadata <- idd_parse_gse70362(primary)
expression <- Biobase::exprs(primary)
feature_annotation <- idd_feature_annotation(primary, "GSE70362")
gene_expression <- idd_collapse_gene_expression(expression, feature_annotation)

run_signed_wgcna <- function(tissue) {
  selected <- metadata$tissue == tissue
  tissue_metadata <- metadata[selected, , drop = FALSE]
  tissue_expression <- gene_expression[, selected, drop = FALSE]
  variances <- apply(tissue_expression, 1L, stats::var)
  feature_count <- min(as.integer(config$wgcna$top_variance_features), length(variances))
  selected_genes <- names(sort(variances, decreasing = TRUE))[seq_len(feature_count)]
  dat_expr <- t(tissue_expression[selected_genes, , drop = FALSE])
  good <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  idd_assert(good$allOK, "%s signed WGCNA input failed goodSamplesGenes", tissue)

  powers <- as.numeric(config$wgcna$powers)
  soft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = "signed",
    verbose = 0
  )
  fit <- soft$fitIndices
  eligible <- fit$Power[fit$SFT.R.sq >= as.numeric(config$wgcna$scale_free_fit_min)]
  power <- if (length(eligible) > 0L) min(eligible) else fit$Power[which.max(fit$SFT.R.sq)]
  network <- WGCNA::blockwiseModules(
    dat_expr,
    power = power,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = as.integer(config$wgcna$min_module_size),
    mergeCutHeight = as.numeric(config$wgcna$merge_cut_height),
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    maxBlockSize = 2500L,
    useCorOptionsThroughout = FALSE,
    verbose = 0,
    randomSeed = IDD_SEED
  )
  module_labels <- as.character(network$colors)
  module_colors <- WGCNA::labels2colors(network$colors)
  eigengenes <- WGCNA::orderMEs(network$MEs)
  trait <- as.numeric(tissue_metadata$group == "high")
  correlation <- stats::cor(eigengenes, trait, use = "pairwise.complete.obs")
  correlation_p <- WGCNA::corPvalueStudent(correlation, nrow(dat_expr))
  module_fit <- idd_limma_two_group(
    t(eigengenes),
    tissue_metadata,
    block = tissue_metadata$donor,
    batch = TRUE
  )$table
  module_ids <- sub("^ME", "", rownames(correlation))
  module_summary <- data.frame(
    tissue = tissue,
    module = module_ids,
    correlation = as.numeric(correlation),
    correlation_p_value = as.numeric(correlation_p),
    adjusted_log2fc = module_fit[rownames(correlation), "logFC"],
    adjusted_p_value = module_fit[rownames(correlation), "P.Value"],
    size = as.integer(table(module_labels)[module_ids]),
    soft_threshold_power = power,
    network_type = "signed",
    tom_type = "signed",
    stringsAsFactors = FALSE
  )
  module_summary$non_grey_bh_fdr <- NA_real_
  interpretable <- module_summary$module != "0"
  module_summary$non_grey_bh_fdr[interpretable] <- stats::p.adjust(
    module_summary$adjusted_p_value[interpretable],
    method = "BH"
  )
  module_summary$interpretation <- ifelse(
    !interpretable,
    "grey_not_biological_module",
    ifelse(module_summary$non_grey_bh_fdr < 0.05, "fdr_significant", "exploratory")
  )

  membership <- stats::cor(dat_expr, eigengenes, use = "pairwise.complete.obs")
  gene_significance <- stats::cor(dat_expr, trait, use = "pairwise.complete.obs")
  eigengene_column <- setNames(
    seq_along(colnames(eigengenes)),
    sub("^ME", "", colnames(eigengenes))
  )
  gene_table <- data.frame(
    tissue = tissue,
    gene_symbol = colnames(dat_expr),
    module = module_labels,
    module_color = module_colors,
    module_membership = vapply(
      seq_along(module_colors),
      function(index) membership[index, eigengene_column[[module_labels[[index]]]]],
      numeric(1L)
    ),
    gene_significance = as.numeric(gene_significance),
    network_type = "signed",
    stringsAsFactors = FALSE
  )
  list(modules = module_summary, genes = gene_table, soft = data.frame(tissue = tissue, fit))
}

results <- lapply(c("AF", "NP"), run_signed_wgcna)
modules <- do.call(rbind, lapply(results, `[[`, "modules"))
genes <- do.call(rbind, lapply(results, `[[`, "genes"))
soft <- do.call(rbind, lapply(results, `[[`, "soft"))

idd_atomic_tsv(modules, file.path(output_dir, "wgcna_signed_modules.tsv"))
idd_atomic_tsv(genes, file.path(output_dir, "wgcna_signed_genes.tsv"))
idd_atomic_tsv(soft, file.path(output_dir, "wgcna_signed_soft_threshold.tsv"))

af_non_grey <- modules[modules$tissue == "AF" & modules$module != "0", , drop = FALSE]
np_non_grey <- modules[modules$tissue == "NP" & modules$module != "0", , drop = FALSE]
summary <- list(
  seed = IDD_SEED,
  network_type = "signed",
  tom_type = "signed",
  af_soft_threshold_power = unique(modules$soft_threshold_power[modules$tissue == "AF"]),
  np_soft_threshold_power = unique(modules$soft_threshold_power[modules$tissue == "NP"]),
  af_non_grey_module_count = nrow(af_non_grey),
  np_non_grey_module_count = nrow(np_non_grey),
  af_fdr_significant_module_count = sum(af_non_grey$non_grey_bh_fdr < 0.05),
  np_fdr_significant_module_count = sum(np_non_grey$non_grey_bh_fdr < 0.05),
  af_top_module = af_non_grey$module[which.min(af_non_grey$adjusted_p_value)],
  af_top_module_p = min(af_non_grey$adjusted_p_value),
  af_top_module_fdr = af_non_grey$non_grey_bh_fdr[which.min(af_non_grey$adjusted_p_value)]
)
idd_atomic_json(
  summary,
  file.path(output_dir, "wgcna_signed_summary.json"),
  digits = 16L
)
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = 16L), "\n")
