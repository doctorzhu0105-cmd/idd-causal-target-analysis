#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
if (length(args) == 2L && identical(args[[1L]], "--repo-root")) repo_root <- args[[2L]]
source(file.path(repo_root, "R", "idd_utils.R"))
idd_require(c("Biobase", "GEOquery", "limma", "WGCNA", "jsonlite", "AnnotationDbi", "org.Hs.eg.db", "digest"))
set.seed(IDD_SEED)
options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()
cor <- WGCNA::cor

config <- jsonlite::fromJSON(file.path(repo_root, "config", "idd_analysis.json"), simplifyVector = TRUE)
output_dir <- file.path(repo_root, "artifacts", "bulk")
log_dir <- file.path(repo_root, "logs", "idd")
figure_dir <- file.path(repo_root, "figures", "idd")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

primary <- idd_read_geo(repo_root, "GSE70362")
metadata <- idd_parse_gse70362(primary)
expression <- Biobase::exprs(primary)
feature_annotation <- idd_feature_annotation(primary, "GSE70362")
gene_expression <- idd_collapse_gene_expression(expression, feature_annotation)
druggable <- utils::read.delim(file.path(repo_root, "data", "processed", "mr", "druggable_genes.tsv"), check.names = FALSE)
druggable_symbols <- unique(druggable$hgnc_names[nzchar(druggable$hgnc_names)])

run_wgcna <- function(tissue) {
  selected <- metadata$tissue == tissue
  tissue_metadata <- metadata[selected, , drop = FALSE]
  tissue_expression <- gene_expression[, selected, drop = FALSE]
  variances <- apply(tissue_expression, 1L, stats::var)
  feature_count <- min(as.integer(config$wgcna$top_variance_features), length(variances))
  selected_genes <- names(sort(variances, decreasing = TRUE))[seq_len(feature_count)]
  dat_expr <- t(tissue_expression[selected_genes, , drop = FALSE])
  good <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  idd_assert(good$allOK, "%s WGCNA input failed goodSamplesGenes", tissue)

  powers <- as.numeric(config$wgcna$powers)
  soft <- WGCNA::pickSoftThreshold(
    dat_expr,
    powerVector = powers,
    networkType = config$wgcna$network_type,
    verbose = 0
  )
  fit <- soft$fitIndices
  eligible <- fit$Power[fit$SFT.R.sq >= as.numeric(config$wgcna$scale_free_fit_min)]
  power <- if (length(eligible) > 0L) min(eligible) else fit$Power[which.max(fit$SFT.R.sq)]
  network <- WGCNA::blockwiseModules(
    dat_expr,
    power = power,
    TOMType = config$wgcna$network_type,
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
  pvalue <- WGCNA::corPvalueStudent(correlation, nrow(dat_expr))
  module_fit <- idd_limma_two_group(
    t(eigengenes),
    tissue_metadata,
    block = tissue_metadata$donor,
    batch = TRUE
  )$table
  module_summary <- data.frame(
    tissue = tissue,
    module = sub("^ME", "", rownames(correlation)),
    correlation = as.numeric(correlation),
    correlation_p_value = as.numeric(pvalue),
    adjusted_log2fc = module_fit[rownames(correlation), "logFC"],
    adjusted_p_value = module_fit[rownames(correlation), "P.Value"],
    adjusted_fdr = module_fit[rownames(correlation), "adj.P.Val"],
    size = as.integer(table(module_labels)[sub("^ME", "", rownames(correlation))]),
    stringsAsFactors = FALSE
  )
  membership <- stats::cor(dat_expr, eigengenes, use = "pairwise.complete.obs")
  gene_significance <- stats::cor(dat_expr, trait, use = "pairwise.complete.obs")
  eigengene_column <- setNames(seq_along(colnames(eigengenes)), sub("^ME", "", colnames(eigengenes)))
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
    is_druggable = colnames(dat_expr) %in% druggable_symbols,
    stringsAsFactors = FALSE
  )
  list(
    summary = module_summary,
    genes = gene_table,
    soft = data.frame(tissue = tissue, fit, stringsAsFactors = FALSE),
    power = power
  )
}

deg_tables <- list()
wgcna_results <- list()
for (tissue in c("AF", "NP")) {
  selected <- metadata$tissue == tissue
  result <- idd_limma_two_group(
    expression[, selected, drop = FALSE],
    metadata[selected, , drop = FALSE],
    block = metadata$donor[selected],
    batch = TRUE
  )
  table <- merge(result$table, feature_annotation, by = "feature_id", all.x = TRUE, sort = FALSE)
  table$dataset <- "GSE70362"
  table$tissue <- tissue
  table$is_druggable <- table$gene_symbol %in% druggable_symbols
  deg_tables[[tissue]] <- table
  wgcna_results[[tissue]] <- run_wgcna(tissue)
}

deg <- do.call(rbind, deg_tables)
module_summary <- do.call(rbind, lapply(wgcna_results, `[[`, "summary"))
module_genes <- do.call(rbind, lapply(wgcna_results, `[[`, "genes"))
soft_threshold <- do.call(rbind, lapply(wgcna_results, `[[`, "soft"))

associated_modules <- module_summary$module[
  module_summary$tissue == "AF" & module_summary$module != "0" & module_summary$adjusted_fdr < 0.05
]
if (length(associated_modules) == 0L) {
  best <- module_summary[module_summary$tissue == "AF" & module_summary$module != "grey", , drop = FALSE]
  associated_modules <- best$module[which.max(abs(best$correlation))]
  module_selection_rule <- "strongest_non_grey_AF_module_because_no_module_passed_BH_FDR_0.05"
} else {
  module_selection_rule <- "AF_modules_passing_donor_blocked_batch_adjusted_BH_FDR_0.05"
}

candidate_pool <- merge(
  deg[deg$tissue == "AF" & deg$pass_primary & !is.na(deg$gene_symbol), ],
  module_genes[module_genes$tissue == "AF", ],
  by = "gene_symbol",
  suffixes = c("_deg", "_wgcna")
)
candidate_pool$in_selected_module <- candidate_pool$module %in% associated_modules
candidate_pool <- merge(candidate_pool, druggable, by.x = "gene_symbol", by.y = "hgnc_names")
candidate_pool$bulk_priority_score <-
  -log10(pmax(candidate_pool$adj.P.Val, .Machine$double.xmin)) +
  abs(candidate_pool$module_membership) +
  abs(candidate_pool$gene_significance) +
  as.numeric(candidate_pool$in_selected_module)
candidate_pool <- candidate_pool[order(-candidate_pool$bulk_priority_score, candidate_pool$gene_symbol), ]

idd_atomic_tsv(deg, file.path(output_dir, "gse70362_deg.tsv"))
idd_atomic_tsv(module_summary, file.path(output_dir, "gse70362_wgcna_modules.tsv"))
idd_atomic_tsv(module_genes, file.path(output_dir, "gse70362_wgcna_genes.tsv"))
idd_atomic_tsv(soft_threshold, file.path(output_dir, "gse70362_wgcna_soft_threshold.tsv"))
idd_atomic_tsv(candidate_pool, file.path(output_dir, "druggable_bulk_candidates.tsv"))
idd_atomic_tsv(metadata, file.path(output_dir, "gse70362_sample_metadata.tsv"))

summary <- list(
  seed = IDD_SEED,
  primary_dataset = "GSE70362",
  deg_counts = lapply(c("AF", "NP"), function(tissue) {
    table <- deg[deg$tissue == tissue, ]
    list(tissue = tissue, fdr = sum(table$pass_fdr), fdr_effect = sum(table$pass_primary))
  }),
  wgcna_power = lapply(names(wgcna_results), function(tissue) list(tissue = tissue, power = wgcna_results[[tissue]]$power)),
  selected_af_modules = associated_modules,
  module_selection_rule = module_selection_rule,
  druggable_bulk_candidate_count = nrow(candidate_pool),
  repeat_combat_performed = FALSE
)
idd_atomic_json(summary, file.path(log_dir, "bulk_summary.json"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
