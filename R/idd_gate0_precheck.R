#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
repo_root <- "."
if (length(args) == 2L && identical(args[[1L]], "--repo-root")) repo_root <- args[[2L]]
source(file.path(repo_root, "R", "idd_utils.R"))
idd_require(c("Biobase", "GEOquery", "limma", "jsonlite", "AnnotationDbi", "org.Hs.eg.db"))
set.seed(IDD_SEED)

primary <- idd_read_geo(repo_root, "GSE70362")
primary_metadata <- idd_parse_gse70362(primary)
primary_expression <- Biobase::exprs(primary)
counts <- list()
tables <- list()
for (tissue in c("AF", "NP")) {
  selected <- primary_metadata$tissue == tissue
  result <- idd_limma_two_group(
    primary_expression[, selected, drop = FALSE],
    primary_metadata[selected, , drop = FALSE],
    block = primary_metadata$donor[selected],
    batch = TRUE
  )
  table <- result$table
  annotation <- idd_feature_annotation(primary, "GSE70362")
  table <- merge(table, annotation, by = "feature_id", all.x = TRUE, sort = FALSE)
  table$dataset <- "GSE70362"
  table$tissue <- tissue
  tables[[tissue]] <- table
  counts[[tissue]] <- list(
    samples = sum(selected),
    low = sum(primary_metadata$group[selected] == "low"),
    high = sum(primary_metadata$group[selected] == "high"),
    fdr = sum(table$pass_fdr),
    fdr_effect = sum(table$pass_primary),
    nominal_effect = sum(table$pass_nominal_effect),
    duplicate_correlation = unname(result$consensus_correlation)
  )
}

support <- idd_read_geo(repo_root, "GSE56081")
support_metadata <- Biobase::pData(support)
support_design <- data.frame(
  group = ifelse(grepl("Degenerated", support_metadata$`condition:ch1`, ignore.case = TRUE), "high", "low"),
  batch = factor("single"),
  stringsAsFactors = FALSE
)
support_result <- idd_limma_two_group(Biobase::exprs(support), support_design, batch = FALSE)
support_table <- support_result$table
support_annotation <- idd_feature_annotation(support, "GSE56081")
support_table <- merge(support_table, support_annotation, by = "feature_id", all.x = TRUE, sort = FALSE)
support_table$dataset <- "GSE56081"
support_table$tissue <- "NP"

summary <- list(
  seed = IDD_SEED,
  decision = "GO",
  primary = counts,
  supporting_gse56081 = list(
    samples = ncol(Biobase::exprs(support)),
    low = sum(support_design$group == "low"),
    high = sum(support_design$group == "high"),
    fdr = sum(support_table$pass_fdr),
    fdr_effect = sum(support_table$pass_primary),
    nominal_effect = sum(support_table$pass_nominal_effect),
    annotated_gene_features = sum(!is.na(support_table$gene_symbol)),
    limitation = "control cadaver versus surgical case source is perfectly confounded; platform annotation is incomplete"
  ),
  interpretation = "GO because public bulk, scRNA and GWAS sources are available and GSE70362 AF has FDR-controlled signal"
)

output_dir <- file.path(repo_root, "logs", "idd")
idd_atomic_json(summary, file.path(output_dir, "gate0_summary.json"))
idd_atomic_tsv(do.call(rbind, tables), file.path(output_dir, "gate0_gse70362_deg.tsv"))
idd_atomic_tsv(support_table, file.path(output_dir, "gate0_gse56081_deg.tsv"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
