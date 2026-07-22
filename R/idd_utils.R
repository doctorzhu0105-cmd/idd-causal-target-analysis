IDD_SEED <- 20260716L

idd_abort <- function(format, ...) {
  stop(sprintf(format, ...), call. = FALSE)
}

idd_assert <- function(condition, format, ...) {
  if (!isTRUE(condition)) idd_abort(format, ...)
  invisible(TRUE)
}

idd_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  idd_assert(length(missing) == 0L, "Missing R packages: %s", paste(missing, collapse = ", "))
  invisible(TRUE)
}

idd_read_geo <- function(repo_root, accession) {
  path <- file.path(repo_root, "data", "raw", "idd", accession, paste0(accession, "_series_matrix.txt.gz"))
  idd_assert(file.exists(path), "Missing GEO matrix: %s", path)
  object <- GEOquery::getGEO(filename = path)
  idd_assert(inherits(object, "ExpressionSet"), "%s did not parse as an ExpressionSet", accession)
  object
}

idd_extract_value <- function(values, prefix) {
  trimws(sub(paste0("^", prefix, "\\s*"), "", as.character(values), ignore.case = TRUE))
}

idd_atomic_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  utils::write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  if (file.exists(path)) file.remove(path)
  idd_assert(file.rename(temporary, path), "Could not move temporary file to %s", path)
  invisible(path)
}

idd_atomic_json <- function(data, path, digits = 4L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  jsonlite::write_json(
    data,
    temporary,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = digits
  )
  write("", file = temporary, append = TRUE)
  if (file.exists(path)) file.remove(path)
  idd_assert(file.rename(temporary, path), "Could not move temporary file to %s", path)
  invisible(path)
}

idd_sha256 <- function(path) {
  idd_require("digest")
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

idd_parse_gse70362 <- function(object) {
  metadata <- Biobase::pData(object)
  result <- data.frame(
    sample_id = metadata$geo_accession,
    donor = idd_extract_value(metadata$`individual:ch1`, "donor"),
    tissue = ifelse(grepl("Annulus", metadata$`tissue:ch1`, ignore.case = TRUE), "AF", "NP"),
    grade = as.character(metadata$`thompson grade:ch1`),
    batch = factor(metadata$`batch:ch1`),
    stringsAsFactors = FALSE
  )
  result$group <- ifelse(result$grade %in% c("I", "I-II", "II"), "low", "high")
  idd_assert(identical(as.integer(table(result$tissue)), c(24L, 24L)), "Unexpected GSE70362 tissue counts")
  result
}

idd_feature_annotation <- function(object, accession) {
  annotation <- Biobase::fData(object)
  feature_id <- rownames(annotation)
  if (accession == "GSE70362") {
    entrez <- trimws(as.character(annotation$ENTREZ_GENE_ID))
    symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(entrez[nzchar(entrez)]),
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
    gene_symbol <- unname(symbol[entrez])
  } else {
    candidates <- intersect(c("Gene Symbol", "GENE_SYMBOL", "GENE_SYMBOL.1"), names(annotation))
    gene_symbol <- if (length(candidates) > 0L) as.character(annotation[[candidates[[1L]]]]) else rep(NA_character_, nrow(annotation))
    entrez_candidates <- intersect(c("ENTREZ_GENE_ID", "ENTREZGENEID"), names(annotation))
    entrez <- if (length(entrez_candidates) > 0L) as.character(annotation[[entrez_candidates[[1L]]]]) else rep(NA_character_, nrow(annotation))
  }
  gene_symbol <- trimws(sub(" ///.*$", "", gene_symbol))
  gene_symbol[!nzchar(gene_symbol)] <- NA_character_
  data.frame(feature_id = feature_id, gene_symbol = gene_symbol, entrez_id = entrez, stringsAsFactors = FALSE)
}

idd_limma_two_group <- function(expression, metadata, block = NULL, batch = FALSE) {
  group <- factor(metadata$group, levels = c("low", "high"))
  design_data <- data.frame(group = group, batch = metadata$batch)
  design <- if (batch) stats::model.matrix(~ batch + group, design_data) else stats::model.matrix(~ group, design_data)
  coefficient <- "grouphigh"
  idd_assert(qr(design)$rank == ncol(design), "Limma design is rank deficient")
  if (is.null(block)) {
    fit <- limma::lmFit(expression, design)
    consensus <- NA_real_
  } else {
    correlation <- limma::duplicateCorrelation(expression, design, block = block)
    consensus <- correlation$consensus.correlation
    fit <- limma::lmFit(expression, design, block = block, correlation = consensus)
  }
  fit <- limma::eBayes(fit)
  table <- limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  table$feature_id <- rownames(table)
  table$pass_fdr <- table$adj.P.Val < 0.05
  table$pass_primary <- table$pass_fdr & abs(table$logFC) >= 0.5
  table$pass_nominal_effect <- table$P.Value < 0.05 & abs(table$logFC) >= 0.5
  list(table = table, consensus_correlation = consensus, design = design)
}

idd_collapse_gene_expression <- function(expression, feature_annotation) {
  keep <- !is.na(feature_annotation$gene_symbol) & nzchar(feature_annotation$gene_symbol)
  expression <- expression[keep, , drop = FALSE]
  genes <- feature_annotation$gene_symbol[keep]
  variances <- apply(expression, 1L, stats::var)
  order_index <- order(genes, -variances, rownames(expression))
  selected <- order_index[!duplicated(genes[order_index])]
  collapsed <- expression[selected, , drop = FALSE]
  rownames(collapsed) <- genes[selected]
  collapsed
}
