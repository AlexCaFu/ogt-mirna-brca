# =====================================================================
# R/03_geo_foldchange.R
# Differential expression in drug-resistant breast cancer models (GEO),
# using limma. For each series, enzyme genes (mRNA series) or candidate
# microRNAs (miRNA series) are tested for resistant-vs-sensitive change.
# Output: results/03_geo_foldchanges.csv
#
# Requires a sample-group file per series at
#   results/geo_<GSE>_samples.csv   (columns: geo_accession, group_guess)
# with group_guess in {sensitive, resistant}. A template is written on the
# first run for any series that lacks one; fill it in and re-run.
# =====================================================================
source("R/config.R")
source("R/00_setup.R")

## helpers (recovered from the analysis session) ----------------------
collapse_targets <- function(feat, lfc, adjp, targets, is_mirna, n_res, n_sen, testable) {
  keys <- if (is_mirna) lapply(feat, canon_tokens) else lapply(feat, gene_tokens)
  out <- lapply(targets, function(f) {
    if (is_mirna) idx <- vapply(keys, function(k) canon_mir(f) %in% k, logical(1))
    else {
      alt <- tolower(unique(c(f, names(GENE_ALIASES)[GENE_ALIASES == f], unname(GENE_ALIASES[f]))))
      idx <- vapply(keys, function(k) any(alt %in% k), logical(1)) }
    if (!any(idx)) return(data.frame(target = f, log2FC = NA_real_, adj_p = NA_real_,
                                     n_probe_hits = 0L, n_resistant = n_res, n_sensitive = n_sen, testable = testable))
    i <- which(idx)[which.max(abs(lfc[idx]))]
    data.frame(target = f, log2FC = round(lfc[i], 3),
               adj_p = if (testable) signif(adjp[i], 3) else NA_real_,
               n_probe_hits = sum(idx), n_resistant = n_res, n_sensitive = n_sen, testable = testable)
  })
  bind_rows(out)
}

fc_two_group <- function(eset, grp, fmap, targets, is_mirna) {
  ex <- to_log2(exprs(eset))
  keep <- grp %in% c("resistant", "sensitive")
  ex <- ex[, keep, drop = FALSE]; g <- factor(grp[keep], levels = c("sensitive", "resistant"))
  n_res <- sum(g == "resistant"); n_sen <- sum(g == "sensitive")
  if (n_res < 1 || n_sen < 1) stop("need both groups present")
  testable <- n_res >= 2 && n_sen >= 2
  feat <- fmap[rownames(ex)]
  if (testable) {
    fit <- eBayes(lmFit(ex, model.matrix(~ g)))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    collapse_targets(feat[match(rownames(tt), rownames(ex))], tt$logFC, tt$adj.P.Val,
                     targets, is_mirna, n_res, n_sen, testable)
  } else {
    lfc <- rowMeans(ex[, g == "resistant", drop = FALSE]) - rowMeans(ex[, g == "sensitive", drop = FALSE])
    collapse_targets(feat, unname(lfc), rep(NA_real_, length(lfc)), targets, is_mirna, n_res, n_sen, testable)
  }
}

get_group_vector <- function(nm, eset) {
  f <- sprintf("results/geo_%s_samples.csv", nm)
  if (!file.exists(f)) {
    tmpl <- data.frame(geo_accession = colnames(eset), group_guess = NA_character_)
    readr::write_csv(tmpl, f)
    message("[", nm, "] wrote sample template ", f, " -- fill group_guess (sensitive/resistant) and re-run.")
    return(NULL)
  }
  samp <- read.csv(f, check.names = FALSE)
  grp <- setNames(samp$group_guess, samp$geo_accession); grp[colnames(eset)]
}

## download + process each series -------------------------------------
mir_targets  <- unique(CANDIDATES$mirna)
gene_targets <- ENZYMES
res <- list()

for (nm in names(GEO_SERIES)) {
  s <- GEO_SERIES[[nm]]
  message("\n[", nm, "] ", s$description)
  eset <- tryCatch(getGEO(nm, GSEMatrix = TRUE, getGPL = TRUE)[[1]],
                   error = function(e) { message("download failed: ", conditionMessage(e)); NULL })
  if (is.null(eset)) next
  is_mir <- s$measures == "miRNA"
  grp <- get_group_vector(nm, eset)
  if (is.null(grp) || all(is.na(grp))) next
  fmap <- feature_map(eset, if (is_mir) "mirna" else "gene")
  targets <- if (is_mir) mir_targets else gene_targets
  out <- tryCatch(fc_two_group(eset, grp, fmap, targets, is_mir),
                  error = function(e) { message("fc failed: ", conditionMessage(e)); NULL })
  if (!is.null(out)) { out$series <- nm; out$measure <- s$measures; res[[nm]] <- out }
}

if (length(res)) {
  geo <- bind_rows(res)
  readr::write_csv(geo, "results/03_geo_foldchanges.csv")
  message("\nWrote results/03_geo_foldchanges.csv"); print(geo)
} else {
  message("\nNo GEO fold changes computed. Fill in the results/geo_*_samples.csv templates and re-run.")
}
