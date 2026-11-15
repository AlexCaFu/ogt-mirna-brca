# =====================================================================
# R/04_geo_foldchange.R    LAYER 4b: COMPUTE fold changes  (CORRECTED)
# All three series are two-group comparisons, resistant minus sensitive.
#
# Convention: log2FC = log2(resistant) - log2(sensitive). Positive = higher in
# the resistant state. miR-200c-3p is expected NEGATIVE (lost in resistance).
#
# Design facts confirmed from the sample metadata (R/03_geo_inspect.R output):
#  - GSE26459 (mRNA, GPL570): 12 sensitive (B7) vs 12 resistant (G11OH-T),
#    balanced treatment factor. Ordinary two-group moderated t-test.
#  - GSE28267 (miRNA, GPL7722): COMMON-REFERENCE two-color array. Each array is
#    log2(cell sample / Exiqon reference pool). The resistant-vs-sensitive fold
#    change is the two-group difference (LY2 minus MCF-7); the common reference
#    cancels. 2 MCF-7 vs 4 LY2. miRBase 11.0 names (e.g. hsa-miR-200c).
#  - GSE66607 (miRNA, GPL16384): 1 parent vs 1 resistant, NO replicates ->
#    fold change only, no p-value (confirmatory direction).
# =====================================================================
source("R/config.R")
suppressMessages({library(GEOquery); library(limma); library(dplyr)})
geo_objects <- readRDS("results/geo_objects.rds")

## Group assignment confirmed from the metadata: the keyword guess is correct
## for all three series, so we accept it. (Override here if ever needed.)
GROUPS <- list(
  GSE26459 = list(accept_guess = TRUE),
  GSE28267 = list(accept_guess = TRUE),
  GSE66607 = list(accept_guess = TRUE)
)

## ---- canonical microRNA name (version-independent) ----
canon_mir <- function(x) {
  x <- tolower(x)
  x <- gsub("hsa-", "", x)
  x <- gsub("^mir-?", "", x)     # drop leading mir / miR
  x <- gsub("-[35]p$", "", x)    # drop arm suffix
  x <- gsub("\\*$", "", x)       # drop star
  x <- gsub("-[0-9]+$", "", x)   # drop locus suffix (-1, -2)
  trimws(x)
}
canon_tokens <- function(s) {
  toks <- unlist(strsplit(as.character(s), "[ ,;/|]+"))
  toks <- toks[nzchar(toks)]
  unique(canon_mir(toks))
}
gene_tokens <- function(s) {
  toks <- unlist(strsplit(as.character(s), "\\s*///\\s*|[;,]"))
  tolower(trimws(toks[nzchar(toks)]))
}

feature_map <- function(eset, kind = c("gene", "mirna")) {
  kind <- match.arg(kind); fd <- fData(eset); cn <- colnames(fd)
  pick <- function(patterns) {
    for (p in patterns) { h <- grep(p, cn, ignore.case = TRUE, value = TRUE); if (length(h)) return(h[1]) }
    NA_character_
  }
  col <- if (kind == "gene") pick(c("^gene.?symbol$", "gene.?symbol", "^symbol$"))
         else pick(c("mirna_id", "mirna", "mir_id", "transcript.?id\\(", "transcript.?id", "^name$"))
  if (is.na(col)) return(setNames(rownames(fd), rownames(fd)))
  setNames(as.character(fd[[col]]), rownames(fd))
}

to_log2 <- function(ex) if (max(ex, na.rm = TRUE) > 100) log2(ex + 1) else ex

collapse_targets <- function(feat, lfc, adjp, targets, is_mirna, n_res, n_sen, testable) {
  keys <- if (is_mirna) lapply(feat, canon_tokens) else lapply(feat, gene_tokens)
  out <- lapply(targets, function(f) {
    if (is_mirna) idx <- vapply(keys, function(k) canon_mir(f) %in% k, logical(1))
    else {
      alt <- tolower(unique(c(f, names(GENE_ALIASES)[GENE_ALIASES == f], unname(GENE_ALIASES[f]))))
      idx <- vapply(keys, function(k) any(alt %in% k), logical(1))
    }
    if (!any(idx)) return(data.frame(target = f, log2FC = NA_real_, adj_p = NA_real_,
                                     n_probe_hits = 0L, n_resistant = n_res, n_sensitive = n_sen, testable = testable))
    i <- which(idx)[which.max(abs(lfc[idx]))]
    data.frame(target = f, log2FC = round(lfc[i], 3),
               adj_p = if (testable) signif(adjp[i], 3) else NA_real_,
               n_probe_hits = sum(idx), n_resistant = n_res, n_sensitive = n_sen, testable = testable)
  })
  bind_rows(out)
}

## Universal two-group fold change (resistant - sensitive), robust to 1 vs 1.
fc_two_group <- function(eset, grp, fmap, targets, is_mirna) {
  ex <- to_log2(exprs(eset))
  keep <- grp %in% c("resistant", "sensitive")
  ex <- ex[, keep, drop = FALSE]; g <- factor(grp[keep], levels = c("sensitive", "resistant"))
  n_res <- sum(g == "resistant"); n_sen <- sum(g == "sensitive")
  if (n_res < 1 || n_sen < 1) stop("need both groups present")
  testable <- n_res >= 2 && n_sen >= 2
  feat <- fmap[rownames(ex)]
  if (testable) {
    fit <- eBayes(lmFit(ex, model.matrix(~ g)))          # coef 2 = resistant - sensitive
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    collapse_targets(feat[match(rownames(tt), rownames(ex))], tt$logFC, tt$adj.P.Val,
                     targets, is_mirna, n_res, n_sen, testable)
  } else {
    lfc <- rowMeans(ex[, g == "resistant", drop = FALSE]) - rowMeans(ex[, g == "sensitive", drop = FALSE])
    collapse_targets(feat, unname(lfc), rep(NA_real_, length(lfc)),
                     targets, is_mirna, n_res, n_sen, testable)
  }
}

get_group_vector <- function(nm, eset) {
  samp <- read.csv(sprintf("results/geo_%s_samples.csv", nm), check.names = FALSE)
  grp <- setNames(samp$group_guess, samp$geo_accession)
  grp[colnames(eset)]
}

## ---- run ----
res <- list()
for (nm in names(GEO_SERIES)) {
  s <- GEO_SERIES[[nm]]; eset <- geo_objects[[nm]]
  if (is.null(eset)) { message("[", nm, "] not downloaded; skipped."); next }
  is_mir  <- s$measures == "miRNA"
  fmap    <- feature_map(eset, if (is_mir) "mirna" else "gene")
  targets <- if (is_mir) CANDIDATES$mirna else GENES
  grp     <- get_group_vector(nm, eset)
  nr <- sum(grp == "resistant", na.rm = TRUE); ns <- sum(grp == "sensitive", na.rm = TRUE)
  message(sprintf("\n[%s] %s | %d resistant vs %d sensitive%s",
                  nm, s$measures, nr, ns, if (nr < 2 || ns < 2) "  (no replication -> fold change only)" else ""))
  out <- tryCatch(fc_two_group(eset, grp, fmap, targets, is_mir),
                  error = function(e) { message("  ", conditionMessage(e)); NULL })
  if (!is.null(out)) { out$series <- nm; out$measures <- s$measures; res[[nm]] <- out }
}
geo_fc <- bind_rows(res)
readr::write_csv(geo_fc, "results/04_geo_foldchange.csv")
message("\nLayer 4b complete -> results/04_geo_foldchange.csv")
print(as.data.frame(geo_fc))
message("\nSANITY CHECK: miR-200c-3p should be strongly NEGATIVE in GSE28267 (LY2 resistant vs MCF-7).")
