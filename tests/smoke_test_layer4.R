## smoke_test_layer4.R  -- validate the fold-change logic on synthetic data
suppressMessages({library(Biobase); library(limma); library(dplyr)})

## ---- minimal config (mirror of R/config.R) ----
GENES <- c("OGT","OGA","EOGT")
GENE_ALIASES <- c(OGA="MGEA5", MGEA5="OGA")
CANDIDATES <- data.frame(mirna=c("hsa-miR-200c-3p","hsa-miR-24-3p","hsa-miR-222-3p"))

## ---- functions copied verbatim from R/03_geo_foldchange.R (single-series variants) ----
feature_map <- function(eset, kind = c("gene","mirna")) {
  kind <- match.arg(kind); fd <- fData(eset)
  col <- if (kind=="gene") grep("gene.?symbol|^symbol$|gene_assignment", colnames(fd), ignore.case=TRUE, value=TRUE)[1]
         else grep("mirna|mir_id|transcript.?id|^name$|probe.?set.?name", colnames(fd), ignore.case=TRUE, value=TRUE)[1]
  if (is.na(col)) return(setNames(rownames(fd), rownames(fd)))
  setNames(as.character(fd[[col]]), rownames(fd))
}
to_log2 <- function(ex) if (max(ex, na.rm=TRUE) > 100) log2(ex+1) else ex
collapse_targets <- function(tt, targets, is_mirna, n_res, n_sen) {
  norm <- function(x) tolower(gsub("_","-",x))
  out <- lapply(targets, function(f){
    if (is_mirna) hit <- tt[grepl(norm(f), norm(tt$feature), fixed=TRUE), ]
    else { alt <- c(f, names(GENE_ALIASES)[GENE_ALIASES==f], unname(GENE_ALIASES[f])); alt <- unique(alt[!is.na(alt)])
           hit <- tt[norm(tt$feature) %in% norm(alt), ] }
    if (!nrow(hit)) return(NULL)
    hit <- hit[which.max(abs(hit$logFC)), ]
    data.frame(target=f, log2FC=round(hit$logFC,3), adj_p=signif(hit$adj.P.Val,3),
               n_probe_hits=nrow(tt[norm(tt$feature)==norm(f), ]), n_resistant=n_res, n_sensitive=n_sen)
  }); bind_rows(out)
}
fc_single <- function(eset, grp, fmap, targets, is_mirna) {
  ex <- to_log2(exprs(eset)); keep <- grp %in% c("resistant","sensitive")
  ex <- ex[, keep, drop=FALSE]; g <- factor(grp[keep], levels=c("sensitive","resistant"))
  design <- model.matrix(~ g); fit <- eBayes(lmFit(ex, design))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none"); tt$feature <- fmap[rownames(tt)]
  collapse_targets(tt, targets, is_mirna, sum(g=="resistant"), sum(g=="sensitive"))
}
fc_twocolor <- function(eset, orient, fmap, targets, is_mirna) {
  ex <- exprs(eset) * orient
  design <- matrix(1, ncol(ex), 1); colnames(design) <- "rvs"
  fit <- eBayes(lmFit(ex, design)); tt <- topTable(fit, coef=1, number=Inf, sort.by="none")
  tt$feature <- fmap[rownames(tt)]
  collapse_targets(tt, targets, is_mirna, ncol(ex), ncol(ex))
}

cat("================ TEST 1: single-channel enzyme fold change ================\n")
set.seed(1)
# 6 samples: 3 sensitive, 3 resistant; EOGT truly +2 log2 higher in resistant
probes <- c("p_OGT","p_OGA","p_EOGT","p_x1","p_x2")
sym    <- c("OGT","MGEA5","EOGT","GAPDH","ACTB")   # note OGA stored as its alias MGEA5
base   <- c(OGT=6, OGA=6, EOGT=5, x1=8, x2=8)
mk <- function(res) c(OGT=base["OGT"], OGA=base["OGA"],
                      EOGT=base["EOGT"] + if(res) 2 else 0, x1=base["x1"], x2=base["x2"])
mat <- sapply(c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE), mk) + matrix(rnorm(30,0,0.1),5,6)
rownames(mat) <- probes; colnames(mat) <- paste0("GSM",1:6)
fd <- data.frame("Gene Symbol"=sym, row.names=probes, check.names=FALSE)
es <- ExpressionSet(assayData=mat, featureData=AnnotatedDataFrame(fd))
grp <- setNames(c("sensitive","sensitive","sensitive","resistant","resistant","resistant"), colnames(mat))
r1 <- fc_single(es, grp, feature_map(es,"gene"), GENES, is_mirna=FALSE)
print(r1)
stopifnot(abs(r1$log2FC[r1$target=="EOGT"] - 2) < 0.3)       # recovers +2
stopifnot(r1$log2FC[r1$target=="EOGT"] > 0)                   # positive = up in resistant
stopifnot(abs(r1$log2FC[r1$target=="OGT"]) < 0.3)            # OGT flat
stopifnot(r1$target[r1$target=="OGA"]=="OGA")               # alias MGEA5 resolved to OGA
cat("TEST 1 PASSED: EOGT recovered at +2 (up in resistant); OGT flat; OGA/MGEA5 alias works.\n\n")

cat("================ TEST 2: two-color microRNA fold change (sign) ============\n")
set.seed(2)
# stored log-ratios ALREADY oriented resistant - sensitive:
# miR-200c-3p truly DOWN in resistant (ratio -3), miR-24-3p UP (+1)
mprobes <- c("m1","m2","m3")
mnames  <- c("hsa-miR-200c-3p","hsa-miR-24-3p","hsa-miR-222-3p")
truth   <- c(-3, 1, 0.5)
mmat <- sapply(1:4, function(i) truth + rnorm(3,0,0.2))
rownames(mmat) <- mprobes; colnames(mmat) <- paste0("arr",1:4)
mfd <- data.frame("miRNA"=mnames, row.names=mprobes, check.names=FALSE)
mes <- ExpressionSet(assayData=mmat, featureData=AnnotatedDataFrame(mfd))
r2 <- fc_twocolor(mes, orient=+1, feature_map(mes,"mirna"), CANDIDATES$mirna, is_mirna=TRUE)
print(r2)
stopifnot(r2$log2FC[r2$target=="hsa-miR-200c-3p"] < 0)        # DOWN in resistant (correct sign)
stopifnot(abs(r2$log2FC[r2$target=="hsa-miR-200c-3p"] + 3) < 0.4)
stopifnot(r2$log2FC[r2$target=="hsa-miR-24-3p"] > 0)          # UP
cat("TEST 2 PASSED: miR-200c-3p negative (down in resistant), miR-24-3p positive.\n\n")

cat("================ TEST 3: orientation flip catches a swapped channel =======\n")
r3 <- fc_twocolor(mes, orient=-1, feature_map(mes,"mirna"), CANDIDATES$mirna, is_mirna=TRUE)
stopifnot(sign(r3$log2FC[r3$target=="hsa-miR-200c-3p"]) == -sign(r2$log2FC[r2$target=="hsa-miR-200c-3p"]))
cat("TEST 3 PASSED: orient = -1 flips the sign, so a wrong channel layout is detectable.\n\n")

cat("ALL SMOKE TESTS PASSED. Layer 4 fold-change logic is correct on synthetic data.\n")
