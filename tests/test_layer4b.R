## test_layer4b.R  -- validate corrected Layer 4b on synthetic mimics of the
## three ACTUAL platforms (gene array; common-reference two-color miRBase 11.0;
## single-channel 1-vs-1 miRNA array).
suppressMessages({library(Biobase); library(limma); library(dplyr)})

GENES <- c("OGT","OGA","EOGT"); GENE_ALIASES <- c(OGA="MGEA5", MGEA5="OGA")
CANDIDATES <- data.frame(mirna = c("hsa-miR-200c-3p","hsa-miR-24-3p","hsa-miR-222-3p","hsa-miR-195-5p"))

## ---- functions copied verbatim from R/03_geo_foldchange.R ----
canon_mir <- function(x){x<-tolower(x);x<-gsub("hsa-","",x);x<-gsub("^mir-?","",x);x<-gsub("-[35]p$","",x);x<-gsub("\\*$","",x);x<-gsub("-[0-9]+$","",x);trimws(x)}
canon_tokens <- function(s){t<-unlist(strsplit(as.character(s),"[ ,;/|]+"));t<-t[nzchar(t)];unique(canon_mir(t))}
gene_tokens <- function(s){t<-unlist(strsplit(as.character(s),"\\s*///\\s*|[;,]"));tolower(trimws(t[nzchar(t)]))}
feature_map <- function(eset, kind=c("gene","mirna")){kind<-match.arg(kind);fd<-fData(eset);cn<-colnames(fd)
  pick<-function(ps){for(p in ps){h<-grep(p,cn,ignore.case=TRUE,value=TRUE);if(length(h))return(h[1])};NA_character_}
  col<-if(kind=="gene")pick(c("^gene.?symbol$","gene.?symbol","^symbol$")) else pick(c("mirna_id","mirna","mir_id","transcript.?id\\(","transcript.?id","^name$"))
  if(is.na(col))return(setNames(rownames(fd),rownames(fd)));setNames(as.character(fd[[col]]),rownames(fd))}
to_log2 <- function(ex) if(max(ex,na.rm=TRUE)>100) log2(ex+1) else ex
collapse_targets <- function(feat,lfc,adjp,targets,is_mirna,n_res,n_sen,testable){
  keys<-if(is_mirna)lapply(feat,canon_tokens) else lapply(feat,gene_tokens)
  bind_rows(lapply(targets,function(f){
    if(is_mirna)idx<-vapply(keys,function(k)canon_mir(f)%in%k,logical(1))
    else{alt<-tolower(unique(c(f,names(GENE_ALIASES)[GENE_ALIASES==f],unname(GENE_ALIASES[f]))));idx<-vapply(keys,function(k)any(alt%in%k),logical(1))}
    if(!any(idx))return(data.frame(target=f,log2FC=NA_real_,adj_p=NA_real_,n_probe_hits=0L,n_resistant=n_res,n_sensitive=n_sen,testable=testable))
    i<-which(idx)[which.max(abs(lfc[idx]))]
    data.frame(target=f,log2FC=round(lfc[i],3),adj_p=if(testable)signif(adjp[i],3) else NA_real_,n_probe_hits=sum(idx),n_resistant=n_res,n_sensitive=n_sen,testable=testable)}))}
fc_two_group <- function(eset,grp,fmap,targets,is_mirna){
  ex<-to_log2(exprs(eset));keep<-grp%in%c("resistant","sensitive");ex<-ex[,keep,drop=FALSE];g<-factor(grp[keep],levels=c("sensitive","resistant"))
  n_res<-sum(g=="resistant");n_sen<-sum(g=="sensitive");if(n_res<1||n_sen<1)stop("need both groups");testable<-n_res>=2&&n_sen>=2;feat<-fmap[rownames(ex)]
  if(testable){fit<-eBayes(lmFit(ex,model.matrix(~g)));tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
    collapse_targets(feat[match(rownames(tt),rownames(ex))],tt$logFC,tt$adj.P.Val,targets,is_mirna,n_res,n_sen,testable)
  }else{lfc<-rowMeans(ex[,g=="resistant",drop=FALSE])-rowMeans(ex[,g=="sensitive",drop=FALSE]);collapse_targets(feat,unname(lfc),rep(NA_real_,length(lfc)),targets,is_mirna,n_res,n_sen,testable)}}

cat("== TEST canon_mir: no false collisions ==\n")
stopifnot(canon_mir("hsa-miR-24-3p") == "24")
stopifnot(canon_mir("hsa-miR-240")   == "240")     # must NOT equal 24
stopifnot(canon_mir("hsa-miR-200c")  == "200c")    # miRBase 11.0 guide name
stopifnot(canon_mir("hsa-miR-200c-3p")== "200c")   # newer name -> same key
stopifnot(canon_mir("hsa-miR-181b-5p")== "181b")
cat("  PASSED: 24 != 240; miR-200c (old) and miR-200c-3p (new) both canonicalize to 200c.\n\n")

cat("== TEST A: gene array (GPL570-like), 12 sensitive vs 12 resistant ==\n")
set.seed(1)
probes<-c("p1","p2","p3","p4","p5"); sym<-c("OGT","MGEA5 /// LOC1","EOGT","GAPDH","ACTB")
base<-c(6,6,5,8,8)
mk<-function(res) c(base[1],base[2],base[3]+if(res)0.66 else 0,base[4],base[5])
grpv<-rep(c("sensitive","resistant"),each=12)
mat<-sapply(grpv=="resistant",mk)+matrix(rnorm(5*24,0,0.1),5,24)
rownames(mat)<-probes;colnames(mat)<-paste0("g",1:24)
fd<-data.frame("Gene Symbol"=sym,row.names=probes,check.names=FALSE)
esA<-ExpressionSet(assayData=mat,featureData=AnnotatedDataFrame(fd))
grp<-setNames(grpv,colnames(mat))
rA<-fc_two_group(esA,grp,feature_map(esA,"gene"),GENES,FALSE);print(rA)
stopifnot(abs(rA$log2FC[rA$target=="EOGT"]-0.66)<0.15, rA$log2FC[rA$target=="EOGT"]>0)
stopifnot(!is.na(rA$log2FC[rA$target=="OGA"]))          # MGEA5 with "///" resolved
stopifnot(rA$testable[1]==TRUE)
cat("  PASSED: EOGT ~ +0.66 (up in resistant); OGA via 'MGEA5 /// ...'; p-values present.\n\n")

cat("== TEST B: common-reference two-color miRBase 11.0 (GPL7722-like) ==\n")
set.seed(2)
mp<-c("m1","m2","m3","m4"); mold<-c("hsa-miR-200c","hsa-miR-24","hsa-miR-222","hsa-miR-195")
ref<-c("200c"=0,"24"=0,"222"=0,"195"=0)
mcf7<-c(3, 1, 0.5, -0.5)     # log2(MCF7/ref)
ly2 <-c(-2.85, 2, 1.0,  0.0) # log2(LY2/ref)  -> LY2-MCF7 for 200c = -5.85
cols<-cbind(mcf7,mcf7, ly2,ly2,ly2,ly2)+matrix(rnorm(4*6,0,0.15),4,6)
colnames(cols)<-paste0("arr",1:6);rownames(cols)<-mp
mfd<-data.frame("miRNA_ID_LIST"=mold,row.names=mp,check.names=FALSE)
esB<-ExpressionSet(assayData=cols,featureData=AnnotatedDataFrame(mfd))
grpB<-setNames(c("sensitive","sensitive","resistant","resistant","resistant","resistant"),colnames(cols))
rB<-fc_two_group(esB,grpB,feature_map(esB,"mirna"),CANDIDATES$mirna,TRUE);print(rB)
stopifnot(rB$log2FC[rB$target=="hsa-miR-200c-3p"]<0)                       # correct sign
stopifnot(abs(rB$log2FC[rB$target=="hsa-miR-200c-3p"]-(ly2[1]-mcf7[1]))<0.4) # ~ -5.85
stopifnot(rB$n_probe_hits[rB$target=="hsa-miR-200c-3p"]==1)                # matched despite old name
stopifnot(rB$testable[1]==TRUE)
cat("  PASSED: miR-200c-3p ~ -5.85 (down in resistant), matched via old name 'hsa-miR-200c', two-group on log-ratios.\n\n")

cat("== TEST C: single-channel 1-vs-1 (GPL16384-like), fold change only ==\n")
set.seed(3)
cp<-c("c1","c2","c3"); cnm<-c("hsa-miR-200c-3p","hsa-miR-24-3p","hsa-miR-222-3p")
sens<-c(10, 5, 6); resi<-c(4.5, 6, 6.2)   # 200c down by -5.5 in resistant
cmat<-cbind(sens,resi);rownames(cmat)<-cp;colnames(cmat)<-c("parent","resistant_cell")
cfd<-data.frame("Transcript ID(Array Design)"=cnm,row.names=cp,check.names=FALSE)
esC<-ExpressionSet(assayData=cmat,featureData=AnnotatedDataFrame(cfd))
grpC<-setNames(c("sensitive","resistant"),colnames(cmat))
rC<-fc_two_group(esC,grpC,feature_map(esC,"mirna"),CANDIDATES$mirna,TRUE);print(rC)
stopifnot(rC$testable[1]==FALSE)                                          # no replication
stopifnot(is.na(rC$adj_p[rC$target=="hsa-miR-200c-3p"]))                  # no p-value
stopifnot(abs(rC$log2FC[rC$target=="hsa-miR-200c-3p"]-(resi[1]-sens[1]))<0.01) # exact difference
stopifnot(rC$log2FC[rC$target=="hsa-miR-200c-3p"]<0)
cat("  PASSED: 1-vs-1 gives fold change only (adj_p = NA), miR-200c-3p negative.\n\n")

cat("ALL LAYER 4b TESTS PASSED on synthetic mimics of the three real platforms.\n")
