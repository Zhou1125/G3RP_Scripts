suppressPackageStartupMessages(library(qvalue))

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=1 || !(args[1] %in% c("M1","M1_PEER5"))) {
  stop("Usage: Rscript 22_calculate_gtex_fdr.R M1|M1_PEER5")
}

MODEL <- args[1]

P <- "."

INPUT <- file.path(
  P,"05_fastqtl/21_full",MODEL,
  "permute",paste0(MODEL,".permute.tsv.gz")
)

OUT <- file.path(P,"05_fastqtl/22_fdr",MODEL)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

FDR <- 0.05

cat("===== STEP22 GTEX QVALUE =====\n")
cat("MODEL =",MODEL,"\n")
cat("INPUT =",INPUT,"\n")

stopifnot(file.exists(INPUT))

D <- read.delim(
  gzfile(INPUT),
  sep="\t",
  header=TRUE,
  check.names=FALSE,
  stringsAsFactors=FALSE,
  na.strings=c("NA","NaN","nan","."),
  quote=""
)

required <- c(
  "gene_id","num_var","beta_shape1",
  "beta_shape2","variant_id","pval_nominal",
  "pval_perm","pval_beta"
)

if (!all(required %in% colnames(D))) {
  stop("Missing required FastQTL output columns")
}

if (nrow(D)!=20093) {
  stop(paste("Unexpected gene count:",nrow(D)))
}

if (anyDuplicated(D$gene_id)) {
  stop("Duplicate gene IDs in permutation output")
}

for (nm in c(
  "pval_beta","pval_perm","pval_nominal",
  "beta_shape1","beta_shape2"
)) {
  D[[nm]] <- as.numeric(D[[nm]])
}

badp <- is.finite(D$pval_beta) &
  (D$pval_beta<0 | D$pval_beta>1)

if (any(badp)) {
  stop("Invalid beta-adjusted p-values")
}

valid <- is.finite(D$pval_beta)

E <- D[valid,,drop=FALSE]

if (nrow(E)==0) {
  stop("No estimable genes")
}

if (any(
  !is.finite(E$beta_shape1) |
  !is.finite(E$beta_shape2) |
  E$beta_shape1<=0 |
  E$beta_shape2<=0
)) {
  stop("Invalid beta-distribution parameters")
}

cat("EXPRESSION_GENES =",nrow(D),"\n")
cat("TESTED_GENES =",nrow(E),"\n")
cat("NON_ESTIMABLE_GENES =",sum(!valid),"\n")

# GTEx-style Storey q-values
Q <- qvalue(E$pval_beta)

E$qval <- Q$qvalues
E$BH_FDR_beta <- p.adjust(E$pval_beta,method="BH")

stopifnot(all(is.finite(E$qval)))

sig <- E$qval <= FDR

cat("PI0 =",Q$pi0,"\n")
cat("EGENES_FDR005 =",sum(sig),"\n")
cat("MIN_QVALUE =",min(E$qval),"\n")

# GTEx gene-specific nominal p-value thresholds
if (any(sig) && any(!sig)) {

  lower <- max(E$pval_beta[sig])
  upper <- min(E$pval_beta[!sig])

  if (lower>upper) {
    stop("Inconsistent q-value significance boundary")
  }

  cutoff <- (lower+upper)/2
  status <- "TWO_SIDED_BOUNDARY"

} else if (all(sig)) {

  cutoff <- 1
  status <- "ALL_GENES_SIGNIFICANT"

} else {

  cutoff <- NA_real_
  status <- "NO_SIGNIFICANT_EGENES"

}

if (is.finite(cutoff)) {

  E$pval_nominal_threshold <- qbeta(
    cutoff,
    shape1=E$beta_shape1,
    shape2=E$beta_shape2
  )

} else {

  E$pval_nominal_threshold <- NA_real_

}

E$is_eGene <- sig

# Preserve all 20,093 genes, including non-estimable genes
D$qval <- NA_real_
D$BH_FDR_beta <- NA_real_
D$pval_nominal_threshold <- NA_real_
D$is_eGene <- FALSE

idx <- which(valid)

D$qval[idx] <- E$qval
D$BH_FDR_beta[idx] <- E$BH_FDR_beta
D$pval_nominal_threshold[idx] <- E$pval_nominal_threshold
D$is_eGene[idx] <- E$is_eGene

egenes <- E[E$is_eGene,,drop=FALSE]

egenes <- egenes[
  order(egenes$qval,egenes$pval_beta),
  ,drop=FALSE
]

write_gz <- function(x,path) {
  con <- gzfile(path,"wt")
  write.table(
    x,con,sep="\t",
    quote=FALSE,row.names=FALSE,
    na="NA"
  )
  close(con)
}

write_gz(
  D,
  file.path(OUT,"all_genes_FDR.tsv.gz")
)

write_gz(
  egenes,
  file.path(OUT,"egenes_FDR005.tsv.gz")
)

write.table(
  head(egenes,20),
  file.path(OUT,"top20_egenes.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

summary <- data.frame(
  Model=MODEL,
  ExpressionGenes=nrow(D),
  TestedGenes=nrow(E),
  NonEstimableGenes=sum(!valid),
  SignificantEGenes=sum(sig),
  FDR=FDR,
  PI0=Q$pi0,
  MinQvalue=min(E$qval),
  BetaP_GlobalCutoff=cutoff,
  ThresholdStatus=status
)

write.table(
  summary,
  file.path(OUT,"FDR_summary.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat("GLOBAL_BETA_P_CUTOFF =",cutoff,"\n")
cat("THRESHOLD_STATUS =",status,"\n")
cat("OUTPUT =",OUT,"\n")
cat("STEP22_",MODEL,"_SUCCESS\n",sep="")
