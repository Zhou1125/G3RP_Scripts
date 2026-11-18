suppressPackageStartupMessages(library(edgeR))

P <- "."
OUT <- file.path(P, "03_expression/16_gtex_style")

cat("===== STEP16B TMM + INT =====\n")

counts <- read.delim(
  gzfile(file.path(P, "03_expression/gene_counts_58.tsv.gz")),
  check.names=FALSE
)

selected <- read.delim(
  file.path(OUT, "GTEx_selected_genes.tsv"),
  check.names=FALSE
)

samples <- read.delim(
  file.path(P, "03_expression/samples_58.tsv"),
  check.names=FALSE
)

ids <- as.character(samples$SampleID)

stopifnot(
  length(ids)==58,
  identical(colnames(counts)[-1], ids),
  !anyDuplicated(counts$Geneid),
  !anyDuplicated(selected$gene_id)
)

M <- as.matrix(counts[, -1, drop=FALSE])
storage.mode(M) <- "numeric"
rownames(M) <- as.character(counts$Geneid)

stopifnot(all(is.finite(M)), all(M>=0))

# Compute TMM from original full gene-count matrix
dge <- DGEList(counts=M)
dge <- calcNormFactors(dge, method="TMM")

cpm_matrix <- cpm(
  dge,
  log=FALSE,
  normalized.lib.sizes=TRUE
)

X <- cpm_matrix[
  as.character(selected$gene_id),
  ,drop=FALSE
]

# Remove zero-variance phenotypes
v <- apply(X, 1, var)
keep <- is.finite(v) & v > 1e-15

selected <- selected[keep, ,drop=FALSE]
X <- X[keep, ,drop=FALSE]

if (nrow(X)==0) stop("No variable genes")

# Inverse normal transformation of each gene
INT <- function(x) {
  qnorm((rank(x,ties.method="average")-0.5)/length(x))
}

Y <- t(apply(X, 1, INT))

rownames(Y) <- rownames(X)
colnames(Y) <- ids

stopifnot(
  identical(rownames(Y), as.character(selected$gene_id)),
  all(is.finite(Y))
)

# Save expression matrix
con <- gzfile(
  file.path(OUT, "GTEx_TMM_INT_matrix.tsv.gz"), "wt"
)

write.table(
  data.frame(gene_id=rownames(Y), Y, check.names=FALSE),
  con, sep="\t", quote=FALSE, row.names=FALSE
)

close(con)

# BED: chromosome, 0-based start, end, gene_id, samples
bed <- data.frame(
  chr=as.character(selected$gtf_chr),
  start=as.integer(selected$gtf_tss)-1L,
  end=as.integer(selected$gtf_tss),
  phenotype_id=as.character(selected$gene_id),
  check.names=FALSE
)

bed <- cbind(
  bed,
  as.data.frame(Y, check.names=FALSE)
)

colnames(bed)[1] <- "#chr"

stopifnot(
  all(bed$start>=0),
  all(bed$end-bed$start==1),
  identical(colnames(bed)[5:ncol(bed)],ids),
  !anyDuplicated(bed$phenotype_id)
)

write.table(
  bed,
  file.path(OUT,"GTEx_style_TMM_INT_58.bed"),
  sep="\t",quote=FALSE,row.names=FALSE
)

# Save TMM factors
factors <- data.frame(
  SampleID=ids,
  SourceGroup=samples$SourceGroup,
  LibrarySize=dge$samples$lib.size,
  TMMFactor=dge$samples$norm.factors,
  EffectiveLibrarySize=
    dge$samples$lib.size*dge$samples$norm.factors
)

write.table(
  factors,
  file.path(OUT,"TMM_factors.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE
)

cat("TMM_FACTOR_RANGE =",
    range(dge$samples$norm.factors),"\n")
cat("FINAL_BED_GENES =",nrow(bed),"\n")
cat("BED_SAMPLES =",length(ids),"\n")
cat("STEP16B_SUCCESS\n")
