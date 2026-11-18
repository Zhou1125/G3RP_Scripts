suppressPackageStartupMessages(library(peer))

P <- "."
O <- file.path(P, "04_covariates/18_PEER/PEER5")

dir.create(O, recursive=TRUE, showWarnings=FALSE)

cat("===== GTEx PEER5 ANALYSIS =====\n")

# Expression: genes x samples
bed <- read.delim(
  gzfile(file.path(
    P, "03_expression/16_gtex_style/GTEx_style_TMM_INT_58.bed.gz"
  )),
  sep="\t",
  check.names=FALSE,
  comment.char="",
  quote=""
)

ids <- colnames(bed)[5:ncol(bed)]

stopifnot(
  length(ids)==58,
  nrow(bed)==20093,
  !anyDuplicated(ids)
)

Y <- t(as.matrix(bed[,5:ncol(bed),drop=FALSE]))
storage.mode(Y) <- "double"

rownames(Y) <- ids

stopifnot(all(is.finite(Y)))

# Known covariates: samples x 9
known <- read.delim(
  file.path(
    P, "04_covariates/17_GTEx/PEER_observed_covariates.tsv"
  ),
  sep="\t",
  row.names=1,
  check.names=FALSE
)

stopifnot(setequal(rownames(known),ids))

C <- as.matrix(known[ids,,drop=FALSE])
storage.mode(C) <- "double"

stopifnot(ncol(C)==9,all(is.finite(C)))

cat("SAMPLES =",nrow(Y),"\n")
cat("GENES =",ncol(Y),"\n")
cat("KNOWN_COVARIATES =",ncol(C),"\n")

# Fit PEER
set.seed(20261031)

model <- PEER()

PEER_setNk(model,5)
PEER_setPhenoMean(model,Y)
PEER_setPriorAlpha(model,0.001,0.01)
PEER_setPriorEps(model,0.1,10)
PEER_setNmax_iterations(model,1000)
PEER_setCovariates(model,C)

cat("===== FITTING PEER =====\n")
print(system.time(PEER_update(model)))

X <- as.matrix(PEER_getX(model))

cat("RETURNED_ROWS =",nrow(X),"\n")
cat("RETURNED_COLUMNS =",ncol(X),"\n")

# Official GTEx ordering:
# known covariates first, inferred factors last
if(nrow(X)!=58 || ncol(X)!=14) {
  stop("Unexpected PEER dimensions")
}

latent <- X[,10:14,drop=FALSE]

# Check this PEER build's returned-column convention before selecting factors.
if (!isTRUE(all.equal(unname(X[,1:9,drop=FALSE]), unname(C),
                     tolerance=1e-6, check.attributes=FALSE))) {
  stop("PEER known-covariate columns differ; verify the installed PEER build")
}

rownames(latent) <- ids
colnames(latent) <- paste0("PEER_",1:5)

stopifnot(all(is.finite(latent)))

latent <- scale(latent)

if(any(!is.finite(latent))) {
  stop("Invalid PEER factor values")
}

# Combined covariate matrix
combined <- cbind(C,latent)

D <- cbind(Intercept=1,combined)

rank <- qr(D)$rank
condition <- kappa(D,exact=TRUE)
residual_df <- 58-rank-1

cat("DESIGN_RANK =",rank,"\n")
cat("CONDITION_NUMBER =",condition,"\n")
cat("SNP_REGRESSION_DF =",residual_df,"\n")

if(rank!=ncol(D)) stop("Rank deficient")
if(!is.finite(condition) || condition>1e6) {
  stop("Severe multicollinearity")
}

# Write FastQTL-style covariates: rows=covariates
write_cov <- function(M,path) {
  Z <- t(M)
  write.table(
    data.frame(id=rownames(Z),Z,check.names=FALSE),
    path,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
  )
}

write_cov(
  latent,
  file.path(O,"PEER5_factors.tsv")
)

write_cov(
  combined,
  file.path(O,"M1_PEER5.covariates.tsv")
)

cat("PEER_FACTORS = 5\n")
cat("STEP18_PEER_SUCCESS\n")
