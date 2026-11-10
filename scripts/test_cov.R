#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

# 1. 数据读入
X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS("./data/MC_indices_n100_B10.rds")

X_raw <- X_all[MC_indices[[ParamsRowIndex]], , drop = FALSE]
subject_median <- apply(X_raw, 1, median)

stopifnot(all(is.finite(subject_median)), all(subject_median > 0))

X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)

Sigma <- tryCatch({
  suppressWarnings(
    matrix(
      as.numeric(
        corpcor::cov.shrink(X, verbose = F)
      ),
      nrow = ncol(X)
    )
  )
}, warning = function(w) {
}, error = function(e) {
  stop(
    "SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0"
  )
}, finally = {
})

dir.create("./results/compare_cov", recursive = TRUE, showWarnings = FALSE)
vec_sigma <- as.vector(Sigma)

capture.output(
  {print(summary(vec_sigma))
  cat("Variance:", var(vec_sigma), "\n")
  cat("SD:", sd(vec_sigma), "\n")},
  file = sprintf("./results/compare_cov/Sigma_summary_%d.txt", ParamsRowIndex)
)