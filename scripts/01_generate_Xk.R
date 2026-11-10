#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

library(knockoff)
library(corpcor)

source("./R/knockoff_functions.R")

X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS("./data/MC_indices_n100_B10.rds")

X_raw <- X_all[MC_indices[[ParamsRowIndex]], , drop = FALSE]
subject_median <- apply(X_raw, 1, median)

stopifnot(all(is.finite(subject_median)), all(subject_median > 0))

X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)

Sigma_sample <- cov(X)

Sigma_shrink <- matrix(
    as.numeric(
        corpcor::cov.shrink(
            X,
            verbose = FALSE
        )
    ),
    nrow = ncol(X)
)

cat("\ndiag sample Sigma summary:\n")
print(summary(diag(Sigma_sample)))

cat("\ndiag shrink Sigma summary:\n")
print(summary(diag(Sigma_shrink)))

upper_idx <- upper.tri(Sigma_sample)

cat("\nSample covariance off-diagonal summary:\n")
print(
    summary(
        Sigma_sample[upper_idx]
    )
)

cat("\nShrink covariance off-diagonal summary:\n")
print(
    summary(
        Sigma_shrink[upper_idx]
    )
)

set.seed(ParamsRowIndex)
system.time({
    Xk <- create.second_order.fast(
        X,
        method = "asdp"
    )
}) -> time

saveRDS(Xk, sprintf("./data/intermediate/Xk/Xk_b%s.rds", ParamsRowIndex))

cat("Finished.\n")
cat("Elapsed time:", time["elapsed"], "seconds\n")
cat("Dimensions:", paste(dim(Xk), collapse = " x "), "\n")