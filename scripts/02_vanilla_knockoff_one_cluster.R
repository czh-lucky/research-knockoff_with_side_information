#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}
# choosing parameter
settingName <- "./results/vanilla_knockoff/one_cluster_simulation"
all_res <- data.frame()
alphalist <- seq(0.3, 0.01, -0.01)
verbose <- FALSE

####################################
## Libraries and sources
####################################
# suppressMessages(library("adaptMT"))
suppressMessages(library("splines"))
suppressMessages(library("knockoff"))
# suppressMessages(library("SNPknock"))
suppressMessages(library("dplyr"))
suppressMessages(library("corpcor"))
suppressMessages(library("glmnet"))
suppressMessages(library("MASS"))
suppressMessages(library("gam"))
suppressMessages(library("randomForest"))
# suppressMessages(library("tidyverse"))
suppressMessages(library("readr"))
suppressMessages(library("mgcv"))
# suppressMessages(library("flare"))
# suppressMessages(library("expm"))
# suppressMessages(library("rdist"))

source_gitfile <- function(filename){
  source(sprintf("%s.R", filename))
}

file_vec <- c("./R/utils/all_other_methods","./R/utils/adaptive_knockoff","./R/utils/filter_EM",
              "./R/utils/filter_gam","./R/utils/filter_glm","./R/utils/filter_randomForest", 
              "./R/utils/All_q_est_functions", "./R/utils/accumulation_test_functions")
getfile <- sapply(file_vec,source_gitfile)

source("./R/knockoff_functions.R")

# 1. 数据读入
X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS("./data/MC_indices_n100_B10.rds")
coords_2d <- readRDS("./data/coords.rds")
n <- length(MC_indices[[ParamsRowIndex]])
p <- ncol(X_all)

# 2. 生成模拟的信号区域
# signal center
center <- c(coord_X = 31, coord_Y = 74)

# signal radius
radius <- 5

# distance from each brain voxel to center
dist_to_center <- sqrt(
  (coords_2d[, "coord_X"] - center["coord_X"])^2 + (coords_2d[, "coord_Y"] - center["coord_Y"])^2
)

# True signal support
signal_idx <- which(dist_to_center <= radius)

# 3. simulating signal
set.seed(24601)
amp <- 3.5
beta0 <- amp * (1:p %in% signal_idx) * sign(rnorm(p, mean = 0, sd = 1)) / sqrt(n)
y.sample = function(X){
  rbinom(nrow(X), 1, plogis(drop(X %*% beta0)))
} # plogis(X %*% beta0) = exp(X %*% beta0) / (1 + exp(X %*% beta0))


####################################
## Generating data
####################################
X_raw <- X_all[MC_indices[[ParamsRowIndex]], , drop = FALSE]
subject_median <- apply(X_raw, 1, median)
X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)

####################################
## PCA diagnostics
####################################
if (verbose){
    pc <- prcomp(
        X,
        center = FALSE,
        scale. = FALSE
    )

    prop_var <- pc$sdev^2 / sum(pc$sdev^2)

    cat("\n")
    cat("===== PCA diagnostics =====\n")

    cat("Top 10 PC variance proportions:\n")
    print(prop_var[1:10])

    cat(
        "PC1 variance proportion =",
        prop_var[1],
        "\n"
    )

    cat(
        "Cumulative variance first 5 PCs =",
        sum(prop_var[1:5]),
        "\n"
    )

    cat(
        "Cumulative variance first 10 PCs =",
        sum(prop_var[1:10]),
        "\n"
    )

    cat(
        "Cumulative variance first 20 PCs =",
        sum(prop_var[1:20]),
        "\n"
    )
}

####################################
## PCA diagnostics_end
####################################

set.seed(ParamsRowIndex)
Y <- y.sample(X)
k <- sum(beta0 != 0)

Xk <- readRDS(
  sprintf("./data/intermediate/Xk_median/Xk_b%s.rds", ParamsRowIndex)
)

stopifnot(
    all(dim(X) == dim(Xk))
)

## Check X and Xk
if (verbose){
    x_sd <- apply(X, 2, sd)
    xk_sd <- apply(Xk, 2, sd)

    paired_cor <- vapply(
        1:p,
        function(j) cor(X[, j], Xk[, j]),
        numeric(1)
    )

    cat("\n")
    cat("===== X / Xk diagnostics =====\n")

    cat(
        "Dimensions:",
        paste(dim(X), collapse = " x "),
        "\n"
    )

    cat(
        "X mean range:",
        paste(round(range(colMeans(X)), 4), collapse = " to "),
        "\n"
    )

    cat(
        "X sd range:",
        paste(round(range(x_sd), 4), collapse = " to "),
        "\n"
    )

    cat(
        "Xk mean range:",
        paste(round(range(colMeans(Xk)), 4), collapse = " to "),
        "\n"
    )

    cat(
        "Xk sd range:",
        paste(round(range(xk_sd), 4), collapse = " to "),
        "\n"
    )

    cat("cor(X_j, Xk_j) summary:\n")
    print(summary(paired_cor))

    cat("\n1 - cor(X_j, Xk_j) summary:\n")
    print(
        summary(
            1 - paired_cor
        )
    )

    x_diff_sd <- vapply(
        1:p,
        function(j) sd(X[, j] - Xk[, j]),
        numeric(1)
    )

    cat("\nsd(X_j - Xk_j) summary:\n")
    print(
        summary(x_diff_sd)
    )
}

####################################
## Vanilla  knockoff
####################################
# Xk <- create.second_order.fast(X)
mdl = cv.glmnet(cbind(X,Xk),Y,alpha=1,family = "binomial")
cvlambda = mdl$lambda.min
lambda_idx <- which.min(
  abs(mdl$glmnet.fit$lambda - mdl$lambda.min)
)
beta = mdl$glmnet.fit$beta[, lambda_idx]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W = T-T_tilde

fdp = c()
power = c()
for (i in alphalist ){
  tau = knockoff.threshold(W,fdr = i,offset = 1)
  rej = which(W>=tau)
  fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
  power = c(power,sum(beta0[rej]!=0)/max(k,1))
}

new_res <- data.frame(power = power, fdp = fdp, method = "vanilla", level = alphalist)
all_res <- rbind(all_res, new_res)

# write out
dir.create("./results/vanilla_knockoff", recursive = TRUE, showWarnings = FALSE)
out_dir <- sprintf("%s_%d.csv", settingName, ParamsRowIndex)
write_csv(all_res, out_dir)