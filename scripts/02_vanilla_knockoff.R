#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

id <- as.integer(args[2])
if (is.na(id) == 1){
  id <- 1
}

amp <- as.numeric(args[3])
if (is.na(amp) == 1){
  amp <- 20
}

sub_num <- as.integer(args[4])
if (is.na(sub_num) == 1){
  sub_num <- 400
}

batch_num <- as.integer(args[5])
if (is.na(batch_num) == 1){
  batch_num <- 10
}

# choosing parameter
settingName <- "./results/vanilla_knockoff/one_cluster_simulation"
all_res <- data.frame()
alphalist <- seq(0.3, 0.01, -0.01)
verbose <- FALSE
good_enough <- FALSE

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
# suppressMessages(library("gam"))
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

file_vec <- c("./utils/utils_side_information/all_other_methods","./utils/utils_side_information/adaptive_knockoff","./utils/utils_side_information/filter_EM",
              "./utils/utils_side_information/filter_gam","./utils/utils_side_information/filter_glm","./utils/utils_side_information/filter_randomForest", 
              "./utils/utils_side_information/All_q_est_functions", "./utils/utils_side_information/accumulation_test_functions")
getfile <- sapply(file_vec,source_gitfile)

source("./utils/knockoff_functions.R")
source("./utils/target_generator.R")
source("./utils/W_overlay.R")
source("./utils/plot_helpers.R")

# parameter setting
# id <- 9
# amp <- 20
# sub_num <- 581
# batch_num <- 1
MC_file <- sprintf("./data/MC_indices_n%d_B%d.rds", sub_num, batch_num)

# 1. 数据读入
X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS(MC_file)
coords_2d <- readRDS("./data/coords.rds")
exp_grid <- read.csv("config/experiments.csv")
n <- length(MC_indices[[ParamsRowIndex]])
p <- ncol(X_all)


# generate signal
signal_idx <- get_signal_index_from_grid(coords_2d, exp_grid, id)
preview_target_on_slice(coords_2d, exp_grid, id)
num <- length(signal_idx)

cat("numbers of signal: ", num, "\n")

# 3. simulating signal
set.seed(24601)
beta0 <- amp * (1:p %in% signal_idx) * sign(rnorm(p, mean = 0, sd = 1)) / sqrt(n)
y.sample = function(X){
  rbinom(nrow(X), 1, plogis(X %*% beta0))
} # plogis(X %*% beta0) = exp(X %*% beta0) / (1 + exp(X %*% beta0))


####################################
## Generating data
####################################
X_raw <- X_all[MC_indices[[ParamsRowIndex]], , drop = FALSE]
subject_median <- apply(X_raw, 1, median)
X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)

## Signal strength diagnostics
eta <- drop(X %*% beta0)
prob <- plogis(eta)

# 绘制图像
plot(
  eta, prob, pch = 16, cex = .5, col = "steelblue", xlab = "eta", ylab = "P(Y = 1)",
  main = "Logistic probability"
)

hist(
  prob, breaks = seq(0, 1, by = 0.05), col = "skyblue", 
  main = "Probability distribution", xlab = "P(Y = 1)"
)

if (verbose){
  cat("\nLinear predictor sd:", sd(eta), "\n")
  cat("\nProbability summary:\n")
  print(summary(plogis(eta)))
  
  eta_cor <- drop(cor(X, eta))
  
  cat("\nTrue signal correlation with eta:\n")
  print(summary(abs(eta_cor[signal_idx])))
  
  cat("\nNull variable correlation with eta:\n")
  print(summary(abs(eta_cor[-signal_idx])))
}

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

file_path <- sprintf("./data/intermediate/Xk_median/n%d_B%d", sub_num, batch_num)
if (file.exists(file_path)){
  Xk <- readRDS(
    sprintf("./data/intermediate/Xk_median/n%d_B%d/Xk_b%02d.rds", sub_num, batch_num, ParamsRowIndex)
  )
}else{stop("Xk doesn't exit")}

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
# mdl_confine = cv.glmnet(cbind(X,Xk),Y,alpha=1,family = "binomial", lambda = 10^seq(-3, -2, length.out = 50))
mdl = cv.glmnet(cbind(X,Xk),Y,alpha=1,family = "binomial")
cvlambda = mdl$lambda.min
lambda_idx <- which.min(
  abs(mdl$glmnet.fit$lambda - mdl$lambda.min)
)
beta = mdl$glmnet.fit$beta[, lambda_idx]
T = abs(beta[1:p])
cat("T>0")
print(sum(T > 0))
T_tilde = abs(beta[(p+1):(2*p)])
# T_max = pmax(T,T_tilde)
W = T-T_tilde

####################################
## Vanilla diagnostics
####################################
cat("\n===== Vanilla diagnostics =====\n")

plot(mdl)

cat("Y distribution:\n")
print(table(Y))

cat("\nNonzero Lasso coefficients:",
    sum(beta != 0), "\n")

cat("\nTrue signals with nonzero Lasso coefficients:",
    sum(beta[signal_idx] != 0), "\n")

cat("\nW summary:\n")
print(summary(W))

cat("\nW > 0:", sum(W > 0), "\n")
cat("W < 0:", sum(W < 0), "\n")
cat("W = 0:", sum(W == 0), "\n")

cat("\nTrue signals with W > 0:",
    sum(W[signal_idx] > 0), "\n")

cat("True signals with W != 0:",
    sum(W[signal_idx] != 0), "\n")

cat("True signals with W < 0:",
    sum(W[signal_idx] < 0), "\n")

cat("\nThreshold at FDR = 0.3:",
    knockoff.threshold(W, fdr = 0.3, offset = 1),
    "\n")

W_all_overlay(W, signal_idx)
W_pos_overlay(W, signal_idx)
W_neg_overlay(W, signal_idx)
plot_W_ordered(W)



if (good_enough){
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
}