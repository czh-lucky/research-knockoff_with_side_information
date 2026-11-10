#!/usr/bin/env Rscript
# Run from the project root:
# Rscript scripts/02_vanilla_knockoff_batch.R <b> <id> <amp> <n> <B>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L){
  stop("Usage: Rscript scripts/02_vanilla_knockoff_batch.R <b> <id> <amp> <n> <B>")
}
values <- suppressWarnings(as.numeric(args))
integer_args <- values[c(1, 2, 4, 5)]
if (any(!is.finite(values)) || any(integer_args < 1) ||
    any(integer_args != floor(integer_args)) ||
    any(integer_args > .Machine$integer.max)){
  stop("b, id, n and B must be positive integers; amp must be finite.")
}
ParamsRowIndex <- as.integer(values[1])
id <- as.integer(values[2])
amp <- values[3]
sub_num <- as.integer(values[4])
batch_num <- as.integer(values[5])
if (ParamsRowIndex > batch_num) stop("b cannot exceed B.")

# choosing parameter
all_res <- data.frame()
alphalist <- seq(0.3, 0.01, -0.01)
verbose <- FALSE
good_enough <- FALSE  # TRUE: also save FDP / Power; n = 100 always skips this.

MC_file <- sprintf("./data/MC_indices_n%d_B%d.rds", sub_num, batch_num)
Xk_file <- sprintf("./data/intermediate/Xk_median/n%d_B%d/Xk_b%02d.rds",
                   sub_num, batch_num, ParamsRowIndex)
input_files <- c("./data/X_all_581x4711.rds", MC_file, Xk_file,
                 "./data/coords.rds", "./config/experiments.csv",
                 "./data/MNI152_T1_2mm_brain_mask.nii.gz")
if (any(!file.exists(input_files))){
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

amp_tag <- gsub(".", "p", as.character(amp), fixed = TRUE)
out_dir <- file.path("./results/vanilla_knockoff",
                     sprintf("n%d_B%d", sub_num, batch_num),
                     sprintf("id%02d_amp%s", id, amp_tag),
                     sprintf("b%02d", ParamsRowIndex))
result_file <- file.path(out_dir, "result.rds")
if (file.exists(result_file)) stop("Result already exists: ", result_file)

# Each plotting call gets its own PNG; close the device even if plotting fails.
save_png <- function(filename, expr){
  png(file.path(out_dir, filename), width = 1400, height = 1000, res = 150)
  on.exit(dev.off())
  force(expr)
  invisible(NULL)
}

####################################
## Libraries and sources
####################################
# Only packages used in this experiment; plotting helpers load oro.nifti.
suppressMessages(library("knockoff"))
suppressMessages(library("glmnet"))
source("./utils/target_generator.R")
source("./utils/W_overlay.R")
source("./utils/plot_helpers.R")

# 1. 数据读入
X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS(MC_file)
coords_2d <- readRDS("./data/coords.rds")
exp_grid <- read.csv("config/experiments.csv")
if (length(MC_indices) != batch_num) stop("MC_indices length does not match B.")
if (!id %in% exp_grid$run_id) stop("id not found in config/experiments.csv.")
sample_indices <- MC_indices[[ParamsRowIndex]]
if (length(sample_indices) != sub_num || anyNA(sample_indices) ||
    any(!is.finite(sample_indices)) || any(sample_indices != floor(sample_indices)) ||
    any(sample_indices < 1 | sample_indices > nrow(X_all)) ||
    anyDuplicated(sample_indices)){
  stop("Invalid subject indices or sample size in MC_indices.")
}
if (nrow(coords_2d) != ncol(X_all)) stop("coords and X_all do not match.")
n <- length(MC_indices[[ParamsRowIndex]])
p <- ncol(X_all)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output:", out_dir, "\nXk:", Xk_file, "\n")

# generate signal
signal_idx <- get_signal_index_from_grid(coords_2d, exp_grid, id)
save_png("target.png", preview_target_on_slice(coords_2d, exp_grid, id))
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
stopifnot(all(is.finite(subject_median)), all(subject_median > 0))
X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)
if (any(!is.finite(X))) stop("X contains non-finite values after scaling.")

## Signal strength diagnostics
eta <- drop(X %*% beta0)
prob <- plogis(eta)

# 绘制图像
save_png("logistic_probability.png", plot(
  eta, prob, pch = 16, cex = .5, col = "steelblue", xlab = "eta", ylab = "P(Y = 1)",
  main = "Logistic probability"
))

save_png("probability_distribution.png", hist(
  prob, breaks = seq(0, 1, by = 0.05), col = "skyblue", 
  main = "Probability distribution", xlab = "P(Y = 1)"
))

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
if (length(unique(Y)) < 2L) stop("Y has only one class; cannot fit binomial Lasso.")
k <- sum(beta0 != 0)

Xk <- readRDS(Xk_file)

stopifnot(
    identical(dim(X), dim(Xk)), all(is.finite(Xk))
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

save_png("lasso_cv.png", plot(mdl))

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

save_png("W_all_overlay.png", W_all_overlay(W, signal_idx))
save_png("W_pos_overlay.png", W_pos_overlay(W, signal_idx))
save_png("W_neg_overlay.png", W_neg_overlay(W, signal_idx))
save_png("W_ordered.png", plot_W_ordered(W))



if (good_enough && sub_num != 100L){
  fdp = c()
  power = c()
  for (i in alphalist ){
    tau = knockoff.threshold(W,fdr = i,offset = 1)
    rej = which(W>=tau)
    fdp = c(fdp,sum(beta0[rej]==0)/max(length(rej),1))
    power = c(power,sum(beta0[rej]!=0)/max(k,1))
  }
  
  # fdp is the realized false discovery proportion in this MC replicate.
  new_res <- data.frame(b = ParamsRowIndex, id = id, amp = amp, n = n, B = batch_num,
                        power = power, fdp = fdp, method = "vanilla", level = alphalist)
  all_res <- rbind(all_res, new_res)
  
  # write out
  write.csv(all_res, file.path(out_dir, "fdp_power.csv"), row.names = FALSE)
}

####################################
## Save experiment results
####################################
params <- list(ParamsRowIndex = ParamsRowIndex, id = id, amp = amp,
               sub_num = sub_num, batch_num = batch_num, beta_seed = 24601L,
               response_seed = ParamsRowIndex, alphalist = alphalist,
               compute_fdp_power = good_enough && sub_num != 100L)
diagnostics <- data.frame(
  b = ParamsRowIndex, id = id, amp = amp, n = n, B = batch_num, p = p,
  signals = k, Y_zero = sum(Y == 0), Y_one = sum(Y == 1), lambda_min = cvlambda,
  lasso_nonzero = sum(beta != 0), W_positive = sum(W > 0),
  W_negative = sum(W < 0), W_zero = sum(W == 0),
  signal_W_positive = sum(W[signal_idx] > 0),
  signal_W_negative = sum(W[signal_idx] < 0)
)
write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)
write.csv(data.frame(feature = seq_len(p), signal = beta0 != 0,
                     beta0 = beta0, T = as.numeric(T),
                     T_tilde = as.numeric(T_tilde), W = as.numeric(W)),
          file.path(out_dir, "W.csv"), row.names = FALSE)

# Record the inputs by path; do not duplicate X_all / Xk in every result.
# Written last: its presence marks a completed experiment.
saveRDS(list(params = params, MC_file = MC_file, Xk_file = Xk_file,
             sample_indices = sample_indices, signal_idx = signal_idx,
             beta0 = beta0, Y = Y, beta = beta, T = T, T_tilde = T_tilde,
             W = W, cvlambda = cvlambda, mdl = mdl,
             diagnostics = diagnostics, fdp_power = all_res,
             session_info = sessionInfo()), result_file)
cat("Finished. Saved to:", out_dir, "\n")
