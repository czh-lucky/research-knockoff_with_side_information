#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- if (length(args) == 0) 1L else as.integer(args[1])

library(knockoff)
library(corpcor)
library(glmnet)

source("./R/knockoff_functions.R")

####################################
## 1. Load data
####################################

X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS("./data/MC_indices_n100_B10.rds")
coords_2d <- readRDS("./data/coords.rds")

idx <- MC_indices[[ParamsRowIndex]]

n <- length(idx)
p <- ncol(X_all)

stopifnot(nrow(coords_2d) == p)

# The 100 subjects used for variable selection
X_raw <- X_all[idx, , drop = FALSE]

# The remaining 481 subjects for covariance estimation
X_ref_raw <- X_all[-idx, , drop = FALSE]

####################################
## 2. Preprocessing
####################################

# Subject-wise median normalization
X_norm <- sweep(
    X_raw, 1, apply(X_raw, 1, median), "/"
)

X_ref_norm <- sweep(
    X_ref_raw, 1, apply(X_ref_raw, 1, median), "/"
)

# Use the reference sample to determine scaling parameters
voxel_mean <- colMeans(X_ref_norm)
voxel_sd <- apply(X_ref_norm, 2, sd)

stopifnot(all(is.finite(voxel_sd)), all(voxel_sd > 0))

X <- scale(
    X_norm,
    center = voxel_mean,
    scale = voxel_sd
)

X_ref <- scale(
    X_ref_norm,
    center = voxel_mean,
    scale = voxel_sd
)

####################################
## 3. Generate beta and Y
####################################

center <- c(coord_X = 31, coord_Y = 74)
radius <- 5

dist_to_center <- sqrt(
    (coords_2d[, "coord_X"] - center["coord_X"])^2 +
    (coords_2d[, "coord_Y"] - center["coord_Y"])^2
)

signal_idx <- which(dist_to_center <= radius)

set.seed(24601)

amp <- 3.5

beta0 <- amp *
    (1:p %in% signal_idx) *
    sign(rnorm(p)) / sqrt(n)

set.seed(ParamsRowIndex)

Y <- rbinom(
    n,
    1,
    plogis(drop(X %*% beta0))
)

stopifnot(length(unique(Y)) == 2)

####################################
## 4. Estimate covariance
####################################

Sigma_100 <- matrix(
    as.numeric(corpcor::cov.shrink(X, verbose = FALSE)),
    nrow = p
)

Sigma_481 <- matrix(
    as.numeric(corpcor::cov.shrink(X_ref, verbose = FALSE)),
    nrow = p
)

Sigma_list <- list(
    sample_100 = Sigma_100,
    reference_481 = Sigma_481
)

####################################
## 5. Generate knockoffs and run Vanilla
####################################

# Same CV folds for both methods
set.seed(10000 + ParamsRowIndex)
foldid <- sample(rep(1:10, length.out = n))

# Signal strength
eta <- drop(X %*% beta0)
cat("Linear predictor sd:", sd(eta), "\n")

# X-only Lasso
mdl_X <- cv.glmnet(
    X, Y,
    family = "binomial",
    alpha = 1,
    foldid = foldid
)

beta_X <- as.numeric(coef(mdl_X, s = "lambda.min"))[-1]
selected_X <- which(beta_X != 0)

cat("X-only selected:", length(selected_X), "\n")
cat("X-only true signals:",
    sum(beta0[selected_X] != 0), "\n")

all_res <- data.frame()

dir.create(
    "./data/intermediate/Xk_cov_compare",
    recursive = TRUE,
    showWarnings = FALSE
)

for (method in names(Sigma_list)) {

    cat("\n===== ", method, " =====\n", sep = "")

    Sigma <- Sigma_list[[method]]

    time <- system.time({

        # Optimized ASDP
        s <- create.solve_asdp.fast(
            Sigma,
            verbose = FALSE
        )

        # Generate knockoffs
        set.seed(ParamsRowIndex)

        Xk <- knockoff::create.gaussian(
            X,
            mu = rep(0, p),
            Sigma = Sigma,
            diag_s = s
        )
    })

    # Knockoff correlations
    paired_cor <- vapply(
        1:p,
        function(j) cor(X[, j], Xk[, j]),
        numeric(1)
    )

    cat("Median s:", median(s), "\n")
    cat("Median cor(Xj, Xkj):", median(paired_cor), "\n")
    cat("Elapsed:", time["elapsed"], "seconds\n")

    # Logistic Lasso
    mdl <- cv.glmnet(
        cbind(X, Xk),
        Y,
        alpha = 1,
        family = "binomial",
        foldid = foldid
    )

    beta <- as.numeric(
        coef(mdl, s = "lambda.min")
    )[-1]

    W <- abs(beta[1:p]) -
         abs(beta[(p + 1):(2 * p)])

    cat("W > 0:", sum(W > 0), "\n")
    cat("True signals W > 0:",
        sum(W[signal_idx] > 0), "\n")

    # FDR and Power
    for (q in seq(0.3, 0.01, by = -0.01)) {

        tau <- knockoff.threshold(
            W,
            fdr = q,
            offset = 1
        )

        rej <- which(W >= tau)

        fdp <- sum(beta0[rej] == 0) /
               max(length(rej), 1)

        power <- sum(beta0[rej] != 0) /
                 length(signal_idx)

        new_res <- data.frame(
            method = method,
            level = q,
            power = power,
            fdp = fdp,
            selected = length(rej),
            median_cor = median(paired_cor),
            median_s = median(s)
        )

        all_res <- rbind(all_res, new_res)
    }

    # Save knockoffs for further analysis
    saveRDS(
        Xk,
        sprintf(
            "./data/intermediate/Xk_cov_compare/%s_b%d.rds",
            method,
            ParamsRowIndex
        )
    )
}

####################################
## 6. Save comparison
####################################

dir.create(
    "./results/covariance_comparison",
    recursive = TRUE,
    showWarnings = FALSE
)

write.csv(
    all_res,
    sprintf(
        "./results/covariance_comparison/b%d.csv",
        ParamsRowIndex
    ),
    row.names = FALSE
)

cat("\nFinished.\n")