library(corpcor)

MC_indices <- readRDS("./data/MC_indices_n100_B10.rds")
X_all <- readRDS("./data/X_all_581x4711.rds")

for (b in c(9, 10)) {
    X_raw <- X_all[MC_indices[[b]], , drop = FALSE]

    X <- scale(sweep(
        X_raw, 1, apply(X_raw, 1, median), "/"
    ))

    lambda <- corpcor::estimate.lambda(
        X, verbose = FALSE
    )

    cat("Dataset", b, ": lambda =", lambda, "\n")
}