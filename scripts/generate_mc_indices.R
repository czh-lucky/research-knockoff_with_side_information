X_all <- readRDS("X_all_581x4711.rds")

B <- 10
n <- 100

MC_indices <- vector("list", B)

for (b in seq_len(B)){
    set.seed(b)
    MC_indices[[b]] <- sample(
        seq_len(nrow(X_all)),
        size = n,
        replace = FALSE
    )
}

saveRDS(
    MC_indices,
    "MC_indices_n100_B10.rds"
)