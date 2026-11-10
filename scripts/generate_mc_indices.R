#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L){
  stop("Usage: Rscript scripts/generate_mc_indices.R <n><B>")
}

n <- suppressWarnings(as.integer(args[1]))
B <- suppressWarnings(as.integer(args[2]))

if (anyNA(c(n, B)) || n < 1 || B < 1){
  stop("n and B must be positive integers.")
}

X_all <- readRDS("./data/X_all_581x4711.rds")

if (n > nrow(X_all)){
  stop("n cannot exceed the number of subjects.")
}

MC_indices <- vector("list", B)

for (b in seq_len(B)){
    set.seed(b)
    MC_indices[[b]] <- sample(
        seq_len(nrow(X_all)),
        size = n,
        replace = FALSE
    )
}

output_file <- sprintf("./data/MC_indices_n%d_B%d.rds", n, B)
if (file.exists(output_file)){
  stop("Output already exists: ", output_file)
}

saveRDS(MC_indices, output_file)

cat("Finished!\n")
cat("n=", n, "\n")
cat("B=", B, "\n")
cat("Saved to:", output_file, "\n")