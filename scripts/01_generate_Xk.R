#!/usr/bin/env Rscript
## Start of problem independent section
args <- commandArgs(trailingOnly = TRUE)
ParamsRowIndex <- as.integer(args[1])
if(is.na(ParamsRowIndex)==1){
  ParamsRowIndex = 1 
}

input_file <- args[2]
if (is.na(input_file)){
  input_file <- "data/MC_indices_n100_B10.rds"
}

# 从文件名中获取n, B
filename <- basename(input_file)
pattern <- "^MC_indices_n([0-9]+)_B([0-9]+)\\.rds$"
matches <- regmatches(filename, regexec(pattern, filename))[[1]]
if (length(matches) != 3L){
  stop("Invalide filename format: ", filename)
}

n <- as.integer(matches[2])
B <- as.integer(matches[3])

# 加载包和函数
library(knockoff)
library(corpcor)

source("./utils/knockoff_functions.R")

# 读取数据和抽样
X_all <- readRDS("./data/X_all_581x4711.rds")
MC_indices <- readRDS(input_file)

# 创建路径并检查路径
output_dir <-sprintf("./data/intermediate/Xk_median/n%d_B%d", n, B)
output_file <- file.path(output_dir, sprintf("Xk_b%02d.rds", ParamsRowIndex))
if (file.exists(output_file)){
  stop("Output already exists: ", output_file)
}

# 生成knockoff变量
X_raw <- X_all[MC_indices[[ParamsRowIndex]], , drop = FALSE]
subject_median <- apply(X_raw, 1, median)

stopifnot(all(is.finite(subject_median)), all(subject_median > 0))

X_norm <- sweep(X_raw, 1, subject_median, "/")
X <- scale(X_norm)

set.seed(ParamsRowIndex)
system.time({
    Xk <- create.second_order.fast(
        X,
        method = "asdp"
    )
}) -> time

# 输出
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(Xk, output_file)

stopifnot(identical(dim(X), dim(Xk)))

## Compare original and knockoff correlations
pair_cor <- vapply(
    seq_len(ncol(X)), function(j) cor(X[, j], Xk[, j]), numeric(1)
)
cat("\nCorrelation between X_j and Xk_j:\n")
print(summary(pair_cor))
cat("\nProportion with correlation > 0.99:", mean(pair_cor > 0.99), "\n")

cat("Finished.\n")
cat("Elapsed time:", time["elapsed"], "seconds\n")
cat("Dimensions:", paste(dim(Xk), collapse = " x "), "\n")