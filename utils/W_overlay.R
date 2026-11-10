#!/usr/bin/env Rscript

library(oro.nifti)
source("utils/nifti_functions.R")

brain_2d_all <- readRDS("data/X_all_581x4711.rds")
brain_slice_real <- vector_to_slice(brain_2d_all[1,])

W_all_overlay <- function(W, signal_idx, brain_slice = brain_slice_real){
  W_pos_vec <- as.numeric(W > 0)
  W_pos_slice <- vector_to_slice(W_pos_vec)
  W_pos_slice[W_pos_slice == 0] <- NA
  
  W_neg_vec <- as.numeric(W < 0)
  W_neg_slice <- vector_to_slice(W_neg_vec)
  W_neg_slice[W_neg_slice == 0] <- NA
  
  target_vec <- rep(NA, length(W))
  target_vec[signal_idx] <- 1
  target_slice <- vector_to_slice(target_vec)
  
  # draw pictures
  x <- seq_len(nrow(brain_slice))
  y <- seq_len(ncol(brain_slice))
  
  image(
    x, y, brain_slice, col = gray.colors(256), asp = 1, axes = FALSE,
    xlab = "", ylab = "", main = "Target + Positive W + Negative W"
  )
  
  image(
    x, y, target_slice, col = rgb(1, 0, 0, alpha =.35), zlim = c(0,1),
    add = TRUE
  )
  
  if (any(is.finite(W_pos_slice))){
    image(
      x, y, W_pos_slice, col = rgb(0, 0, 1, alpha =.5), zlim = c(0,1), add = TRUE
    )
  }
  
  if (any(is.finite(W_neg_slice))){
    image(
      x, y, W_neg_slice, col = rgb(0, 1, 0, alpha = 0.5), zlim = c(0,1), add = TRUE
    )
  }
}

W_pos_overlay <- function(W, signal_idx, brain_slice = brain_slice_real){
  W_vec <- as.numeric(W > 0)
  W_slice <- vector_to_slice(W_vec)
  W_slice[W_slice == 0] <- NA
  target_vec <- rep(NA, length(W))
  target_vec[signal_idx] <- 1
  target_slice <- vector_to_slice(target_vec)
  
  # draw pictures
  x <- seq_len(nrow(brain_slice))
  y <- seq_len(ncol(brain_slice))
  
  image(
    x, y, brain_slice, col = gray.colors(256), asp = 1, axes = FALSE,
    xlab = "", ylab = "", main = "Target + Positive W"
  )
  
  image(
    x, y, target_slice, col = rgb(1, 0, 0, alpha =.35), zlim = c(0,1),
    add = TRUE
  )
  
  if (any(is.finite(W_slice))){
    image(
      x, y, W_slice, col = rgb(0, 0, 1, alpha =.5), zlim = c(0,1), add = TRUE
    )
  }
}

W_neg_overlay <- function(W, signal_idx, brain_slice = brain_slice_real){
  W_vec <- as.numeric(W < 0)
  W_slice <- vector_to_slice(W_vec)
  W_slice[W_slice == 0] <- NA
  target_vec <- rep(NA, length(W))
  target_vec[signal_idx] <- 1
  target_slice <- vector_to_slice(target_vec)
  
  # draw pictures
  x <- seq_len(nrow(brain_slice))
  y <- seq_len(ncol(brain_slice))
  
  image(
    x, y, brain_slice, col = gray.colors(256), asp = 1, axes = FALSE,
    xlab = "", ylab = "", main = "Target + Negative W"
  )
  
  image(
    x, y, target_slice, col = rgb(1, 0, 0, alpha =.35), zlim = c(0,1),
    add = TRUE
  )
  
  if (any(is.finite(W_slice))){
    image(
      x, y, W_slice, col = rgb(0, 1, 0, alpha = 0.5), zlim = c(0,1), add = TRUE
    )
  }
}