#!/usr/bin/env Rscript

#' @param coords_2d
#' @param exp_grid
#' @param current_run_id

library(oro.nifti)
source("utils/nifti_functions.R")

brain_2d_all <- readRDS("data/X_all_581x4711.rds")
brain_2d <- vector_to_nifti(brain_2d_all[1, ])

get_signal_index_from_grid <- function(coords_2d, exp_grid, current_run_id){
  specs <- exp_grid[exp_grid$run_id ==current_run_id, ]
  
  if (nrow(specs) == 0){
    stop(sprintf("Run ID %d not found in experiment grid!", current_run_id))
  }
  
  p <- nrow(coords_2d)
  combined_signal <- rep(FALSE, p)
  
  for (i in seq_len(nrow(specs))){
    cx <- specs$center_x[i]
    cy <- specs$center_y[i]
    r <- specs$radius[i]
    
    dist_sq <- (coords_2d[, "coord_X"] - cx)^2 + (coords_2d[, "coord_Y"] - cy)^2
    region <- dist_sq <= r^2
    
    combined_signal <- combined_signal | region
  }
  
  return(which(combined_signal))
}


preview_target_on_slice <- function(coords_2d, exp_grid, current_run_id, brain_2d_img = brain_2d){
  
  signal_idx <- get_signal_index_from_grid(coords_2d, exp_grid, current_run_id)
  
  p <- nrow(coords_2d)
  x_vec <- rep(NA, p)
  x_vec[signal_idx] <- 1
  
  target_nifti <- vector_to_nifti(x_vec)
  
  specs <- exp_grid[exp_grid$run_id == current_run_id, ]
  
  oro.nifti::overlay(
    brain_2d_img,
    target_nifti,
    z = 46, plot.type = "single", NA.y = TRUE, col.y = rgb(1, 0, 0, alpha = 0.35),
    main = sprintf("Run ID: %d | Voxels: %d",
                   current_run_id, length(signal_idx))
  )
  
}



