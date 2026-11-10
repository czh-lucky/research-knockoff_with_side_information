library(oro.nifti)

mask <- readNIfTI("data/MNI152_T1_2mm_brain_mask.nii.gz")

k <- 46
mask_2d <- mask[ , , k] > 0


vector_to_slice <- function(x, mask2d = mask_2d) {
  stopifnot(length(x) == sum(mask2d))
  
  result <- matrix(
    NA,
    nrow = nrow(mask2d),
    ncol = ncol(mask2d)
  )
  
  result[mask2d] <- x
  
  return(result)
}


vector_to_nifti <- function(x, mask_img = mask, k = 46) {
  brain_array <- array(0, dim = dim(mask_img))
  
  mask2d <- mask_img[ , , k] > 0
  
  brain_array[ , , k] <- vector_to_slice(x, mask2d)
  
  brain <- as.nifti(brain_array, value = mask_img)
  
  return(brain)
}