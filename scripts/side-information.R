library(oro.nifti)

mask_file <- "../MNI152_T1_2mm_brain_mask.nii.gz"

mask <- readNIfTI(mask_file)
mask_2d <- mask[, , 46] > 0

coords_2d <- which(mask_2d, arr.ind = TRUE)
colnames(coords_2d) <- c("coord_X", "coord_Y")

saveRDS(
    coords_2d, "coords.rds"
)