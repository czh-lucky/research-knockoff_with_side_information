library(oro.nifti)

img_dir <- "data/mni_images"
mask_file <- "data/MNI152_T1_2mm_brain_mask.nii.gz"

mask <- readNIfTI(mask_file)
mask_2d <- mask[,,46]
voxel_idx <- which(mask_2d > 0)

img_files <- list.files(
    img_dir,
    pattern = "_MNI_T1_2mm\\.nii\\.gz",
    full.names = TRUE
)

X_all <- matrix(
    NA_real_,
    nrow = length(img_files),
    ncol = length(voxel_idx)
)
system.time({
    for (i in seq_along(img_files)) {
        img <- readNIfTI(img_files[i])
        img_2d <- img[,,46]
        X_all[i, ] <- img_2d[mask_2d > 0]
    }
})
saveRDS(X_all, "X_all_581x4711.rds")

