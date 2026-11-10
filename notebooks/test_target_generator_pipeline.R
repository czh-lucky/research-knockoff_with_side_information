library(oro.nifti)
source("utils/target_generator.R")

# 1. 载入基础数据
coords_2d <- readRDS("./data/coords.rds")
exp_grid <- read.csv("config/experiments.csv")
mask <- readNIfTI("data/MNI152_T1_2mm_brain_mask.nii.gz")

id <- 7

# 2. 测试函数：获取第 1 个实验的靶点信息
job_1_res <- get_signal_index_from_grid(coords_2d, exp_grid, current_run_id = id)
cat("Job 1 Voxel Count:", length(job_1_res), "\n")

# 3. 测试 overlay 预览图
# 构造 3D NIfTI 并在脑切片上渲染
preview_target_on_slice(coords_2d, exp_grid, id) 