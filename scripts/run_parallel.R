library(parallel)

# 开启5个并行进程
n_cores <- 5
cl <- makeCluster(n_cores)

# 动态获取当下环境的 Rscript 路径
rscript_path <- file.path(R.home("bin"),"Rscript")

cat(sprintf("[%s] 已启动 %d 个并行进程，开始运行前 10 组模拟(1~10)...\n", Sys.time(), n_cores))

# 运行 1 到 10 组
parLapply(cl, 1:10, function(i, exec_path){
  
  log_file <- sprintf("./Xk/b%s.log", i)

  system2(
    command = exec_path,
    args = c("generate_Xk.R", as.character(i)),
    stdout = log_file,
    stderr = log_file
  )
}, exec_path = rscript_path)

# 释放资源
stopCluster(cl)

cat("Finished.\n")
# cat(sprintf("[%s] 前 10 组模拟已全部计算完成！", Sys.time()))