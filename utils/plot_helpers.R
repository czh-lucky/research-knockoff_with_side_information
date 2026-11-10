plot_W_ordered <- function(
    W,
    main = "Vanilla Knockoff",
    positive_col = "#56B4E9",
    negative_col = "#17365D",
    show_threshold = FALSE,
    fdr = 0.1,
    offset = 1
) {
  
  # 1. 只保留有限的非零 W
  nonzero_idx <- which(is.finite(W) & W != 0)
  
  W_nonzero <- W[nonzero_idx]
  
  # 如果没有非零 W
  if (length(W_nonzero) == 0) {
    plot.new()
    title(main = main)
    text(0.5, 0.5, "No nonzero W values")
    
    return(invisible(NULL))
  }
  
  # 2. 按照 |W| 从小到大排序
  ord <- order(abs(W_nonzero))
  
  W_sorted <- W_nonzero[ord]
  
  # 保留原始体素索引
  ordered_idx <- nonzero_idx[ord]
  
  # 3. 根据 W 的正负设置颜色
  W_col <- ifelse(
    W_sorted > 0,
    positive_col,
    negative_col
  )
  
  # 4. 设置横坐标
  x_order <- seq_along(W_sorted)
  
  # 5. 设置纵轴范围
  yr <- range(c(0, W_sorted))
  padding <- 0.05 * diff(yr)
  
  # 6. 创建空白绘图区域
  plot(
    x_order,
    W_sorted,
    type = "n",
    xlab = "Ordering (increasing |W|)",
    ylab = "W",
    main = main,
    ylim = yr + c(-padding, padding)
  )
  
  # 7. 绘制 W 的竖线
  segments(
    x0 = x_order,
    y0 = 0,
    x1 = x_order,
    y1 = W_sorted,
    col = W_col,
    lwd = 1.5
  )
  
  # 8. W = 0 参考线
  abline(
    h = 0,
    col = "black",
    lwd = 0.8
  )
  
  # 9. 图例
  legend(
    "topleft",
    legend = c("W > 0", "W < 0"),
    col = c(positive_col, negative_col),
    lwd = 2,
    bty = "n"
  )
  
  # 10. 可选：添加 Knockoff+ 阈值
  tau <- NULL
  
  if (show_threshold) {
    
    tau <- knockoff::knockoff.threshold(
      W,
      fdr = fdr,
      offset = offset
    )
    
    if (is.finite(tau) && tau <= max(yr)) {
      
      abline(
        h = tau,
        col = "red",
        lty = 2,
        lwd = 1.5
      )
      
    } else if (is.infinite(tau)) {
      
      mtext(
        "No discoveries: threshold = Inf",
        side = 3,
        col = "red",
        cex = 0.8
      )
    }
  }
  
  # 11. 返回排序结果，方便后续分析
  invisible(list(
    W_sorted = W_sorted,
    ordered_idx = ordered_idx,
    colors = W_col,
    threshold = tau
  ))
}
