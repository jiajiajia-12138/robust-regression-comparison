# =====================================================================
# 图1：离群点对 OLS 与稳健回归的破坏作用演示
#
# 设计与论文图1一致：95 个正常点来自 y = 1 + x + N(0,1)，
# 另 5 个离群点的 y 被人为抬高约 10 个单位。
# 直观预期：OLS 拟合线被明显拉偏，LAD / Huber 基本贴合真实线。
#
# 用法：Rscript R/03_figure1_outlier_demo.R
# =====================================================================

set.seed(123)
suppressPackageStartupMessages({
  library(MASS); library(quantreg)
})

n0 <- 95; n1 <- 5
x0 <- runif(n0, 0, 5); y0 <- 1 + x0 + rnorm(n0)
x1 <- runif(n1, 0, 5); y1 <- 1 + x1 + rnorm(n1) + 10   # 离群点
x <- c(x0, x1); y <- c(y0, y1)

fits <- rbind(
  OLS   = coef(lm(y ~ x)),
  LAD   = coef(rq(y ~ x, tau = 0.5)),
  Huber = coef(rlm(y ~ x)))
cat("===== 各方法估计（真实 b0 = b1 = 1）=====\n")
print(round(fits, 3))

dir.create("figures", showWarnings = FALSE)
font_cn <- if (.Platform$OS.type == "windows") "Microsoft YaHei" else "sans"
png("figures/fig1_outlier_demo.png", width = 1800, height = 1300,
    res = 220, type = "cairo", family = font_cn)
plot(x, y, pch = 19, col = c(rep("gray40", n0), rep("red", n1)),
     xlab = "x", ylab = "y",
     main = "5个离群点（红色）对OLS与稳健回归的影响")
abline(a = 1, b = 1, col = "black", lwd = 2)                    # 真实线
cols <- c(OLS = "blue", LAD = "green3", Huber = "orange")
for (m in rownames(fits)) abline(fits[m, 1], fits[m, 2], col = cols[m], lwd = 2, lty = 2)
legend("topleft", legend = c("真实线 y=1+x", names(cols)),
       col = c("black", cols), lwd = 2, lty = c(1, 2, 2, 2), bty = "n")
dev.off()
cat("\nfigures/fig1_outlier_demo.png saved\n")
