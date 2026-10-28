# =====================================================================
# 实证分析：R 内置 cars 数据（50 辆汽车车速 speed 与刹车距离 dist）
#
# 内容：
#   1) 五种方法回归系数估计（对应论文表2）
#   2) OLS 标准化残差诊断（对应论文 5.3 节）
#   3) 污染实验：第 50 号观测 (25, 85) 的 dist 增加 60（对应论文表3）
#   4) 图2：散点图 + 五种方法拟合直线
#
# 用法：Rscript R/02_empirical_cars.R
# =====================================================================

set.seed(1)
options(digits = 6)
suppressPackageStartupMessages({
  library(MASS); library(quantreg); library(robustbase)
})
data(cars)

f5 <- function(d) rbind(
  OLS   = coef(lm(dist ~ speed, data = d)),
  LAD   = coef(rq(dist ~ speed, data = d, tau = 0.5)),
  Huber = coef(rlm(dist ~ speed, data = d)),
  MM    = coef(lmrob(dist ~ speed, data = d)),
  LTS   = coef(ltsReg(dist ~ speed, data = d)))

cat("===== 表2 cars 五方法估计 =====\n")
fit0 <- f5(cars); print(round(fit0, 4))

cat("\n===== OLS 标准化残差绝对值大于 2 的观测 =====\n")
rstd <- rstandard(lm(dist ~ speed, data = cars))
print(data.frame(obs = which(abs(rstd) > 2),
                 speed = cars$speed[abs(rstd) > 2],
                 dist  = cars$dist[abs(rstd) > 2],
                 rstd  = round(rstd[abs(rstd) > 2], 2)))

cat("\n===== 表3 污染实验：第50号观测(25,85)的dist增加60 =====\n")
cars2 <- cars
cars2$dist[50] <- cars2$dist[50] + 60
fit1 <- f5(cars2); print(round(fit1, 4))

cat("\n===== OLS 污染前后变化（截距 -17.58 -> -22.85, 斜率 3.93 -> 4.35）=====\n")
print(round(cbind(before = fit0["OLS", ], after = fit1["OLS", ],
                  delta = fit1["OLS", ] - fit0["OLS", ]), 4))

# ---- 图2 ----------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
font_cn <- if (.Platform$OS.type == "windows") "Microsoft YaHei" else "sans"
png("figures/fig2_cars_fits.png", width = 1800, height = 1300,
    res = 220, type = "cairo", family = font_cn)
plot(dist ~ speed, data = cars, pch = 19, col = "gray40",
     xlab = "车速 (mph)", ylab = "刹车距离 (英尺)",
     main = "cars数据：散点图与五种方法的拟合直线")
cols <- c(OLS = "blue", LAD = "green3", Huber = "orange", MM = "purple", LTS = "brown")
for (m in rownames(fit0)) abline(fit0[m, 1], fit0[m, 2], col = cols[m], lwd = 2)
legend("bottomright", legend = names(cols), col = cols, lwd = 2, bty = "n")
dev.off()
cat("\nfigures/fig2_cars_fits.png saved\n")
