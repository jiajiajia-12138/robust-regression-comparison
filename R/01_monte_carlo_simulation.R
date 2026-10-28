# =====================================================================
# 蒙特卡洛模拟：OLS 与 LAD / Huber / MM / LTS 四类稳健估计的比较
#
# 设计与课程论文一致：
#   y = b0 + b1 * x + eps,  真实参数 b0 = b1 = 1
#   x  ~ Uniform(0, 5)，n = 100，重复 M = 1000
#   6 种误差情景：正态 / 拉普拉斯 / t(3) / t(5) / 混合正态 / 柯西
#   评价指标：beta0、beta1 估计的 bias 与 MSE（重点看 beta1 的 MSE）
#
# 用法：
#   Rscript R/01_monte_carlo_simulation.R          # 默认 M = 1000
#   Rscript R/01_monte_carlo_simulation.R 100      # 指定 M，快速试跑
# =====================================================================

args <- commandArgs(trailingOnly = TRUE)
M <- if (length(args) >= 1) as.integer(args[1]) else 1000
n <- 100

set.seed(2024)

suppressPackageStartupMessages({
  library(MASS)        # rlm (Huber M 估计)
  library(quantreg)    # rq (LAD = tau 0.5 分位数回归)
  library(robustbase)  # lmrob (MM) / ltsReg (LTS)
})

# ---- 六种误差情景 -------------------------------------------------
gen_eps <- list(
  normal   = function(m) rnorm(m),
  laplace  = function(m) sample(c(-1, 1), m, TRUE) * rexp(m),   # Laplace(0,1)
  t3       = function(m) rt(m, df = 3),
  t5       = function(m) rt(m, df = 5),
  mixture  = function(m) ifelse(runif(m) < 0.95, rnorm(m), rnorm(m, 5, 10)),  # 0.95N(0,1)+0.05N(5,100)
  cauchy   = function(m) rcauchy(m)
)
scn_names <- c("正态", "拉普拉斯", "t(3)", "t(5)", "混合正态", "柯西")

# ---- 五种估计方法（返回长度 2 的 c(b0, b1)）-----------------------
fit5 <- function(x, y) {
  c(OLS   = coef(lm(y ~ x)),
    LAD   = coef(rq(y ~ x, tau = 0.5)),
    Huber = coef(rlm(y ~ x)),
    MM    = coef(lmrob(y ~ x)),
    LTS   = coef(ltsReg(y ~ x)))
}
mth_names <- c("OLS", "LAD", "Huber", "MM", "LTS")

# ---- 主循环 -------------------------------------------------------
est <- array(NA_real_, dim = c(M, 5, 2, 6),
             dimnames = list(rep = NULL, method = mth_names,
                             coef = c("b0", "b1"), scenario = scn_names))
t0 <- proc.time()[3]
for (s in seq_along(gen_eps)) {
  for (m in seq_len(M)) {
    x <- runif(n, 0, 5)
    y <- 1 + x + gen_eps[[s]](n)
    cf <- suppressWarnings(tryCatch(fit5(x, y), error = function(e) rep(NA_real_, 10)))
    est[m, , , s] <- matrix(cf, nrow = 5, ncol = 2, byrow = TRUE)  # 每行 = 一个方法的 (b0, b1)
  }
}
cat(sprintf("完成 %d 次重复 x %d 种情景，用时 %.1f 秒；失败拟合 %d 个\n",
            M, length(gen_eps), proc.time()[3] - t0,
            sum(is.na(est[, "OLS", 1, ])) ))

# ---- bias 与 MSE --------------------------------------------------
res <- data.frame()
for (s in scn_names) for (mt in mth_names) for (cf in c("b0", "b1")) {
  v <- est[, mt, cf, s]
  v <- v[!is.na(v)]
  res <- rbind(res, data.frame(
    scenario = s, method = mt, coef = cf,
    bias = mean(v) - 1,
    mse  = mean((v - 1)^2)
  ))
}

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/mse_simulation.csv", row.names = FALSE)

mse1 <- res[res$coef == "b1", c("scenario", "method", "mse")]
tab1 <- reshape(mse1, idvar = "method", timevar = "scenario", direction = "wide")
tab1 <- tab1[match(mth_names, tab1$method), ]
names(tab1) <- c("method", scn_names)
m <- as.matrix(tab1[, -1]); rownames(m) <- tab1$method
cat("\n===== beta1 的 MSE（方法 x 情景）=====\n")
print(round(m, 5))
cat("\n论文报告值（Table 1, M=1000）：正态下 OLS 最小 0.00484；柯西下 OLS 79.31，\n",
    "稳健方法 0.013~0.018；混合正态下 MM 最优 0.00580。因随机数种子不同，\n",
    "重新运行结果会有细微差异。\n", sep = "")
