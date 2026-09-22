############################################################
# Sequential learning (theta = +1), plot only l_n
# Baseline (non-private) vs randomized release (epsilon = 1, 3)
############################################################
library(scales)  # 用于 alpha() 设置透明度
# -------- Non-private simulator (returns l_n series) --------
simulate_seq_learning <- function(N, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  l <- numeric(N + 1); l[1] <- 0
  for (n in 1:N) {
    ln <- l[n]
    t  <- - (sigma^2) * ln / 2
    s  <- rnorm(1, mean = +1, sd = sigma)
    a  <- if (s > t) +1L else -1L
    
    z1 <- (t - 1) / sigma
    z2 <- (t + 1) / sigma
    if (a == +1L) {
      pnum <- 1 - pnorm(z1)   # P(a=+1 | theta=+1)
      pden <- 1 - pnorm(z2)   # P(a=+1 | theta=-1)
    } else {
      pnum <- pnorm(z1)       # P(a=-1 | theta=+1)
      pden <- pnorm(z2)       # P(a=-1 | theta=-1)
    }
    l[n + 1] <- ln + log(pnum / pden)
  }
  l[1:N]  # return l_n for n=1..N
}

# ---- Closed-form P(x=-1 | l, mu) for randomized actions ----
p_x_minus1 <- function(l, sigma, eps, mu, a = 0.5) {
  t <- - (sigma^2) * l / 2
  Z  <- pnorm((t - mu) / sigma)
  B1 <- exp(-eps * t) * exp(mu * eps + 0.5 * eps^2 * sigma^2) *
    pnorm((t - mu - eps * sigma^2) / sigma)
  B2 <- exp(+eps * t) * exp(-mu * eps + 0.5 * eps^2 * sigma^2) *
    (1 - pnorm((t - mu + eps * sigma^2) / sigma))
  Z - a * B1 + a * B2
}

# -------- Randomized-release simulator (returns l_n) --------
simulate_seq_learning_rand <- function(N, sigma, eps, a = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  l <- numeric(N + 1); l[1] <- 0
  
  for (n in 1:N) {
    ln <- l[n]
    t  <- - (sigma^2) * ln / 2
    s  <- rnorm(1, mean = +1, sd = sigma)
    
    # baseline action
    a_det <- if (s > t) +1L else -1L
    
    # randomized released action x_n
    p_flip <- a * exp(-eps * abs(s - t))
    x_n <- if (runif(1) < p_flip) -a_det else a_det
    
    # Bayes update using released-action probabilities
    p_num <- if (x_n == -1L)
      p_x_minus1(l = ln, sigma = sigma, eps = eps, mu = +1, a = a)
    else
      1 - p_x_minus1(l = ln, sigma = sigma, eps = eps, mu = +1, a = a)
    p_den <- if (x_n == -1L)
      p_x_minus1(l = ln, sigma = sigma, eps = eps, mu = -1, a = a)
    else
      1 - p_x_minus1(l = ln, sigma = sigma, eps = eps, mu = -1, a = a)
    
    l[n + 1] <- ln + log(p_num / p_den)
  }
  l[1:N]
}

############################################################
# Run and plot (only l_n)
############################################################
N <- 100
sigma <- 3

ln_nonpriv <- simulate_seq_learning(N, sigma, seed = 5)
ln_eps1    <- simulate_seq_learning_rand(N, sigma, eps = 0.5, seed = 4)
ln_eps3    <- simulate_seq_learning_rand(N, sigma, eps = 1, seed = 4)

# 自动放大纵轴范围（留 10% 边距）
y_all <- c(ln_nonpriv, ln_eps1, ln_eps3)
yr <- range(y_all)
pad <- 0.10 * diff(yr)
ylim_use <- c(yr[1] - pad, yr[2] + pad)

# 只画 l_n
plot(1:N, ln_nonpriv, type = "l", lwd = 2,
     xlab = "Agent index n", ylab = expression(l[n]),
     ylim = ylim_use)
grid()
lines(1:N, ln_eps1, lwd = 2, lty = 2)
lines(1:N, ln_eps3, lwd = 2, lty = 3)

legend("bottomright",
       legend = c("non-private", expression(epsilon==1), expression(epsilon==3)),
       lty = c(1,2,3), lwd = 2, bty = "n")







# ---- log-safe 版本：P(x=-1 | l, sigma, eps, mu) ----
############################################################
# Plot l_n under epsilon ~ Uniform[0,3]
############################################################

# ---- log-safe: P(x=-1 | l, sigma, eps, mu) ----
p_x_minus1 <- function(l, sigma, eps, mu, a = 0.5) {
  t  <- - (sigma^2) * l / 2
  Z  <- pnorm((t - mu) / sigma)
  
  z1     <- (t - mu - eps * sigma^2) / sigma
  logB1  <- -eps * t + mu * eps + 0.5 * (eps^2) * (sigma^2) +
    pnorm(z1, log.p = TRUE)
  B1     <- exp(logB1)
  
  z2     <- (t - mu + eps * sigma^2) / sigma
  logB2  <-  eps * t - mu * eps + 0.5 * (eps^2) * (sigma^2) +
    pnorm(z2, lower.tail = FALSE, log.p = TRUE)
  B2     <- exp(logB2)
  
  Z - a * B1 + a * B2
}

# ---- E_{eps ~ U[lo,hi]}[P(x=-1 | ...)] with integrate->Simpson fallback ----
p_x_minus1_eps_uniform_ab <- function(l, sigma, mu, a = 0.5, lo = 0, hi = 3,
                                      rel.tol = 1e-7, n_grid = 801L) {
  stopifnot(hi > lo)
  f <- function(eps) p_x_minus1(l = l, sigma = sigma, eps = eps, mu = mu, a = a)
  
  ans <- try(integrate(f, lower = lo, upper = hi,
                       rel.tol = rel.tol, subdivisions = 2000L,
                       stop.on.error = FALSE), silent = TRUE)
  if (!inherits(ans, "try-error") && is.finite(ans$value)) {
    return( (1 / (hi - lo)) * ans$value )
  }
  
  # Simpson 复合公式（n_grid 取奇数）
  if (n_grid %% 2 == 0) n_grid <- n_grid + 1L
  xs <- seq(lo, hi, length.out = n_grid)
  fx <- vapply(xs, f, numeric(1))
  h  <- (hi - lo) / (n_grid - 1L)
  S  <- fx[1] + fx[n_grid] +
    4 * sum(fx[seq(2, n_grid - 1L, by = 2)]) +
    2 * sum(fx[seq(3, n_grid - 2L, by = 2)])
  (1 / (hi - lo)) * (h / 3) * S
}

# ---- Simulator: a_n -> (eps_n ~ U[lo,hi]) flip -> x_n; update l using marginalized probs ----
simulate_seq_learning_epsU_sample_ab <- function(N, sigma, a = 0.5,
                                                 lo , hi , seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  l <- numeric(N + 1); l[1] <- 0
  s <- numeric(N); a_det <- integer(N); x_rel <- integer(N)
  eps_n <- numeric(N); pflip <- numeric(N)
  
  for (n in 1:N) {
    ln <- l[n]
    t  <- - (sigma^2) * ln / 2
    
    # private signal (theta=+1) & baseline action
    s[n]     <- rnorm(1, mean = +1, sd = sigma)
    a_det[n] <- if (s[n] > t) +1L else -1L
    
    # draw eps_n ~ U[lo,hi], flip a_n to x_n with p_flip
    eps_n[n] <- runif(1, lo, hi)
    pflip[n] <- a * exp(-eps_n[n] * abs(s[n] - t))
    x_rel[n] <- if (runif(1) < pflip[n]) -a_det[n] else a_det[n]
    
    # marginalized (public) probabilities
    p_m1_t1  <- p_x_minus1_eps_uniform_ab(l = ln, sigma = sigma, mu = +1, a = a, lo = lo, hi = hi)
    p_m1_tm1 <- p_x_minus1_eps_uniform_ab(l = ln, sigma = sigma, mu = -1, a = a, lo = lo, hi = hi)
    
    if (x_rel[n] == -1L) {
      pnum <- p_m1_t1;      pden <- p_m1_tm1
    } else {
      pnum <- 1 - p_m1_t1;  pden <- 1 - p_m1_tm1
    }
    l[n + 1] <- ln + log(pnum / pden)
  }
  
  data.frame(
    n        = 1:N,
    x_n      = x_rel,
    a_n      = a_det,
    eps_n    = eps_n,
    p_flip_n = pflip,
    l_n      = l[1:N],
    l_next   = l[2:(N + 1)]
  )
}

# # 假设你之前已经有：
# # - simulate_seq_learning_rand(N, sigma, eps, seed=...)
# # - simulate_seq_learning_epsU_sample_ab(N, sigma, lo, hi, seed=...)
# #10 17
# N <- 100; sigma <- 3; seed <- 10
# 
# df_fixed <- simulate_seq_learning_rand(N, sigma, eps = 1, seed = seed)
# df_U03   <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = 0.5, hi = 1, seed = seed)
# 
# # 自动放大纵轴 10% 边距
# y_all <- c(df_fixed, df_U03$l_n)
# yr <- range(y_all); pad <- 0.10 * diff(yr); if (pad == 0) pad <- 1
# 
# plot(1:N, df_fixed, type = "l", lwd = 2,
#      xlab = "Agent index n", ylab = expression(l[n]),
#      ylim = c(yr[1] - pad, yr[2] + pad))
# grid()
# lines(df_U03$n, df_U03$l_n, lwd = 2, lty = 2)
# 
# legend("bottomright",
#        legend = c(expression(epsilon==1), expression(epsilon %~% U(0,1))),
#        lty = c(1,2), lwd = 2, bty = "n")
# 
# 
# N <- 200; sigma <- 3; seed <- 7
# 
# # 三个区间的随机 epsilon
# df_U05_1 <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = 0., hi = 0.5,   seed = seed)
# df_U1_15 <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = 0.5,   hi = 1,   seed = seed+1)
# df_U15_2 <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = 1,   hi = 1.5,   seed = seed+2)
# ln_nonpriv <- simulate_seq_learning(N, sigma, seed = 5)  # returns a numeric vector
# 
# # 自动放大纵轴 10% 边距
# y_all <- c(df_U05_1$l_n, df_U1_15$l_n, df_U15_2$l_n, ln_nonpriv)
# yr <- range(y_all); pad <- 0.10 * diff(yr); if (pad == 0) pad <- 1
# 
# # 先画第一条曲线
# plot(df_U05_1$n, df_U05_1$l_n, type = "l", lwd = 2, col = "red",
#      xlab = "n", ylab = expression(l[n]),
#      ylim = c(yr[1] - pad, yr[2] + pad))
# grid()
# 
# # 添加另外三条曲线
# lines(df_U1_15$n,  df_U1_15$l_n,  lwd = 2, lty = 2, col = "blue")
# lines(df_U15_2$n, df_U15_2$l_n, lwd = 2, lty = 3, col = "darkgreen")
# lines(1:N, ln_nonpriv, lwd = 2, lty = 4, col = "black")  # Non-private line
# 
# # 添加图例
# legend("bottomright",
#        legend = c(expression(Low),
#                   expression(Medium),
#                   expression(High),
#                   "Non-private"),
#        lty = c(1,2,3,4), lwd = 2,
#        col = c("red","blue","darkgreen","black"),
#        bty = "n")
# 
# 
# 
# N <- 2000; sigma <- 3; seed <- 100
# nsim <- 5  # 每种情况模拟次数
# 
# # 自动放大纵轴前先收集所有数据
# all_vals <- c()
# 
# # 创建画布（先画空图）
# plot(1, type = "n",
#      xlim = c(1, N),
#      ylim = c(-3.5, 5.5),  # 会根据你实际情况调整，下面有更自动的做法
#      xlab = "n", ylab = expression(l[n]))
# grid()
# 
# # 模拟 + 绘图函数
# plot_eps_range <- function(lo, hi, base_seed, col, lty) {
#   for (i in 0:(nsim - 1)) {
#     df <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = lo, hi = hi, seed = base_seed + i)
#     lines(df$n, df$l_n, lwd = 1.5, col = col, lty = lty)
#   }
# }
# 
# # 1. Low (红色, 实线)
# plot_eps_range(0.0, 0.5, seed, col = "red", lty = 1)
# 
# # 2. Medium (蓝色, 虚线)
# plot_eps_range(0.5, 1.0, seed + 10, col = "blue", lty = 2)
# 
# # 3. High (深绿, 点虚线)
# plot_eps_range(1.0, 1.5, seed + 20, col = "darkgreen", lty = 3)
# 
# # 4. Non-private (黑色, 长虚线)
# for (i in 0:(nsim - 1)) {
#   ln <- simulate_seq_learning(N, sigma, seed = 100 + i)
#   lines(1:N, ln, lwd = 1.5, col = "black", lty = 4)
# }
# 
# # 添加图例
# legend("bottomright",
#        legend = c(expression(Low),
#                   expression(Medium),
#                   expression(High),
#                   "Non-private"),
#        lty = c(1,2,3,4), lwd = 2,
#        col = c("red","blue","darkgreen","black"),
#        bty = "n")
# 
# 
# 
# 
# 
# 
# N <- 2000; sigma <- 2; seed <- 100
# nsim <- 5  # 每种情况模拟次数
# 
# # 创建空图
# plot(1, type = "n", xlim = c(1, N), ylim = c(-3.5, 10.5),
#      xlab = "n", ylab = expression(l[n]))
# grid()
# 
# # 找出“最中间的那一条曲线”（整体最接近平均）
# get_center_index <- function(mat) {
#   row_means <- rowMeans(mat)
#   avg <- mean(row_means)
#   which.min(abs(row_means - avg))  # 离平均最近的那一条
# }
# 
# # 通用函数：绘制某组 epsilon 区间的曲线
# plot_eps_group <- function(lo, hi, base_seed, col_base, lty) {
#   all_lns <- matrix(NA, nrow = nsim, ncol = N)
#   
#   for (i in 1:nsim) {
#     df <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = lo, hi = hi, seed = base_seed + i)
#     all_lns[i, ] <- df$l_n
#   }
#   
#   center_idx <- get_center_index(all_lns)
#   
#   for (i in 1:nsim) {
#     col_use <- if (i == center_idx) col_base else alpha(col_base, 0.3)
#     lw_use  <- if (i == center_idx) 2.5 else 1
#     lines(1:N, all_lns[i, ], col = col_use, lwd = lw_use, lty = lty)
#   }
# }
# 
# # Non-private 特例
# plot_nonpriv_group <- function(base_seed, col_base, lty) {
#   all_lns <- matrix(NA, nrow = nsim, ncol = N)
#   for (i in 1:nsim) {
#     all_lns[i, ] <- simulate_seq_learning(N, sigma, seed = base_seed + i)
#   }
#   
#   center_idx <- get_center_index(all_lns)
#   
#   for (i in 1:nsim) {
#     col_use <- if (i == center_idx) col_base else alpha(col_base, 0.3)
#     lw_use  <- if (i == center_idx) 2.5 else 1
#     lines(1:N, all_lns[i, ], col = col_use, lwd = lw_use, lty = lty)
#   }
# }
# 
# # 绘图
# plot_eps_group(0.0, 0.5, seed,       col_base = "red",       lty = 1)
# plot_eps_group(0.5, 1.0, seed + 10,  col_base = "blue",      lty = 2)
# plot_eps_group(1.0, 1.5, seed + 20,  col_base = "darkgreen", lty = 3)
# plot_nonpriv_group(seed + 30,       col_base = "black",     lty = 4)
# 
# # 图例
# legend("bottomright",
#        legend = c(expression(Low),
#                   expression(Medium),
#                   expression(High),
#                   "Non-private"),
#        lty = c(1,2,3,4), lwd = 2,
#        col = c("red","blue","darkgreen","black"),
#        bty = "n")






N <- 2000; sigma <- 3; seed <- 6
nsim <- 5  # number of simulations per setting

# 空图
plot(1, type = "n", xlim = c(1, N), ylim = c(-3.5, 7.5),
     xlab = expression(italic(n)), ylab = expression(l[n]))
grid()

# 找“中间曲线”的函数
get_center_index <- function(mat) {
  row_means <- rowMeans(mat)
  avg <- mean(row_means)
  which.min(abs(row_means - avg))
}

# 通用函数: uniform epsilon
plot_uniform_group <- function(lo, hi, base_seed, col_base, lty) {
  all_lns <- matrix(NA, nrow = nsim, ncol = N)
  for (i in 1:nsim) {
    df <- simulate_seq_learning_epsU_sample_ab(N, sigma, lo = lo, hi = hi, seed = base_seed + i)
    all_lns[i, ] <- df$l_n
  }
  center_idx <- get_center_index(all_lns)
  for (i in 1:nsim) {
    col_use <- if (i == center_idx) col_base else alpha(col_base, 0.3)
    lw_use  <- if (i == center_idx) 2.5 else 1
    lines(1:N, all_lns[i, ], col = col_use, lwd = lw_use, lty = lty)
  }
}

# 通用函数: fixed epsilon
plot_fixed_group <- function(eps, base_seed, col_base, lty) {
  all_lns <- matrix(NA, nrow = nsim, ncol = N)
  for (i in 1:nsim) {
    all_lns[i, ] <- simulate_seq_learning_rand(N, sigma, eps = eps, seed = base_seed + i)
  }
  center_idx <- get_center_index(all_lns)
  for (i in 1:nsim) {
    col_use <- if (i == center_idx) col_base else alpha(col_base, 0.3)
    lw_use  <- if (i == center_idx) 2.5 else 1
    lines(1:N, all_lns[i, ], col = col_use, lwd = lw_use, lty = lty)
  }
}

# Non-private
plot_nonpriv_group <- function(base_seed, col_base, lty) {
  all_lns <- matrix(NA, nrow = nsim, ncol = N)
  for (i in 1:nsim) {
    all_lns[i, ] <- simulate_seq_learning(N, sigma, seed = base_seed + i)
  }
  center_idx <- get_center_index(all_lns)
  for (i in 1:nsim) {
    col_use <- if (i == center_idx) col_base else alpha(col_base, 0.3)
    lw_use  <- if (i == center_idx) 2.5 else 1
    lines(1:N, all_lns[i, ], col = col_use, lwd = lw_use, lty = lty)
  }
}

# ==== 绘图五种情况 ====
plot_uniform_group(0.0, 1.0, seed,       col_base = "purple",   lty = 1)  # Uniform[0,1]
plot_fixed_group(0.1,     seed + 10, col_base = "red",      lty = 2)      # eps=0.5
plot_fixed_group(0.5,     seed + 30, col_base = "blue",     lty = 3)      # eps=1.0
plot_fixed_group(1,     seed + 40, col_base = "darkgreen",lty = 4)      # eps=1.5
plot_nonpriv_group(seed + 50,        col_base = "black",    lty = 5)      # Non-private

# 图例
legend("bottomright",
       legend = c("Uniform[0,1]",
                  expression(Low),
                  expression(Medium),
                  expression(High),
                  "Non-private"),
       lty = c(1,2,3,4,5), lwd = 2,
       col = c("purple","red","blue","darkgreen","black"),
       bty = "n")
