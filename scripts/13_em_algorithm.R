# 13_em_algorithm.R — 4-режимная копула-HMM, алгоритм Баума–Велша (EM)
#
# Состояния:
#   1 = Frank    — симметричная умеренная положительная зависимость (норма)
#   2 = Gumbel   — верхний хвост (нефтяной бум, рубль укрепляется совместно)
#   3 = Student  — тяжёлые хвосты, кризис (ρ < 0)
#   4 = Gumbel°  — 90°-поворот, отрицательная зависимость (санкции)
#
# E-шаг : фильтр Гамильтона + сглаживатель Кима
#          → ξ_{t|T}(s)  и  ξ_{t,t+1|T}(s,s')
# M-шаг : обновляем P (Baum-Welch) + θ_s (взвешенный MLE)

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr); library(tidyr); library(copula); library(ggplot2)

# задайте рабочую директорию проекта перед запуском
PARENT  <- ".."
OUT_DIR <- "results/13_em_algorithm"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source("scripts/00_helpers.R")

# ДАННЫЕ

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

u_all <- as.data.frame(pobs(cbind(lr$log_OIL, lr$log_FX)))
colnames(u_all) <- c("u_oil", "u_fx")
u_all$Date <- lr$Date
T_obs <- nrow(u_all)
u_mat <- as.matrix(u_all[, c("u_oil","u_fx")])
u_rot <- cbind(1 - u_mat[,1], u_mat[,2])   # для Gumbel° (90°-поворот)

BREAK_DATES <- as.Date(c("2020-01-18","2020-06-20","2022-02-26",
                          "2022-06-11","2023-09-09"))

K           <- 4
STATE_NAMES <- c("Frank","Gumbel","Student","Gumbel°")

cat(sprintf("Данные: N=%d  (%s — %s)\n", T_obs,
            format(min(lr$Date)), format(max(lr$Date))))

# ПЛОТНОСТИ КОПУЛ

compute_densities <- function(u_mat, u_rot, params) {
  n <- nrow(u_mat)
  safe_d <- function(cop, U) tryCatch(
    pmax(dCopula(U, cop), 1e-300),
    error = function(e) rep(1e-300, n)
  )
  f1 <- safe_d(frankCopula(param = params$frank, dim = 2), u_mat)
  f2 <- safe_d(gumbelCopula(param = params$gumbel, dim = 2), u_mat)
  f3 <- safe_d(tCopula(param = params$rho, df = params$df, dim = 2), u_mat)
  f4 <- safe_d(gumbelCopula(param = params$gumbel_rot, dim = 2), u_rot)
  cbind(f1, f2, f3, f4)
}

# 3. ФИЛЬТР ГАМИЛЬТОНА (K состояний)

hamilton_filter <- function(f_mat, P, pi0 = NULL) {
  K <- ncol(f_mat); T <- nrow(f_mat)
  if (is.null(pi0)) {
    ev  <- eigen(t(P))$vectors[, 1]
    pi0 <- pmax(Re(ev) / sum(Re(ev)), 0)
    pi0 <- pi0 / sum(pi0)
  }
  xi_filt <- matrix(0, T, K)
  xi_pred <- matrix(0, T, K)
  LL_vec  <- numeric(T)
  p_pred  <- pi0

  for (t in 1:T) {
    xi_pred[t,] <- p_pred
    joint <- f_mat[t,] * p_pred
    denom <- sum(joint)
    if (denom <= 0 || is.nan(denom)) {
      xi_filt[t,] <- p_pred
      LL_vec[t]   <- log(.Machine$double.eps)
    } else {
      xi_filt[t,] <- joint / denom
      LL_vec[t]   <- log(denom)
    }
    p_pred <- pmax(as.vector(t(P) %*% xi_filt[t,]), 0)
    p_pred <- p_pred / sum(p_pred)
  }
  list(xi_filt = xi_filt, xi_pred = xi_pred, LL = sum(LL_vec))
}

# 4. СГЛАЖИВАТЕЛЬ КИМА + СОВМЕСТНЫЕ АПОСТЕРИОРНЫЕ ξ_{t,t+1|T}

kim_smoother <- function(xi_filt, xi_pred, P) {
  K <- ncol(xi_filt); T <- nrow(xi_filt)

  xi_smooth <- matrix(0, T, K)
  xi_smooth[T,] <- xi_filt[T,]

  # xi2[t, j, k] = P(S_t=j, S_{t+1}=k | data) для M-шага
  xi2 <- array(0, c(T - 1, K, K))

  for (t in (T - 1):1) {
    ratio <- xi_smooth[t + 1,] /
             pmax(xi_pred[t + 1,], .Machine$double.eps)
    for (j in 1:K) {
      xi_smooth[t, j] <- xi_filt[t, j] * sum(P[j,] * ratio)
      xi2[t, j,]      <- xi_filt[t, j] * P[j,] * ratio
    }
    s <- sum(xi_smooth[t,])
    if (s > 0) {
      xi_smooth[t,] <- xi_smooth[t,] / s
      xi2[t,,]      <- xi2[t,,]      / s
    }
  }
  list(xi_smooth = xi_smooth, xi2 = xi2)
}

# 5. M-ШАГ: ОБНОВЛЕНИЕ P И θ

update_P <- function(xi2) {
  # Baum-Welch: p̂_{jk} = Σ_t ξ_{t,t+1}(j,k) / Σ_t ξ_{t}(j)
  K <- dim(xi2)[2]
  P_new <- matrix(0, K, K)
  for (j in 1:K) {
    denom <- sum(xi2[, j,])
    P_new[j,] <- if (denom > 0) colSums(xi2[, j,]) / denom else rep(1/K, K)
  }
  pmax(P_new, 1e-6) / rowSums(pmax(P_new, 1e-6))
}

wll <- function(d, w) -sum(w * log(pmax(d, 1e-300)))

update_copula_params <- function(u_mat, u_rot, w_list, p_old) {
  p_new <- p_old

  # Frank (θ > 0)
  if (sum(w_list[[1]]) > 2) {
    opt <- tryCatch(
      optimize(function(th) wll(dCopula(u_mat,frankCopula(th,2)), w_list[[1]]),
               interval = c(0.01, 25), tol = 1e-5),
      error = function(e) NULL)
    if (!is.null(opt)) p_new$frank <- opt$minimum
  }

  # Gumbel (θ ≥ 1)
  if (sum(w_list[[2]]) > 2) {
    opt <- tryCatch(
      optimize(function(th) wll(dCopula(u_mat,gumbelCopula(th,2)), w_list[[2]]),
               interval = c(1.001, 12), tol = 1e-5),
      error = function(e) NULL)
    if (!is.null(opt)) p_new$gumbel <- opt$minimum
  }

  # Student (ρ, df)
  if (sum(w_list[[3]]) > 2) {
    opt <- tryCatch(
      optim(c(p_old$rho, p_old$df),
            function(par) wll(dCopula(u_mat,
              tCopula(par[1], df=par[2], dim=2)), w_list[[3]]),
            method = "L-BFGS-B",
            lower = c(-0.999, 0.3), upper = c(0.999, 20)),
      error = function(e) NULL)
    if (!is.null(opt) && opt$convergence == 0) {
      p_new$rho <- opt$par[1]
      p_new$df  <- opt$par[2]
    }
  }

  # Gumbel° (на повёрнутых данных)
  if (sum(w_list[[4]]) > 2) {
    opt <- tryCatch(
      optimize(function(th) wll(dCopula(u_rot,gumbelCopula(th,2)), w_list[[4]]),
               interval = c(1.001, 12), tol = 1e-5),
      error = function(e) NULL)
    if (!is.null(opt)) p_new$gumbel_rot <- opt$minimum
  }

  p_new
}

# em-АЛГОРИТМ

cat("\n=== EM-АЛГОРИТМ (4-режимная копула-HMM) ===\n")

params <- list(frank = 2.5066, gumbel = 1.50,
               rho = -0.6347, df = 0.5551, gumbel_rot = 1.0869)

P_em <- matrix(c(
  0.90, 0.05, 0.03, 0.02,
  0.05, 0.88, 0.05, 0.02,
  0.03, 0.03, 0.90, 0.04,
  0.01, 0.02, 0.02, 0.95
), K, K, byrow = TRUE)

MAX_ITER <- 150
TOL      <- 1e-5
LL_prev  <- -Inf
LL_hist  <- numeric(MAX_ITER)
n_iter   <- MAX_ITER

for (iter in 1:MAX_ITER) {
  f_mat_em  <- compute_densities(u_mat, u_rot, params)
  res_filt  <- hamilton_filter(f_mat_em, P_em)
  res_smth  <- kim_smoother(res_filt$xi_filt, res_filt$xi_pred, P_em)

  LL_curr       <- res_filt$LL
  LL_hist[iter] <- LL_curr

  cat(sprintf("Iter %3d | LL = %10.4f | ΔLL = %+.6f | params: F=%.3f G=%.3f G°=%.3f\n",
              iter, LL_curr, LL_curr - LL_prev,
              params$frank, params$gumbel, params$gumbel_rot))

  if (abs(LL_curr - LL_prev) < TOL && iter > 5) {
    cat(sprintf("Сходимость на итерации %d (ΔLL < %.1e)\n", iter, TOL))
    n_iter <- iter; break
  }

  P_em    <- update_P(res_smth$xi2)
  params  <- update_copula_params(u_mat, u_rot,
                                  lapply(1:K, function(k) res_smth$xi_smooth[,k]),
                                  params)
  LL_prev <- LL_curr
}

# Финальный прогон
f_mat_em  <- compute_densities(u_mat, u_rot, params)
res_filt  <- hamilton_filter(f_mat_em, P_em)
res_smth  <- kim_smoother(res_filt$xi_filt, res_filt$xi_pred, P_em)
xi_smooth <- res_smth$xi_smooth

# AIC / BIC (12 переходных параметров + 5 параметров копул = 17)
N_PAR <- K * (K - 1) + 5
AIC_em <- 2 * N_PAR - 2 * res_filt$LL
BIC_em <- N_PAR * log(T_obs) - 2 * res_filt$LL

rownames(P_em) <- colnames(P_em) <- STATE_NAMES

cat(sprintf("\n=== ИТОГОВЫЕ ПАРАМЕТРЫ EM ===\n"))
cat(sprintf("LogLik = %.4f  AIC = %.4f  BIC = %.4f\n",
            res_filt$LL, AIC_em, BIC_em))
cat(sprintf("Параметры копул:\n"))
cat(sprintf("  Frank:    θ = %.4f\n", params$frank))
cat(sprintf("  Gumbel:   θ = %.4f   (τ ≈ %.3f)\n",
            params$gumbel, 1 - 1/params$gumbel))
cat(sprintf("  Student:  ρ = %.4f, df = %.4f\n", params$rho, params$df))
cat(sprintf("  Gumbel°:  θ = %.4f   (τ ≈ %.3f)\n",
            params$gumbel_rot, 1 - 1/params$gumbel_rot))

cat("\nМатрица переходов P_EM:\n")
print(round(P_em, 4))

ev_stat <- eigen(t(P_em))$vectors[, 1]
pi_stat <- pmax(Re(ev_stat) / sum(Re(ev_stat)), 0)
pi_stat <- pi_stat / sum(pi_stat)
dur_em  <- 1 / (1 - diag(P_em))

cat("\nСтационарные вероятности и средние длительности:\n")
for (k in 1:K) {
  cat(sprintf("  %-8s: π=%.3f  (%.1f%%)  ср.длит.=%.1f нед.\n",
              STATE_NAMES[k], pi_stat[k], pi_stat[k]*100, dur_em[k]))
}

# КРИВАЯ СХОДИМОСТИ

ll_df   <- data.frame(iter = seq_len(n_iter), LL = LL_hist[seq_len(n_iter)])
p_conv  <- ggplot(ll_df, aes(x = iter, y = LL)) +
  geom_line(color = "#1565C0", linewidth = 1) +
  geom_point(color = "#1565C0", size = 2) +
  labs(title = "Кривая сходимости EM-алгоритма (4-режимная HMM)",
       x = "Итерация", y = "Log-правдоподобие") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "13a_em_convergence.png"),
       p_conv, width = 8, height = 4, dpi = 150)

# АПОСТЕРИОРНЫЕ ВЕРОЯТНОСТИ

probs_df <- data.frame(
  Date           = u_all$Date,
  Frank_smooth   = xi_smooth[,1],
  Gumbel_smooth  = xi_smooth[,2],
  Student_smooth = xi_smooth[,3],
  GumbelR_smooth = xi_smooth[,4],
  Frank_filt     = res_filt$xi_filt[,1],
  Gumbel_filt    = res_filt$xi_filt[,2],
  Student_filt   = res_filt$xi_filt[,3],
  GumbelR_filt   = res_filt$xi_filt[,4],
  Regime         = STATE_NAMES[apply(xi_smooth, 1, which.max)]
)
write.csv(probs_df,
          file.path(OUT_DIR, "13b_regime_probs.csv"), row.names = FALSE)

regime_colors <- c("Frank"="#1565C0","Gumbel"="#2E7D32",
                   "Student"="#B71C1C","Gumbel°"="#FF8F00")

prob_long <- probs_df %>%
  select(Date, Frank_smooth, Gumbel_smooth, Student_smooth, GumbelR_smooth) %>%
  pivot_longer(-Date, names_to = "Regime", values_to = "Prob") %>%
  mutate(Regime = recode(Regime,
    Frank_smooth   = "Frank",
    Gumbel_smooth  = "Gumbel",
    Student_smooth = "Student",
    GumbelR_smooth = "Gumbel°"
  ))

# Стековая диаграмма
p_stack <- ggplot(prob_long, aes(x = Date, y = Prob, fill = Regime)) +
  geom_area(position = "stack", alpha = 0.85) +
  geom_vline(xintercept = BREAK_DATES, linetype = "dashed",
             color = "white", linewidth = 0.7) +
  scale_fill_manual(values = regime_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "4-режимная копула-HMM: апостериорные вероятности (EM)",
       subtitle = "Baum–Welch, совместная оценка P и θ",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "13c_stacked_probs.png"),
       p_stack, width = 13, height = 5, dpi = 150)

# Отдельные линии (overlay)
p_lines <- ggplot(prob_long, aes(x = Date, y = Prob, color = Regime)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  geom_vline(xintercept = BREAK_DATES, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  scale_color_manual(values = regime_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title    = "Апостериорные вероятности режимов (сглаженные)",
       subtitle = "EM-оценка, 4 состояния",
       x = NULL, y = "P(режим | данные)", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "13d_regime_lines.png"),
       p_lines, width = 13, height = 5, dpi = 150)

# 9. СРАВНЕНИЕ 3-STATE (ФИКС.) vs 4-STATE (EM)

# Из скрипта 12: 3-state с фиксированными параметрами
AIC_3 <- 2*6  - 2*(-30.13/2)     # грубая оценка — пересчитаем точно
# На самом деле пересчитаем с теми же данными
cop3_frank   <- frankCopula(2.5066, 2)
cop3_student <- tCopula(-0.6347, df=0.5551, 2)
cop3_gumbel  <- gumbelCopula(1.0869, 2)
f3 <- cbind(
  pmax(dCopula(u_mat, cop3_frank),   1e-300),
  pmax(dCopula(u_mat, cop3_student), 1e-300),
  pmax(dCopula(u_rot, cop3_gumbel),  1e-300)
)
# Загружаем P_hat из скрипта 12 если есть, иначе используем диагональную
p12_file <- "results/12_regime_switching/12b_transition_matrix.csv"
P3 <- if (file.exists(p12_file)) {
  as.matrix(read.csv(p12_file, row.names=1))
} else {
  matrix(c(0.949,0.028,0.023, 0.024,0.879,0.097, 0.003,0.000,0.997),
         3, 3, byrow=TRUE)
}
res3 <- hamilton_filter(f3, P3)
AIC_3_real <- 2*6  - 2*res3$LL
BIC_3_real <- 6*log(T_obs) - 2*res3$LL

cat(sprintf("\n=== СРАВНЕНИЕ МОДЕЛЕЙ ===\n"))
cat(sprintf("3-state (фикс. θ, скрипт 12): LogLik=%8.4f  AIC=%8.4f  BIC=%8.4f\n",
            res3$LL, AIC_3_real, BIC_3_real))
cat(sprintf("4-state EM   (скрипт 13):     LogLik=%8.4f  AIC=%8.4f  BIC=%8.4f\n",
            res_filt$LL, AIC_em, BIC_em))
cat(sprintf("ΔAIC = %.4f  ΔBIC = %.4f  (%s)\n",
            AIC_em - AIC_3_real, BIC_em - BIC_3_real,
            ifelse(AIC_em < AIC_3_real,
                   "4-state лучше по AIC", "3-state лучше по AIC")))

# Likelihood ratio test: 4-state vs 3-state (вложение не строгое, но ориентировочно)
LR_stat <- 2 * (res_filt$LL - res3$LL)
df_lr   <- N_PAR - 6
cat(sprintf("LR-статистика = %.4f  df=%d  p≈%.4f\n",
            LR_stat, df_lr, 1 - pchisq(LR_stat, df_lr)))

# СОХРАНЕНИЕ РЕЗУЛЬТАТОВ

write.csv(as.data.frame(round(P_em, 6)),
          file.path(OUT_DIR, "13e_transition_matrix.csv"))
write.csv(data.frame(
  Model   = c("3-state fixed", "4-state EM"),
  LogLik  = round(c(res3$LL, res_filt$LL), 4),
  N_par   = c(6, N_PAR),
  AIC     = round(c(AIC_3_real, AIC_em), 4),
  BIC     = round(c(BIC_3_real, BIC_em), 4)
), file.path(OUT_DIR, "13f_model_comparison.csv"), row.names = FALSE)

saveRDS(
  list(params = params, P = P_em, pi_stat = pi_stat,
       probs  = probs_df, LL = res_filt$LL,
       xi_filt = res_filt$xi_filt, xi_smooth = xi_smooth,
       state_names = STATE_NAMES),
  file.path(OUT_DIR, "13_em_results.rds")
)

cat("\nГрафики сохранены в", OUT_DIR, "\n")
cat("=== Этап 13 завершён ===\n")
