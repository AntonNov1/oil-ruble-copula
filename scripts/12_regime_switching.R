# 12_regime_switching.R — Копула-базированная модель переключения режимов
#
# Модель: скрытая марковская цепь (HMM) с тремя режимами зависимости
#   S_t ∈ {1=Frank, 2=Student, 3=Gumbel°}
#
# Ядро: фильтр Гамильтона (1989) + сглаживатель Кима (1994)
#
# Структура:
#   1. Матрица переходов P (3×3), параметризованная через softmax
#   2. Плотности копул f_s(u_t) при фиксированных θ из скрипта 09
#   3. Фильтрованные вероятности: P(S_t = s | u_1,...,u_t)
#   4. Сглаженные вероятности:   P(S_t = s | u_1,...,u_T)
#   5. Декодирование режимов (Viterbi)
#
# Параметры:
#   Вариант A: только P (6 свободных параметров) при фиксированных θ
#   Вариант B: P + θ совместно (8 параметров)

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(tidyr)
library(copula)
library(ggplot2)

# задайте рабочую директорию проекта перед запуском

PARENT  <- ".."
OUT_DIR <- "results/12_regime_switching"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ДАННЫЕ

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

# Псевдонаблюдения на полной выборке
u_all <- pobs(cbind(lr$log_OIL, lr$log_FX))
colnames(u_all) <- c("u_oil","u_fx")
u_all <- as.data.frame(u_all)
u_all$Date <- lr$Date

T_obs <- nrow(u_all)
cat(sprintf("Полная выборка: N=%d (%s — %s)\n",
            T_obs, min(lr$Date), max(lr$Date)))

# Известные даты структурных сдвигов (из скрипта 09)
BREAK_DATES <- as.Date(c(
  "2020-01-18",   # P1→P2
  "2020-06-20",   # P2→P3
  "2022-02-26",   # P3→P4
  "2022-06-11",   # P4→P5
  "2023-09-09"    # P5→P6
))

# 2. ПАРАМЕТРЫ КОПУЛ (фиксированные из скрипта 09)

THETA_FRANK   <- 2.5066
RHO_STUDENT   <- -0.6347
DF_STUDENT    <- 0.5551
THETA_GUMBEL  <- 1.0869

cop_frank   <- frankCopula(param = THETA_FRANK,  dim = 2)
cop_student <- tCopula(param = RHO_STUDENT, df = DF_STUDENT, dim = 2)
cop_gumbel  <- gumbelCopula(param = THETA_GUMBEL, dim = 2)

# ВЫЧИСЛЕНИЕ ПЛОТНОСТЕЙ КОПУЛ

cat("Вычисляю плотности копул...\n")

u_mat   <- as.matrix(u_all[, c("u_oil","u_fx")])
u_rot   <- cbind(1 - u_mat[,1], u_mat[,2])   # 90°-поворот для Gumbel°

f_frank   <- tryCatch(pmax(dCopula(u_mat, cop_frank),   .Machine$double.eps),
                      error = function(e) rep(1e-10, T_obs))
f_student <- tryCatch(pmax(dCopula(u_mat, cop_student), .Machine$double.eps),
                      error = function(e) rep(1e-10, T_obs))
f_gumbel  <- tryCatch(pmax(dCopula(u_rot, cop_gumbel),  .Machine$double.eps),
                      error = function(e) rep(1e-10, T_obs))

f_mat <- cbind(f_frank, f_student, f_gumbel)  # T × 3
cat(sprintf("Плотности: min=%.6f  max=%.4f\n", min(f_mat), max(f_mat)))

# 4. ФИЛЬТР ГАМИЛЬТОНА (3 состояния)

hamilton_filter <- function(f_mat, P, pi0 = NULL) {
  K <- ncol(f_mat)
  T <- nrow(f_mat)

  if (is.null(pi0)) {
    # Стационарное распределение
    ev  <- eigen(t(P))$vectors[, 1]
    pi0 <- Re(ev) / sum(Re(ev))
    pi0 <- pmax(pi0, 0);  pi0 <- pi0 / sum(pi0)
  }

  xi_filt <- matrix(0, T, K)   # фильтрованные вероятности
  xi_pred <- matrix(0, T, K)   # предсказанные
  LL_vec  <- numeric(T)

  p_pred <- pi0

  for (t in 1:T) {
    xi_pred[t,] <- p_pred
    eta          <- f_mat[t,]
    joint        <- eta * p_pred
    denom        <- sum(joint)
    if (denom <= 0 || is.nan(denom)) {
      xi_filt[t,] <- p_pred
      LL_vec[t]   <- log(.Machine$double.eps)
    } else {
      xi_filt[t,] <- joint / denom
      LL_vec[t]   <- log(denom)
    }
    p_pred <- as.vector(t(P) %*% xi_filt[t,])
    p_pred <- pmax(p_pred, 0);  p_pred <- p_pred / sum(p_pred)
  }

  list(xi_filt = xi_filt, xi_pred = xi_pred, LL = sum(LL_vec), LL_vec = LL_vec)
}

# 5. СГЛАЖИВАТЕЛЬ КИМА (3 состояния)

kim_smoother <- function(xi_filt, xi_pred, P) {
  K <- ncol(xi_filt)
  T <- nrow(xi_filt)

  xi_smooth <- matrix(0, T, K)
  xi_smooth[T,] <- xi_filt[T,]

  for (t in (T - 1):1) {
    for (j in 1:K) {
      num <- P[j,] * xi_smooth[t+1,]
      den <- pmax(xi_pred[t+1,], .Machine$double.eps)
      xi_smooth[t,j] <- xi_filt[t,j] * sum(num / den)
    }
    s <- sum(xi_smooth[t,])
    if (s > 0) xi_smooth[t,] <- xi_smooth[t,] / s
  }

  xi_smooth
}

# 6. ОЦЕНКА МАТРИЦЫ ПЕРЕХОДОВ (Вариант A: θ фиксированы)

cat("\n=== ОЦЕНКА МАТРИЦЫ ПЕРЕХОДОВ ===\n")

# Параметризация: softmax для каждой строки P
# arr = [a11_logit, a12_logit, a21_logit, a22_logit, a31_logit, a32_logit]
# p_ii = 1 / (1 + exp(a_i1) + exp(a_i2)), p_ij = exp(a_ij) / (...)

softmax3 <- function(x) {
  e <- exp(x - max(x))
  e / sum(e)
}

arr_to_P <- function(arr) {
  # arr длиной 6: два свободных логита на строку
  P <- matrix(0, 3, 3)
  P[1,] <- softmax3(c(arr[1], arr[2], 0))
  P[2,] <- softmax3(c(arr[3], arr[4], 0))
  P[3,] <- softmax3(c(arr[5], arr[6], 0))
  P
}

neg_loglik_P <- function(arr) {
  P <- tryCatch(arr_to_P(arr), error = function(e) NULL)
  if (is.null(P) || any(P < 0) || any(P > 1)) return(1e9)
  res <- hamilton_filter(f_mat, P)
  -res$LL
}

# Начальное приближение: диагональная P (высокая персистентность)
P_init <- matrix(c(0.90, 0.05, 0.05,
                   0.05, 0.90, 0.05,
                   0.05, 0.05, 0.90), 3, 3, byrow = TRUE)

init_arr <- function(P) {
  row_to_logit <- function(p) {
    ref <- p[3]
    c(log(p[1]/ref), log(p[2]/ref))
  }
  c(row_to_logit(P[1,]),
    row_to_logit(P[2,]),
    row_to_logit(P[3,]))
}

x0 <- init_arr(P_init)

cat("Оптимизация матрицы переходов (Nelder-Mead)...\n")
opt_res <- optim(
  x0, neg_loglik_P,
  method  = "Nelder-Mead",
  control = list(maxit = 10000, reltol = 1e-8)
)

cat(sprintf("Статус: %s  |  LogLik = %.4f  |  AIC = %.4f\n",
            opt_res$message,
            -opt_res$value,
            2 * opt_res$value + 2 * 6))

P_hat <- arr_to_P(opt_res$par)
rownames(P_hat) <- colnames(P_hat) <- c("Frank","Student","Gumbel°")
cat("\nОценённая матрица переходов P:\n")
print(round(P_hat, 4))

# Стационарное распределение
ev_P   <- eigen(t(P_hat))$vectors[, 1]
pi_stat <- Re(ev_P) / sum(Re(ev_P))
pi_stat <- pmax(pi_stat, 0);  pi_stat <- pi_stat / sum(pi_stat)
cat(sprintf("\nСтационарные вероятности: Frank=%.3f  Student=%.3f  Gumbel°=%.3f\n",
            pi_stat[1], pi_stat[2], pi_stat[3]))

# Средняя продолжительность режимов
dur <- 1 / (1 - diag(P_hat))
cat(sprintf("Средняя дл. режимов (недели): Frank=%.1f  Student=%.1f  Gumbel°=%.1f\n",
            dur[1], dur[2], dur[3]))

# ФИЛЬТРОВАННЫЕ И СГЛАЖЕННЫЕ ВЕРОЯТНОСТИ

cat("\nФильтр Гамильтона + сглаживатель Кима...\n")

res_filter <- hamilton_filter(f_mat, P_hat)
xi_smooth  <- kim_smoother(res_filter$xi_filt, res_filter$xi_pred, P_hat)

cat(sprintf("LogLik (финал) = %.4f  AIC = %.4f  BIC = %.4f\n",
            res_filter$LL,
            2 * 6 - 2 * res_filter$LL,
            6 * log(T_obs) - 2 * res_filter$LL))

# Сохраняем временны́е ряды вероятностей
probs_df <- data.frame(
  Date          = u_all$Date,
  Frank_filt    = res_filter$xi_filt[,1],
  Student_filt  = res_filter$xi_filt[,2],
  Gumbel_filt   = res_filter$xi_filt[,3],
  Frank_smooth  = xi_smooth[,1],
  Student_smooth= xi_smooth[,2],
  Gumbel_smooth = xi_smooth[,3],
  Regime_filt   = c("Frank","Student","Gumbel°")[apply(res_filter$xi_filt, 1, which.max)],
  Regime_smooth = c("Frank","Student","Gumbel°")[apply(xi_smooth, 1, which.max)]
)
write.csv(probs_df, file.path(OUT_DIR, "12a_regime_probs.csv"), row.names = FALSE)

# 8. ДЕКОДИРОВАНИЕ VITERBI (наиболее вероятная последовательность режимов)

viterbi <- function(f_mat, P, pi0) {
  K <- ncol(f_mat);  T <- nrow(f_mat)
  log_P   <- log(pmax(P, .Machine$double.eps))
  log_f   <- log(pmax(f_mat, .Machine$double.eps))
  log_pi0 <- log(pmax(pi0, .Machine$double.eps))

  delta <- matrix(-Inf, T, K)
  psi   <- matrix(0L,   T, K)
  delta[1,] <- log_pi0 + log_f[1,]

  for (t in 2:T) {
    for (j in 1:K) {
      scores    <- delta[t-1,] + log_P[,j]
      psi[t,j]  <- which.max(scores)
      delta[t,j]<- max(scores) + log_f[t,j]
    }
  }

  path <- integer(T)
  path[T] <- which.max(delta[T,])
  for (t in (T-1):1) path[t] <- psi[t+1, path[t+1]]
  path
}

vit_path <- viterbi(f_mat, P_hat, pi_stat)
probs_df$Regime_viterbi <- c("Frank","Student","Gumbel°")[vit_path]

cat("\n=== ДЕКОДИРОВАННЫЕ РЕЖИМЫ (Viterbi — первые/последние) ===\n")
print(table(probs_df$Regime_viterbi))

# ВИЗУАЛИЗАЦИИ

cat("\n=== ВИЗУАЛИЗАЦИИ ===\n")

regime_colors <- c("Frank"="#1565C0","Student"="#B71C1C","Gumbel°"="#FF8F00")

# --- 9.1 Сглаженные вероятности ---
prob_long <- probs_df %>%
  select(Date, Frank_smooth, Student_smooth, Gumbel_smooth) %>%
  pivot_longer(-Date, names_to = "Regime", values_to = "Prob") %>%
  mutate(Regime = gsub("_smooth","", Regime),
         Regime = gsub("Gumbel","Gumbel°", Regime))

p_smooth <- ggplot(prob_long, aes(x = Date, y = Prob, color = Regime, fill = Regime)) +
  geom_area(alpha = 0.3, position = "identity") +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = BREAK_DATES, linetype = "dashed",
             color = "black", linewidth = 0.5, alpha = 0.6) +
  scale_color_manual(values = regime_colors) +
  scale_fill_manual(values  = regime_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
  labs(
    title    = "Апостериорные вероятности режимов (сглаженные)",
    subtitle = "Фильтр Гамильтона + сглаживатель Кима | Вертикальные линии = структурные сдвиги",
    x = NULL, y = "P(режим | данные)", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "12a_smoothed_probs.png"), p_smooth,
       width = 13, height = 5, dpi = 150)

# --- 9.2 Стековая диаграмма (stacked area) ---
p_stack <- ggplot(prob_long, aes(x = Date, y = Prob, fill = Regime)) +
  geom_area(position = "stack", alpha = 0.85) +
  geom_vline(xintercept = BREAK_DATES, linetype = "dashed",
             color = "white", linewidth = 0.7) +
  scale_fill_manual(values = regime_colors) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Структура режимов во времени (сглаженные вероятности)",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "12b_stacked_probs.png"), p_stack,
       width = 13, height = 4.5, dpi = 150)

# --- 9.3 Сравнение: Viterbi vs сглаженные vs детерминированные границы ---
# Детерминированное назначение на основе структурных сдвигов из скрипта 09
period_regime <- data.frame(
  Start  = c(as.Date("2017-08-25"), BREAK_DATES),
  End    = c(BREAK_DATES - 1, as.Date("2025-12-31")),
  Regime = c("Frank","Frank","Frank","Student","Gumbel°","Gumbel°")
)

probs_df$Regime_determ <- sapply(probs_df$Date, function(d) {
  idx <- which(period_regime$Start <= d & period_regime$End >= d)
  if (length(idx) == 0) NA else period_regime$Regime[idx[1]]
})

# Согласованность Viterbi с детерминированными границами
agree <- mean(probs_df$Regime_viterbi == probs_df$Regime_determ, na.rm = TRUE)
cat(sprintf("  Согласованность Viterbi с детерм. границами: %.1f%%\n", agree*100))

# --- 9.4 Фильтрованные vs сглаженные вероятности для Frank ---
comp_df <- probs_df %>%
  select(Date, Frank_filt, Frank_smooth) %>%
  pivot_longer(-Date, names_to = "Type", values_to = "Prob") %>%
  mutate(Type = ifelse(Type == "Frank_filt",
                       "Фильтрованные", "Сглаженные"))

p_comp <- ggplot(comp_df, aes(x = Date, y = Prob, color = Type)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = BREAK_DATES, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c("Фильтрованные"="steelblue","Сглаженные"="darkblue")) +
  scale_x_date(date_labels = "%Y") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title    = "Frank-режим: фильтрованные vs сглаженные вероятности",
       subtitle = "Сглаживание устраняет запаздывание фильтра",
       x = NULL, y = "P(Frank | данные)", color = NULL) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "12c_frank_filt_vs_smooth.png"),
       p_comp, width = 12, height = 4, dpi = 150)

# СВОДНАЯ СТАТИСТИКА

cat("\n", strrep("=", 60), "\n")
cat("СВОДКА МОДЕЛИ ПЕРЕКЛЮЧЕНИЯ РЕЖИМОВ\n")
cat(strrep("=", 60), "\n\n")

cat("Матрица переходов P:\n")
print(round(P_hat, 4))

cat(sprintf("\nLogLik = %.4f  AIC = %.4f  BIC = %.4f\n",
            res_filter$LL,
            2*6 - 2*res_filter$LL,
            6*log(T_obs) - 2*res_filter$LL))

cat("\nСтационарные вероятности (теоретические):\n")
cat(sprintf("  Frank   = %.3f (%.1f%%)  ср. продолж. %.1f нед.\n",
            pi_stat[1], pi_stat[1]*100, dur[1]))
cat(sprintf("  Student = %.3f (%.1f%%)  ср. продолж. %.1f нед.\n",
            pi_stat[2], pi_stat[2]*100, dur[2]))
cat(sprintf("  Gumbel° = %.3f (%.1f%%)  ср. продолж. %.1f нед.\n",
            pi_stat[3], pi_stat[3]*100, dur[3]))

cat("\nФактическое распределение наблюдений (сглаж. вероятности):\n")
frac_smooth <- colMeans(xi_smooth)
cat(sprintf("  Frank   = %.3f (%.1f%%)\n", frac_smooth[1], frac_smooth[1]*100))
cat(sprintf("  Student = %.3f (%.1f%%)\n", frac_smooth[2], frac_smooth[2]*100))
cat(sprintf("  Gumbel° = %.3f (%.1f%%)\n", frac_smooth[3], frac_smooth[3]*100))

cat(sprintf("\nСогласованность Viterbi с детерм. границами: %.1f%%\n", agree*100))

# Сохраняем P_hat
write.csv(as.data.frame(round(P_hat, 6)),
          file.path(OUT_DIR, "12b_transition_matrix.csv"))

cat("\nГрафики сохранены в", OUT_DIR, "\n")
cat("=== Этап 12 завершён ===\n")
