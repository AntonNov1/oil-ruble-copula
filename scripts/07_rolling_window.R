# 07_rolling_window.R — Оптимизация размера скользящего окна
# ЦЕЛЬ: найти оптимальный размер окна W* для своевременного обнаружения
#       структурного сдвига при минимальном числе ложных срабатываний.
#
# Метод: для каждого кандидата W из WINDOW_CANDIDATES:
#   1. Применяем rolling KS-тест с шагом STEP
#   2. Считаем: (a) число обнаруженных известных сдвигов,
#              (b) число ложных срабатываний (сигналов вне сдвигов),
#              (c) задержку обнаружения (лаг)
#   3. Строим Pareto-фронт: точность vs. задержка

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(ggplot2)
library(copula)
library(zoo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/07_rolling_window"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

WINDOW_CANDIDATES <- c(26, 39, 52, 65, 78, 104)   # кандидаты (недель)
STEP              <- 4                              # шаг (недель)
ALPHA             <- 0.95                           # уровень доверия KS-теста
TOLERANCE         <- 13                             # недель — допуск для "поймал сдвиг"
DECISION_LEVEL    <- 0.95

# Известные даты структурных сдвигов 
KNOWN_BREAKS <- as.Date(c(
  "2020-01-18",   # COVID обвал нефти
  "2020-06-20",   # Восстановление
  "2022-02-26",   # Начало СВО
  "2022-06-11",   # Адаптация к санкциям
  "2023-09-09"    # Новый режим
))

# ЗАГРУЗКА

returns_full <- load_log_returns(data_dir = PARENT)
n_obs <- nrow(returns_full)

cat(sprintf("Данные: N=%d  (%s — %s)\n", n_obs,
            min(returns_full$Date), max(returns_full$Date)))
cat(sprintf("Известных сдвигов: %d\n", length(KNOWN_BREAKS)))

# ROLLING KS-СТАТИСТИКА

# Вычисляет KS-статистику для одного окна (быстрая версия для 2D)
rolling_ks_one_window <- function(oil, fx) {
  n  <- length(oil)
  u1 <- ecdf(fx)(fx)
  u2 <- ecdf(oil)(oil)

  grid   <- seq(1/n, 1, length.out = min(n, 20))   # грубая сетка для скорости
  z1 <- rep(grid, each = length(grid))
  z2 <- rep(grid, times = length(grid))

  KS_vec <- numeric(n - 1)
  for (t in 1:(n - 1)) {
    u01 <- u1[1:t];       u02 <- u2[1:t]
    u11 <- u1[(t+1):n];   u12 <- u2[(t+1):n]
    EC0 <- numeric(length(z1))
    EC1 <- numeric(length(z1))
    for (l in seq_along(z1)) {
      EC0[l] <- sum(u01 <= z1[l] & u02 <= z2[l]) / t
      EC1[l] <- sum(u11 <= z1[l] & u12 <= z2[l]) / (n - t)
    }
    KS_vec[t] <- max(abs(EC0 - EC1) * sqrt(t * (n - t)) / n)
  }
  list(max_KS = max(KS_vec), break_idx = which.max(KS_vec))
}

# ОСНОВНОЙ ЦИКЛ

summary_results <- data.frame()

for (W in WINDOW_CANDIDATES) {
  cat(sprintf("\n=== Окно W = %d недель ===\n", W))

  signal_dates <- c()
  idx_seq <- seq(W, n_obs - 1, by = STEP)
  Cv <- critical.value.ks.test(W, ALPHA)

  for (t_end in idx_seq) {
    t_start  <- t_end - W + 1
    window_df <- returns_full[t_start:t_end, ]
    oil  <- window_df$log_OIL
    fx   <- window_df$log_FX

    res <- tryCatch(rolling_ks_one_window(oil, fx), error = function(e) NULL)
    if (is.null(res)) next

    if (res$max_KS >= Cv) {
      break_date <- window_df$Date[res$break_idx]
      signal_dates <- c(signal_dates, as.numeric(break_date))
    }
  }

  signal_dates <- as.Date(signal_dates)

  # Считаем TP, FP, FN
  tp     <- 0
  delays <- c()
  matched <- rep(FALSE, length(KNOWN_BREAKS))

  for (sig in signal_dates) {
    diffs <- abs(as.numeric(sig - KNOWN_BREAKS))
    min_d <- min(diffs)
    min_i <- which.min(diffs)
    if (min_d <= TOLERANCE * 7 && !matched[min_i]) {  # 7 дней/нед
      tp <- tp + 1
      matched[min_i] <- TRUE
      delays <- c(delays, min_d / 7)
    }
  }

  fp <- length(signal_dates) - tp
  fn <- sum(!matched)

  row <- data.frame(
    W              = W,
    n_signals      = length(signal_dates),
    TP             = tp,
    FP             = fp,
    FN             = fn,
    Recall         = round(tp / length(KNOWN_BREAKS), 3),
    Precision      = round(ifelse(length(signal_dates) > 0, tp / length(signal_dates), 0), 3),
    Avg_delay_wk   = round(ifelse(length(delays) > 0, mean(delays), NA), 1)
  )
  summary_results <- rbind(summary_results, row)
  cat(sprintf("  TP=%d  FP=%d  FN=%d  Recall=%.2f  Precision=%.2f  Задержка=%.1f нед\n",
              tp, fp, fn, row$Recall, row$Precision, ifelse(is.na(row$Avg_delay_wk), 0, row$Avg_delay_wk)))
}

write.csv(summary_results, file.path(OUT_DIR, "07a_window_comparison.csv"), row.names = FALSE)

cat("\n=== ИТОГОВАЯ ТАБЛИЦА ===\n")
print(summary_results)

# ОПТИМАЛЬНОЕ ОКНО: максимальный F1-score с минимальной задержкой

summary_results <- summary_results %>%
  mutate(
    F1    = round(ifelse((Precision + Recall) > 0,
                         2 * Precision * Recall / (Precision + Recall), 0), 3),
    Score = F1 - 0.01 * replace_na(Avg_delay_wk, 0)   # штраф за задержку
  )

best_W <- summary_results$W[which.max(summary_results$Score)]
cat(sprintf("\n★ Оптимальный размер окна: W* = %d недель\n", best_W))

# ГРАФИКИ

p_f1 <- ggplot(summary_results, aes(x = W)) +
  geom_line(aes(y = Recall,    color = "Recall"),    linewidth = 1) +
  geom_line(aes(y = Precision, color = "Precision"), linewidth = 1) +
  geom_line(aes(y = F1,        color = "F1"),        linewidth = 1.5, linetype = "dashed") +
  geom_vline(xintercept = best_W, color = "gold", linewidth = 1.2) +
  annotate("text", x = best_W + 2, y = 0.9, label = paste0("W*=", best_W), color = "goldenrod") +
  scale_color_manual(values = c("Recall" = "steelblue", "Precision" = "tomato", "F1" = "darkgreen")) +
  scale_x_continuous(breaks = WINDOW_CANDIDATES) +
  labs(title = "Качество детекции структурных сдвигов vs. размер окна",
       x = "Размер окна W (недели)", y = "Значение метрики", color = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "07a_f1_vs_window.png"), p_f1, width = 8, height = 5, dpi = 150)

p_delay <- ggplot(summary_results, aes(x = W, y = Avg_delay_wk)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = best_W, color = "gold", linewidth = 1.2) +
  scale_x_continuous(breaks = WINDOW_CANDIDATES) +
  labs(title = "Средняя задержка обнаружения vs. размер окна",
       x = "Размер окна W (недели)", y = "Задержка (недели)") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "07b_delay_vs_window.png"), p_delay, width = 8, height = 4, dpi = 150)

# Pareto: Recall vs. Задержка
p_pareto <- ggplot(summary_results, aes(x = Avg_delay_wk, y = Recall, label = paste0("W=", W))) +
  geom_point(aes(color = as.factor(W)), size = 5) +
  geom_text(nudge_y = 0.02, size = 3.5) +
  geom_point(data = filter(summary_results, W == best_W),
             size = 8, shape = 1, stroke = 2, color = "gold") +
  labs(title = "Pareto-фронт: Recall vs. Задержка",
       subtitle = "Золотой круг = оптимальное окно W*",
       x = "Средняя задержка (недели)", y = "Recall (доля найденных сдвигов)",
       color = "W (недели)") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "07c_pareto_front.png"), p_pareto, width = 7, height = 5, dpi = 150)

cat(sprintf("\nОптимальный размер окна: W* = %d недель (F1=%.3f)\n",
            best_W, max(summary_results$F1)))
