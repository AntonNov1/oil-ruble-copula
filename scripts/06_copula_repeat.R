# 06_copula_repeat.R — Повторяемость копул (обоснование модели со сменой режимов)
# ЦЕЛЬ: показать, что одни и те же типы копул повторяются в схожих
#       экономических условиях, что обосновывает марковскую модель с режимами.
#
# Анализ:
#   1. Скользящее окно: лучшая копула в каждый момент времени
#   2. Матрица переходов между типами копул (Markov chain)
#   3. Стационарные вероятности режимов
#   4. Кластеризация периодов по параметрам копулы

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(copula)
library(ggplot2)
library(tidyr)
library(zoo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/06_repeatability"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

WINDOW_SIZE  <- 52     # ширина скользящего окна (недели)
STEP_SIZE    <- 1      # шаг скользящего окна (недели)
COP_METHOD   <- "ml"
COP_TYPES    <- c("Frank", "Gumbel", "Clayton", "Gaussian", "Student")

# ЗАГРУЗКА ДАННЫХ

returns_full <- load_log_returns(data_dir = PARENT)
best_models  <- load_best_models(data_dir = PARENT)
n_obs <- nrow(returns_full)

cat(sprintf("Данные: %s — %s (N=%d)\n",
            min(returns_full$Date), max(returns_full$Date), n_obs))

# 1. СКОЛЬЗЯЩЕЕ ОКНО: ЛУЧШАЯ КОПУЛА

cat(sprintf("\n=== СКОЛЬЗЯЩЕЕ ОКНО (W=%d недель) ===\n", WINDOW_SIZE))

rolling_results <- data.frame()
dates_seq <- seq(WINDOW_SIZE, n_obs, by = STEP_SIZE)

for (t_end in dates_seq) {
  t_start <- t_end - WINDOW_SIZE + 1
  window_df <- returns_full[t_start:t_end, ]
  center_date <- window_df$Date[ceiling(WINDOW_SIZE / 2)]
  end_date    <- window_df$Date[WINDOW_SIZE]

  # Псевдонаблюдения
  u_mat <- tryCatch(
    pobs(cbind(window_df$log_OIL, window_df$log_FX)),
    error = function(e) NULL
  )
  if (is.null(u_mat)) next

  # Подбор лучшей копулы по BIC
  best_bic  <- Inf
  best_type <- NA_character_
  best_param <- NA_real_

  for (ct in COP_TYPES) {
    res <- fit_copula_bivariate(u_mat, ct, method = COP_METHOD)
    if (!res$ok || is.na(res$BIC)) next
    if (res$BIC < best_bic) {
      best_bic   <- res$BIC
      best_type  <- ct
      best_param <- res$params[1]
    }
  }

  # Если лучший BIC > 0 → независимость
  if (is.na(best_type) || best_bic > 0) {
    best_type  <- "independence"
    best_param <- 0
    best_bic   <- 0
  }

  rolling_results <- rbind(rolling_results, data.frame(
    center_date = center_date,
    end_date    = end_date,
    best_copula = best_type,
    param1      = best_param,
    BIC         = best_bic
  ))

  progress_bar(which(dates_seq == t_end), length(dates_seq), "  Rolling:")
}

write.csv(rolling_results, file.path(OUT_DIR, "06a_rolling_copula.csv"), row.names = FALSE)

# ВИЗУАЛИЗАЦИЯ СКОЛЬЗЯЩЕЙ КОПУЛЫ

# Цвета для типов копул
cop_colors <- c(
  "Frank"        = "#2196F3",
  "Gumbel"       = "#FF9800",
  "Clayton"      = "#4CAF50",
  "Gaussian"     = "#9C27B0",
  "Student"      = "#F44336",
  "independence" = "#9E9E9E"
)

# График типа копулы во времени
p_rolling <- ggplot(rolling_results, aes(x = center_date, y = best_copula, color = best_copula)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = cop_colors, drop = FALSE) +
  # Отметки идентифицированных периодов из best_models
  geom_vline(xintercept = as.numeric(best_models$Start), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  labs(title = sprintf("Лучшая копула в скользящем окне (W=%d недель)", WINDOW_SIZE),
       x = NULL, y = "Тип копулы", color = "Копула") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "06a_rolling_copula_type.png"), p_rolling, width = 12, height = 5, dpi = 150)

# Частота появления каждого типа
freq_table <- rolling_results %>%
  count(best_copula, name = "Count") %>%
  mutate(Pct = round(Count / sum(Count) * 100, 1)) %>%
  arrange(desc(Count))
write.csv(freq_table, file.path(OUT_DIR, "06a_copula_frequency.csv"), row.names = FALSE)

cat("\n=== Частота копул в скользящем окне ===\n")
print(freq_table)

p_freq <- ggplot(freq_table, aes(x = reorder(best_copula, -Count), y = Count, fill = best_copula)) +
  geom_col() +
  geom_text(aes(label = paste0(Pct, "%")), vjust = -0.3, size = 4) +
  scale_fill_manual(values = cop_colors, guide = "none") +
  labs(title = "Частота типов копул в скользящем окне",
       x = "Тип копулы", y = "Число окон") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "06b_copula_frequency.png"), p_freq, width = 7, height = 5, dpi = 150)

# МАТРИЦА ПЕРЕХОДОВ МЕЖДУ КОПУЛАМИ 

cat("\n=== МАТРИЦА ПЕРЕХОДОВ ===\n")

cop_seq   <- rolling_results$best_copula
all_types <- sort(unique(cop_seq))
n_types   <- length(all_types)

trans_mat <- matrix(0, nrow = n_types, ncol = n_types,
                    dimnames = list(From = all_types, To = all_types))

for (k in seq_len(length(cop_seq) - 1)) {
  from <- cop_seq[k]
  to   <- cop_seq[k + 1]
  trans_mat[from, to] <- trans_mat[from, to] + 1
}

# Нормируем по строкам → вероятности переходов
trans_prob <- sweep(trans_mat, 1, rowSums(trans_mat), FUN = "/")
trans_prob[is.nan(trans_prob)] <- 0

cat("Матрица переходных вероятностей:\n")
print(round(trans_prob, 3))

write.csv(as.data.frame(trans_prob), file.path(OUT_DIR, "06c_transition_matrix.csv"))

# Визуализация матрицы переходов
trans_df <- as.data.frame(trans_prob) %>%
  rownames_to_column("From") %>%
  pivot_longer(-From, names_to = "To", values_to = "Prob")

p_trans <- ggplot(trans_df, aes(x = To, y = From, fill = Prob)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Prob, 2)), size = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Вероятность") +
  labs(title = "Матрица переходных вероятностей между типами копул",
       x = "Следующий тип", y = "Текущий тип") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "06c_transition_matrix.png"), p_trans, width = 7, height = 6, dpi = 150)

# СТАЦИОНАРНЫЕ ВЕРОЯТНОСТИ 

cat("\n=== СТАЦИОНАРНЫЕ ВЕРОЯТНОСТИ ===\n")

steady_state_prob <- function(P, n_iter = 10000) {
  p <- rep(1 / nrow(P), nrow(P))
  for (i in seq_len(n_iter)) p <- p %*% P
  p[1, ]
}

if (n_types > 1 && sum(trans_prob) > 0) {
  ss <- steady_state_prob(trans_prob)
  ss_df <- data.frame(Copula = all_types, Stationary_Prob = round(ss, 4))
  cat("Стационарное распределение:\n")
  print(ss_df)
  write.csv(ss_df, file.path(OUT_DIR, "06d_stationary_probs.csv"), row.names = FALSE)
}

# СВЯЗЬ ПАРАМЕТРА КОПУЛЫ С ЦЕНОЙ НЕФТИ 

cat("\n=== СВЯЗЬ ПАРАМЕТРА И ВОЛАТИЛЬНОСТИ НЕФТИ ===\n")

# Добавляем скользящую волатильность нефти
oil_sd <- rollapply(returns_full$log_OIL, width = WINDOW_SIZE,
                    FUN = sd, fill = NA, align = "right")
returns_full$oil_vol <- oil_sd

# Объединяем с rolling_results по end_date
vol_df <- returns_full %>% select(Date, oil_vol) %>% filter(!is.na(oil_vol))

rolling_results <- rolling_results %>%
  left_join(vol_df, by = c("end_date" = "Date"))

p_vol_cop <- ggplot(
  filter(rolling_results, !is.na(param1) & !is.na(oil_vol) & best_copula != "independence"),
  aes(x = oil_vol, y = param1, color = best_copula)
) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  facet_wrap(~best_copula, scales = "free_y") +
  scale_color_manual(values = cop_colors, guide = "none") +
  labs(title = "Параметр копулы vs. волатильность нефти",
       x = "Скользящее SD(Oil returns)", y = "Параметр θ копулы") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "06e_param_vs_vol.png"), p_vol_cop, width = 10, height = 6, dpi = 150)
