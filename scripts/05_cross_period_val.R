# 05_cross_period_val.R — Кросс-периодная валидация структурных сдвигов
# ЦЕЛЬ: доказать наличие структурных сдвигов двумя способами:
#
# Метод A (BIC-таблица): для каждой пары (копула_i, период_j)
#   - Re-fit копулы типа i на данных периода j
#   - Если лучшая копула для периода j != лучшей для периода i -> сдвиг
#   - Строим тепловую карту BIC[i,j]
#
# Метод B (фиксированные параметры): для каждой пары (период_i, период_j)
#   - Берём параметры θ_i (оцененные на данных периода i)
#   - Вычисляем LL(θ_i | данные_j) — out-of-sample log-правдоподобие
#   - Диагональ LL[i,i] >> недиагональных LL[i,j] -> структурные сдвиги
#
# Выход: матрицы BIC и LL, тепловые карты;

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(copula)
library(ggplot2)
library(tidyr)
library(tibble)
library(zoo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/05_cross_period"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ЗАГРУЗКА ДАННЫХ

returns_full <- load_log_returns(data_dir = PARENT)
best_models  <- load_best_models(data_dir = PARENT)

# Периоды (используем те же, что в best_models)
periods <- best_models %>%
  select(Period_ID, Start, End, Copula, Params, Lag_USD, Lag_Oil)
n_periods <- nrow(periods)

cat(sprintf("Загружено %d периодов для кросс-валидации.\n", n_periods))
print(periods[, c("Period_ID", "Start", "End", "Copula")])

# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

# Получить псевдонаблюдения lag=0 для периода [start, end]
get_pobs_bivar <- function(start_date, end_date) {
  df <- preprocess(returns_full, start_date, end_date)
  if (nrow(df) < 5) return(NULL)
  mat <- pobs(cbind(df$log_OIL, df$log_FX))
  colnames(mat) <- c("u_oil", "u_fx")
  mat
}

# Подгонка бивариатной копулы с возвратом объекта + LL + BIC
fit_biv <- function(u_mat, cop_type) {
  n <- nrow(u_mat)
  tryCatch({
    cop_obj <- make_copula(cop_type, dim = 2)
    fit     <- fitCopula(cop_obj, data = u_mat, method = "ml")
    k       <- length(fit@copula@parameters)
    L       <- fit@loglik
    list(fit = fit, copula = fit@copula, loglik = L,
         AIC = count_AIC(L, k), BIC = count_BIC(L, k, n),
         params = fit@copula@parameters, k = k, ok = TRUE)
  }, error = function(e) {
    list(fit = NULL, copula = NULL, loglik = NA, AIC = NA, BIC = NA, params = NA, k = NA, ok = FALSE)
  })
}

# Оценить LL фиксированной копулы на новых данных
fixed_loglik <- function(copula_obj, u_new) {
  if (is.null(copula_obj)) return(NA_real_)
  tryCatch({
    dens <- dCopula(as.matrix(u_new), copula_obj)
    dens <- pmax(dens, .Machine$double.eps)
    sum(log(dens))
  }, error = function(e) NA_real_)
}

# МЕТОД A: BIC-матрица 

cat("\n=== МЕТОД A: BIC-МАТРИЦА (re-fit) ===\n")

COP_TYPES <- c("Frank", "Gumbel", "Clayton", "Gaussian", "Student")

# Загружаем псевдонаблюдения для всех периодов
pobs_list <- lapply(seq_len(n_periods), function(i) {
  u <- get_pobs_bivar(periods$Start[i], periods$End[i])
  cat(sprintf("  Период %d: N=%s\n", i, ifelse(is.null(u), "NULL", nrow(u))))
  u
})

# Матрица BIC[cop_type, period]
bic_matrix <- matrix(NA, nrow = length(COP_TYPES), ncol = n_periods,
                     dimnames = list(COP_TYPES, paste0("P", periods$Period_ID)))

for (ci in seq_along(COP_TYPES)) {
  for (pi in seq_len(n_periods)) {
    u <- pobs_list[[pi]]
    if (is.null(u)) next
    res <- fit_biv(u, COP_TYPES[ci])
    bic_matrix[ci, pi] <- if (res$ok) res$BIC else NA
    cat(sprintf("  [A] %s × Period %d: BIC = %.3f\n",
                COP_TYPES[ci], pi, ifelse(is.na(bic_matrix[ci, pi]), NA, bic_matrix[ci, pi])))
  }
}

write.csv(as.data.frame(bic_matrix), file.path(OUT_DIR, "05a_bic_matrix.csv"))

# Строим тепловую карту BIC[копула, период]
bic_df <- as.data.frame(bic_matrix) %>%
  rownames_to_column("Copula") %>%
  pivot_longer(-Copula, names_to = "Period", values_to = "BIC")

# Отмечаем минимум по столбцу (лучшая копула для каждого периода)
best_per_period <- bic_df %>%
  group_by(Period) %>%
  filter(!is.na(BIC)) %>%
  slice_min(order_by = BIC, n = 1, with_ties = FALSE) %>%
  mutate(is_best = TRUE)

bic_df <- bic_df %>%
  left_join(best_per_period %>% select(Copula, Period, is_best), by = c("Copula", "Period")) %>%
  mutate(is_best = replace_na(is_best, FALSE))

# Добавляем период_labels
period_labels <- paste0("P", periods$Period_ID, "\n(", format(periods$Start, "%Y-%m"), ")")
names(period_labels) <- paste0("P", periods$Period_ID)
bic_df$Period_label <- period_labels[bic_df$Period]

p_bic <- ggplot(bic_df, aes(x = Period_label, y = Copula, fill = BIC)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(data = filter(bic_df, !is.na(BIC)),
            aes(label = round(BIC, 2), fontface = ifelse(is_best, "bold", "plain")),
            size = 3.2, color = "black") +
  geom_tile(data = filter(bic_df, is_best), fill = NA, color = "gold", linewidth = 1.8) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                        midpoint = 0, name = "BIC", na.value = "grey90") +
  labs(title = "Метод A: Таблица BIC (re-fit каждой копулы на каждом периоде)",
       subtitle = "Золотая рамка = лучшая копула для периода; синий = хорошая подгонка",
       x = "Период", y = "Тип копулы") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 8), panel.grid = element_blank())

ggsave(file.path(OUT_DIR, "05a_bic_heatmap.png"), p_bic, width = 10, height = 5, dpi = 150)
cat("График 05a сохранён.\n")

# МЕТОД B: Log-likelihood с ФИКСИРОВАННЫМИ параметрами (из периода i)

cat("\n=== МЕТОД B: LL-матрица (фиксированные параметры из периода i) ===\n")

# Для каждого периода i подгоняем его лучшую копулу (bivariate lag=0)
# и запоминаем объект копулы с фиксированными параметрами
best_cop_fitted <- lapply(seq_len(n_periods), function(i) {
  u <- pobs_list[[i]]
  if (is.null(u)) return(NULL)
  cop_type <- as.character(periods$Copula[i])
  if (tolower(cop_type) == "independence") return(list(copula = indepCopula(dim=2), type="independence"))
  res <- fit_biv(u, cop_type)
  if (res$ok) res else NULL
})

# Матрица LL[source_period, target_period]
# LL[i,j] = sum(log dCopula(pobs_j; θ_i))  — параметры из периода i, данные из периода j
ll_matrix <- matrix(NA, nrow = n_periods, ncol = n_periods,
                    dimnames = list(paste0("Params_P", periods$Period_ID),
                                    paste0("Data_P",   periods$Period_ID)))

for (i in seq_len(n_periods)) {
  if (is.null(best_cop_fitted[[i]])) next
  cop_i <- best_cop_fitted[[i]]$copula

  for (j in seq_len(n_periods)) {
    u_j <- pobs_list[[j]]
    if (is.null(u_j)) next
    ll_matrix[i, j] <- fixed_loglik(cop_i, u_j)
    cat(sprintf("  [B] θ_P%d → Данные_P%d: LL = %.3f\n", i, j,
                ifelse(is.na(ll_matrix[i, j]), NA, ll_matrix[i, j])))
  }
}

write.csv(as.data.frame(ll_matrix), file.path(OUT_DIR, "05b_fixed_loglik_matrix.csv"))

# Нормализуем: LL[i,j] / n_j  (на единицу наблюдения)
ll_norm <- ll_matrix
for (j in seq_len(n_periods)) {
  u_j <- pobs_list[[j]]
  if (!is.null(u_j)) ll_norm[, j] <- ll_matrix[, j] / nrow(u_j)
}
write.csv(as.data.frame(ll_norm), file.path(OUT_DIR, "05b_normalized_loglik_matrix.csv"))

# Тепловая карта нормализованного LL
ll_df <- as.data.frame(ll_norm) %>%
  rownames_to_column("Source") %>%
  pivot_longer(-Source, names_to = "Target", values_to = "LL_per_obs")

diag_df <- data.frame(
  Source = paste0("Params_P", periods$Period_ID),
  Target = paste0("Data_P",   periods$Period_ID),
  is_diag = TRUE
)
ll_df <- ll_df %>%
  left_join(diag_df, by = c("Source", "Target")) %>%
  mutate(is_diag = replace_na(is_diag, FALSE))

# Красивые метки
src_labels <- setNames(paste0("P", periods$Period_ID, "\n(", as.character(periods$Copula), ")"),
                        paste0("Params_P", periods$Period_ID))
tgt_labels <- setNames(paste0("P", periods$Period_ID, "\n", format(periods$Start, "%Y-%m")),
                        paste0("Data_P",   periods$Period_ID))
ll_df$Source_label <- src_labels[ll_df$Source]
ll_df$Target_label <- tgt_labels[ll_df$Target]

p_ll <- ggplot(ll_df, aes(x = Target_label, y = Source_label, fill = LL_per_obs)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(data = filter(ll_df, !is.na(LL_per_obs)),
            aes(label = round(LL_per_obs, 3), fontface = ifelse(is_diag, "bold", "plain")),
            size = 3, color = "black") +
  geom_tile(data = filter(ll_df, is_diag), fill = NA, color = "gold", linewidth = 2) +
  scale_fill_gradient2(low = "tomato", mid = "lightyellow", high = "steelblue",
                        midpoint = 0, name = "LL/n", na.value = "grey90") +
  labs(title = "Метод B: LL/n с фиксированными параметрами из периода i, данные из периода j",
       subtitle = "Диагональ (золотая рамка) = in-sample. Чем ниже off-diagonal → тем сильнее структурный сдвиг.",
       x = "Данные (период j)", y = "Параметры θ_i (из периода i)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(file.path(OUT_DIR, "05b_loglik_heatmap.png"), p_ll, width = 10, height = 7, dpi = 150)
cat("График 05b сохранён.\n")

# ИТОГОВАЯ ТАБЛИЦА: Деградация LL при неверном периоде

degradation <- data.frame()
for (i in seq_len(n_periods)) {
  ll_in  <- ll_norm[i, i]    # in-sample (диагональ)
  ll_out <- ll_norm[i, -i]   # out-of-sample (остальные периоды)
  avg_out <- mean(ll_out, na.rm = TRUE)
  degradation <- rbind(degradation, data.frame(
    Period    = i,
    Copula    = as.character(periods$Copula[i]),
    LL_in     = round(ll_in, 4),
    LL_out_avg = round(avg_out, 4),
    Degr_pct  = round((avg_out - ll_in) / abs(ll_in) * 100, 1)
  ))
}
write.csv(degradation, file.path(OUT_DIR, "05c_degradation_table.csv"), row.names = FALSE)

cat("\n=== ДЕГРАДАЦИЯ LL ПРИ НЕВЕРНОМ ПЕРИОДЕ ===\n")
cat("(Отрицательный % = out-of-sample хуже in-sample -> структурный сдвиг подтверждён)\n")
print(degradation)