# 08_forecasting.R — Прогнозы с доверительными интервалами
#                    
# ЦЕЛЬ:
#   1. Backtesting: для каждой недели в тесте — прогнозируем USD/RUB
#      через условную копулу C(U_fx | U_oil = v), используя режим текущего периода
#   2. Сценарный анализ: при заданном изменении нефти (−20%, 0%, +20%)
#      строим условное распределение курса рубля
#   3. VaR / CVaR из копульной модели

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
OUT_DIR <- "results/08_forecasting"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

N_SAMPLES   <- 10000    # сэмплов для симуляции условного распределения
OIL_SCENARIOS <- c(-0.20, -0.10, 0.00, +0.10, +0.20)   # сценарии лог-доходности нефти
VAR_LEVEL   <- 0.95     # уровень VaR
BACKTEST_START_FRAC <- 0.3   # начало бэктеста от конца выборки

# ЗАГРУЗКА

returns_full <- load_log_returns(data_dir = PARENT)
weekly_prices <- load_weekly_prices(data_dir = PARENT)
best_models   <- load_best_models(data_dir = PARENT)

n_obs <- nrow(returns_full)
cat(sprintf("Данные: N=%d  (%s — %s)\n", n_obs,
            min(returns_full$Date), max(returns_full$Date)))

# ФУНКЦИИ ПРОГНОЗА

# Определяет режим (строку best_models) для заданной даты
get_regime <- function(date) {
  row <- best_models %>%
    filter(Start <= date, End >= date) %>%
    head(1)
  if (nrow(row) == 0) {
    row <- best_models %>% filter(End < date) %>% tail(1)
  }
  row
}

# Преобразует лог-доходность в u = ecdf(history)(x), зажатый в (0,1)
to_u <- function(x, history) {
  u <- ecdf(history)(x)
  pmin(pmax(u, 1e-5), 1 - 1e-5)
}

# Условный квантиль FX | OIL через cCopula (обращение)
conditional_quantile_fx <- function(cop_obj, u_oil, prob, tol = 1e-6) {
  u_oil <- pmin(pmax(u_oil, 1e-5), 1 - 1e-5)
  prob  <- pmin(pmax(prob,  1e-5), 1 - 1e-5)

  f_inv <- function(u_fx) {
    u_fx <- pmin(pmax(u_fx, 1e-5), 1 - 1e-5)
    tryCatch(
      cCopula(cbind(u_oil, u_fx), cop_obj, indices = 2, inverse = FALSE) - prob,
      error = function(e) NA_real_
    )
  }

  tryCatch({
    lo <- 1e-5; hi <- 1 - 1e-5
    flo <- f_inv(lo); fhi <- f_inv(hi)
    if (!is.finite(flo) || !is.finite(fhi)) return(NA_real_)
    if (sign(flo) == sign(fhi)) return(pmin(pmax(prob, 1e-5), 1 - 1e-5))
    uniroot(f_inv, interval = c(lo, hi), tol = tol)$root
  }, error = function(e) NA_real_)
}

# Полный условный прогноз: для заданного u_oil и режима
# Возвращает квантили q05, q50, q95 лог-доходности FX
conditional_forecast <- function(cop_obj, cop_type, u_oil, history_fx) {
  if (is.null(cop_obj) || tolower(cop_type) == "independence") {
    return(list(q05 = quantile(history_fx, 0.05),
                q50 = quantile(history_fx, 0.50),
                q95 = quantile(history_fx, 0.95)))
  }

  # Обращение через uniroot для трёх квантилей
  u_05 <- conditional_quantile_fx(cop_obj, u_oil, 0.05)
  u_50 <- conditional_quantile_fx(cop_obj, u_oil, 0.50)
  u_95 <- conditional_quantile_fx(cop_obj, u_oil, 0.95)

  list(
    q05 = ifelse(is.na(u_05), quantile(history_fx, 0.05), quantile(history_fx, u_05, type = 8)),
    q50 = ifelse(is.na(u_50), quantile(history_fx, 0.50), quantile(history_fx, u_50, type = 8)),
    q95 = ifelse(is.na(u_95), quantile(history_fx, 0.95), quantile(history_fx, u_95, type = 8))
  )
}

# 1. БЭКТЕСТ: прогноз на каждую неделю тестового периода

cat("\n=== БЭКТЕСТ ===\n")

backtest_start_idx <- max(30, n_obs - round(n_obs * BACKTEST_START_FRAC))
cat(sprintf("Бэктест с наблюдения %d (%s) по %d (%s)\n",
            backtest_start_idx, returns_full$Date[backtest_start_idx],
            n_obs - 1, returns_full$Date[n_obs - 1]))

backtest_results <- data.frame()

for (i in backtest_start_idx:(n_obs - 1)) {
  curr_date <- returns_full$Date[i]
  next_date <- returns_full$Date[i + 1]

  regime <- get_regime(curr_date)
  if (nrow(regime) == 0) next

  cop_type <- as.character(regime$Copula)

  # История внутри режима (с начала периода до curr_date)
  history <- returns_full %>%
    filter(Date >= regime$Start, Date <= curr_date)
  if (nrow(history) < 15) next

  # Псевдонаблюдения истории (для ECDF)
  hist_oil <- history$log_OIL
  hist_fx  <- history$log_FX

  # Наблюдаемое изменение нефти на следующей неделе
  next_oil_ret <- returns_full$log_OIL[i + 1]
  u_oil <- to_u(next_oil_ret, hist_oil)

  # Строим копулу по истории
  cop_obj <- tryCatch({
    if (tolower(cop_type) == "independence") return(NULL)
    u_mat <- pobs(cbind(hist_oil, hist_fx))
    cop   <- make_copula(cop_type, dim = 2)
    fit   <- fitCopula(cop, data = u_mat, method = "ml")
    fit@copula
  }, error = function(e) NULL)

  # Прогноз
  preds <- conditional_forecast(cop_obj, cop_type, u_oil, hist_fx)

  # Истинная лог-доходность FX
  true_fx <- returns_full$log_FX[i + 1]

  backtest_results <- rbind(backtest_results, data.frame(
    Date       = next_date,
    Period_ID  = regime$Period_ID,
    Copula     = cop_type,
    Oil_next   = next_oil_ret,
    FX_true    = true_fx,
    FX_q05     = preds$q05,
    FX_q50     = preds$q50,
    FX_q95     = preds$q95,
    Hit_90     = (true_fx >= preds$q05) & (true_fx <= preds$q95)
  ))

  progress_bar(i - backtest_start_idx + 1, n_obs - backtest_start_idx, "  Бэктест:")
}

write.csv(backtest_results, file.path(OUT_DIR, "08a_backtest_results.csv"), row.names = FALSE)

# Метрики бэктеста
cat("\n=== МЕТРИКИ БЭКТЕСТА ===\n")
bt_metrics <- backtest_results %>%
  summarise(
    N        = n(),
    Coverage = round(mean(Hit_90, na.rm = TRUE), 3),
    RMSE_50  = round(sqrt(mean((FX_true - FX_q50)^2, na.rm = TRUE)), 5),
    MAE_50   = round(mean(abs(FX_true - FX_q50), na.rm = TRUE), 5)
  )
cat(sprintf("Наблюдений: %d | Покрытие 90%% ИД: %.1f%% (цель: 90%%)\n",
            bt_metrics$N, bt_metrics$Coverage * 100))
cat(sprintf("RMSE медианы: %.5f | MAE медианы: %.5f\n", bt_metrics$RMSE_50, bt_metrics$MAE_50))
write.csv(bt_metrics, file.path(OUT_DIR, "08b_backtest_metrics.csv"), row.names = FALSE)

# График бэктеста
p_bt <- ggplot(backtest_results, aes(x = Date)) +
  geom_ribbon(aes(ymin = FX_q05, ymax = FX_q95, fill = "90% ИД"), alpha = 0.25) +
  geom_line(aes(y = FX_true, color = "Факт"),   linewidth = 0.8) +
  geom_line(aes(y = FX_q50,  color = "Медиана"), linetype = "dashed", linewidth = 0.7) +
  geom_point(data = filter(backtest_results, !Hit_90),
             aes(x = Date, y = FX_true), color = "red", size = 2, shape = 4) +
  scale_fill_manual(values = c("90% ИД" = "steelblue")) +
  scale_color_manual(values = c("Факт" = "black", "Медиана" = "red")) +
  labs(title = "Бэктест: условный прогноз лог-доходности USD/RUB",
       subtitle = sprintf("Покрытие 90%% ИД = %.1f%% (цель 90%%). Красные × = промахи.",
                          bt_metrics$Coverage * 100),
       x = NULL, y = "Лог-доходность FX (USD/RUB)", fill = NULL, color = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "08a_backtest_plot.png"), p_bt, width = 12, height = 5, dpi = 150)

# 2. СЦЕНАРНЫЙ АНАЛИЗ: прогноз при разных ценах нефти

cat("\n=== СЦЕНАРНЫЙ АНАЛИЗ ===\n")

# Берём последний период (текущий режим)
last_regime <- best_models %>% slice_tail(n = 1)
last_history <- returns_full %>% filter(Date >= last_regime$Start)
cat(sprintf("Текущий режим: %s | Копула: %s | N_hist=%d\n",
            last_regime$Start, last_regime$Copula, nrow(last_history)))

hist_oil <- last_history$log_OIL
hist_fx  <- last_history$log_FX

# Строим копулу для текущего режима
cop_type <- as.character(last_regime$Copula)
u_mat_last <- pobs(cbind(hist_oil, hist_fx))

cop_obj_last <- tryCatch({
  if (tolower(cop_type) == "independence") return(NULL)
  cop <- make_copula(cop_type, dim = 2)
  fit <- fitCopula(cop, data = u_mat_last, method = "ml")
  fit@copula
}, error = function(e) {
  cat("Ошибка подгонки копулы для текущего режима:", conditionMessage(e), "\n")
  NULL
})

# Для каждого сценария нефти
scenario_results <- data.frame()

for (oil_scen in OIL_SCENARIOS) {
  u_oil_scen <- to_u(oil_scen, hist_oil)

  preds <- conditional_forecast(cop_obj_last, cop_type, u_oil_scen, hist_fx)

  # Последняя известная цена USD/RUB
  last_price <- tail(weekly_prices$UsdRub, 1)
  last_date  <- tail(weekly_prices$Date, 1)

  scenario_results <- rbind(scenario_results, data.frame(
    Oil_scenario_pct  = round(oil_scen * 100, 0),
    Oil_logret        = oil_scen,
    FX_q05            = preds$q05,
    FX_q50            = preds$q50,
    FX_q95            = preds$q95,
    USDRUB_q05        = round(last_price * exp(preds$q95), 2),   # инвертируем: меньший q05 FX = больший USD/RUB
    USDRUB_q50        = round(last_price * exp(preds$q50), 2),
    USDRUB_q95        = round(last_price * exp(preds$q05), 2)
  ))

  cat(sprintf("  Нефть %+.0f%%: USD/RUB q05=%.1f  медиана=%.1f  q95=%.1f\n",
              oil_scen * 100,
              last_price * exp(preds$q95),
              last_price * exp(preds$q50),
              last_price * exp(preds$q05)))
}

write.csv(scenario_results, file.path(OUT_DIR, "08c_scenario_analysis.csv"), row.names = FALSE)

# График сценарного анализа
p_scen <- ggplot(scenario_results,
                 aes(x = Oil_scenario_pct, y = USDRUB_q50,
                     ymin = USDRUB_q05, ymax = USDRUB_q95)) +
  geom_ribbon(fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  geom_hline(yintercept = last_price, linetype = "dashed", color = "grey50") +
  annotate("text", x = min(OIL_SCENARIOS) * 100, y = last_price + 1,
           label = sprintf("Текущий: %.1f", last_price), hjust = 0, color = "grey40") +
  labs(title = sprintf("Условный прогноз USD/RUB при сценариях нефти\n(режим: %s, последняя дата: %s)",
                        cop_type, format(last_date, "%Y-%m-%d")),
       x = "Изменение цены нефти, %",
       y = "USD/RUB (прогноз на следующую неделю)",
       caption = "Лента = 90% доверительный интервал") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "08c_scenario_plot.png"), p_scen, width = 8, height = 5, dpi = 150)

# 3. VaR И CVaR ИЗ КОПУЛЬНОЙ МОДЕЛИ

cat("\n=== VaR / CVaR ===\n")

# VaR(α) для лог-доходности FX при нейтральном сценарии нефти (0%)
u_oil_neutral <- to_u(0, hist_oil)
u_var <- conditional_quantile_fx(cop_obj_last, u_oil_neutral, 1 - VAR_LEVEL)
VaR_val <- if (!is.na(u_var)) quantile(hist_fx, u_var, type = 8) else quantile(hist_fx, 1 - VAR_LEVEL)

# CVaR: среднее за хвостом
CVaR_val <- mean(hist_fx[hist_fx > VaR_val], na.rm = TRUE)

cat(sprintf("Текущий режим: %s\n", cop_type))
cat(sprintf("VaR(%.0f%%):   %.5f  (USD/RUB: %.1f → %.1f)\n",
            VAR_LEVEL * 100, VaR_val,
            last_price, last_price * exp(VaR_val)))
cat(sprintf("CVaR(%.0f%%): %.5f  (USD/RUB: %.1f → %.1f)\n",
            VAR_LEVEL * 100, CVaR_val,
            last_price, last_price * exp(CVaR_val)))

var_df <- data.frame(
  Metric = c(sprintf("VaR(%.0f%%)", VAR_LEVEL*100), sprintf("CVaR(%.0f%%)", VAR_LEVEL*100)),
  LogRet = round(c(VaR_val, CVaR_val), 5),
  USDRUB = round(c(last_price * exp(VaR_val), last_price * exp(CVaR_val)), 2)
)
write.csv(var_df, file.path(OUT_DIR, "08d_var_cvar.csv"), row.names = FALSE)

cat("\nРезультаты VaR/CVaR:\n")
print(var_df)
