# 01_data_prep.R — Загрузка данных, EDA, тесты на стационарность
# Данные: нефть Brent (Oil, USD/барр.) и курс USD/RUB (UsdRub)
# Источник: FINAL_Daily_Prices.csv, FINAL_Weekly_Prices.csv
# Выход: описательная статистика, графики, тесты; 

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(ggplot2)
library(tseries)   # adf.test
library(moments)   # skewness, kurtosis
library(zoo)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")   

PARENT  <- ".."             
OUT_DIR <- "results/01_data_prep"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ЗАГРУЗКА

cat("=== ЗАГРУЗКА ДАННЫХ ===\n")

# Дневные цены
daily <- read.csv(file.path(PARENT, "FINAL_Daily_Prices.csv"), stringsAsFactors = FALSE)
names(daily)[1] <- "Date"
daily$Date <- as.Date(daily$Date)
daily <- daily %>% arrange(Date)
cat(sprintf("Дневные данные: %s — %s (%d наблюдений)\n",
            min(daily$Date), max(daily$Date), nrow(daily)))

# Недельные цены
weekly <- read.csv(file.path(PARENT, "FINAL_Weekly_Prices.csv"), stringsAsFactors = FALSE)
names(weekly)[1] <- "Date"
weekly$Date <- as.Date(weekly$Date)
weekly <- weekly %>% arrange(Date)
cat(sprintf("Недельные данные: %s — %s (%d наблюдений)\n",
            min(weekly$Date), max(weekly$Date), nrow(weekly)))

# Недельные лог-доходности
lr <- read.csv(file.path(PARENT, "FINAL_Log_Returns.csv"), stringsAsFactors = FALSE)
names(lr)[1] <- "Date"
lr$Date <- as.Date(lr$Date)
lr <- lr %>% arrange(Date)
cat(sprintf("Лог-доходности: %s — %s (%d наблюдений)\n",
            min(lr$Date), max(lr$Date), nrow(lr)))

# ОПИСАТЕЛЬНАЯ СТАТИСТИКА

cat("\n=== ОПИСАТЕЛЬНАЯ СТАТИСТИКА ===\n")

stats_table <- function(x, name) {
  data.frame(
    Ряд        = name,
    N          = length(x),
    Среднее    = round(mean(x, na.rm = TRUE), 5),
    Медиана    = round(median(x, na.rm = TRUE), 5),
    Мин        = round(min(x, na.rm = TRUE), 5),
    Макс       = round(max(x, na.rm = TRUE), 5),
    СКО        = round(sd(x, na.rm = TRUE), 5),
    Асим       = round(skewness(x, na.rm = TRUE), 4),
    Эксцесс    = round(kurtosis(x, na.rm = TRUE) - 3, 4)
  )
}

tab <- rbind(
  stats_table(weekly$Oil,    "Oil (уровень)"),
  stats_table(weekly$UsdRub, "UsdRub (уровень)"),
  stats_table(lr$Oil,        "LogRet Oil"),
  stats_table(lr$UsdRub,     "LogRet UsdRub")
)

print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT_DIR, "descriptive_stats.csv"), row.names = FALSE)

# КОРРЕЛЯЦИЯ

cat("\n=== КОРРЕЛЯЦИЯ ЛОГ-ДОХОДНОСТЕЙ ===\n")
pearson <- cor(lr$Oil, lr$UsdRub, method = "pearson")
kendall <- cor(lr$Oil, lr$UsdRub, method = "kendall")
cat(sprintf("  Пирсон (линейная):   %.4f\n", pearson))
cat(sprintf("  Кендалл (ранговая):  %.4f\n", kendall))

# ТЕСТ АДФ НА СТАЦИОНАРНОСТЬ

cat("\n=== ТЕСТ АДФ (ДИКИ-ФУЛЛЕРА) ===\n")
run_adf <- function(x, name) {
  t <- adf.test(x, alternative = "stationary")
  cat(sprintf("  %-30s  стат = %6.3f  p = %.4f  %s\n",
              name, t$statistic, t$p.value,
              ifelse(t$p.value < 0.05, "[стационарен]", "[нестационарен]")))
}
run_adf(weekly$Oil,    "Oil (уровень)")
run_adf(weekly$UsdRub, "UsdRub (уровень)")
run_adf(lr$Oil,        "LogRet Oil")
run_adf(lr$UsdRub,     "LogRet UsdRub")

# 5. ТЕСТ НА НОРМАЛЬНОСТЬ (КОЛМОГОРОВ–СМИРНОВ)

cat("\n=== ТЕСТ КС НА НОРМАЛЬНОСТЬ (лог-доходности) ===\n")
ks_oil <- ks.test(scale(lr$Oil),    "pnorm")
ks_rub <- ks.test(scale(lr$UsdRub), "pnorm")
cat(sprintf("  Oil:    стат = %.4f  p = %.4f\n", ks_oil$statistic, ks_oil$p.value))
cat(sprintf("  UsdRub: стат = %.4f  p = %.4f\n", ks_rub$statistic, ks_rub$p.value))

# ГРАФИКИ

cat("\n=== ПОСТРОЕНИЕ ГРАФИКОВ ===\n")

# 6.1. Динамика цен
df_prices <- tidyr::pivot_longer(weekly, cols = c("Oil", "UsdRub"),
                                 names_to = "Series", values_to = "Value")
p_prices <- ggplot(df_prices, aes(x = Date, y = Value, color = Series)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~Series, scales = "free_y", ncol = 1) +
  labs(title = "Еженедельные цены: Brent и USD/RUB",
       x = NULL, y = NULL) +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "01_weekly_prices.png"), p_prices, width = 10, height = 5, dpi = 150)

# 6.2. Лог-доходности
df_lr <- tidyr::pivot_longer(lr, cols = c("Oil", "UsdRub"),
                             names_to = "Series", values_to = "Return")
p_lr <- ggplot(df_lr, aes(x = Date, y = Return, color = Series)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~Series, scales = "free_y", ncol = 1) +
  labs(title = "Недельные лог-доходности",
       x = NULL, y = "Лог-доходность") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "02_log_returns.png"), p_lr, width = 10, height = 5, dpi = 150)

# 6.3. Диаграмма рассеяния лог-доходностей
p_scatter <- ggplot(lr, aes(x = Oil, y = UsdRub)) +
  geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(title = "Связь: лог-доходности Brent и USD/RUB",
       x = "LogRet(Oil)", y = "LogRet(UsdRub)") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "03_scatter_returns.png"), p_scatter, width = 6, height = 5, dpi = 150)

# 6.4. Псевдонаблюдения 
u_oil <- rank(lr$Oil) / (nrow(lr) + 1)
u_rub <- rank(lr$UsdRub) / (nrow(lr) + 1)
p_pobs <- ggplot(data.frame(u_oil, u_rub), aes(x = u_oil, y = u_rub)) +
  geom_point(alpha = 0.3, size = 1.2, color = "steelblue") +
  labs(title = "Псевдонаблюдения (копульная плоскость [0,1]²)",
       x = "U(Oil)", y = "U(UsdRub)") +
  theme_minimal()
ggsave(file.path(OUT_DIR, "04_pseudo_observations.png"), p_pobs, width = 6, height = 5, dpi = 150)
