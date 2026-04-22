# 00_helpers.R — Вспомогательные функции, используемые во всех скриптах

library(copula)
library(dplyr)
library(tidyr)
library(zoo)

# ПУТИ К ДАННЫМ

get_project_root <- function() {
  dirname(dirname(sys.frame(1)$ofile))
}

DATA_DIR    <- file.path(dirname(getwd()), "")   # заполнится в каждом скрипте
RESULTS_DIR <- file.path(dirname(getwd()), "results")

# ЗАГРУЗКА ДАННЫХ

#' Загружает еженедельные лог-доходности из FINAL_Log_Returns.csv
#' Возвращает data.frame с колонками: Date, log_OIL, log_FX
#' log_FX = -log(UsdRub) => положительный, если рубль укрепляется
load_log_returns <- function(data_dir = ".") {
  f <- file.path(data_dir, "FINAL_Log_Returns.csv")
  if (!file.exists(f)) stop(paste("Файл не найден:", f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  names(df)[1] <- "Date"
  df$Date <- as.Date(df$Date)
  df <- df %>%
    rename(OIL_ret = Oil, FX_ret = UsdRub) %>%
    mutate(
      log_OIL = OIL_ret,
      log_FX  = -FX_ret    # знак "–": рост USDRUB = ослабление рубля => берём с минусом
    ) %>%
    filter(!is.na(log_OIL), !is.na(log_FX)) %>%
    arrange(Date)
  df
}

#' Загружает еженедельные цены из FINAL_Weekly_Prices.csv
load_weekly_prices <- function(data_dir = ".") {
  f <- file.path(data_dir, "FINAL_Weekly_Prices.csv")
  if (!file.exists(f)) stop(paste("Файл не найден:", f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  names(df)[1] <- "Date"
  df$Date <- as.Date(df$Date)
  df <- df %>% rename(Oil = Oil, UsdRub = UsdRub) %>% arrange(Date)
  df
}

#' Загружает результаты лучших моделей из 2_Best_Models.csv
load_best_models <- function(data_dir = ".") {
  f <- file.path(data_dir, "2_Best_Models.csv")
  if (!file.exists(f)) stop(paste("Файл не найден:", f))
  df <- read.csv(f, stringsAsFactors = FALSE)
  # Нормализуем колонки для разных форматов файла
  if (!"Period_ID" %in% names(df)) df$Period_ID <- seq_len(nrow(df))
  if ("Lag_USD"   %in% names(df) && !"Lag_FX"  %in% names(df)) df$Lag_FX  <- df$Lag_USD
  if ("Lag_Oil"   %in% names(df) && !"Lag_Oil" %in% names(df)) {}  # уже есть
  df$Start <- as.Date(df$Start)
  df$End   <- as.Date(df$End)
  df
}

# ПРЕДОБРАБОТКА

#' Фильтрует данные по диапазону дат
preprocess <- function(df, start_date, end_date) {
  df %>% filter(Date >= as.Date(start_date), Date <= as.Date(end_date))
}

#' Трансформирует данные в псевдонаблюдения [0,1]^2 (бивариатный случай, lag=0)
uniforming_bivariate <- function(df) {
  mat <- pobs(as.matrix(df[, c("log_OIL", "log_FX")]))
  as.data.frame(mat) %>% rename(u_oil = log_OIL, u_fx = log_FX)
}

# КРИТИЧЕСКИЕ ЗНАЧЕНИЯ KS-ТЕСТА НА СТРУКТУРНЫЙ СДВИГ
# Источник: стандартные значения теста Колмогорова–Смирнова,
#           масштабированные на размер выборки

critical.value.ks.test <- function(N, confidence_level = 0.95) {
  const <- switch(as.character(round(confidence_level, 2)),
    "0.9"  = 1.224,
    "0.95" = 1.358,
    "0.99" = 1.628,
    stop(paste("Неподдерживаемый уровень доверия:", confidence_level))
  )
  const / sqrt(N)
}

# ИНФОРМАЦИОННЫЕ КРИТЕРИИ

count_AIC <- function(L, k)      2 * k - 2 * L
count_BIC <- function(L, k, n)   k * log(n) - 2 * L

# ПОДБОР И ОЦЕНКА КОПУЛ

COPULA_TYPES <- c("Frank", "Gumbel", "Clayton", "Gaussian", "Student")

#' Создаёт объект копулы по типу
make_copula <- function(copula_type, dim = 2) {
  switch(copula_type,
    "Frank"    = frankCopula(dim = dim),
    "Gumbel"   = gumbelCopula(dim = dim),
    "Clayton"  = claytonCopula(dim = dim),
    "Gaussian" = normalCopula(dim = dim),
    "Student"  = tCopula(dim = dim),
    "Joe"      = joeCopula(dim = dim),
    stop(paste("Неизвестный тип копулы:", copula_type))
  )
}

#' Подбирает бивариатную копулу (lag=0) методом максимального правдоподобия.
#' Возвращает список: type, fit, copula, loglik, AIC, BIC, params, n
fit_copula_bivariate <- function(u_mat, copula_type, method = "ml") {
  n <- nrow(u_mat)
  u_mat <- as.matrix(u_mat)

  result <- tryCatch({
    cop_obj <- make_copula(copula_type, dim = 2)
    fit     <- fitCopula(cop_obj, data = u_mat, method = method)
    k       <- length(fit@copula@parameters)
    L       <- fit@loglik
    list(
      type   = copula_type,
      fit    = fit,
      copula = fit@copula,
      loglik = L,
      AIC    = count_AIC(L, k),
      BIC    = count_BIC(L, k, n),
      params = fit@copula@parameters,
      k      = k,
      n      = n,
      ok     = TRUE
    )
  }, error = function(e) {
    list(type = copula_type, fit = NULL, copula = NULL,
         loglik = NA, AIC = NA, BIC = NA, params = NA, k = NA, n = n, ok = FALSE)
  })
  result
}

#' Перебирает все типы копул и возвращает лучшую по BIC
fit_best_copula_bivariate <- function(u_mat, types = COPULA_TYPES, method = "ml") {
  results <- lapply(types, function(ct) fit_copula_bivariate(u_mat, ct, method))
  bic_vals <- sapply(results, function(r) ifelse(is.na(r$BIC), Inf, r$BIC))
  best <- results[[which.min(bic_vals)]]
  list(best = best, all = results)
}

#' Оценивает логарифмическое правдоподобие заданной копулы (с фиксированными параметрами)
#' на новых данных u_mat
eval_loglik <- function(copula_obj, u_mat) {
  if (is.null(copula_obj)) return(NA_real_)
  tryCatch({
    dens <- dCopula(as.matrix(u_mat), copula_obj)
    dens <- pmax(dens, .Machine$double.eps)
    sum(log(dens))
  }, error = function(e) NA_real_)
}

#' Создаёт объект копулы с заданными параметрами из строки best_models
copula_from_model_row <- function(row) {
  fam    <- row$Copula
  params <- as.numeric(unlist(strsplit(gsub(";", " ", as.character(row$Params)), "\\s+")))
  params <- params[!is.na(params)]
  tryCatch({
    switch(fam,
      "Frank"    = frankCopula(param = params[1], dim = 2),
      "Gumbel"   = gumbelCopula(param = params[1], dim = 2),
      "Clayton"  = claytonCopula(param = params[1], dim = 2),
      "Gaussian" = normalCopula(param = params[1], dim = 2),
      "Student"  = tCopula(param = params[1],
                           df = if (length(params) >= 2) params[2] else 4,
                           dim = 2),
      NULL
    )
  }, error = function(e) NULL)
}

# ПРОГРЕСС-БАР 

progress_bar <- function(i, total, label = "") {
  pct <- round(i / total * 100)
  cat(sprintf("\r%s %3d%%", label, pct))
  if (i == total) cat("\n")
}
