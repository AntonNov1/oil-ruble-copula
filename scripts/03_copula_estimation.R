# 03_copula_estimation.R — Оценка копул по однородным периодам
# Метод: полупараметрический 
# Копулы: Frank, Gumbel, Clayton, Gaussian, Student (с лагами до max_lag)
# Выбор лучшей: по BIC
# Выход: CSV со всеми моделями и лучшими по BIC

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(copula)
library(zoo)
library(R.utils)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/03_copula_estimation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

MAX_LAG   <- 5      # максимальный лаг для обоих рядов
TIMEOUT   <- 120    # секунд на одну копулу
COP_METHOD <- "ml"  # метод оценки: "ml" или "itau"

periods_file <- if (file.exists("results/02_structural_breaks/02c_final_homogeneous_periods.csv")) {
  "results/02_structural_breaks/02c_final_homogeneous_periods.csv"
} else {
  file.path(PARENT, "2_final_periods_30.csv")
}
cat("Читаем периоды из:", periods_file, "\n")

homogeneous_periods <- read.csv(periods_file, stringsAsFactors = FALSE)
# Поддержка разных форматов файлов
if ("StartDate" %in% names(homogeneous_periods)) {
  homogeneous_periods <- homogeneous_periods %>%
    rename(start_date = StartDate, end_date = EndDate)
} else if ("start_date" %in% names(homogeneous_periods)) {
  # уже в нужном формате
} else if ("Start" %in% names(homogeneous_periods)) {
  homogeneous_periods <- homogeneous_periods %>% rename(start_date = Start, end_date = End)
}
homogeneous_periods$start_date <- as.Date(homogeneous_periods$start_date)
homogeneous_periods$end_date   <- as.Date(homogeneous_periods$end_date)

cat(sprintf("Загружено %d периодов.\n", nrow(homogeneous_periods)))

returns_full <- load_log_returns(data_dir = PARENT)

# ФУНКЦИИ

# Создаёт матрицу псевдонаблюдений с лагами (dim = (1+lag_fx)+(1+lag_oil))
make_pobs_with_lags <- function(period_df, lag_fx, lag_oil) {
  n <- nrow(period_df)
  oil <- period_df$log_OIL
  fx  <- period_df$log_FX

  shift_vec <- function(v, k) c(rep(NA, k), v[1:(length(v) - k)])

  # Ряды: [fx_0, fx_1,...fx_lag_fx, oil_0, oil_1,...oil_lag_oil]
  mat_raw <- data.frame(fx_0 = fx, oil_0 = oil)
  for (i in seq_len(lag_fx))  mat_raw[[paste0("fx_",  i)]] <- shift_vec(fx,  i)
  for (i in seq_len(lag_oil)) mat_raw[[paste0("oil_", i)]] <- shift_vec(oil, i)

  mat_raw <- mat_raw[complete.cases(mat_raw), ]
  if (nrow(mat_raw) < 10) return(NULL)
  as.matrix(pobs(mat_raw))
}

# Строка результата копулы
make_result_row <- function(period_id, start, end, cop_type, lag_fx, lag_oil,
                             loglik, k, n, params, se_params) {
  data.frame(
    Period_ID    = period_id,
    Start        = as.character(start),
    End          = as.character(end),
    Lag_FX       = lag_fx,
    Lag_Oil      = lag_oil,
    Dim          = 2 + lag_fx + lag_oil,
    Copula       = cop_type,
    Params       = paste(round(params, 4), collapse = "; "),
    SE_Params    = paste(round(se_params, 4), collapse = "; "),
    LogLik       = round(loglik, 4),
    AIC          = round(count_AIC(loglik, k), 4),
    BIC          = round(count_BIC(loglik, k, n), 4),
    stringsAsFactors = FALSE
  )
}

# Основная функция оценки одной копулы с лагами
estimate_one <- function(u_mat, cop_type, period_id, start, end, lag_fx, lag_oil) {
  n <- nrow(u_mat)
  dim_cop <- ncol(u_mat)

  tryCatch({
    cop_fit <- NA
    withTimeout({
      cop_obj <- switch(cop_type,
        "Frank"    = frankCopula(dim = dim_cop),
        "Gumbel"   = gumbelCopula(dim = dim_cop),
        "Clayton"  = claytonCopula(dim = dim_cop),
        "Gaussian" = normalCopula(dim = dim_cop, dispstr = "un"),
        "Student"  = tCopula(dim = dim_cop, dispstr = "un"),
        stop("Unknown copula")
      )
      cop_fit <- fitCopula(cop_obj, data = u_mat, method = COP_METHOD)
    }, timeout = TIMEOUT)

    if (is.logical(cop_fit) && is.na(cop_fit)) return(NULL)

    s     <- summary(cop_fit)
    coefs <- s$coefficients
    params    <- coefs[, 1]
    se_params <- coefs[, 2]
    L     <- cop_fit@loglik
    k     <- length(params)

    make_result_row(period_id, start, end, cop_type, lag_fx, lag_oil, L, k, n, params, se_params)

  }, error = function(e) {
    make_result_row(period_id, start, end, cop_type, lag_fx, lag_oil, NA, NA, n, NA, NA)
  })
}

# ОСНОВНОЙ ЦИКЛ

all_results  <- list()
COPULA_TYPES <- c("Frank", "Gumbel", "Clayton", "Gaussian", "Student")

for (i in seq_len(nrow(homogeneous_periods))) {
  p_start <- homogeneous_periods$start_date[i]
  p_end   <- homogeneous_periods$end_date[i]
  period_df <- preprocess(returns_full, p_start, p_end)
  n_obs <- nrow(period_df)

  cat(sprintf("\n[Период %d/%d] %s — %s (N=%d)\n",
              i, nrow(homogeneous_periods), p_start, p_end, n_obs))

  if (n_obs < 15) {
    cat("  Слишком мало наблюдений, пропускаем.\n"); next
  }

  for (lag_fx in 0:MAX_LAG) {
    for (lag_oil in 0:MAX_LAG) {
      u_mat <- make_pobs_with_lags(period_df, lag_fx, lag_oil)
      if (is.null(u_mat)) next

      for (cop_type in COPULA_TYPES) {
        cat(sprintf("  %s lag_fx=%d lag_oil=%d... ", cop_type, lag_fx, lag_oil))
        row <- estimate_one(u_mat, cop_type, i, p_start, p_end, lag_fx, lag_oil)
        if (!is.null(row)) {
          all_results <- c(all_results, list(row))
          cat(sprintf("BIC=%.3f\n", row$BIC))
        } else {
          cat("FAILED\n")
        }
      }
    }
  }
}

# СОХРАНЕНИЕ РЕЗУЛЬТАТОВ

all_df <- do.call(rbind, all_results)

if (!is.null(all_df) && nrow(all_df) > 0) {
  write.csv(all_df, file.path(OUT_DIR, "03a_all_copulas.csv"), row.names = FALSE)

  # Лучшая по BIC для каждого периода
  best_df <- all_df %>%
    filter(!is.na(BIC)) %>%
    group_by(Period_ID) %>%
    slice_min(order_by = BIC, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Если лучший BIC > 0 → независимость (нет значимой связи)
  best_df <- best_df %>%
    mutate(
      Copula = ifelse(BIC > 0, "independence", Copula),
      Params = ifelse(BIC > 0, NA_character_, Params)
    )

  write.csv(best_df, file.path(OUT_DIR, "03b_best_copulas_by_BIC.csv"), row.names = FALSE)

  cat("\n=== ЛУЧШИЕ КОПУЛЫ ПО ПЕРИОДАМ ===\n")
  print(best_df[, c("Period_ID", "Start", "End", "Copula", "Lag_FX", "Lag_Oil", "Params", "BIC")])
} else {
  cat("Нет результатов для сохранения.\n")
}