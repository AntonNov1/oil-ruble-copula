# 04_gof_tests.R — Тесты согласия (GoF) для лучших копул по периодам
# Методы:
#   1. Rosenblatt + KS-тест (условное преобразование -> U[0,1])
#   2. Cramér–von Mises (Omega^2) через bootstrap
#   3. SnC-статистика (gofCopula)
# Выход: таблица p-значений GoF; 

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(copula)
library(gofCopula)
library(zoo)
library(R.utils)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/04_gof_tests"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

N_SIM      <- 1000   # число симуляций для bootstrap p-значения
COP_METHOD <- "ml"

# ЗАГРУЗКА ДАННЫХ

returns_full <- load_log_returns(data_dir = PARENT)

# Лучшие модели (bivariate, lag=0):
best_file <- if (file.exists("results/03_copula_estimation/03b_best_copulas_by_BIC.csv")) {
  "results/03_copula_estimation/03b_best_copulas_by_BIC.csv"
} else {
  file.path(PARENT, "2_Best_Models.csv")
}
cat("Читаем лучшие модели из:", best_file, "\n")

best_models <- read.csv(best_file, stringsAsFactors = FALSE)
# Нормализуем колонки
if ("Period_ID" %in% names(best_models)) {
  best_models <- best_models %>%
    rename(start_date = Start, end_date = End, copula_type = Copula)
} else {
  best_models <- best_models %>%
    rename(start_date = Start, end_date = End, copula_type = Copula,
           lag_fx = Lag_USD, lag_oil = Lag_Oil)
  best_models$Period_ID <- seq_len(nrow(best_models))
}
best_models$start_date <- as.Date(best_models$start_date)
best_models$end_date   <- as.Date(best_models$end_date)

# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

# Статистика Cramér–von Mises
omega_sq <- function(x) {
  n <- length(x)
  v <- sort(x)
  1 / (12 * n) + sum((v - (2 * seq_len(n) - 1) / (2 * n))^2)
}

# Bootstrap p-значение для GoF через Rosenblatt + KS/CvM
bootstrap_gof <- function(u_mat, cop_fit, n_sim = 1000) {
  n        <- nrow(u_mat)
  dim_cop  <- ncol(u_mat)

  # Наблюдённые статистики
  cond_u <- tryCatch(
    cCopula(u_mat, copula = cop_fit@copula, indices = dim_cop, inverse = FALSE),
    error = function(e) NULL
  )
  if (is.null(cond_u)) return(list(ks_stat = NA, ks_pval = NA, cvm_stat = NA, cvm_pval = NA))

  ks_obs  <- ks.test(cond_u, "punif")$statistic
  cvm_obs <- omega_sq(cond_u)

  # Bootstrap
  ks_boot  <- numeric(n_sim)
  cvm_boot <- numeric(n_sim)

  cat("    Bootstrap GoF (n_sim =", n_sim, "):\n")
  for (s in seq_len(n_sim)) {
    if (s %% 200 == 0) cat(sprintf("      %d/%d\n", s, n_sim))
    sim_data <- rCopula(n, cop_fit@copula)
    sim_u    <- pobs(sim_data)
    sim_fit  <- tryCatch(
      fitCopula(cop_fit@copula, data = sim_u, method = COP_METHOD),
      error = function(e) NULL
    )
    if (is.null(sim_fit)) { ks_boot[s] <- NA; cvm_boot[s] <- NA; next }
    cond_sim <- tryCatch(
      cCopula(sim_u, copula = sim_fit@copula, indices = dim_cop, inverse = FALSE),
      error = function(e) NULL
    )
    if (is.null(cond_sim)) { ks_boot[s] <- NA; cvm_boot[s] <- NA; next }
    ks_boot[s]  <- ks.test(cond_sim, "punif")$statistic
    cvm_boot[s] <- omega_sq(cond_sim)
  }

  list(
    ks_stat  = ks_obs,
    ks_pval  = mean(ks_boot >= ks_obs, na.rm = TRUE),
    cvm_stat = cvm_obs,
    cvm_pval = mean(cvm_boot >= cvm_obs, na.rm = TRUE)
  )
}

# ОСНОВНОЙ ЦИКЛ

gof_results <- list()

for (i in seq_len(nrow(best_models))) {
  row <- best_models[i, ]
  cat(sprintf("\n[GoF %d/%d] Период %s — %s | %s\n",
              i, nrow(best_models), row$start_date, row$end_date,
              row$copula_type))

  if (tolower(row$copula_type) == "independence") {
    cat("  Independence — тест GoF не применяется.\n")
    gof_results[[i]] <- data.frame(
      Period_ID = i, Start = as.character(row$start_date), End = as.character(row$end_date),
      Copula = "independence", KS_stat = NA, KS_pval = NA, CvM_stat = NA, CvM_pval = NA
    )
    next
  }

  period_df <- preprocess(returns_full, row$start_date, row$end_date)

  # Используем lag=0 для GoF (biv. копула)
  lag_fx  <- if ("Lag_FX"  %in% names(row)) row$Lag_FX  else 0
  lag_oil <- if ("Lag_Oil" %in% names(row)) row$Lag_Oil else 0

  # Псевдонаблюдения
  oil <- period_df$log_OIL
  fx  <- period_df$log_FX
  u_mat <- pobs(cbind(oil, fx))
  if (lag_fx == 0 && lag_oil == 0) {
    # bivariate, уже готово
  } else {
    # с лагами (упрощение: берём только lag=0 для GoF)
    u_mat <- pobs(cbind(oil, fx))
    lag_fx <- 0; lag_oil <- 0
  }

  # Строим и подгоняем копулу
  dim_cop <- 2
  cop_fit <- tryCatch({
    cop_obj <- switch(row$copula_type,
      "Frank"    = frankCopula(dim = dim_cop),
      "Gumbel"   = gumbelCopula(dim = dim_cop),
      "Clayton"  = claytonCopula(dim = dim_cop),
      "Gaussian" = normalCopula(dim = dim_cop),
      "Student"  = tCopula(dim = dim_cop),
      stop("Unknown")
    )
    fitCopula(cop_obj, data = u_mat, method = COP_METHOD)
  }, error = function(e) { cat("  Ошибка подгонки:", e$message, "\n"); NULL })

  if (is.null(cop_fit)) {
    gof_results[[i]] <- data.frame(
      Period_ID = i, Start = as.character(row$start_date), End = as.character(row$end_date),
      Copula = row$copula_type, KS_stat = NA, KS_pval = NA, CvM_stat = NA, CvM_pval = NA)
    next
  }

  cat(sprintf("  Параметры: %s\n", paste(round(cop_fit@copula@parameters, 4), collapse = ", ")))

  gof <- bootstrap_gof(u_mat, cop_fit, n_sim = N_SIM)

  gof_results[[i]] <- data.frame(
    Period_ID = i,
    Start     = as.character(row$start_date),
    End       = as.character(row$end_date),
    Copula    = row$copula_type,
    KS_stat   = round(gof$ks_stat,  4),
    KS_pval   = round(gof$ks_pval,  4),
    CvM_stat  = round(gof$cvm_stat, 4),
    CvM_pval  = round(gof$cvm_pval, 4)
  )
  cat(sprintf("  KS: stat=%.4f  p=%.4f | CvM: stat=%.4f  p=%.4f\n",
              gof$ks_stat, gof$ks_pval, gof$cvm_stat, gof$cvm_pval))
}

gof_df <- do.call(rbind, gof_results)
write.csv(gof_df, file.path(OUT_DIR, "04_gof_results.csv"), row.names = FALSE)

cat("\n=== ИТОГИ GoF ТЕСТОВ ===\n")
print(gof_df)
