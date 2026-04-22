# 02_structural_breaks.R — Поиск структурных сдвигов в зависимости нефть–рубль
# Метод: KS-статистика на эмпирических копулах 
# Алгоритм: рекурсивное разбиение периода (Split) + объединение смежных подпериодов (Join)
# Входные данные: FINAL_Log_Returns.csv
# Выход: CSV с однородными периодами, KS-графики;

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(zoo)
library(ggplot2)
library(R.utils)      # withTimeout

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")

PARENT  <- ".."
OUT_DIR <- "results/02_structural_breaks"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ПАРАМЕТРЫ

DECISION_LEVEL     <- 0.99   # уровень доверия для теста: 0.90, 0.95, 0.99
SMALLEST_SUBPERIOD <- 30     # минимальный размер подпериода (недель)
DISTR_TYPE         <- "emp"  # "emp" = эмпирическое, "ghyp" = GH-распределение

# ЗАГРУЗКА ДАННЫХ

returns_full <- load_log_returns(data_dir = PARENT)
cat(sprintf("Загружено %d наблюдений: %s — %s\n",
            nrow(returns_full), min(returns_full$Date), max(returns_full$Date)))

START_DATE <- min(returns_full$Date)
END_DATE   <- max(returns_full$Date)

# ГЛОБАЛЬНЫЕ СТРУКТУРЫ ДАННЫХ

splitting_results <- list(
  to_test = data.frame(start_date = as.Date(character()), end_date = as.Date(character())),
  homogeneous_periods = data.frame(start_date = as.Date(character()), end_date = as.Date(character())),
  splitting_iterations = data.frame(
    iteration = numeric(), start_date = as.Date(character()), end_date = as.Date(character()),
    N = numeric(), possible_break_date = as.Date(character()),
    K_S = numeric(), Critical_90 = numeric(), Critical_95 = numeric(), Critical_99 = numeric(),
    is_break_90 = logical(), is_break_95 = logical(), is_break_99 = logical()
  )
)
splitting_results$to_test <- rbind(splitting_results$to_test,
  data.frame(start_date = START_DATE, end_date = END_DATE))

joining_results <- list(
  joining_iterations = data.frame(
    iteration = numeric(), start_date1 = as.Date(character()), end_date1 = as.Date(character()),
    start_date2 = as.Date(character()), end_date2 = as.Date(character()),
    joined_possible_break_date = as.Date(character()), decision = character(),
    joined_KS = numeric(), joined_Critical_90 = numeric(),
    joined_Critical_95 = numeric(), joined_Critical_99 = numeric(),
    is_break_90 = logical(), is_break_95 = logical(), is_break_99 = logical()
  ),
  homogeneous_periods = data.frame(start_date = as.Date(character()), end_date = as.Date(character()))
)

plot_counter <- 0

# KS-ТЕСТ НА СТРУКТУРНЫЙ СДВИГ

ks_copula_test <- function(period_df, period_start, period_end, save_plot = TRUE) {
  N <- nrow(period_df)
  cat(sprintf("  KS-тест: %s — %s (N=%d)\n", period_start, period_end, N))

  init_oil <- period_df$log_OIL
  init_fx  <- period_df$log_FX

  # Псевдонаблюдения
  u1 <- ecdf(init_fx)(init_fx)
  u2 <- ecdf(init_oil)(init_oil)

  # Сетка для эмпирической копулы
  grid <- seq(1/N, 1, length.out = N)
  z1 <- rep(grid, each = N)
  z2 <- rep(grid, times = N)

  KS_vec <- numeric(N - 1)
  max_KS  <- 0
  max_t   <- 1

  pb <- txtProgressBar(min = 0, max = N - 1, style = 3)
  for (t in 1:(N - 1)) {
    setTxtProgressBar(pb, t)
    u01 <- u1[1:t];       u02 <- u2[1:t]
    u11 <- u1[(t+1):N];   u12 <- u2[(t+1):N]

    EC0 <- numeric(N^2)
    EC1 <- numeric(N^2)
    for (l in 1:(N^2)) {
      EC0[l] <- sum(u01 <= z1[l] & u02 <= z2[l]) / t
      EC1[l] <- sum(u11 <= z1[l] & u12 <= z2[l]) / (N - t)
    }
    KS_t     <- max(abs(EC0 - EC1) * sqrt(t * (N - t)) / N)
    KS_vec[t] <- KS_t
    if (KS_t > max_KS) { max_KS <- KS_t; max_t <- t }
  }
  close(pb)

  Cv90 <- critical.value.ks.test(N, 0.90)
  Cv95 <- critical.value.ks.test(N, 0.95)
  Cv99 <- critical.value.ks.test(N, 0.99)

  # Сохранение графика
  if (save_plot) {
    plot_counter <<- plot_counter + 1
    Time  <- period_df$Date[1:(N - 1)]
    y_max <- max(max(KS_vec), Cv99) * 1.1

    png(file.path(OUT_DIR, sprintf("%02d_KS_%s_to_%s.png", plot_counter, period_start, period_end)),
        width = 1200, height = 700, res = 130)
    plot(Time, KS_vec, type = "l", lwd = 1.5,
         main  = sprintf("KS-статистика: %s — %s", period_start, period_end),
         xlab  = "Дата", ylab = "KS-статистика",
         ylim  = c(0, y_max), xaxt = "n")
    ticks <- pretty(Time, n = 6)
    axis(1, at = ticks, labels = format(ticks, "%Y-%m"), cex.axis = 0.85)
    abline(h = Cv90, col = "darkgreen", lty = 2, lwd = 1.2)
    abline(h = Cv95, col = "blue",      lty = 2, lwd = 1.2)
    abline(h = Cv99, col = "red",       lty = 2, lwd = 1.5)
    break_date <- Time[max_t]
    abline(v = break_date, col = "purple", lty = 1, lwd = 1.5)
    text(break_date, max_KS,
         labels = sprintf("Max KS=%.3f\n%s", max_KS, format(break_date, "%Y-%m-%d")),
         pos = 4, cex = 0.8)
    legend("topright", bty = "n", cex = 0.8,
           legend = c("KS", sprintf("90%% (%.3f)", Cv90),
                      sprintf("95%% (%.3f)", Cv95), sprintf("99%% (%.3f)", Cv99)),
           col = c("black", "darkgreen", "blue", "red"), lty = c(1, 2, 2, 2))
    dev.off()
  }

  data.frame(
    possible_break_date = period_df$Date[max_t],
    KS = max_KS,
    Critical_90 = Cv90, Critical_95 = Cv95, Critical_99 = Cv99,
    is_break_90 = max_KS >= Cv90,
    is_break_95 = max_KS >= Cv95,
    is_break_99 = max_KS >= Cv99
  )
}

# АЛГОРИТМ РАЗБИЕНИЯ (SPLIT)

find_all_breaks <- function() {
  iteration <- 0
  while (nrow(splitting_results$to_test) > 0) {
    iteration <- iteration + 1
    period <- splitting_results$to_test[1, ]
    splitting_results$to_test <<- splitting_results$to_test[-1, , drop = FALSE]
    p_start <- period$start_date
    p_end   <- period$end_date

    period_df <- preprocess(returns_full, p_start, p_end)

    if (nrow(period_df) <= SMALLEST_SUBPERIOD) {
      splitting_results$homogeneous_periods <<- rbind(
        splitting_results$homogeneous_periods,
        data.frame(start_date = p_start, end_date = p_end)
      )
      cat(sprintf("[Split iter %d] Мало данных (%d <= %d), добавляем как однородный.\n",
                  iteration, nrow(period_df), SMALLEST_SUBPERIOD))
      next
    }

    res <- ks_copula_test(period_df, p_start, p_end)
    is_break <- switch(as.character(DECISION_LEVEL),
      "0.99" = res$is_break_99, "0.95" = res$is_break_95, "0.9" = res$is_break_90)

    new_row <- cbind(data.frame(iteration = iteration, start_date = p_start,
                                end_date = p_end, N = nrow(period_df)), res)
    splitting_results$splitting_iterations <<- rbind(splitting_results$splitting_iterations, new_row)

    if (!is_break) {
      splitting_results$homogeneous_periods <<- rbind(
        splitting_results$homogeneous_periods,
        data.frame(start_date = p_start, end_date = p_end)
      )
      cat(sprintf("[Split iter %d] Сдвига нет. Период %s — %s однороден.\n",
                  iteration, p_start, p_end))
    } else {
      bd <- res$possible_break_date
      cat(sprintf("[Split iter %d] Сдвиг обнаружен %s (KS=%.4f). Разбиваем.\n",
                  iteration, bd, res$KS))
      splitting_results$to_test <<- rbind(splitting_results$to_test,
        data.frame(start_date = p_start, end_date = bd))
      splitting_results$to_test <<- rbind(splitting_results$to_test,
        data.frame(start_date = bd + 1,  end_date = p_end))
    }
  }
  splitting_results$homogeneous_periods <<- splitting_results$homogeneous_periods[
    order(splitting_results$homogeneous_periods$start_date), ]
}

# АЛГОРИТМ ОБЪЕДИНЕНИЯ (JOIN)

join_periods <- function() {
  hom <- splitting_results$homogeneous_periods
  continue <- TRUE
  iteration <- 0

  while (continue && nrow(hom) > 1) {
    best_j  <- 0
    best_KS <- Inf
    best_row <- NULL

    for (j in 1:(nrow(hom) - 1)) {
      iteration <- iteration + 1
      joint_df <- preprocess(returns_full, hom[j, 1], hom[j + 1, 2])

      if (nrow(joint_df) <= SMALLEST_SUBPERIOD) {
        best_j <- j; best_KS <- -Inf
        best_row <- data.frame(
          iteration = iteration,
          start_date1 = hom[j, 1], end_date1 = hom[j, 2],
          start_date2 = hom[j+1, 1], end_date2 = hom[j+1, 2],
          decision = "JOIN_BY_SIZE",
          joined_KS = NA, joined_Critical_99 = NA,
          is_break_99 = FALSE
        )
        next
      }

      res <- ks_copula_test(joint_df, hom[j, 1], hom[j + 1, 2])
      is_break <- switch(as.character(DECISION_LEVEL),
        "0.99" = res$is_break_99, "0.95" = res$is_break_95, "0.9" = res$is_break_90)

      if (!is_break && res$KS < best_KS) {
        best_KS <- res$KS
        best_j  <- j
        best_row <- data.frame(
          iteration = iteration,
          start_date1 = hom[j, 1], end_date1 = hom[j, 2],
          start_date2 = hom[j+1, 1], end_date2 = hom[j+1, 2],
          decision = "JOIN_BY_KS",
          joined_KS = res$KS, joined_Critical_99 = res$Critical_99,
          is_break_99 = res$is_break_99
        )
      }
    }

    if (best_j > 0) {
      cat(sprintf("[Join] Объединяем периоды %d и %d\n", best_j, best_j + 1))
      hom[best_j, 2] <- hom[best_j + 1, 2]
      hom <- hom[-(best_j + 1), , drop = FALSE]
      joining_results$joining_iterations <<- rbind(joining_results$joining_iterations, best_row)
    } else {
      continue <- FALSE
    }
  }
  joining_results$homogeneous_periods <<- hom[order(hom$start_date), ]
}

cat("\n=== ЭТАП РАЗБИЕНИЯ (SPLIT) ===\n")
t0 <- Sys.time()
find_all_breaks()
cat(sprintf("\nРазбиение завершено за %.1f мин\n", as.numeric(Sys.time() - t0, units = "mins")))

cat(sprintf("Подпериодов после Split: %d\n", nrow(splitting_results$homogeneous_periods)))
print(splitting_results$homogeneous_periods)

write.csv(splitting_results$homogeneous_periods,
          file.path(OUT_DIR, "02a_periods_after_split.csv"), row.names = FALSE)
write.csv(splitting_results$splitting_iterations,
          file.path(OUT_DIR, "02b_splitting_iterations.csv"), row.names = FALSE)

cat("\n=== ЭТАП ОБЪЕДИНЕНИЯ (JOIN) ===\n")
t0 <- Sys.time()
join_periods()
cat(sprintf("\nОбъединение завершено за %.1f мин\n", as.numeric(Sys.time() - t0, units = "mins")))

cat(sprintf("Финальных однородных периодов: %d\n", nrow(joining_results$homogeneous_periods)))
print(joining_results$homogeneous_periods)

write.csv(joining_results$homogeneous_periods,
          file.path(OUT_DIR, "02c_final_homogeneous_periods.csv"), row.names = FALSE)
write.csv(joining_results$joining_iterations,
          file.path(OUT_DIR, "02d_joining_iterations.csv"), row.names = FALSE)
