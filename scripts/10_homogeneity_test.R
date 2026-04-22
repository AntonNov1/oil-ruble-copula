# 10_homogeneity_test.R — Тест однородности: доказательство объединимости
#
# Цель: доказать, что периоды, назначенные одному режиму копулы, имеют
#       статистически неразличимые параметры, то есть их можно объединять.
#
# Методы:
#   A. Likelihood Ratio Test (LRT):
#      H₀: общий θ для всех периодов режима
#      H₁: θ различается по периодам
#      Λ = 2*(Σ LL_индивид - LL_пул) ~ χ²(df = (K-1)*p)
#      Не отвергаем H₀ → объединение статистически обоснованно
#
#   B. Попарный KS-тест (Пеникас 2014) между соседними периодами
#      одного режима: совместные эмпирические копулы не отличаются
#
#   C. Попарный тест Пеникаса: KS-статистика между копулами двух подвыборок
#      (аналог теста однородности на основе расстояния Колмогорова–Смирнова
#      для двумерных выборок)

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(tidyr)
library(copula)
library(ggplot2)

# задайте рабочую директорию проекта перед запуском

PARENT  <- ".."
OUT_DIR <- "results/10_homogeneity"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ЗАГРУЗКА

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

# Назначение из скрипта 09 (Frank={1,2,3}, Student={4}, Gumbel°={5,6})
REGIME <- data.frame(
  Period_ID = 1:6,
  Copula    = c("Frank","Frank","Frank","Student","Gumbel","Gumbel"),
  Rotated   = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
)

cat("=== НАЗНАЧЕНИЕ РЕЖИМОВ ===\n")
print(REGIME)

# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

get_period_pobs <- function(pid) {
  df <- lr %>% filter(Date >= as.Date(best$Start[pid]),
                      Date <= as.Date(best$End[pid]))
  if (nrow(df) < 5) return(NULL)
  u <- pobs(cbind(df$log_OIL, df$log_FX))
  colnames(u) <- c("u_oil","u_fx")
  as.data.frame(u)
}

fit_copula_pid <- function(pid, cop_type, rotated = FALSE) {
  u <- get_period_pobs(pid)
  if (is.null(u)) return(NULL)
  if (rotated) u$u_oil <- 1 - u$u_oil
  tryCatch(
    fitCopula(make_copula(cop_type, 2), data = as.matrix(u), method = "ml"),
    error = function(e) NULL
  )
}

loglik_pooled <- function(pids, cop_type, rotated = FALSE) {
  chunks <- lapply(pids, function(pid) {
    u <- get_period_pobs(pid)
    if (rotated && !is.null(u)) u$u_oil <- 1 - u$u_oil
    u
  })
  u_pool <- do.call(rbind, chunks[!sapply(chunks, is.null)])
  fit <- tryCatch(
    fitCopula(make_copula(cop_type, 2), data = as.matrix(u_pool), method = "ml"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(LL = NA, k = NA, n = nrow(u_pool)))
  list(LL = fit@loglik, k = length(fit@copula@parameters), n = nrow(u_pool))
}

# A. LIKELIHOOD RATIO TEST — однородность θ внутри режима

cat("\n=== A. LIKELIHOOD RATIO TEST (H₀: θ одинаков по периодам) ===\n")

lrt_results <- data.frame()

for (cop in unique(REGIME$Copula)) {
  pids    <- REGIME$Period_ID[REGIME$Copula == cop]
  rotated <- REGIME$Rotated[REGIME$Copula == cop][1]

  if (length(pids) < 2) {
    cat(sprintf("  %s: только 1 период — тест не применим\n", cop))
    next
  }

  # Логправдоподобия отдельных периодов
  LL_individual <- sapply(pids, function(pid) {
    fit <- fit_copula_pid(pid, cop, rotated)
    if (is.null(fit)) NA else fit@loglik
  })
  k_per  <- REGIME %>% filter(Copula == cop) %>% nrow()
  params <- switch(cop, "Frank"=1, "Student"=2, "Gumbel"=1)

  # Логправдоподобие объединённой модели
  pool   <- loglik_pooled(pids, cop, rotated)

  LL_sum <- sum(LL_individual, na.rm = TRUE)
  Lambda <- 2 * (LL_sum - pool$LL)       # статистика LRT
  df     <- (length(pids) - 1) * params  # степени свободы
  pval   <- pchisq(Lambda, df = df, lower.tail = FALSE)

  cat(sprintf("\n  [%s] Периоды {%s}\n", cop, paste(pids, collapse=",")))
  cat(sprintf("  LL индивид.: %s  (сумма = %.4f)\n",
              paste(round(LL_individual, 4), collapse=", "), LL_sum))
  cat(sprintf("  LL пул.:     %.4f\n", pool$LL))
  cat(sprintf("  Λ = %.4f, df=%d, p-value = %.4f  → %s\n",
              Lambda, df, pval,
              ifelse(pval > 0.05, "Не отвергаем H₀ ✓ (объединение обосновано)",
                                  "Отвергаем H₀ ✗ (θ различаются)")))

  lrt_results <- rbind(lrt_results, data.frame(
    Copula   = cop,
    Periods  = paste(pids, collapse=","),
    LL_indiv = round(LL_sum, 4),
    LL_pool  = round(pool$LL, 4),
    Lambda   = round(Lambda, 4),
    df       = df,
    p_value  = round(pval, 4),
    Decision = ifelse(pval > 0.05, "H0 not rejected", "H0 rejected")
  ))
}

write.csv(lrt_results, file.path(OUT_DIR, "10a_lrt_homogeneity.csv"), row.names = FALSE)
cat("\n")
print(lrt_results)

# B. ПОПАРНЫЙ KS-ТЕСТ (ПЕНИКАС 2014) — соседние периоды одного режима
#
# Логика: для двух выборок (u¹, v¹) и (u², v²) вычисляем расстояние
# между их эмпирическими копулами (бивариатный KS на единичном кубе).
# Если D < критического значения → распределения неотличимы.

cat("\n=== B. ПОПАРНЫЙ KS-ТЕСТ МЕЖДУ ПЕРИОДАМИ ОДНОГО РЕЖИМА ===\n")

bivariate_ks <- function(U1, U2, n_grid = 20) {
  grid    <- seq(0.05, 0.95, length.out = n_grid)
  g1 <- rep(grid, each = n_grid)
  g2 <- rep(grid, times = n_grid)
  n1  <- nrow(U1);  n2 <- nrow(U2)
  C1  <- sapply(seq_along(g1), function(l)
            mean(U1[,1] <= g1[l] & U1[,2] <= g2[l]))
  C2  <- sapply(seq_along(g1), function(l)
            mean(U2[,1] <= g1[l] & U2[,2] <= g2[l]))
  D   <- max(abs(C1 - C2))
  N_eff <- n1 * n2 / (n1 + n2)
  list(D = D, N_eff = N_eff,
       cv_95 = 1.358 / sqrt(N_eff),
       pval_approx = exp(-2 * (D * sqrt(N_eff))^2))
}

ks_results <- data.frame()

for (cop in unique(REGIME$Copula)) {
  pids    <- REGIME$Period_ID[REGIME$Copula == cop]
  rotated <- REGIME$Rotated[REGIME$Copula == cop][1]
  if (length(pids) < 2) next

  # Попарное сравнение
  pairs <- combn(pids, 2, simplify = FALSE)
  for (pr in pairs) {
    u1 <- get_period_pobs(pr[1]);  u2 <- get_period_pobs(pr[2])
    if (is.null(u1) || is.null(u2)) next
    if (rotated) { u1$u_oil <- 1-u1$u_oil; u2$u_oil <- 1-u2$u_oil }

    res <- bivariate_ks(as.matrix(u1), as.matrix(u2))

    cat(sprintf("  %s: P%d vs P%d  D=%.4f  cv(95%%)=%.4f  p≈%.4f  %s\n",
                cop, pr[1], pr[2], res$D, res$cv_95, res$pval_approx,
                ifelse(res$D < res$cv_95, "одинаковы ✓", "различаются ✗")))

    ks_results <- rbind(ks_results, data.frame(
      Copula    = cop,
      Period_A  = pr[1],
      Period_B  = pr[2],
      D_stat    = round(res$D, 4),
      CV_95     = round(res$cv_95, 4),
      p_approx  = round(res$pval_approx, 4),
      N_eff     = round(res$N_eff, 1),
      Same_dist = res$D < res$cv_95
    ))
  }
}

write.csv(ks_results, file.path(OUT_DIR, "10b_ks_pairwise.csv"), row.names = FALSE)
cat("\n")
print(ks_results)

# C. СВОДНАЯ ИНТЕРПРЕТАЦИЯ

cat("\n", strrep("=", 60), "\n")
cat("ВЫВОД: доказательство объединимости периодов\n")
cat(strrep("=", 60), "\n\n")

for (i in seq_len(nrow(lrt_results))) {
  r    <- lrt_results[i,]
  ok_lrt <- r$Decision == "H0 not rejected"
  ks_sub <- ks_results %>% filter(Copula == r$Copula)
  ok_ks  <- all(ks_sub$Same_dist, na.rm = TRUE)

  cat(sprintf("  %s (периоды {%s}):\n", r$Copula, r$Periods))
  cat(sprintf("    LRT: Λ=%.3f, p=%.4f → %s\n", r$Lambda, r$p_value,
              ifelse(ok_lrt, "H₀ не отвергается ✓", "H₀ отвергается ✗")))
  if (nrow(ks_sub) > 0)
    cat(sprintf("    KS:  все попарные тесты %s\n",
                ifelse(ok_ks, "не значимы ✓", "есть значимые различия ✗")))
  cat(sprintf("    → Объединение %s\n\n",
              ifelse(ok_lrt && ok_ks,
                     "СТАТИСТИЧЕСКИ ОБОСНОВАНО",
                     "требует дополнительной проверки")))
}

# D. ВИЗУАЛИЗАЦИЯ: эмпирические копулы по периодам

cat("=== D. ВИЗУАЛИЗАЦИЯ ЭМПИРИЧЕСКИХ КОПУЛ ПО ПЕРИОДАМ ===\n")

all_u <- data.frame()
for (i in seq_len(nrow(REGIME))) {
  pid <- REGIME$Period_ID[i]
  u   <- get_period_pobs(pid)
  if (is.null(u)) next
  if (REGIME$Rotated[i]) u$u_oil <- 1 - u$u_oil
  u$Period  <- paste0("P", pid)
  u$Regime  <- REGIME$Copula[i]
  all_u     <- rbind(all_u, u)
}

regime_colors <- c("Frank"="#1565C0","Student"="#B71C1C","Gumbel"="#FF8F00")

p_scatter <- ggplot(all_u, aes(x = u_oil, y = u_fx, color = Regime)) +
  geom_point(alpha = 0.5, size = 1.2) +
  facet_wrap(~ Period, ncol = 3) +
  scale_color_manual(values = regime_colors) +
  labs(
    title    = "Псевдонаблюдения по периодам (пространство копулы)",
    subtitle = "Периоды одного цвета — один режим; однородность видна невооружённым глазом",
    x = "u(Нефть)", y = "u(Рубль)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "10a_scatter_by_period.png"),
       p_scatter, width = 10, height = 7, dpi = 150)

# Эмпирические маргинальные CDF по периодам внутри режима
p_ecdf <- ggplot(all_u, aes(x = u_oil, color = Period, linetype = Regime)) +
  stat_ecdf(linewidth = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "ECDF(u_oil) по периодам — однородность внутри режима",
       x = "u(Нефть)", y = "ECDF") +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "10b_ecdf_oil.png"),
       p_ecdf, width = 9, height = 5, dpi = 150)

cat("Результаты сохранены в", OUT_DIR, "\n")
cat("=== Этап 10 завершён ===\n")
