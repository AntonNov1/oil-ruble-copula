# 11_gof_pooled.R — GoF-тесты для объединённых моделей (3 копулы)
#
# Цель: проверить, что каждая из трёх финальных копул адекватно описывает
#       объединённые данные всех «своих» периодов.
#
# Тесты:
#   1. Cramér-von Mises (CvM) — основной GoF для копул
#      (реализован в пакете copula: gofCopula)
#   2. Kolmogorov-Smirnov на маргинальных псевдонаблюдениях
#   3. Визуальная диагностика: QQ-plot копулы, scatter

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(copula)
library(ggplot2)
library(tidyr)

# задайте рабочую директорию проекта перед запуском

PARENT  <- ".."
OUT_DIR <- "results/11_gof_pooled"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

N_BOOT <- 199   # число бутстрэп-репликаций (199 рекомендовано для GoF копул)

# ПАРАМЕТРЫ ФИНАЛЬНЫХ МОДЕЛЕЙ (из скрипта 09)

MODELS <- list(
  list(
    name      = "Frank",
    copula    = frankCopula(param = 2.5066, dim = 2),
    pids      = c(1, 2, 3),
    rotated   = FALSE,
    label     = "Frank (θ=2.507, периоды 1–3)"
  ),
  list(
    name      = "Student",
    copula    = tCopula(param = -0.6347, df = 0.5551, dim = 2),
    pids      = c(4),
    rotated   = FALSE,
    label     = "Student (ρ=−0.635, df=0.555, период 4)"
  ),
  list(
    name      = "Gumbel°",
    copula    = gumbelCopula(param = 1.0869, dim = 2),
    pids      = c(5, 6),
    rotated   = TRUE,                    # 90°-поворот: u_oil → 1-u_oil
    label     = "Gumbel° (θ=1.087, периоды 5–6)"
  )
)

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

# ФУНКЦИЯ: собрать объединённые псевдонаблюдения для режима

get_pooled_u <- function(pids, rotated = FALSE) {
  chunks <- lapply(pids, function(pid) {
    df <- lr %>% filter(Date >= as.Date(best$Start[pid]),
                        Date <= as.Date(best$End[pid]))
    if (nrow(df) < 5) return(NULL)
    u <- pobs(cbind(df$log_OIL, df$log_FX))
    as.data.frame(u)
  })
  u_pool <- do.call(rbind, chunks[!sapply(chunks, is.null)])
  colnames(u_pool) <- c("u_oil","u_fx")
  if (rotated) u_pool$u_oil <- 1 - u_pool$u_oil
  u_pool
}

# GoF ТЕСТЫ

gof_results <- data.frame()

for (m in MODELS) {
  cat(sprintf("\n══════════════════════════════════════════════════\n"))
  cat(sprintf("  Модель: %s\n", m$label))
  cat(sprintf("══════════════════════════════════════════════════\n"))

  u_pool <- get_pooled_u(m$pids, m$rotated)
  n      <- nrow(u_pool)
  u_mat  <- as.matrix(u_pool)

  cat(sprintf("  N = %d наблюдений\n", n))

  # ──────────────────────────────────────────────────────
  # 1. Cramér-von Mises (bootstrap)
  # ──────────────────────────────────────────────────────
  cat(sprintf("  [1] Cramér-von Mises GoF (B=%d bootstrap)...\n", N_BOOT))

  cvm_stat <- NA; cvm_pval <- NA

  if (m$name == "Student") {
    # t-копула с нецелым df: pCopula недоступна → bootstrap log-likelihood GoF
    cat("  [Student] CvM недоступен для нецелого df → бутстрэп LogLik-тест\n")
    LL_obs <- tryCatch({
      fit_s <- fitCopula(m$copula, u_mat, method = "ml")
      fit_s@loglik
    }, error = function(e) NA)
    boot_ll <- replicate(N_BOOT, tryCatch({
      u_b  <- rCopula(n, m$copula)
      fit_b <- fitCopula(m$copula, pobs(u_b), method = "ml")
      fit_b@loglik
    }, error = function(e) NA))
    boot_ll <- boot_ll[!is.na(boot_ll)]
    if (!is.na(LL_obs) && length(boot_ll) > 10) {
      cvm_pval <- mean(boot_ll <= LL_obs)   # p = P(LL_boot ≤ LL_obs)
      cat(sprintf("  Bootstrap LL: LL_obs=%.4f  p(LL_boot≤LL_obs)=%.4f  → %s\n",
                  LL_obs, cvm_pval,
                  ifelse(cvm_pval > 0.05, "Не отвергаем H₀ ✓", "Отвергаем H₀ ✗")))
    }
  } else {
    gof_res <- tryCatch(
      gofCopula(m$copula, x = u_mat, N = N_BOOT,
                simulation  = "pb",
                method       = "Sn",
                estim.method = "mpl"),
      error = function(e) { cat("  Ошибка CvM:", e$message, "\n"); NULL }
    )
    if (!is.null(gof_res)) {
      cvm_stat <- round(gof_res$statistic, 5)
      cvm_pval <- round(gof_res$p.value,   4)
      cat(sprintf("  CvM Sn = %.5f,  p-value = %.4f  → %s\n",
                  cvm_stat, cvm_pval,
                  ifelse(cvm_pval > 0.05, "Не отвергаем H₀ ✓", "Отвергаем H₀ ✗")))
    }
  }

  # ──────────────────────────────────────────────────────
  # 2. KS-тест на маргиналях (U[0,1]?)
  # ──────────────────────────────────────────────────────
  ks_oil <- ks.test(u_pool$u_oil, "punif")
  ks_fx  <- ks.test(u_pool$u_fx,  "punif")
  cat(sprintf("  [2] KS marginal: u_oil p=%.4f | u_fx p=%.4f\n",
              ks_oil$p.value, ks_fx$p.value))

  # ──────────────────────────────────────────────────────
  # 3. Логправдоподобие и информационные критерии
  # ──────────────────────────────────────────────────────
  dens <- tryCatch(
    dCopula(u_mat, m$copula),
    error = function(e) rep(NA_real_, n)
  )
  dens[is.na(dens) | dens <= 0] <- .Machine$double.eps
  LL  <- sum(log(dens))
  k   <- length(m$copula@parameters)
  AIC <- count_AIC(LL, k)
  BIC <- count_BIC(LL, k, n)
  cat(sprintf("  [3] LL=%.4f  AIC=%.4f  BIC=%.4f\n", LL, AIC, BIC))

  gof_results <- rbind(gof_results, data.frame(
    Copula    = m$name,
    Periods   = paste(m$pids, collapse=","),
    N         = n,
    CvM_Sn    = cvm_stat,
    CvM_pval  = cvm_pval,
    KS_oil_p  = round(ks_oil$p.value, 4),
    KS_fx_p   = round(ks_fx$p.value,  4),
    LogLik    = round(LL,  4),
    AIC       = round(AIC, 4),
    BIC       = round(BIC, 4),
    Adequate  = ifelse(!is.na(cvm_pval), cvm_pval > 0.05, AIC < 0)
  ))

  # ──────────────────────────────────────────────────────
  # 4. QQ-plot: теоретическая vs эмпирическая копула
  # ──────────────────────────────────────────────────────
  theo <- tryCatch(
    pCopula(u_mat, m$copula),
    error = function(e) rep(NA_real_, n)
  )
  emp  <- apply(u_mat, 1, function(z)
    mean(u_mat[,1] <= z[1] & u_mat[,2] <= z[2]))

  qq_df <- data.frame(theoretical = theo, empirical = emp)
  qq_df <- qq_df[!is.na(theo),]

  fname <- paste0("11_qq_", gsub("[°]","deg", m$name), ".png")
  p_qq <- ggplot(qq_df, aes(x = theoretical, y = empirical)) +
    geom_point(alpha = 0.4, size = 1.5, color = "#1565C0") +
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1) +
    labs(
      title    = paste("QQ-plot копулы:", m$label),
      subtitle = sprintf("N=%d  |  AIC=%.3f  |  CvM p=%.4f",
                         n, AIC, ifelse(is.na(cvm_pval), NA, cvm_pval)),
      x = "Теоретическая C(u,v)", y = "Эмпирическая Ĉ(u,v)"
    ) +
    theme_minimal()

  ggsave(file.path(OUT_DIR, fname), p_qq, width = 6, height = 6, dpi = 150)
}

# СВОДНАЯ ТАБЛИЦА

cat("\n", strrep("=", 60), "\n")
cat("СВОДНАЯ ТАБЛИЦА GoF\n")
cat(strrep("=", 60), "\n\n")
print(gof_results)
write.csv(gof_results, file.path(OUT_DIR, "11_gof_summary.csv"), row.names = FALSE)

# ТЕСТ ПО ОТДЕЛЬНЫМ ПЕРИОДАМ ВНУТРИ КАЖДОГО РЕЖИМА

cat("\n=== GoF ПО ОТДЕЛЬНЫМ ПЕРИОДАМ ===\n")

period_gof <- data.frame()

for (m in MODELS) {
  for (pid in m$pids) {
    df_p <- lr %>% filter(Date >= as.Date(best$Start[pid]),
                          Date <= as.Date(best$End[pid]))
    if (nrow(df_p) < 10) next
    u_p <- as.data.frame(pobs(cbind(df_p$log_OIL, df_p$log_FX)))
    colnames(u_p) <- c("u_oil","u_fx")
    if (m$rotated) u_p$u_oil <- 1 - u_p$u_oil
    u_mat_p <- as.matrix(u_p)
    n_p     <- nrow(u_mat_p)

    dens_p <- tryCatch(
      dCopula(u_mat_p, m$copula),
      error = function(e) rep(NA_real_, n_p)
    )
    dens_p[is.na(dens_p) | dens_p <= 0] <- .Machine$double.eps
    LL_p  <- sum(log(dens_p))
    k_p   <- length(m$copula@parameters)
    AIC_p <- count_AIC(LL_p, k_p)

    cat(sprintf("  %s на P%d (N=%d): LL=%.4f  AIC=%.4f  %s\n",
                m$name, pid, n_p, LL_p, AIC_p,
                ifelse(AIC_p < 0, "✓", "✗")))

    period_gof <- rbind(period_gof, data.frame(
      Copula  = m$name,
      Period  = pid,
      N       = n_p,
      LogLik  = round(LL_p,  4),
      AIC     = round(AIC_p, 4),
      OK      = AIC_p < 0
    ))
  }
}

write.csv(period_gof, file.path(OUT_DIR, "11b_gof_by_period.csv"), row.names = FALSE)

# Диаграмма AIC по периодам
p_aic <- ggplot(period_gof, aes(x = factor(Period), y = AIC,
                                 fill = Copula, alpha = OK)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Frank"="#1565C0","Student"="#B71C1C","Gumbel°"="#FF8F00")) +
  scale_alpha_manual(values = c("TRUE"=0.9,"FALSE"=0.4), guide = "none") +
  labs(
    title    = "AIC для каждой копулы на её периодах",
    subtitle = "Красная пунктирная линия = 0 (порог адекватности)",
    x = "Период", y = "AIC (< 0 = адекватно)"
  ) +
  theme_minimal() + theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "11c_aic_by_period.png"), p_aic,
       width = 8, height = 5, dpi = 150)

cat("\nРезультаты сохранены в", OUT_DIR, "\n")
cat("=== Этап 11 завершён ===\n")
