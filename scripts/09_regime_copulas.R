# 09_regime_copulas.R — Выбор 2–3 копул, покрывающих все периоды 
#
# Алгоритм:
#   1. Загружаем кросс-периодную таблицу AIC 
#   2. Вычисляем Kendall's tau для каждого периода (знак зависимости)
#   3. Для каждого типа копулы определяем периоды с AIC < 0
#   4. Жадный алгоритм минимального покрытия -> 3 копулы
#   5. Объединяем назначенные периоды -> переоцениваем параметры
#      Для периодов с tau < 0 используем повёрнутую копулу Гумбеля (90°)
#   6. Финальная сводка и визуализации

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr)
library(tidyr)
library(copula)
library(ggplot2)

# задайте рабочую директорию проекта перед запуском

PARENT  <- ".."
OUT_DIR <- "results/09_regime_copulas"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source("scripts/00_helpers.R")

# ЗАГРУЗКА

cat("=== ЗАГРУЗКА ===\n")
cross <- read.csv(file.path(PARENT, "Copula_Cross_Period_Results.csv"),
                  stringsAsFactors = FALSE)
cross$Start_Date <- as.Date(cross$Start_Date)
cross$End_Date   <- as.Date(cross$End_Date)

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

n_periods <- length(unique(cross$Period_ID))
cat(sprintf("Периодов: %d\n", n_periods))

# 2. KENDALL'S TAU И ЗНАК ЗАВИСИМОСТИ ПО ПЕРИОДАМ

cat("\n=== KENDALL'S TAU ПО ПЕРИОДАМ ===\n")

tau_df <- data.frame()
for (pid in sort(unique(cross$Period_ID))) {
  df <- lr %>% filter(Date >= best$Start[pid], Date <= best$End[pid])
  tau_val <- cor(df$log_OIL, df$log_FX, method = "kendall")
  tau_df <- rbind(tau_df, data.frame(
    Period_ID = pid,
    Start     = as.character(best$Start[pid]),
    End       = as.character(best$End[pid]),
    N         = nrow(df),
    tau       = round(tau_val, 4),
    dep_sign  = ifelse(tau_val > 0, "positive", "negative")
  ))
}
print(tau_df)
write.csv(tau_df, file.path(OUT_DIR, "09_tau_by_period.csv"), row.names = FALSE)

cat("\nНАПОМИНАНИЕ О СОВМЕСТИМОСТИ КОПУЛ С ЗНАКОМ TAU:\n")
cat("  Gumbel  (θ≥1):     tau > 0   верхний хвост\n")
cat("  Clayton (θ≥0):     tau > 0   нижний хвост\n")
cat("  Gumbel° (90°-rot): tau < 0   обращённая зависимость (нефть↑ → рубль↓)\n")
cat("  Frank   (θ∈R):     любой tau симметричная\n")
cat("  Student (ρ∈−1,1):  любой tau тяжёлые хвосты\n")

# 3. ТАБЛИЦА ПОКРЫТИЯ (AIC < 0) — ПО ИСХОДНОЙ КРОСС-ПЕРИОДНОЙ МАТРИЦЕ

cat("\n=== ПОКРЫТИЕ (AIC < 0) ===\n")

coverage <- cross %>%
  group_by(Copula) %>%
  summarise(
    covered_ids = list(sort(Period_ID[AIC < 0])),
    n_covered   = sum(AIC < 0),
    best_aic    = min(AIC),
    .groups = "drop"
  ) %>%
  arrange(desc(n_covered), best_aic)

for (i in seq_len(nrow(coverage))) {
  ids <- paste(unlist(coverage$covered_ids[i]), collapse = ",")
  cat(sprintf("  %-10s: %d/%d периодов {%s}  (лучший AIC=%.4f)\n",
              coverage$Copula[i], coverage$n_covered[i], n_periods, ids,
              coverage$best_aic[i]))
}

# 4. НАЗНАЧЕНИЕ РЕЖИМОВ — РЕШЕНИЕ B: Frank + Student + Gumbel(90°)
#
# Экономическая интерпретация:
#   Frank   {1,2,3}: классическая положительная зависимость нефть–рубль (до 2022)
#   Student {4}    : шок санкций — экстремальные события, отрицательная зависимость
#   Gumbel° {5,6}  : адаптация и новый режим — обращённая зависимость (tau < 0),
#                    моделируем повёрнутой 90° копулой Гумбеля

cat("\n=== НАЗНАЧЕНИЕ РЕЖИМОВ (Frank + Student + Gumbel°) ===\n")

# Для каждого периода — лучшая копула из {Frank, Student, Gumbel}
assign_B <- data.frame()
for (pid in sort(unique(cross$Period_ID))) {
  sub <- cross %>%
    filter(Period_ID == pid, Copula %in% c("Frank", "Student", "Gumbel")) %>%
    arrange(AIC) %>%
    slice_min(AIC, n = 1, with_ties = FALSE)
  assign_B <- rbind(assign_B, data.frame(
    Period_ID  = pid,
    Start      = as.character(best$Start[pid]),
    End        = as.character(best$End[pid]),
    Best_orig  = best$Copula[pid],
    tau        = tau_df$tau[tau_df$Period_ID == pid],
    Regime_B   = sub$Copula,
    AIC_cross  = round(sub$AIC, 4),
    Adequate   = sub$AIC < 0
  ))
}
print(assign_B)
write.csv(assign_B, file.path(OUT_DIR, "09a_assignment_3copulas.csv"), row.names = FALSE)

# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ

get_pooled_pobs <- function(period_ids) {
  chunks <- lapply(period_ids, function(pid) {
    df <- lr %>% filter(Date >= as.Date(best$Start[pid]),
                        Date <= as.Date(best$End[pid]))
    if (nrow(df) < 5) return(NULL)
    u <- pobs(cbind(df$log_OIL, df$log_FX))
    as.data.frame(u)
  })
  chunks <- chunks[!sapply(chunks, is.null)]
  if (length(chunks) == 0) return(NULL)
  u_pool <- do.call(rbind, chunks)
  colnames(u_pool) <- c("u_oil", "u_fx")
  u_pool
}

# Подбор копулы с автоматическим выбором вращения для tau<0
fit_pooled <- function(cop_type, period_ids, label = "") {
  cat(sprintf("\n  [%s] Периоды {%s}:\n", cop_type,
              paste(period_ids, collapse = ",")))
  u_pool <- get_pooled_pobs(period_ids)
  if (is.null(u_pool)) { cat("  Нет данных!\n"); return(NULL) }
  n <- nrow(u_pool)
  tau_pool <- cor(u_pool$u_oil, u_pool$u_fx, method = "kendall")
  cat(sprintf("  N = %d | Kendall's tau = %.4f\n", n, tau_pool))

  # Если Gumbel и tau<0 → используем 90°-повёрнутую Gumbel (rotCopula)
  use_rotated <- FALSE
  if (cop_type == "Gumbel" && tau_pool < 0) {
    cat("  tau < 0 → используем Gumbel° (90°-вращение): обращённая зависимость\n")
    use_rotated <- TRUE
  }

  # Для повёрнутой Gumbel флипаем u_oil → 1-u_oil (90° вращение вручную)
  fit_data <- if (use_rotated) {
    u_flip <- u_pool
    u_flip$u_oil <- 1 - u_pool$u_oil
    as.matrix(u_flip)
  } else {
    as.matrix(u_pool)
  }

  fit <- tryCatch({
    cop_obj <- make_copula(cop_type, dim = 2)
    fitCopula(cop_obj, data = fit_data, method = "ml")
  }, error = function(e) {
    cat(sprintf("  Ошибка ML: %s\n  Пробуем itau...\n", e$message))
    tryCatch(
      fitCopula(make_copula(cop_type, dim = 2), data = fit_data, method = "itau"),
      error = function(e2) { cat("  itau тоже ошибка:", e2$message, "\n"); NULL }
    )
  })

  if (is.null(fit)) return(NULL)

  params <- fit@copula@parameters
  L      <- fit@loglik
  k      <- length(params)

  AIC_val <- count_AIC(L, k)
  BIC_val <- count_BIC(L, k, n)

  param_str <- paste(round(params, 4), collapse = " | ")
  cat(sprintf("  Параметры: %s\n", param_str))
  cat(sprintf("  LL=%.4f  AIC=%.4f  BIC=%.4f  [%s]\n",
              L, AIC_val, BIC_val,
              ifelse(AIC_val < 0, "АДЕКВАТНО", "НЕ АДЕКВАТНО")))
  if (use_rotated) cat("  (Модель: Gumbel°, повёрнут на 90°, обращённая хвостовая зависимость)\n")

  list(copula   = cop_type,
       rotated  = use_rotated,
       pids     = period_ids,
       n        = n,
       tau      = tau_pool,
       params   = params,
       loglik   = L,
       AIC      = AIC_val,
       BIC      = BIC_val)
}

# ПЕРЕОЦЕНКА ПАРАМЕТРОВ НА ОБЪЕДИНЁННЫХ ДАННЫХ

cat("\n=== ПЕРЕОЦЕНКА НА ОБЪЕДИНЁННЫХ ДАННЫХ ===\n")
cat("──────────────────────────────────────────\n")

frank_pids   <- assign_B %>% filter(Regime_B == "Frank")   %>% pull(Period_ID)
student_pids <- assign_B %>% filter(Regime_B == "Student") %>% pull(Period_ID)
gumbel_pids  <- assign_B %>% filter(Regime_B == "Gumbel")  %>% pull(Period_ID)

res_frank   <- fit_pooled("Frank",   frank_pids)
res_student <- fit_pooled("Student", student_pids)
res_gumbel  <- fit_pooled("Gumbel",  gumbel_pids)

# ПАРАМЕТРЫ ПО ОТДЕЛЬНЫМ ПЕРИОДАМ ВНУТРИ РЕЖИМА b

cat("\n=== ПАРАМЕТРЫ ПО ПЕРИОДАМ ВНУТРИ РЕЖИМОВ ===\n")

period_params <- data.frame()
for (pid in sort(unique(cross$Period_ID))) {
  df  <- lr %>% filter(Date >= as.Date(best$Start[pid]),
                       Date <= as.Date(best$End[pid]))
  u   <- pobs(cbind(df$log_OIL, df$log_FX))
  rg  <- assign_B$Regime_B[assign_B$Period_ID == pid]
  tau_p <- tau_df$tau[tau_df$Period_ID == pid]
  use_rot <- (rg == "Gumbel" && tau_p < 0)

  fit_u <- if (use_rot) { u_r <- u; u_r[,1] <- 1 - u[,1]; u_r } else u
  r <- tryCatch(
    fitCopula(make_copula(rg, 2), data = fit_u, method = "ml"),
    error = function(e) NULL
  )

  period_params <- rbind(period_params, data.frame(
    Period_ID = pid,
    Start     = as.character(best$Start[pid]),
    End       = as.character(best$End[pid]),
    Regime    = ifelse(use_rot, "Gumbel°", rg),
    tau_obs   = round(tau_p, 4),
    Param1    = if (!is.null(r)) round(r@copula@parameters[1], 4) else NA,
    Param2    = if (!is.null(r) && length(r@copula@parameters) > 1)
                  round(r@copula@parameters[2], 4) else NA,
    AIC_indiv = if (!is.null(r))
                  round(count_AIC(r@loglik, length(r@copula@parameters)), 4) else NA
  ))
}
cat("\n")
print(period_params)
write.csv(period_params, file.path(OUT_DIR, "09b_period_params.csv"), row.names = FALSE)

# ФИНАЛЬНАЯ СВОДКА

cat("\n", strrep("=", 60), "\n")
cat("★★★  ФИНАЛЬНАЯ СВОДКА: 3 копулы  ★★★\n")
cat(strrep("=", 60), "\n\n")

final_rows <- list()
for (res in list(res_frank, res_student, res_gumbel)) {
  if (is.null(res)) next
  label <- if (res$rotated) "Gumbel° (90°-вращение)" else res$copula
  cat(sprintf("  %s  (периоды {%s}, N=%d)\n",
              label, paste(res$pids, collapse=","), res$n))
  cat(sprintf("    tau_pool = %.4f\n", res$tau))
  cat(sprintf("    θ = %s\n", paste(round(res$params, 4), collapse="; ")))
  cat(sprintf("    LL=%.4f  AIC=%.4f  BIC=%.4f  → %s\n\n",
              res$loglik, res$AIC, res$BIC,
              ifelse(res$AIC < 0, "АДЕКВАТНО ✓", "НЕ АДЕКВАТНО ✗")))

  asgn <- assign_B %>% filter(Regime_B == res$copula)
  final_rows[[length(final_rows)+1]] <- data.frame(
    Copula       = ifelse(res$rotated, "Gumbel°", res$copula),
    Periods      = paste(res$pids, collapse=","),
    Dates        = paste(min(asgn$Start), "–", max(asgn$End)),
    N_pooled     = res$n,
    tau_pool     = round(res$tau, 4),
    Params       = paste(round(res$params, 4), collapse="; "),
    LogLik       = round(res$loglik, 4),
    AIC          = round(res$AIC, 4),
    BIC          = round(res$BIC, 4),
    Adequate     = res$AIC < 0,
    stringsAsFactors = FALSE
  )
}

final_tbl <- do.call(rbind, final_rows)
write.csv(final_tbl, file.path(OUT_DIR, "09c_final_3copulas.csv"), row.names = FALSE)
print(final_tbl[, c("Copula","Periods","Dates","N_pooled","Params","AIC","Adequate")])

# ВИЗУАЛИЗАЦИИ

cat("\n=== ВИЗУАЛИЗАЦИИ ===\n")

# --- 9.1 Тепловая карта AIC ---
all_cop <- c("Frank","Gumbel","Clayton","Student")
hm <- cross %>%
  filter(Copula %in% all_cop) %>%
  left_join(tau_df[, c("Period_ID","tau","dep_sign")], by = "Period_ID") %>%
  mutate(
    adequate = AIC < 0,
    selected = (Copula == "Frank"   & Period_ID %in% frank_pids)   |
               (Copula == "Student" & Period_ID %in% student_pids) |
               (Copula == "Gumbel"  & Period_ID %in% gumbel_pids)
  )

period_lbl <- setNames(
  paste0("P", best$Period_ID, "\n",
         format(best$Start, "%Y-%m"), "\n–\n",
         format(best$End,   "%Y-%m")),
  best$Period_ID)
hm$Period_lbl <- factor(period_lbl[as.character(hm$Period_ID)],
                        levels = period_lbl)

# Метка tau под периодом
tau_label <- setNames(paste0("τ=", tau_df$tau), tau_df$Period_ID)
hm <- hm %>%
  mutate(Period_lbl2 = paste0(Period_lbl, "\nτ=", round(tau, 4)))

p_hm <- ggplot(hm, aes(x = Period_lbl, y = Copula, fill = AIC)) +
  geom_tile(aes(color = adequate), linewidth = 0.8, alpha = 0.95) +
  geom_text(aes(label = sprintf("%.2f", AIC)), size = 3.0, color = "black") +
  geom_tile(data = filter(hm, selected & adequate),
            fill = NA, color = "gold", linewidth = 2) +
  scale_fill_gradient2(low = "#1565C0", mid = "#E3F2FD", high = "#B71C1C",
                       midpoint = 0, name = "AIC") +
  scale_color_manual(values = c("TRUE" = "#1565C0", "FALSE" = "#B71C1C"),
                     name = "AIC<0",
                     guide = guide_legend(override.aes = list(fill = "white"))) +
  labs(
    title    = "Матрица AIC: тип копулы × период",
    subtitle = "Золотая рамка = выбранные (Frank+Student+Gumbel°)\nСинее = AIC<0 (адекватно)",
    x = "Период", y = "Тип копулы"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, lineheight = 0.85))

ggsave(file.path(OUT_DIR, "09a_aic_heatmap.png"), p_hm,
       width = 11, height = 5.5, dpi = 150)

# --- 9.2 Временна́я шкала режимов ---
regime_labels <- c(
  "Frank"   = "Frank (симметричная)",
  "Student" = "Student (тяжёлые хвосты)",
  "Gumbel"  = "Gumbel° (обращённая)"
)
regime_colors <- c("Frank"="#1565C0", "Student"="#B71C1C", "Gumbel"="#FF8F00")

tl <- assign_B %>%
  mutate(Start = as.Date(Start), End = as.Date(End),
         regime_label = regime_labels[Regime_B])

p_tl <- ggplot(tl, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1, fill = Regime_B)) +
  geom_rect(alpha = 0.82) +
  geom_text(
    aes(x     = Start + (End - Start) / 2, y = 0.5,
        label = paste0("P", Period_ID, "\n", Regime_B,
                       "\nτ=", tau, "\nAIC=", AIC_cross)),
    size = 2.7, color = "white", fontface = "bold", lineheight = 0.9
  ) +
  scale_fill_manual(values = regime_colors,
                    labels = regime_labels, name = "Режим") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "6 months") +
  labs(
    title    = "Временна́я шкала режимов зависимости нефть–рубль",
    subtitle = "Frank (синий): классическая пол. зависимость  |  Student (красный): шок санкций  |  Gumbel° (жёлтый): обращённая зависимость",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid  = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(OUT_DIR, "09b_regime_timeline.png"), p_tl,
       width = 12, height = 4.5, dpi = 150)

# --- 9.3 Параметр θ₁ по периодам ---
regime_colors_pp <- c("Frank"="#1565C0","Student"="#B71C1C","Gumbel°"="#FF8F00","Gumbel"="#FF8F00")

p_par <- ggplot(period_params %>% filter(!is.na(Param1)),
                aes(x = as.Date(Start), y = Param1, color = Regime, shape = Regime)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = regime_colors_pp) +
  scale_shape_manual(values = c("Frank"=16,"Student"=17,"Gumbel°"=15,"Gumbel"=15)) +
  labs(
    title    = "Параметр θ₁ по периодам (Решение B: Frank + Student + Gumbel°)",
    subtitle = "Frank: θ>0 = положительная зависимость | Student: ρ<0 = отрицательная | Gumbel°: θ>1 = обращённый хвост",
    x = NULL, y = "θ₁ копулы", color = "Режим", shape = "Режим"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "09c_params_by_period.png"), p_par,
       width = 10, height = 5, dpi = 150)

cat("Графики сохранены в", OUT_DIR, "\n")
cat("=== Этап 09 завершён ===\n")
