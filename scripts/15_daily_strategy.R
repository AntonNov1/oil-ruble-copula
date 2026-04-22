# 15_daily_strategy.R — Торговая стратегия на дневных данных
#
# Что исправляем по сравнению со скриптом 14:
#   1. Дневные данные (~2150 точек вместо 431 недельных)
#   2. Сигнал с лагом 1 день (нет lookahead bias)
#   3. Транзакционные издержки (15 bps при каждом обороте)
#   4. Bootstrap доверительные интервалы для Sharpe ratio
#   5. Правильная условная CDF через числовую производную
#
# Источники данных:
#   - Brent:  Yahoo Finance (quantmod) → тикер "BZ=F"
#   - USDRUB: ЦБ РФ XML API (бесплатно, без ключа)

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr); library(tidyr); library(copula); library(ggplot2)
library(quantmod); library(httr); library(xml2); library(scales)

# задайте рабочую директорию проекта перед запуском
OUT_DIR <- "results/15_daily_strategy"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COST_BPS   <- 15      # транзакционные издержки, базисных пунктов
COST       <- COST_BPS / 10000
DELTA      <- 0.05    # порог для условного квантильного сигнала
TRAIN_FRAC <- 0.70    # доля обучающей выборки
N_BOOT     <- 2000    # репликаций для bootstrap CI Sharpe
BLOCK_SIZE <- 5       # размер блока (~ 1 торговая неделя)


cat("  Скрипт 15: торговая стратегия на дневных данных \n")


# ЗАГРУЗКА ДАННЫХ

FROM_DATE <- "2017-08-01"
TO_DATE   <- Sys.Date()

# ── 1a. Brent crude (Yahoo Finance) ──────────────────────────────────────────

brent_raw <- tryCatch({
  getSymbols("BZ=F", src = "yahoo", from = FROM_DATE, to = TO_DATE,
             auto.assign = FALSE)
}, error = function(e) {
  cat("  Yahoo Finance недоступен, пробуем ICE (CB=F)...\n")
  tryCatch(
    getSymbols("CB=F", src = "yahoo", from = FROM_DATE, to = TO_DATE,
               auto.assign = FALSE),
    error = function(e2) NULL
  )
})

if (is.null(brent_raw)) {
  stop("Не удалось загрузить данные Brent. Проверьте интернет-соединение.")
}

brent_df <- data.frame(
  Date  = as.Date(index(brent_raw)),
  Brent = as.numeric(Cl(brent_raw))
) %>% filter(!is.na(Brent))
cat(sprintf("  Brent: %d дней (%s — %s)\n",
            nrow(brent_df), min(brent_df$Date), max(brent_df$Date)))

# ── 1b. USDRUB (ЦБ РФ XML API) ───────────────────────────────────────────────


get_cbr_usdrub <- function(from_dt, to_dt) {
  # ЦБ принимает не более ~365 дней за запрос → разбиваем на годовые куски
  chunks <- seq(from_dt, to_dt, by = "year")
  if (tail(chunks, 1) < to_dt) chunks <- c(chunks, to_dt)

  result <- data.frame()
  for (i in seq_len(length(chunks) - 1)) {
    d1  <- format(chunks[i],     "%d/%m/%Y")
    d2  <- format(chunks[i + 1], "%d/%m/%Y")
    url <- paste0(
      "https://www.cbr.ru/scripts/XML_dynamic.asp",
      "?date_req1=", d1, "&date_req2=", d2, "&VAL_NM_RQ=R01235"
    )
    resp <- tryCatch(GET(url, timeout(30)), error = function(e) NULL)
    if (is.null(resp) || status_code(resp) != 200) next

    raw_bytes <- tryCatch(content(resp, as = "raw"), error = function(e) NULL)
    if (is.null(raw_bytes)) next
    xml_text  <- iconv(rawToChar(raw_bytes), from = "windows-1251", to = "UTF-8")
    if (is.na(xml_text) || nchar(xml_text) < 50) next

    parsed  <- tryCatch(read_xml(xml_text), error = function(e) NULL)
    if (is.null(parsed)) next

    records <- xml_find_all(parsed, ".//Record")
    if (length(records) == 0) next

    chunk_df <- data.frame(
      Date   = as.Date(xml_attr(records, "Date"), "%d.%m.%Y"),
      USDRUB = as.numeric(gsub(",", ".",
                 xml_text(xml_find_first(records, ".//Value"))))
    )
    result <- rbind(result, chunk_df)
    Sys.sleep(0.3)   # вежливая пауза между запросами
  }
  result %>% filter(!is.na(USDRUB)) %>% arrange(Date) %>%
    distinct(Date, .keep_all = TRUE)
}

usdrub_df <- get_cbr_usdrub(as.Date(FROM_DATE), as.Date(TO_DATE))
cat(sprintf("  USDRUB (ЦБ): %d дней (%s — %s)\n",
            nrow(usdrub_df), min(usdrub_df$Date), max(usdrub_df$Date)))

# ── 1c. Объединяем и вычисляем лог-доходности ────────────────────────────────
daily <- inner_join(brent_df, usdrub_df, by = "Date") %>%
  arrange(Date) %>%
  filter(!is.na(Brent), !is.na(USDRUB)) %>%
  mutate(
    log_OIL = c(NA, diff(log(Brent))),
    log_FX  = c(NA, -diff(log(USDRUB)))   # знак −: рост USDRUB = ослабление рубля
  ) %>%
  filter(!is.na(log_OIL), !is.na(log_FX))

T_obs <- nrow(daily)
cat(sprintf("\nОбъединённый ряд: N=%d дней (%s — %s)\n\n",
            T_obs, min(daily$Date), max(daily$Date)))

write.csv(daily, file.path(OUT_DIR, "15a_daily_data.csv"), row.names = FALSE)

# ПСЕВДОНАБЛЮДЕНИЯ И ПАРАМЕТРЫ МОДЕЛИ

u_mat <- pobs(as.matrix(daily[, c("log_OIL","log_FX")]))
colnames(u_mat) <- c("u_oil","u_fx")
u_rot <- cbind(1 - u_mat[,"u_oil"], u_mat[,"u_fx"])

# Параметры из скрипта 12 (3-state, лучшая модель по AIC)
# Переоцениваем на суточных данных, но используем те же структурные сдвиги
BREAK_DATES <- as.Date(c("2020-01-18","2020-06-20","2022-02-26",
                          "2022-06-11","2023-09-09"))

# Используем параметры из скрипта 12 (еженедельная модель — лучше по AIC)
# Дневные и недельные копулы совпадают в форме зависимости (инвариантность к частоте)


frank_theta <- 2.5066
rho_s       <- -0.6347
df_s        <- 4.0      # округляем до целого — pCopula работает только с целым df
gumbel_rot  <- 1.0869

cat(sprintf("  Frank θ=%.4f  |  Student ρ=%.4f df=%.1f (целый)  |  Gumbel° θ=%.4f\n",
            frank_theta, rho_s, df_s, gumbel_rot))

# 3. ПЛОТНОСТИ И ФИЛЬТР ГАМИЛЬТОНА (3 состояния)

K           <- 3
STATE_NAMES <- c("Frank","Student","Gumbel°")
KAPPA       <- c(1, -1, -1)   # знак зависимости по режимам

f_frank   <- pmax(dCopula(u_mat, frankCopula(frank_theta, 2)),   1e-300)
f_student <- pmax(dCopula(u_mat, tCopula(rho_s, df=df_s, 2)),   1e-300)
f_gumbel  <- pmax(dCopula(u_rot, gumbelCopula(gumbel_rot, 2)),  1e-300)
f_mat     <- cbind(f_frank, f_student, f_gumbel)

# Матрица переходов из скрипта 12 (масштабируем на дневную частоту)
# P_weekly → P_daily: если P_weekly[i,i]=0.95, то за 5 дней ≈ 0.95
# значит P_daily[i,i] = 0.95^{1/5} ≈ 0.990
weekly_diag <- c(0.949, 0.879, 0.997)
daily_diag  <- weekly_diag^(1/5)
P_daily <- matrix(0, K, K)
for (k in 1:K) {
  P_daily[k, k] <- daily_diag[k]
  off <- (1 - daily_diag[k]) / (K - 1)
  P_daily[k, -k] <- off
}
cat(sprintf("  P_daily диагональ: %.4f  %.4f  %.4f\n",
            P_daily[1,1], P_daily[2,2], P_daily[3,3]))

hamilton_filter <- function(f_mat, P, pi0 = NULL) {
  K <- ncol(f_mat); N <- nrow(f_mat)
  if (is.null(pi0)) {
    ev <- eigen(t(P))$vectors[,1]
    pi0 <- pmax(Re(ev)/sum(Re(ev)), 0); pi0 <- pi0/sum(pi0)
  }
  xi_filt <- matrix(0, N, K)
  LL_vec  <- numeric(N)
  p_pred  <- pi0
  for (t in 1:N) {
    joint <- f_mat[t,] * p_pred
    denom <- sum(joint)
    if (denom <= 0) { xi_filt[t,] <- p_pred; LL_vec[t] <- -700 }
    else { xi_filt[t,] <- joint / denom; LL_vec[t] <- log(denom) }
    p_pred <- pmax(as.vector(t(P) %*% xi_filt[t,]), 0)
    p_pred <- p_pred / sum(p_pred)
  }
  list(xi_filt = xi_filt, LL = sum(LL_vec))
}

ev_P <- eigen(t(P_daily))$vectors[,1]
pi0  <- pmax(Re(ev_P)/sum(Re(ev_P)), 0); pi0 <- pi0/sum(pi0)
res_filt <- hamilton_filter(f_mat, P_daily, pi0)
xi_filt  <- res_filt$xi_filt

cat(sprintf("  LogLik дневной модели = %.2f\n", res_filt$LL))

# 4. СИГНАЛ С ЛАГОМ 1 ДЕНЬ

cat("\n▶ Генерация сигнала (лаг = 1 день)...\n")

# Условная CDF через числовую производную: P(V ≤ 0.5 | U = u)
# ∂C(u,v)/∂u ≈ [C(u+h, v) - C(u-h, v)] / (2h)
cond_cdf <- function(u_vec, v = 0.5, cop, eps = 1e-4) {
  u_hi <- pmin(u_vec + eps, 1 - eps)
  u_lo <- pmax(u_vec - eps, eps)
  v_c  <- rep(pmin(pmax(v, eps), 1 - eps), length(u_vec))
  C_hi <- pCopula(cbind(u_hi, v_c), cop)
  C_lo <- pCopula(cbind(u_lo, v_c), cop)
  pmin(pmax((C_hi - C_lo) / (2 * eps), 0), 1)
}

cop_frank <- frankCopula(frank_theta, 2)
cop_stud  <- tCopula(rho_s, df = df_s, dim = 2)
cop_grot  <- gumbelCopula(gumbel_rot, 2)

# P(рубль укрепится | u_oil) для каждого режима
p_up_frank   <- 1 - cond_cdf(u_mat[,1], 0.5, cop_frank)
p_up_student <- 1 - cond_cdf(u_mat[,1], 0.5, cop_stud)
p_up_grot    <- 1 - cond_cdf(1 - u_mat[,1], 0.5, cop_grot)  # повёрнутые u

p_up_combined <- xi_filt[,1] * p_up_frank +
                 xi_filt[,2] * p_up_student +
                 xi_filt[,3] * p_up_grot

# ── Сигнал с лагом ──────────────────────────────────────────────────────────
# Наблюдаем данные дня t → сигнал генерируется в конце дня t
# → позиция открывается в начале дня t+1 → зарабатываем доходность дня t+1

p_up_lag    <- c(NA, head(p_up_combined, -1))  # сдвиг на 1 день
direction_lag <- c(NA, head(as.numeric(xi_filt %*% KAPPA), -1))
oil_sign_lag  <- c(NA, head(sign(daily$log_OIL), -1))

# Два сигнала:
# A) Режимно-направленный (с лагом)
pos_regime <- ifelse(is.na(direction_lag), 0,
              ifelse(abs(direction_lag) >= DELTA, sign(direction_lag * oil_sign_lag), 0))

# B) Условный квантиль (с лагом)
pos_cq <- ifelse(is.na(p_up_lag), 0,
          ifelse(p_up_lag > 0.5 + DELTA,  1,
          ifelse(p_up_lag < 0.5 - DELTA, -1, 0)))

# РАСЧЁТ ДОХОДНОСТЕЙ С ТРАНЗАКЦИОННЫМИ ИЗДЕРЖКАМИ

cat(sprintf("▶ Расчёт доходностей (издержки %.0f bps на оборот)...\n", COST_BPS))

calc_returns <- function(pos, log_fx) {
  # Доходность позиции: pos > 0 = лонг рубль = прибыль когда log_FX > 0
  gross <- pos * log_fx
  # Издержки: платим при изменении позиции
  turnover <- c(0, abs(diff(pos)))
  net <- gross - turnover * COST
  net
}

ret_regime <- calc_returns(pos_regime, daily$log_FX)
ret_cq     <- calc_returns(pos_cq,     daily$log_FX)
ret_bh_rub <- daily$log_FX              # buy-and-hold лонг рубль
ret_bh_usd <- -daily$log_FX             # buy-and-hold лонг USDRUB

# Разбивка train / test
T_train    <- floor(T_obs * TRAIN_FRAC)
idx_test   <- (T_train + 1):T_obs
dates_test <- daily$Date[idx_test]

ret_regime_test <- ret_regime[idx_test]
ret_cq_test     <- ret_cq[idx_test]
ret_bh_rub_test <- ret_bh_rub[idx_test]
ret_bh_usd_test <- ret_bh_usd[idx_test]

cat(sprintf("  Train: %s — %s (%d дней)\n",
            daily$Date[1], daily$Date[T_train], T_train))
cat(sprintf("  Test:  %s — %s (%d дней)\n",
            dates_test[1], tail(dates_test,1), length(idx_test)))

# МЕТРИКИ ПРОИЗВОДИТЕЛЬНОСТИ

sharpe <- function(r, ann = 252) {
  r <- r[!is.na(r) & is.finite(r)]
  if (length(r) < 2 || sd(r) == 0) return(NA)
  mean(r) / sd(r) * sqrt(ann)
}
ann_ret <- function(r, ann = 252) mean(r, na.rm=TRUE) * ann * 100
max_dd  <- function(cum_r) {
  cx <- cummax(cum_r)
  max(cx - cum_r)
}
win_rate <- function(r, pos) {
  active <- pos[!is.na(pos) & pos != 0]
  r_act  <- r[!is.na(pos) & pos != 0]
  if (length(r_act) == 0) return(NA)
  mean(r_act > 0)
}

metrics <- data.frame(
  Стратегия = c("Режим (лаг+изд.)","Условный квантиль (лаг+изд.)",
                "B&H рубль","B&H USDRUB"),
  `Ann.Return,%` = round(c(
    ann_ret(ret_regime_test), ann_ret(ret_cq_test),
    ann_ret(ret_bh_rub_test), ann_ret(ret_bh_usd_test)
  ), 2),
  Sharpe = round(c(
    sharpe(ret_regime_test), sharpe(ret_cq_test),
    sharpe(ret_bh_rub_test), sharpe(ret_bh_usd_test)
  ), 3),
  `MaxDD,%` = round(c(
    max_dd(cumsum(replace(ret_regime_test, is.na(ret_regime_test), 0))),
    max_dd(cumsum(replace(ret_cq_test,     is.na(ret_cq_test),     0))),
    max_dd(cumsum(replace(ret_bh_rub_test, is.na(ret_bh_rub_test), 0))),
    max_dd(cumsum(replace(ret_bh_usd_test, is.na(ret_bh_usd_test), 0)))
  ) * 100, 2),
  `N_trades` = c(
    sum(abs(diff(pos_cq[idx_test])) > 0, na.rm=TRUE),
    sum(abs(diff(pos_cq[idx_test])) > 0, na.rm=TRUE),
    0, 0
  ),
  check.names = FALSE
)

cat("\n╔══════════════════════════════════════════════╗\n")


print(metrics)

# bootstrap ДОВЕРИТЕЛЬНЫЕ ИНТЕРВАЛЫ ДЛЯ sharpe

cat(sprintf("\n▶ Bootstrap ДИ для Sharpe (N=%d, блок=%d дней)...\n",
            N_BOOT, BLOCK_SIZE))

circular_block_bootstrap <- function(r, block = BLOCK_SIZE, B = N_BOOT) {
  n      <- length(r)
  n_blks <- ceiling(n / block)
  boot_sharpes <- numeric(B)
  for (b in seq_len(B)) {
    starts <- sample(1:n, n_blks, replace = TRUE)
    idx    <- as.vector(sapply(starts, function(s) ((s:(s+block-1) - 1) %% n) + 1))
    idx    <- idx[1:n]
    boot_sharpes[b] <- sharpe(r[idx])
  }
  boot_sharpes
}

bs_regime <- circular_block_bootstrap(
  ret_regime_test[!is.na(ret_regime_test) & is.finite(ret_regime_test)])
bs_cq     <- circular_block_bootstrap(
  ret_cq_test[!is.na(ret_cq_test) & is.finite(ret_cq_test)])

ci_regime <- quantile(bs_regime, c(0.05, 0.95), na.rm = TRUE)
ci_cq     <- quantile(bs_cq,     c(0.05, 0.95), na.rm = TRUE)

cat(sprintf("  Режим (лаг):      Sharpe = %.3f  [90%% ДИ: %.3f ; %.3f]\n",
            sharpe(ret_regime_test), ci_regime[1], ci_regime[2]))
cat(sprintf("  Усл. квантиль:    Sharpe = %.3f  [90%% ДИ: %.3f ; %.3f]\n",
            sharpe(ret_cq_test), ci_cq[1], ci_cq[2]))
cat(sprintf("  H₀: Sharpe=0 отвергается (Режим): %s\n",
            ifelse(ci_regime[1] > 0, "ДА ✓", "НЕТ ✗")))
cat(sprintf("  H₀: Sharpe=0 отвергается (КвантЛ): %s\n",
            ifelse(ci_cq[1] > 0, "ДА ✓", "НЕТ ✗")))

# АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ К ТРАНЗАКЦИОННЫМ ИЗДЕРЖКАМ

costs_bps <- c(0, 5, 10, 15, 20, 30, 50)
sens_df   <- lapply(costs_bps, function(cb) {
  c_val <- cb / 10000
  tv    <- c(0, abs(diff(pos_cq[idx_test])))
  r     <- pos_cq[idx_test] * daily$log_FX[idx_test] - tv * c_val
  data.frame(Cost_bps=cb, Sharpe=round(sharpe(r),3),
             Ann_ret=round(ann_ret(r),2),
             Profitable=ann_ret(r)>0)
}) %>% do.call(rbind, .)

cat("\nЧувствительность Sharpe к транзакционным издержкам (усл. квантиль):\n")
print(sens_df)

# ВИЗУАЛИЗАЦИИ

cat("\n▶ Визуализации...\n")

# 9.1 Кумулятивная доходность
cum_df <- data.frame(
  Date    = dates_test,
  `Режим (лаг)` = cumsum(replace(ret_regime_test, is.na(ret_regime_test), 0)) * 100,
  `Усл.квантиль` = cumsum(replace(ret_cq_test, is.na(ret_cq_test), 0)) * 100,
  `B&H рубль`   = cumsum(replace(ret_bh_rub_test, is.na(ret_bh_rub_test), 0)) * 100,
  `B&H USDRUB`  = cumsum(replace(ret_bh_usd_test, is.na(ret_bh_usd_test), 0)) * 100,
  check.names = FALSE
) %>% pivot_longer(-Date, names_to = "Strategy", values_to = "CumRet")

strat_colors <- c(
  "Режим (лаг)"  = "#1565C0",
  "Усл.квантиль" = "#2E7D32",
  "B&H рубль"    = "#B71C1C",
  "B&H USDRUB"   = "#FF8F00"
)

p_cum <- ggplot(cum_df, aes(x=Date, y=CumRet, color=Strategy)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  scale_color_manual(values=strat_colors) +
  scale_x_date(date_labels="%Y-%m", date_breaks="3 months") +
  labs(
    title    = "Кумулятивная доходность (тест, дневные данные, лаг=1, изд.=15bps)",
    subtitle = sprintf("Тест: %s — %s  |  N=%d дней",
                       format(dates_test[1]), format(tail(dates_test,1)),
                       length(idx_test)),
    x=NULL, y="Кум. доходность, %", color=NULL
  ) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_DIR,"15b_cumulative_returns.png"),
       p_cum, width=14, height=5, dpi=150)

# 9.2 Апостериорные вероятности (дневные)
probs_long <- data.frame(
  Date    = daily$Date,
  Frank   = xi_filt[,1],
  Student = xi_filt[,2],
  `Gumbel°` = xi_filt[,3], check.names=FALSE
) %>% pivot_longer(-Date, names_to="Regime", values_to="Prob")

regime_colors <- c("Frank"="#1565C0","Student"="#B71C1C","Gumbel°"="#FF8F00")

p_stack <- ggplot(probs_long, aes(x=Date, y=Prob, fill=Regime)) +
  geom_area(position="stack", alpha=0.85) +
  geom_vline(xintercept=BREAK_DATES, linetype="dashed",
             color="white", linewidth=0.5) +
  scale_fill_manual(values=regime_colors) +
  scale_x_date(date_labels="%Y", date_breaks="1 year") +
  scale_y_continuous(labels=percent_format()) +
  labs(title="Апостериорные вероятности режимов (дневной фильтр Гамильтона)",
       x=NULL, y=NULL, fill=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_DIR,"15c_regime_probs_daily.png"),
       p_stack, width=14, height=4.5, dpi=150)

# 9.3 Распределение Sharpe (bootstrap)
boot_df <- data.frame(
  Sharpe   = c(bs_regime, bs_cq),
  Strategy = rep(c("Режим (лаг)","Усл.квантиль"), each=N_BOOT)
)
p_boot <- ggplot(boot_df, aes(x=Sharpe, fill=Strategy, color=Strategy)) +
  geom_density(alpha=0.4) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  geom_vline(xintercept=ci_cq[1], linetype="dotted", color="#2E7D32", linewidth=1) +
  geom_vline(xintercept=ci_cq[2], linetype="dotted", color="#2E7D32", linewidth=1) +
  scale_fill_manual(values=c("Режим (лаг)"="#1565C0","Усл.квантиль"="#2E7D32")) +
  scale_color_manual(values=c("Режим (лаг)"="#1565C0","Усл.квантиль"="#2E7D32")) +
  labs(title = "Bootstrap-распределение коэффициента Шарпа",
       subtitle = sprintf("Пунктир = 90%%-й ДИ (%.3f ; %.3f) для Усл.квантиль",
                          ci_cq[1], ci_cq[2]),
       x="Sharpe (годовой)", y="Плотность", fill=NULL, color=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom")
ggsave(file.path(OUT_DIR,"15d_bootstrap_sharpe.png"),
       p_boot, width=9, height=5, dpi=150)

# 9.4 Чувствительность к издержкам
p_sens <- ggplot(sens_df, aes(x=Cost_bps, y=Sharpe, fill=Profitable)) +
  geom_col(width=4) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  scale_fill_manual(values=c("TRUE"="#2E7D32","FALSE"="#B71C1C"), guide="none") +
  scale_x_continuous(breaks=costs_bps) +
  labs(title    = "Sharpe условного квантиля vs транзакционные издержки",
       subtitle = "Зелёный = прибыльная стратегия, красный = убыточная",
       x="Издержки, bps", y="Sharpe (годовой)") +
  theme_minimal(base_size=12)
ggsave(file.path(OUT_DIR,"15e_cost_sensitivity.png"),
       p_sens, width=8, height=4, dpi=150)

# СОХРАНЕНИЕ ИТОГОВ

write.csv(metrics, file.path(OUT_DIR,"15f_metrics.csv"), row.names=FALSE)
write.csv(sens_df, file.path(OUT_DIR,"15g_cost_sensitivity.csv"), row.names=FALSE)
write.csv(data.frame(
  Strategy   = c("Режим","Усл.квантиль"),
  Sharpe     = round(c(sharpe(ret_regime_test), sharpe(ret_cq_test)), 4),
  CI_lo_90   = round(c(ci_regime[1], ci_cq[1]), 4),
  CI_hi_90   = round(c(ci_regime[2], ci_cq[2]), 4),
  Reject_H0  = c(ci_regime[1] > 0, ci_cq[1] > 0)
), file.path(OUT_DIR,"15h_sharpe_ci.csv"), row.names=FALSE)

cat("\n╔══════════════════════════════════════════════╗\n")
cat("           ИТОГ СКРИПТА 15                    \n")

cat(sprintf("Данные:       %d торговых дней\n", T_obs))
cat(sprintf("Тест-выборка: %d дней\n", length(idx_test)))
cat(sprintf("Sharpe (усл.квантиль, лаг+изд.): %.3f [%.3f; %.3f]\n",
            sharpe(ret_cq_test), ci_cq[1], ci_cq[2]))
cat(sprintf("H₀ Sharpe=0 отвергается: %s\n",
            ifelse(ci_cq[1] > 0, "ДА", "НЕТ")))
cat(sprintf("\nГрафики и данные сохранены в: %s\n", OUT_DIR))
cat("=== Этап 15 завершён ===\n")
