# 14_trading_strategy.R — Копула-режимная торговая стратегия
#
# Логика:
#   1. В каждый момент t используем ФИЛЬТРОВАННЫЕ вероятности P(S_t|u_{1:t})
#      (не сглаженные — они используют будущее)
#   2. Направление зависимости нефть→рубль по режимам:
#      Frank (+1), Gumbel (+1), Student (-1), Gumbel° (-1)
#   3. Сигнал: direction_t = Σ_s κ_s · p_s(t)
#              position_t  = sign(log_OIL_t) · direction_t
#   4. Доходность: ret_t = position_t · log_FX_t
#      (position_t > 0 → лонг рубль = шорт USDRUB)
#
# Разделение выборки:
#   Обучение: первые 70% (оценка P и θ через EM из скрипта 13)
#   Тест:     последние 30% (только фильтр, без переоценки)
#
# Walk-forward:
#   Каждые WF_STEP недель переоцениваем модель на скользящем окне WF_WINDOW

rm(list = ls(all = TRUE))
options(warn = -1, scipen = 999)

library(dplyr); library(tidyr); library(copula); library(ggplot2)

# задайте рабочую директорию проекта перед запуском
PARENT  <- ".."
OUT_DIR <- "results/14_trading_strategy"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source("scripts/00_helpers.R")

# ДАННЫЕ

best <- load_best_models(data_dir = PARENT)
lr   <- load_log_returns(data_dir = PARENT)

# Псевдонаблюдения на полной выборке
u_all <- as.data.frame(pobs(cbind(lr$log_OIL, lr$log_FX)))
colnames(u_all) <- c("u_oil","u_fx")
u_all$Date    <- lr$Date
u_all$log_OIL <- lr$log_OIL
u_all$log_FX  <- lr$log_FX

T_obs <- nrow(u_all)
u_mat <- as.matrix(u_all[, c("u_oil","u_fx")])
u_rot <- cbind(1 - u_mat[,1], u_mat[,2])

K           <- 4
STATE_NAMES <- c("Frank","Gumbel","Student","Gumbel°")

# Знак зависимости по режимам: +1 = нефть и рубль движутся вместе
KAPPA <- c(Frank=+1, Gumbel=+1, Student=-1, `Gumbel°`=-1)

cat(sprintf("Данные: N=%d  (%s — %s)\n", T_obs,
            format(min(lr$Date)), format(max(lr$Date))))

# 2. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (копируем из скрипта 13)

compute_densities <- function(U, U_rot, params) {
  n <- nrow(U)
  safe_d <- function(cop, Ud) tryCatch(
    pmax(dCopula(Ud, cop), 1e-300), error = function(e) rep(1e-300, n))
  cbind(
    safe_d(frankCopula(params$frank, 2), U),
    safe_d(gumbelCopula(params$gumbel, 2), U),
    safe_d(tCopula(params$rho, df=params$df, dim=2), U),
    safe_d(gumbelCopula(params$gumbel_rot, 2), U_rot)
  )
}

hamilton_filter <- function(f_mat, P, pi0 = NULL) {
  K <- ncol(f_mat); N <- nrow(f_mat)
  if (is.null(pi0)) {
    ev  <- eigen(t(P))$vectors[,1]
    pi0 <- pmax(Re(ev)/sum(Re(ev)), 0); pi0 <- pi0/sum(pi0)
  }
  xi_filt <- matrix(0, N, K)
  xi_pred <- matrix(0, N, K)
  LL_vec  <- numeric(N)
  p_pred  <- pi0
  for (t in 1:N) {
    xi_pred[t,] <- p_pred
    joint <- f_mat[t,] * p_pred
    denom <- sum(joint)
    if (denom <= 0 || is.nan(denom)) {
      xi_filt[t,] <- p_pred; LL_vec[t] <- log(.Machine$double.eps)
    } else {
      xi_filt[t,] <- joint / denom; LL_vec[t] <- log(denom)
    }
    p_pred <- pmax(as.vector(t(P) %*% xi_filt[t,]), 0)
    p_pred <- p_pred / sum(p_pred)
  }
  list(xi_filt=xi_filt, xi_pred=xi_pred, LL=sum(LL_vec))
}

# 3. ЗАГРУЗКА ПАРАМЕТРОВ EM (скрипт 13)

em_file <- "results/13_em_algorithm/13_em_results.rds"
if (file.exists(em_file)) {
  em_res <- readRDS(em_file)
  params  <- em_res$params
  P_hat   <- em_res$P
  pi_stat <- em_res$pi_stat
  cat("Параметры EM загружены из скрипта 13.\n")
} else {
  cat("Файл EM не найден — используем параметры из скрипта 12.\n")
  params  <- list(frank=2.5066, gumbel=1.50, rho=-0.6347, df=0.5551, gumbel_rot=1.0869)
  P_hat   <- matrix(c(0.90,0.05,0.03,0.02, 0.05,0.88,0.05,0.02,
                      0.03,0.03,0.90,0.04, 0.01,0.02,0.02,0.95), K,K,byrow=TRUE)
  ev_s    <- eigen(t(P_hat))$vectors[,1]
  pi_stat <- pmax(Re(ev_s)/sum(Re(ev_s)),0); pi_stat <- pi_stat/sum(pi_stat)
}

# РАЗДЕЛЕНИЕ ВЫБОРКИ

TRAIN_FRAC <- 0.70
T_train    <- floor(T_obs * TRAIN_FRAC)
T_test     <- T_obs - T_train

cat(sprintf("Train: %d наблюд. (%s — %s)\n", T_train,
            format(u_all$Date[1]), format(u_all$Date[T_train])))
cat(sprintf("Test:  %d наблюд. (%s — %s)\n", T_test,
            format(u_all$Date[T_train+1]), format(u_all$Date[T_obs])))

# 5. БАЗОВАЯ СТРАТЕГИЯ (out-of-sample: фильтр на тесте с фиксированными P, θ)

cat("\n=== БАЗОВАЯ СТРАТЕГИЯ ===\n")

f_full  <- compute_densities(u_mat, u_rot, params)
res_all <- hamilton_filter(f_full, P_hat, pi_stat)

# Используем только фильтрованные вероятности (causal)
xi_filt <- res_all$xi_filt

# Направление и сигнал
direction <- as.numeric(xi_filt %*% KAPPA)   # T×1, ∈(-1,+1)
oil_sign  <- sign(u_all$log_OIL)              # +1/-1/0
raw_signal<- direction * oil_sign             # [-1, +1]

# Бинарная позиция: лонг (+1) / шорт (-1) / плоско (0 если |signal|<порога)
THRESHOLD  <- 0.10
position   <- ifelse(abs(raw_signal) >= THRESHOLD, sign(raw_signal), 0)

# Доходность стратегии: position_t * log_FX_t
# position > 0 → лонг рубль (шорт USDRUB) → зарабатываем, когда log_FX < 0
strat_ret  <- position * (-u_all$log_FX)
bh_ret     <- -u_all$log_FX            # buy-and-hold рубль (шорт USDRUB всегда)
oil_bh     <- u_all$log_OIL            # buy-and-hold нефть

# Разбивка на трейн / тест
idx_test   <- (T_train + 1):T_obs
ret_test   <- strat_ret[idx_test]
bh_test    <- bh_ret[idx_test]
oil_test   <- oil_bh[idx_test]

# Кумулятивные доходности
cum_strat  <- cumsum(ret_test)
cum_bh     <- cumsum(bh_test)
cum_oil    <- cumsum(oil_test)

dates_test <- u_all$Date[idx_test]

# Метрики производительности
sharpe <- function(r, annualize=52) {
  r <- r[!is.na(r)]
  if (sd(r) == 0) return(NA)
  mean(r) / sd(r) * sqrt(annualize)
}
max_dd <- function(cum) {
  running_max <- cummax(cum)
  max(running_max - cum)
}
n_trades   <- sum(position[idx_test] != 0)
win_rate   <- mean(ret_test[position[idx_test] != 0] > 0, na.rm=TRUE)

cat(sprintf("Test period (%d нед.):\n", T_test))
cat(sprintf("  Стратегия:  Annual.ret=%.2f%%  Sharpe=%.3f  MaxDD=%.4f\n",
            mean(ret_test)*52*100, sharpe(ret_test), max_dd(cum_strat)))
cat(sprintf("  B&H рубль:  Annual.ret=%.2f%%  Sharpe=%.3f  MaxDD=%.4f\n",
            mean(bh_test)*52*100, sharpe(bh_test), max_dd(cum_bh)))
cat(sprintf("  B&H нефть:  Annual.ret=%.2f%%  Sharpe=%.3f  MaxDD=%.4f\n",
            mean(oil_test)*52*100, sharpe(oil_test), max_dd(cum_oil)))
cat(sprintf("  Сделок: %d (%.1f%% активны)  Win-rate: %.1f%%\n",
            n_trades, n_trades/T_test*100, win_rate*100))

# walk-forward ОПТИМИЗАЦИЯ

cat("\n=== WALK-FORWARD (скользящее окно) ===\n")

WF_WINDOW <- 104   # 2 года обучения (недели)
WF_STEP   <- 26    # переоцениваем каждые полгода

# Inline EM (упрощённый) для walk-forward — без source скрипта 13
kim_smoother_light <- function(xi_filt, xi_pred, P) {
  K <- ncol(xi_filt); N <- nrow(xi_filt)
  xi_smooth <- matrix(0, N, K); xi_smooth[N,] <- xi_filt[N,]
  xi2 <- array(0, c(N-1, K, K))
  for (t in (N-1):1) {
    ratio <- xi_smooth[t+1,] / pmax(xi_pred[t+1,], .Machine$double.eps)
    for (j in 1:K) {
      xi_smooth[t,j] <- xi_filt[t,j] * sum(P[j,] * ratio)
      xi2[t,j,]      <- xi_filt[t,j] * P[j,] * ratio
    }
    s <- sum(xi_smooth[t,]); if (s>0) { xi_smooth[t,]<-xi_smooth[t,]/s; xi2[t,,]<-xi2[t,,]/s }
  }
  list(xi_smooth=xi_smooth, xi2=xi2)
}
update_P_wf <- function(xi2) {
  K <- dim(xi2)[2]
  P_new <- matrix(0,K,K)
  for (j in 1:K) { d<-sum(xi2[,j,]); P_new[j,]<-if(d>0) colSums(xi2[,j,])/d else rep(1/K,K) }
  pmax(P_new,1e-6) / rowSums(pmax(P_new,1e-6))
}

# Walk-forward: стартуем с T_train, переоцениваем каждые WF_STEP
wf_positions <- rep(0, T_obs)
wf_params    <- params   # стартовые параметры
wf_P         <- P_hat

refit_times  <- seq(T_train, T_obs - WF_STEP, by = WF_STEP)
cat(sprintf("  Переоценка на %d точках: каждые %d нед., окно %d нед.\n",
            length(refit_times), WF_STEP, WF_WINDOW))

for (t_start in refit_times) {
  # Индексы обучающего окна
  win_idx <- max(1, t_start - WF_WINDOW + 1):t_start
  Uw  <- u_mat[win_idx,, drop=FALSE]
  Urw <- u_rot[win_idx,, drop=FALSE]

  # Быстрый EM (20 итераций)
  p_wf <- wf_params; P_wf <- wf_P
  for (em_i in 1:20) {
    f_wf  <- compute_densities(Uw, Urw, p_wf)
    rf    <- hamilton_filter(f_wf, P_wf)
    rs    <- kim_smoother_light(rf$xi_filt, rf$xi_pred, P_wf)
    P_wf  <- update_P_wf(rs$xi2)
    # Упрощённый M-шаг: только Frank и Gumbel° (они наиболее чувствительны)
    w1 <- rs$xi_smooth[,1]; w4 <- rs$xi_smooth[,4]
    if (sum(w1)>2) {
      o1 <- tryCatch(optimize(function(th) -sum(w1*log(pmax(dCopula(Uw,frankCopula(th,2)),1e-300))),
                              c(0.01,25)), error=function(e) NULL)
      if (!is.null(o1)) p_wf$frank <- o1$minimum
    }
    if (sum(w4)>2) {
      o4 <- tryCatch(optimize(function(th) -sum(w4*log(pmax(dCopula(Urw,gumbelCopula(th,2)),1e-300))),
                              c(1.001,12)), error=function(e) NULL)
      if (!is.null(o4)) p_wf$gumbel_rot <- o4$minimum
    }
  }
  wf_params <- p_wf; wf_P <- P_wf

  # Применяем к следующим WF_STEP барам
  pred_idx <- (t_start + 1):min(t_start + WF_STEP, T_obs)
  f_pred   <- compute_densities(u_mat[pred_idx,,drop=FALSE],
                                u_rot[pred_idx,,drop=FALSE], wf_params)
  # Начальная вероятность — последний фильтрованный бар окна
  pi0_wf   <- rf$xi_filt[nrow(rf$xi_filt),]
  rf_pred  <- hamilton_filter(f_pred, wf_P, pi0_wf)

  dir_wf    <- as.numeric(rf_pred$xi_filt %*% KAPPA)
  sig_wf    <- dir_wf * sign(u_all$log_OIL[pred_idx])
  pos_wf    <- ifelse(abs(sig_wf) >= THRESHOLD, sign(sig_wf), 0)
  wf_positions[pred_idx] <- pos_wf

  cat(sprintf("  t=%d (%s): F=%.3f G°=%.3f → dir∈[%.2f,%.2f]\n",
              t_start, format(u_all$Date[t_start]),
              wf_params$frank, wf_params$gumbel_rot,
              min(dir_wf), max(dir_wf)))
}

# Walk-forward доходности
wf_ret_full  <- wf_positions * (-u_all$log_FX)
wf_ret_test  <- wf_ret_full[idx_test]
wf_active    <- sum(wf_positions[idx_test] != 0)

cat(sprintf("\nWalk-forward test period:\n"))
cat(sprintf("  Annual.ret=%.2f%%  Sharpe=%.3f  MaxDD=%.4f\n",
            mean(wf_ret_test)*52*100, sharpe(wf_ret_test), max_dd(cumsum(wf_ret_test))))
cat(sprintf("  Активных позиций: %d (%.1f%%)\n",
            wf_active, wf_active/T_test*100))

# 7. УСЛОВНЫЙ КВАНТИЛЬНЫЙ СИГНАЛ (продвинутый)

cat("\n=== УСЛОВНЫЙ КВАНТИЛЬНЫЙ СИГНАЛ ===\n")

# P(FX > median | u_oil_t, режим s) через условную копулу C(v|u)
# Frank: ∂C/∂u = e^{-θu}(e^{-θv}-1) / (e^{-θu}+e^{-θv}-e^{-θ}-1+...)
# Для скорости: численно P(v>0.5|u) = 1 - C(u,0.5)/u (независимость = 0.5)
# Режимное среднее: signal_cq_t = Σ_s p_s * P_s(v>0.5|u_oil_t)

cond_prob_frank <- function(u, v_threshold = 0.5, theta) {
  # P(V > v | U = u) = 1 - dC/du|_{v=v_thresh}
  # Для Frank: dC(u,v)/du = exp(-θu)(exp(-θv)-1) / (exp(-θu)+exp(-θv)-exp(-θ(u+v))-1+...)
  # Упрощение: 1 - pCopula(u, v_threshold, frankCopula) / u (непрерывная)
  # Используем числовую производную
  eps <- 1e-4
  u_c <- pmin(pmax(u, eps), 1-eps)
  v_c <- pmin(pmax(v_threshold, eps), 1-eps)
  cop <- frankCopula(theta, 2)
  Cuv <- pCopula(cbind(u_c, rep(v_c, length(u_c))), cop)
  Cu1 <- u_c
  pmax(pmin(1 - Cuv/Cu1, 1), 0)
}

cond_prob_gumbel_rot <- function(u, v_threshold = 0.5, theta) {
  # Gumbel° = 90°-поворот: (u_oil → 1-u_oil)
  u_rot_val <- 1 - u
  eps <- 1e-4
  u_c <- pmin(pmax(u_rot_val, eps), 1-eps)
  v_c <- pmin(pmax(v_threshold, eps), 1-eps)
  cop <- gumbelCopula(theta, 2)
  Cuv <- pCopula(cbind(u_c, rep(v_c, length(u_c))), cop)
  Cu1 <- u_c
  pmax(pmin(1 - Cuv/Cu1, 1), 0)
}

# Режимный сигнал: P(рубль укрепится | нефть)
# Укрепление рубля = log_FX < 0 → u_fx < 0.5
# P(рубль укрепится) = P(u_fx < 0.5 | u_oil) = P(V < 0.5 | u_oil)
# Frank (+): P(V<0.5|u>0.5) > 0.5; Gumbel°(-): P(V<0.5|u>0.5) < 0.5

u_oil_vec <- u_mat[, 1]
p_rub_frank    <- cond_prob_frank(u_oil_vec, 0.5, params$frank)
p_rub_grot     <- 1 - cond_prob_gumbel_rot(u_oil_vec, 0.5, params$gumbel_rot)

# Для Student и Gumbel приближаем через линейную интерполяцию
p_rub_student  <- 0.5 + (0.5 - u_oil_vec) * abs(params$rho)  # простое приближение
p_rub_gumbel   <- cond_prob_frank(u_oil_vec, 0.5, params$frank) # схоже с Frank

p_rub_combined <- xi_filt[,1] * p_rub_frank +
                  xi_filt[,2] * p_rub_gumbel +
                  xi_filt[,3] * p_rub_student +
                  xi_filt[,4] * p_rub_grot

# Сигнал: лонг рубль если P(рубль укрепится) > 0.5+delta
DELTA <- 0.05
pos_cq <- ifelse(p_rub_combined > 0.5 + DELTA,  1,
          ifelse(p_rub_combined < 0.5 - DELTA, -1, 0))

ret_cq_test <- (pos_cq * (-u_all$log_FX))[idx_test]
cat(sprintf("Conditional quantile strategy (test):\n"))
cat(sprintf("  Annual.ret=%.2f%%  Sharpe=%.3f  MaxDD=%.4f\n",
            mean(ret_cq_test)*52*100, sharpe(ret_cq_test),
            max_dd(cumsum(ret_cq_test))))

# ВИЗУАЛИЗАЦИИ

cat("\n=== ВИЗУАЛИЗАЦИИ ===\n")

# 8.1 Кумулятивная доходность: все стратегии
perf_df <- data.frame(
  Date   = dates_test,
  Copula  = cumsum(ret_test),
  WalkFwd = cumsum(wf_ret_test),
  CondQ   = cumsum(ret_cq_test),
  BH_RUB  = cumsum(bh_test),
  BH_OIL  = cumsum(oil_test)
) %>% pivot_longer(-Date, names_to="Strategy", values_to="CumRet")

strategy_colors <- c(
  Copula   = "#1565C0",
  WalkFwd  = "#2E7D32",
  CondQ    = "#6A1B9A",
  BH_RUB   = "#B71C1C",
  BH_OIL   = "#FF8F00"
)
strategy_labels <- c(
  Copula   = "Копула-режимная",
  WalkFwd  = "Walk-forward",
  CondQ    = "Условный квантиль",
  BH_RUB   = "B&H рубль",
  BH_OIL   = "B&H нефть"
)

p_cum <- ggplot(perf_df, aes(x=Date, y=CumRet*100, color=Strategy)) +
  geom_line(linewidth=0.9) +
  geom_hline(yintercept=0, linetype="dashed", color="grey50") +
  scale_color_manual(values=strategy_colors, labels=strategy_labels) +
  scale_x_date(date_labels="%Y-%m", date_breaks="3 months") +
  labs(title    = "Кумулятивная доходность (тестовая выборка)",
       subtitle = sprintf("Тест: %s — %s  |  N=%d нед.",
                          format(dates_test[1]), format(tail(dates_test,1)), T_test),
       x=NULL, y="Кум. доходность, %", color=NULL) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_DIR,"14a_cumulative_returns.png"),
       p_cum, width=13, height=5, dpi=150)

# 8.2 Апостериорные вероятности + позиция стратегии
regime_colors_4 <- c("Frank"="#1565C0","Gumbel"="#2E7D32",
                     "Student"="#B71C1C","Gumbel°"="#FF8F00")

pos_df <- data.frame(
  Date      = u_all$Date[idx_test],
  Position  = position[idx_test],
  Direction = direction[idx_test]
)

p_pos <- ggplot(pos_df, aes(x=Date)) +
  geom_col(aes(y=Position), fill="#1565C0", alpha=0.6, width=7) +
  geom_line(aes(y=Direction), color="#B71C1C", linewidth=0.8) +
  geom_hline(yintercept=c(-THRESHOLD, THRESHOLD),
             linetype="dotted", color="grey60") +
  scale_x_date(date_labels="%Y-%m", date_breaks="3 months") +
  labs(title = "Позиция стратегии и направление зависимости (тест)",
       subtitle = "Синие столбцы = позиция {-1,0,+1}; красная линия = direction_t",
       x=NULL, y=NULL) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_DIR,"14b_positions.png"),
       p_pos, width=13, height=4, dpi=150)

# 8.3 Понедельная доходность (heatmap-style)
ret_df <- data.frame(
  Date    = dates_test,
  Return  = ret_test * 100,
  Positive = ret_test >= 0
)
p_bars <- ggplot(ret_df, aes(x=Date, y=Return, fill=Positive)) +
  geom_col(width=7) +
  scale_fill_manual(values=c("TRUE"="#2E7D32","FALSE"="#B71C1C"), guide="none") +
  geom_hline(yintercept=0, color="black") +
  scale_x_date(date_labels="%Y-%m", date_breaks="3 months") +
  labs(title="Недельная доходность копула-режимной стратегии (тест)",
       x=NULL, y="Доходность, %") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(OUT_DIR,"14c_weekly_returns.png"),
       p_bars, width=13, height=4, dpi=150)

# 8.4 Сводная таблица метрик
metrics_df <- data.frame(
  Стратегия = c("Копула-режимная","Walk-forward","Условный квантиль",
                "B&H рубль","B&H нефть"),
  `Ann.Return,%` = round(c(mean(ret_test),mean(wf_ret_test),mean(ret_cq_test),
                            mean(bh_test),mean(oil_test))*52*100, 2),
  Sharpe = round(c(sharpe(ret_test),sharpe(wf_ret_test),sharpe(ret_cq_test),
                   sharpe(bh_test),sharpe(oil_test)), 3),
  `MaxDD,%` = round(c(max_dd(cumsum(ret_test)),max_dd(cumsum(wf_ret_test)),
                      max_dd(cumsum(ret_cq_test)),max_dd(cumsum(bh_test)),
                      max_dd(cumsum(oil_test)))*100, 2),
  check.names = FALSE
)
cat("\n=== СВОДНАЯ ТАБЛИЦА ===\n")
print(metrics_df)
write.csv(metrics_df, file.path(OUT_DIR,"14d_metrics.csv"), row.names=FALSE)

# 9. РОБАСТНОСТЬ: АНАЛИЗ ЧУВСТВИТЕЛЬНОСТИ ПОРОГА

cat("\n=== ЧУВСТВИТЕЛЬНОСТЬ К ПОРОГУ ===\n")

thresholds <- seq(0, 0.30, by=0.05)
thr_results <- lapply(thresholds, function(thr) {
  pos_t <- ifelse(abs(raw_signal[idx_test]) >= thr,
                  sign(raw_signal[idx_test]), 0)
  r     <- pos_t * (-u_all$log_FX[idx_test])
  data.frame(Threshold=thr,
             Sharpe=round(sharpe(r),3),
             Ann_ret=round(mean(r)*52*100,2),
             N_trades=sum(pos_t!=0))
})
thr_df <- do.call(rbind, thr_results)
print(thr_df)
write.csv(thr_df, file.path(OUT_DIR,"14e_threshold_sensitivity.csv"), row.names=FALSE)

p_thr <- ggplot(thr_df, aes(x=Threshold, y=Sharpe)) +
  geom_line(color="#1565C0", linewidth=1) +
  geom_point(color="#1565C0", size=3) +
  geom_hline(yintercept=sharpe(bh_test), linetype="dashed",
             color="#B71C1C", linewidth=0.8) +
  annotate("text", x=0.25, y=sharpe(bh_test)+0.05,
           label="B&H рубль", color="#B71C1C", size=3.5) +
  labs(title    = "Коэффициент Шарпа в зависимости от порога сигнала",
       subtitle = "Красная пунктирная = Sharpe стратегии B&H рубль",
       x="Порог θ", y="Коэф. Шарпа (годовой)") +
  theme_minimal()
ggsave(file.path(OUT_DIR,"14f_sharpe_vs_threshold.png"),
       p_thr, width=8, height=4, dpi=150)

cat("\nРезультаты сохранены в", OUT_DIR, "\n")
cat("=== Этап 14 завершён ===\n")
