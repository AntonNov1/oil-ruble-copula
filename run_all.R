# ==============================================================================
# run_all.R — Запуск всего pipeline по порядку
# ==============================================================================
# Запустите этот скрипт из папки ThesisCopula/ в RStudio.
# Каждый этап независим — можно запускать отдельно.
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

run_step <- function(script, label) {
  cat(sprintf("\n%s\n=== %s ===\n%s\n", strrep("=", 60), label, strrep("=", 60)))
  t0 <- Sys.time()
  tryCatch(
    source(file.path("scripts", script)),
    error = function(e) cat(sprintf("ОШИБКА в %s: %s\n", script, e$message))
  )
  cat(sprintf("→ Завершено за %.1f мин\n", as.numeric(Sys.time() - t0, units = "mins")))
}

# Раскомментируйте нужные этапы:

run_step("01_data_prep.R",        "01. Подготовка данных и EDA")
# run_step("02_structural_breaks.R","02. Поиск структурных сдвигов (ДОЛГО)")
# run_step("03_copula_estimation.R","03. Оценка копул по периодам (ДОЛГО)")
run_step("04_gof_tests.R",        "04. GoF-тесты")
run_step("05_cross_period_val.R", "05. Кросс-периодная валидация")
run_step("06_copula_repeat.R",    "06. Повторяемость копул")
run_step("07_rolling_window.R",   "07. Оптимизация размера окна")
run_step("08_forecasting.R",      "08. Прогнозы и сценарный анализ")

cat("\n=== ВСЕ ЭТАПЫ ЗАВЕРШЕНЫ ===\n")
