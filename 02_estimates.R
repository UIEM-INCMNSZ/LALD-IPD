if (!exists("lald_out")) source("01_prepare.R")

local({
  set.seed(42)
  db <- readRDS(lald_cache("M0_db_derived.rds"))

  COVS <- c("age_dx_y", "severity_index", "LLT_harmonized")

  fit_one <- function(d, y, adjust) {
    d$.y <- d[[y]]
    vars <- c(".y", "ert_received", if (adjust) COVS)
    d <- d[stats::complete.cases(d[, vars]), vars]
    if (nrow(d) < 8) return(NULL)
    if (length(unique(d$.y)) < 2 || length(unique(d$ert_received)) < 2) return(NULL)
    keep <- c("ert_received", if (adjust) COVS[sapply(COVS, function(v) length(unique(d[[v]])) > 1)])
    f <- as.formula(paste(".y ~", paste(keep, collapse = " + ")))
    fit <- tryCatch(lald_firth(f, d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    j <- which(fit$names == "ert_received")
    if (!length(j)) return(NULL)
    e <- d[d$ert_received == 1, ]; n <- d[d$ert_received == 0, ]
    list(n = nrow(d), n_ert = nrow(e), n_non = nrow(n),
         ev_ert = sum(e$.y == 1), ev_non = sum(n$.y == 1),
         beta = fit$coef[j], lcl = fit$ci[j, 1], ucl = fit$ci[j, 2], p = fit$p[j],
         covs = paste(setdiff(keep, "ert_received"), collapse = ","), engine = fit$engine)
  }

  OUT <- list()
  spec <- list(
    list("O1",  "Mortality (any cause)",       "mortality",     "CESD", TRUE),
    list("O1",  "Mortality (any cause)",       "mortality",     "WD",   TRUE),
    list("O2",  "ALT normalization",           "ALT_norm",      "CESD", TRUE),
    list("O2x", "AST normalization",           "AST_norm",      "CESD", TRUE),
    list("O5",  "Steatosis regression",        "O5_steatosis",  "CESD", TRUE),
    list("O6",  "Liver transplantation",       "liver_tx",      "CESD", TRUE),
    list("O6",  "Liver transplantation",       "liver_tx",      "WD",   TRUE),
    list("O4",  "Favourable lipid response",   "O4_lipid",      "CESD", FALSE),
    list("O5s", "Steatosis regression (unrestricted)", "O5_steatosis_unrestricted", "CESD", FALSE))

  for (s in spec) {
    id <- s[[1]]; lab <- s[[2]]; v <- s[[3]]; ph <- s[[4]]; primary <- s[[5]]
    d <- db[db$phenotype %in% ph, ]
    for (adj in c(FALSE, TRUE)) {
      r <- fit_one(d, v, adj)
      if (is.null(r)) {
        OUT[[length(OUT) + 1]] <- data.frame(outcome_id = id, outcome = lab, phenotype = ph,
          model = ifelse(adj, "Firth adjusted", "Firth unadjusted"), primary = primary,
          n = NA, n_ert = NA, n_non = NA, ev_ert = NA, ev_non = NA,
          OR = NA, lcl = NA, ucl = NA, p = NA, d = NA, covariates = "", engine = "",
          note = "not estimable in this analysis set")
        next
      }
      OUT[[length(OUT) + 1]] <- data.frame(outcome_id = id, outcome = lab, phenotype = ph,
        model = ifelse(adj, "Firth adjusted", "Firth unadjusted"), primary = primary,
        n = r$n, n_ert = r$n_ert, n_non = r$n_non, ev_ert = r$ev_ert, ev_non = r$ev_non,
        OR = exp(r$beta), lcl = exp(r$lcl), ucl = exp(r$ucl), p = r$p,
        d = d_from_logOR(r$beta), covariates = r$covs, engine = r$engine, note = "")
    }
  }
  R <- do.call(rbind, OUT)

  R$p_FDR <- NA
  for (ph in unique(R$phenotype)) {
    i <- which(R$phenotype == ph & R$primary & R$model == "Firth adjusted" & !is.na(R$p))
    if (length(i)) R$p_FDR[i] <- bh(R$p[i])
  }
  R$estimate_fmt <- ifelse(is.na(R$OR), "-", sprintf("OR %.3f [%.3f, %.3f]", R$OR, R$lcl, R$ucl))
  R$p_fmt     <- ifelse(is.na(R$p), "-", ifelse(R$p < 0.001, "<0.001", sprintf("%.4f", R$p)))
  R$p_FDR_fmt <- ifelse(is.na(R$p_FDR), "", ifelse(R$p_FDR < 0.001, "<0.001", sprintf("%.4f", R$p_FDR)))
  write.csv(R, lald_out("M2_binary_results.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db <- readRDS(lald_cache("M0_db_derived.rds"))
  COVS <- c("age_dx_y", "severity_index", "LLT_harmonized")

  OUTS <- list(c("O2b", "ALT percent change", "ALT_pct", "TRUE"),
               c("O2c", "AST percent change", "AST_pct", "TRUE"),
               c("O3",  "LDL-C percent change", "LDL_pct", "TRUE"),
               c("O3b", "TG percent change",  "TG_pct",  "TRUE"),
               c("O3c", "Total cholesterol percent change", "TC_pct", "FALSE"))

  res <- list()
  for (o in OUTS) {
    id <- o[1]; lab <- o[2]; v <- o[3]; primary <- as.logical(o[4])
    d <- db[db$phenotype %in% "CESD" & !is.na(db[[v]]) & !is.na(db$ert_received), ]
    e <- d[[v]][d$ert_received == 1]; n <- d[[v]][d$ert_received == 0]
    ts <- theil_sen(d$ert_received, d[[v]])
    hl_stat <- function(dat, i) {
      s <- dat[i, ]
      a <- s[[v]][s$ert_received == 1]; b <- s[[v]][s$ert_received == 0]
      if (length(a) < 2 || length(b) < 2) return(NA_real_)
      median(outer(a, b, "-"))
    }
    bc <- tryCatch(bca_ci(d[, c(v, "ert_received")], function(dat, i) {
          s <- dat[i, ]; a <- s[[1]][s[[2]] == 1]; b <- s[[1]][s[[2]] == 0]
          if (length(a) < 2 || length(b) < 2) return(NA_real_); median(outer(a, b, "-")) },
          R = 2000), error = function(x) NULL)
    dd <- d[stats::complete.cases(d[, c(v, "ert_received", COVS)]), ]
    keep <- COVS[sapply(COVS, function(cv) length(unique(dd[[cv]])) > 1)]
    fit <- if (nrow(dd) >= 10)
      lm(as.formula(paste(v, "~ ert_received +", paste(keep, collapse = " + "))), data = dd) else NULL
    cf <- if (!is.null(fit)) summary(fit)$coefficients else NULL
    b  <- if (!is.null(cf) && "ert_received" %in% rownames(cf)) cf["ert_received", ] else NULL
    ci <- if (!is.null(fit)) confint(fit)["ert_received", ] else c(NA, NA)
    g  <- hedges_g(e, n)
    res[[length(res) + 1]] <- data.frame(
      outcome_id = id, outcome = lab, phenotype = "CESD", primary = primary,
      n = nrow(d), n_ert = length(e), n_non = length(n),
      median_ert = median(e), median_non = median(n),
      ts_estimate = ts$slope,
      ts_lcl = if (!is.null(bc)) bc$lcl else NA, ts_ucl = if (!is.null(bc)) bc$ucl else NA,
      ts_ci_type = if (!is.null(bc)) bc$type else "", ts_p = ts$p,
      ols_n = if (!is.null(fit)) nrow(dd) else NA,
      ols_beta = if (!is.null(b)) b[1] else NA,
      ols_lcl = ci[1], ols_ucl = ci[2],
      ols_p = if (!is.null(b)) b[4] else NA,
      ols_r2 = if (!is.null(fit)) summary(fit)$r.squared else NA,
      ols_covs = paste(keep, collapse = ","),
      g = g[1], g_lcl = g[2], g_ucl = g[3])
  }
  R <- do.call(rbind, res)
  R$p_FDR_ols <- NA
  i <- which(R$primary & !is.na(R$ols_p)); if (length(i)) R$p_FDR_ols[i] <- bh(R$ols_p[i])
  R$p_FDR_ts <- NA
  i <- which(R$primary & !is.na(R$ts_p)); if (length(i)) R$p_FDR_ts[i] <- bh(R$ts_p[i])
  write.csv(R, lald_out("M3_continuous_results.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  suppressPackageStartupMessages(library(survival))
  db <- readRDS(lald_cache("M0_db_derived.rds"))

  rmst <- function(fit, tau) {
    s <- summary(fit, times = sort(unique(pmin(fit$time, tau))), extend = TRUE)
    tt <- c(0, s$time[s$time <= tau], tau); ss <- c(1, s$surv[s$time <= tau])
    ss <- c(ss, ss[length(ss)])
    sum(diff(tt) * head(ss, -1))
  }

  res <- list()
  for (ph in c("CESD", "WD")) {
    d <- db[db$phenotype %in% ph & !is.na(db$fu_years) & !is.na(db$mortality) &
            !is.na(db$ert_received) & db$fu_years >= 0, ]
    if (nrow(d) < 10) next
    d$fu_years[d$fu_years == 0] <- 0.01
    S <- Surv(d$fu_years, d$mortality)
    km <- survfit(S ~ ert_received, data = d)
    lr <- tryCatch(survdiff(S ~ ert_received, data = d), error = function(e) NULL)
    cox <- tryCatch(coxph(S ~ ert_received, data = d), error = function(e) NULL)
    covs <- c("age_dx_y", "severity_index", "LLT_harmonized")
    da <- d[stats::complete.cases(d[, c("fu_years", "mortality", "ert_received", covs)]), ]
    coxa <- if (nrow(da) >= 20) tryCatch(
        coxph(Surv(fu_years, mortality) ~ ert_received + age_dx_y + severity_index + LLT_harmonized,
              data = da), error = function(e) NULL) else NULL
    tau <- as.numeric(quantile(d$fu_years, 0.8))
    k1 <- survfit(Surv(fu_years, mortality) ~ 1, data = d[d$ert_received == 1, ])
    k0 <- survfit(Surv(fu_years, mortality) ~ 1, data = d[d$ert_received == 0, ])
    r1 <- tryCatch(rmst(k1, tau), error = function(e) NA)
    r0 <- tryCatch(rmst(k0, tau), error = function(e) NA)
    zph <- if (!is.null(cox)) tryCatch(cox.zph(cox)$table["GLOBAL", "p"], error = function(e) NA) else NA
    degenerate <- sum(d$mortality[d$ert_received == 1]) == 0 ||
                  sum(d$mortality[d$ert_received == 0]) == 0
    if (degenerate) { cox <- NULL; coxa <- NULL }

    res[[length(res) + 1]] <- data.frame(
      phenotype = ph, n = nrow(d), n_ert = sum(d$ert_received == 1), n_non = sum(d$ert_received == 0),
      deaths_ert = sum(d$mortality[d$ert_received == 1]), deaths_non = sum(d$mortality[d$ert_received == 0]),
      median_fu = median(d$fu_years),
      logrank_p = if (!is.null(lr)) pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE) else NA,
      HR = if (!is.null(cox)) exp(coef(cox)[1]) else NA,
      HR_lcl = if (!is.null(cox)) exp(confint(cox)[1, 1]) else NA,
      HR_ucl = if (!is.null(cox)) exp(confint(cox)[1, 2]) else NA,
      HR_p = if (!is.null(cox)) summary(cox)$coefficients[1, 5] else NA,
      HR_adj = if (!is.null(coxa)) exp(coef(coxa)[1]) else NA,
      HR_adj_lcl = if (!is.null(coxa)) exp(confint(coxa)[1, 1]) else NA,
      HR_adj_ucl = if (!is.null(coxa)) exp(confint(coxa)[1, 2]) else NA,
      HR_adj_p = if (!is.null(coxa)) summary(coxa)$coefficients[1, 5] else NA,
      n_adj = if (!is.null(coxa)) nrow(da) else NA,
      ph_assumption_p = zph,
      tau = tau, rmst_ert = r1, rmst_non = r0, rmst_diff = r1 - r0,
      cox_note = if (degenerate)
        "Cox not estimable: zero deaths in one arm; log-rank and RMST are the reported comparisons" else "")

  }
  R <- do.call(rbind, res)
  write.csv(R, lald_out("M4_survival_results.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db  <- readRDS(lald_cache("M0_db_derived.rds"))
  agg <- read.csv(lald_in("aggregate"), check.names = FALSE, na.strings = c("", "NA"))
  agg <- agg[agg$row_kind == "result", ]

  pick <- function(rid, outcome, arm, timepoint = NULL, n_expect = NULL, need = "events") {
    s <- agg[agg$R_ID == rid & agg$outcome == outcome & agg$arm == arm, ]
    if (!is.null(timepoint)) s <- s[grepl(timepoint, s$timepoint, fixed = TRUE), ]
    if (!is.null(n_expect))  s <- s[!is.na(s$n) & s$n == n_expect, ]
    s <- s[!is.na(s[[need]]) & !is.na(s$n), ]
    if (nrow(s) == 0) stop(sprintf("M6: no usable row for %s / %s / %s / need=%s", rid, outcome, arm, need))
    s[1, ]
  }
  ipd_prop <- function(d, v, ert) {
    keep <- !is.na(d[[v]]) & !is.na(d$ert_received) & d$ert_received == ert
    d <- d[keep, ]
    c(events = sum(d[[v]] == 1), n = nrow(d))
  }
  ipd_mean <- function(d, v, ert) {
    keep <- !is.na(d[[v]]) & !is.na(d$ert_received) & d$ert_received == ert
    x <- d[[v]][keep]
    c(mean = mean(x), sd = sd(x), n = length(x))
  }
  CESD <- db[db$phenotype %in% "CESD", ]
  WD   <- db[db$phenotype %in% "WD", ]
  WD$surv12 <- with(WD, ifelse(
    !is.na(mortality) & mortality == 1 & !is.na(age_death_y), as.integer(age_death_y > 1),
    ifelse(!is.na(mortality) & mortality == 0 & !is.na(fu_years) & fu_years >= 1, 1L, NA)))

  FOREST <- list()
  RES <- list()

  a_ipd_e <- ipd_prop(CESD, "ALT_norm", 1); a_ipd_c <- ipd_prop(CESD, "ALT_norm", 0)
  r27e <- pick("R027", "ALT normalization", "ERT", "week 20"); r27c <- pick("R027", "ALT normalization", "Placebo", "week 20")
  r01e <- pick("R001", "ALT normalization", "ERT", "1 year");  r01c <- pick("R001", "ALT normalization", "Untreated", "1 year")
  F1 <- data.frame(
    study = c("IPD (this review)", "ARISE, Burton 2015", "LAL-D Registry, D'Antiga 2025"),
    group = c("IPD", "ARISE / LAL-CL02", "International LAL-D Registry"),
    a = c(a_ipd_e["events"], r27e$events, r01e$events),
    n1 = c(a_ipd_e["n"], r27e$n, r01e$n),
    c = c(a_ipd_c["events"], r27c$events, r01c$events),
    n2 = c(a_ipd_c["n"], r27c$n, r01c$n))
  b_ipd_e <- ipd_prop(CESD, "AST_norm", 1); b_ipd_c <- ipd_prop(CESD, "AST_norm", 0)
  s27e <- pick("R027", "AST normalization", "ERT", "week 20"); s27c <- pick("R027", "AST normalization", "Placebo", "week 20")
  s01e <- pick("R001", "AST normalization", "ERT", "1 year");  s01c <- pick("R001", "AST normalization", "Untreated", "1 year")
  F2 <- data.frame(
    study = c("IPD (this review)", "ARISE, Burton 2015", "LAL-D Registry, D'Antiga 2025"),
    group = c("IPD", "ARISE / LAL-CL02", "International LAL-D Registry"),
    a = c(b_ipd_e["events"], s27e$events, s01e$events),
    n1 = c(b_ipd_e["n"], s27e$n, s01e$n),
    c = c(b_ipd_c["events"], s27c$events, s01c$events),
    n2 = c(b_ipd_c["n"], s27c$n, s01c$n))
  w_e <- ipd_prop(WD, "surv12", 1); w_c <- ipd_prop(WD, "surv12", 0)
  v31 <- pick("R031", "12-month survival", "ERT", n_expect = 19)
  h62 <- pick("R062", "12-month survival", "Untreated", n_expect = 21)
  F5 <- data.frame(
    study = c("IPD (this review)", "VITAL/CL08 vs LAL-1-NH01"),
    group = c("IPD", "VITAL + historical control"),
    a = c(w_e["events"], v31$events), n1 = c(w_e["n"], v31$n),
    c = c(w_c["events"], h62$events), n2 = c(w_c["n"], h62$n))
  m_e <- ipd_prop(WD, "mortality", 1); m_c <- ipd_prop(WD, "mortality", 0)
  f68e <- pick("R068", "mortality", "ERT"); f68c <- pick("R068", "mortality", "Untreated")
  F6 <- data.frame(
    study = c("IPD, Wolman disease (this review)", "French real-world, Lacaille 2026"),
    group = c("IPD", "French real-world cohort"),
    a = c(m_e["events"], f68e$events), n1 = c(m_e["n"], f68e$n),
    c = c(m_c["events"], f68c$events), n2 = c(m_c["n"], f68c$n))

  or_forest <- function(F, title, xlab, file, invert = FALSE) {
    es <- escalc_OR(F$a, F$n1 - F$a, F$c, F$n2 - F$c)
    r  <- rma_re(es$yi, es$vi, "REML")
    F$yi <- es$yi; F$vi <- es$vi
    F$or <- exp(es$yi); F$lcl <- exp(es$yi - 1.96 * sqrt(es$vi)); F$ucl <- exp(es$yi + 1.96 * sqrt(es$vi))
    F$w  <- (1 / (es$vi + r$tau2)) / sum(1 / (es$vi + r$tau2)) * 100
    list(F = F, r = r)
  }
  mk <- function(F, title, xlab, file, id) {
    o <- or_forest(F, title, xlab, file)
    FOREST[[length(FOREST) + 1]] <<- cbind(forest = id, o$F)
    RES[[length(RES) + 1]] <<- data.frame(forest = id, measure = "OR", k = o$r$k,
      estimate = exp(o$r$b), lcl = exp(o$r$ci.lb), ucl = exp(o$r$ci.ub),
      p = o$r$pval, tau2 = o$r$tau2, I2 = o$r$I2, Q = o$r$Q, Qp = o$r$Qp, method = o$r$method)
    o
  }
  o1 <- mk(F1, "CESD: ALT normalization, ERT versus control",
           "odds ratio (>1 favours ERT)", "M6_forest_ALT_norm", "CESD_ALT_normalization")
  o2 <- mk(F2, "CESD: AST normalization, ERT versus control",
           "odds ratio (>1 favours ERT)", "M6_forest_AST_norm", "CESD_AST_normalization")
  o5 <- mk(F5, "Wolman disease: survival to 12 months, ERT versus untreated",
           "odds ratio (>1 favours ERT)", "M6_forest_WD_survival", "WD_12mo_survival")
  o6 <- mk(F6, "Wolman disease: mortality, ERT versus untreated",
           "odds ratio (<1 favours ERT)", "M6_forest_mortality", "WD_mortality")

  md_block <- function(v, r27, lab, id) {
    ie <- ipd_mean(CESD, v, 1); ic <- ipd_mean(CESD, v, 0)
    e <- pick("R027", r27, "ERT", "week 20", need = "mean"); c2 <- pick("R027", r27, "Placebo", "week 20", need = "mean")
    m1 <- c(ie["mean"], e$mean); s1 <- c(ie["sd"], e$sd); n1 <- c(ie["n"], e$n)
    m2 <- c(ic["mean"], c2$mean); s2 <- c(ic["sd"], c2$sd); n2 <- c(ic["n"], c2$n)
    es <- escalc_MD(m1, s1, n1, m2, s2, n2)
    r  <- rma_re(es$yi, es$vi, "REML")
    F <- data.frame(forest = id, study = c("IPD (this review)", "ARISE, Burton 2015"),
                    group = c("IPD", "ARISE / LAL-CL02"),
                    m1 = m1, sd1 = s1, n1 = n1, m2 = m2, sd2 = s2, n2 = n2,
                    yi = es$yi, vi = es$vi,
                    lcl = es$yi - 1.96 * sqrt(es$vi), ucl = es$yi + 1.96 * sqrt(es$vi))
    RES[[length(RES) + 1]] <<- data.frame(forest = id, measure = "MD (percentage points)", k = r$k,
      estimate = r$b, lcl = r$ci.lb, ucl = r$ci.ub, p = r$pval, tau2 = r$tau2,
      I2 = r$I2, Q = r$Q, Qp = r$Qp, method = r$method)
    list(F = F, r = r, lab = lab)
  }
  o3 <- md_block("LDL_pct", "LDL-C percent change", "CESD: LDL-C percent change", "CESD_LDL_pct_MD")
  o4 <- md_block("TG_pct",  "TG percent change",   "CESD: triglyceride percent change", "CESD_TG_pct_MD")

  R <- do.call(rbind, RES)
  write.csv(R, lald_out("M6_aggregate_results.csv"), row.names = FALSE)
  FD <- do.call(rbind, lapply(FOREST, function(x) x[, c("forest","study","group","a","n1","c","n2","or","lcl","ucl","w")]))
  write.csv(FD, lald_out("M6_forest_data.csv"), row.names = FALSE)
  write.csv(rbind(o3$F, o4$F), lald_out("M6_forest_data_MD.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  suppressPackageStartupMessages(library(survival))
  db <- readRDS(lald_cache("M0_db_derived.rds"))

  db$any_tx <- with(db, ifelse(is.na(liver_tx) & is.na(hsct), NA,
                        as.integer((!is.na(liver_tx) & liver_tx == 1) | (!is.na(hsct) & hsct == 1))))
  db$modality <- with(db, ifelse(is.na(ert_received) | is.na(any_tx), NA,
    ifelse(ert_received == 0 & any_tx == 0, "none",
    ifelse(ert_received == 1 & any_tx == 0, "ERT",
    ifelse(ert_received == 0 & any_tx == 1, "transplant", "both")))))
  LEV <- c("none", "ERT", "transplant", "both")
  db$modality <- factor(db$modality, levels = LEV)
  saveRDS(db, lald_cache("M8_db_modality.rds"))

  comp <- list()
  for (ph in c("CESD", "WD")) {
    d <- db[db$phenotype %in% ph, ]
    for (m in LEV) {
      s <- d[d$modality %in% m, ]
      o <- s[!is.na(s$mortality), ]
      comp[[length(comp) + 1]] <- data.frame(
        phenotype = ph, modality = m, n = nrow(s),
        n_mortality_observed = nrow(o), deaths = sum(o$mortality == 1),
        pct_died = if (nrow(o)) round(100 * mean(o$mortality == 1), 1) else NA,
        n_liver_tx = sum(s$liver_tx %in% 1), n_hsct = sum(s$hsct %in% 1),
        median_age_dx = if (any(!is.na(s$age_dx_y))) round(median(s$age_dx_y, na.rm = TRUE), 2) else NA,
        median_fu = if (any(!is.na(s$fu_years))) round(median(s$fu_years, na.rm = TRUE), 1) else NA,
        modelled = m != "both")
    }
    comp[[length(comp) + 1]] <- data.frame(phenotype = ph, modality = "(flags incomplete)",
      n = sum(is.na(d$modality)), n_mortality_observed = NA, deaths = NA, pct_died = NA,
      n_liver_tx = NA, n_hsct = NA, median_age_dx = NA, median_fu = NA, modelled = FALSE)
  }
  COMP <- do.call(rbind, comp)
  write.csv(COMP, lald_out("M8_modality_composition.csv"), row.names = FALSE)

  fit_mod <- function(d, adjust) {
    d <- d[d$modality %in% c("none", "ERT", "transplant") & !is.na(d$mortality), ]
    d$modality <- droplevels(factor(d$modality, levels = c("none", "ERT", "transplant")))
    if (nlevels(d$modality) < 2) return(NULL)
    d$m_ERT <- as.integer(d$modality == "ERT"); d$m_TX <- as.integer(d$modality == "transplant")
    vars <- c("mortality", "m_ERT", "m_TX", if (adjust) c("age_dx_y", "severity_index"))
    d <- d[stats::complete.cases(d[, vars]), vars]
    if (nrow(d) < 12) return(NULL)
    keep <- c("m_ERT", "m_TX")
    keep <- keep[sapply(keep, function(v) length(unique(d[[v]])) > 1)]
    if (adjust) keep <- c(keep, c("age_dx_y", "severity_index")[
      sapply(c("age_dx_y", "severity_index"), function(v) length(unique(d[[v]])) > 1)])
    f <- as.formula(paste("mortality ~", paste(keep, collapse = " + ")))
    fit <- tryCatch(lald_firth(f, d), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    list(fit = fit, n = nrow(d), d = d)
  }
  res <- list()
  for (ph in c("CESD", "WD")) {
    d0 <- db[db$phenotype %in% ph, ]
    for (adj in c(FALSE, TRUE)) {
      r <- fit_mod(d0, adj)
      if (is.null(r)) next
      for (term in c("m_ERT", "m_TX")) {
        j <- match(term, r$fit$names); if (is.na(j)) next
        lvl <- ifelse(term == "m_ERT", "ERT", "transplant")
        s <- r$d[r$d[[term]] == 1, ]
        res[[length(res) + 1]] <- data.frame(
          phenotype = ph, contrast = paste(lvl, "vs none"),
          model = ifelse(adj, "Firth adjusted", "Firth unadjusted"), n_model = r$n,
          n_group = nrow(s), deaths_group = sum(s$mortality == 1),
          OR = exp(r$fit$coef[j]), lcl = exp(r$fit$ci[j, 1]), ucl = exp(r$fit$ci[j, 2]),
          p = r$fit$p[j], engine = r$fit$engine)
      }
      dd <- d0[d0$modality %in% c("ERT", "transplant") & !is.na(d0$mortality), ]
      dd$m_ERT <- as.integer(dd$modality == "ERT")
      vars <- c("mortality", "m_ERT", if (adj) c("age_dx_y", "severity_index"))
      dd <- dd[stats::complete.cases(dd[, vars]), vars]
      if (nrow(dd) >= 12 && length(unique(dd$mortality)) > 1 && length(unique(dd$m_ERT)) > 1) {
        keep <- "m_ERT"
        if (adj) keep <- c(keep, c("age_dx_y", "severity_index")[
          sapply(c("age_dx_y", "severity_index"), function(v) length(unique(dd[[v]])) > 1)])
        fp <- tryCatch(lald_firth(as.formula(paste("mortality ~", paste(keep, collapse = " + "))), dd),
                       error = function(e) NULL)
        if (!is.null(fp)) {
          j <- match("m_ERT", fp$names)
          res[[length(res) + 1]] <- data.frame(
            phenotype = ph, contrast = "ERT vs transplant",
            model = ifelse(adj, "Firth adjusted", "Firth unadjusted"), n_model = nrow(dd),
            n_group = sum(dd$m_ERT == 1), deaths_group = sum(dd$mortality[dd$m_ERT == 1] == 1),
            OR = exp(fp$coef[j]), lcl = exp(fp$ci[j, 1]), ucl = exp(fp$ci[j, 2]),
            p = fp$p[j], engine = fp$engine)
        }
      }
    }
  }
  R <- do.call(rbind, res)
  write.csv(R, lald_out("M8_modality_results.csv"), row.names = FALSE)

  km_rows <- list()
  for (ph in c("CESD", "WD")) {
    d <- db[db$phenotype %in% ph & !is.na(db$fu_years) & !is.na(db$mortality) &
            db$modality %in% c("none", "ERT", "transplant"), ]
    if (nrow(d) < 12) next
    d$fu_years[d$fu_years == 0] <- 0.01
    d$modality <- droplevels(d$modality)
    lr <- tryCatch(survdiff(Surv(fu_years, mortality) ~ modality, data = d), error = function(e) NULL)
    km <- survfit(Surv(fu_years, mortality) ~ modality, data = d)
    tau <- as.numeric(quantile(d$fu_years, 0.8))
    rm1 <- function(g) { s <- d[d$modality == g, ]
      if (nrow(s) < 2) return(NA)
      k <- survfit(Surv(fu_years, mortality) ~ 1, data = s)
      su <- summary(k, times = sort(unique(pmin(k$time, tau))), extend = TRUE)
      tt <- c(0, su$time[su$time <= tau], tau); ss <- c(1, su$surv[su$time <= tau])
      ss <- c(ss, ss[length(ss)]); sum(diff(tt) * head(ss, -1)) }
    for (g in levels(d$modality))
      km_rows[[length(km_rows) + 1]] <- data.frame(phenotype = ph, modality = g,
        n = sum(d$modality == g), deaths = sum(d$mortality[d$modality == g] == 1),
        median_fu = round(median(d$fu_years[d$modality == g]), 2),
        tau = round(tau, 2), rmst = round(rm1(g), 2),
        logrank_p = if (!is.null(lr)) pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE) else NA)
  }
  KM <- do.call(rbind, km_rows)
  write.csv(KM, lald_out("M8_survival_by_modality.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db  <- readRDS(lald_cache("M0_db_derived.rds"))
  sub <- read.csv(lald_in("ipd"), check.names = FALSE, na.strings = c("", "NA"))
  gv <- function(r, c) { v <- r[[c]]; if (is.null(v) || is.na(v)) "" else trimws(as.character(v)) }

  COLS <- list(c("lal_after_liver_tx", "after liver transplant"),
               c("lal_after_hsct", "after HSCT"))
  rows <- list()
  for (cc in COLS) {
    v <- cc[1]; lab <- cc[2]
    if (!(v %in% names(sub))) next
    s <- sub[!is.na(sub[[v]]) & trimws(sub[[v]]) != "", ]
    for (i in seq_len(nrow(s))) {
      r <- s[i, ]
      rows[[length(rows) + 1]] <- data.frame(
        graft = lab, source_R_ID = gv(r, "source_R_ID"), row_in_source = gv(r, "row_in_source"),
        phenotype = gv(r, "phenotype"),
        value_as_stored = gv(r, v),
        value_kind = if (grepl("^[01]$", gv(r, v))) "flag_only" else "measurement",
        baseline_LAL = substr(gv(r, "lal_baseline"), 1, 90),
        liver_tx = gv(r, "liver_tx"), hsct = gv(r, "hsct"),
        ert = gv(r, "ert_received"), mortality = gv(r, "mortality"),
        artefact_level = gv(r, "artefact_level"), field = v)
    }
  }
  E <- do.call(rbind, rows)
  write.csv(E, lald_out("M9_enzyme_after_transplant.csv"), row.names = FALSE)
  EM <- E[E$value_kind == "measurement", ]
  EF <- E[E$value_kind == "flag_only", ]

  olt  <- db[db$liver_tx %in% 1, ]
  hsct <- db[db$hsct %in% 1, ]
  rec  <- olt[!is.na(olt$recurrence_post_liver_tx), ]
  loss <- sub[!is.na(sub$graft_loss_after_hsct), ]
  tt_rec  <- sub$time_to_recurrence_liver_tx_y
  tt_loss <- sub$time_to_graft_loss_hsct_y
  O <- rbind(
    data.frame(graft = "liver transplant", quantity = "recipients", events = NA, n = nrow(olt), detail = ""),
    data.frame(graft = "liver transplant", quantity = "recurrence in the graft",
               events = sum(rec$recurrence_post_liver_tx == 1), n = nrow(rec),
               detail = "of recipients with a recurrence field"),
    data.frame(graft = "liver transplant", quantity = "time to recurrence, years",
               events = NA, n = sum(!is.na(tt_rec)),
               detail = paste(sort(tt_rec[!is.na(tt_rec)]), collapse = ", ")),
    data.frame(graft = "liver transplant", quantity = "also received ERT",
               events = sum(olt$ert_received %in% 1), n = nrow(olt), detail = ""),
    data.frame(graft = "HSCT", quantity = "recipients", events = NA, n = nrow(hsct), detail = ""),
    data.frame(graft = "HSCT", quantity = "graft loss or death",
               events = sum(loss$graft_loss_after_hsct == 1), n = nrow(loss),
               detail = "of recipients with the field recorded"),
    data.frame(graft = "HSCT", quantity = "time to graft loss, years",
               events = NA, n = sum(!is.na(tt_loss)),
               detail = paste(sort(tt_loss[!is.na(tt_loss)]), collapse = ", ")),
    data.frame(graft = "HSCT", quantity = "also received ERT",
               events = sum(hsct$ert_received %in% 1), n = nrow(hsct), detail = ""))
  write.csv(O, lald_out("M9_transplant_outcomes.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db <- if (file.exists(lald_cache("M8_db_modality.rds"))) readRDS(lald_cache("M8_db_modality.rds")) else {
    d <- readRDS(lald_cache("M0_db_derived.rds"))
    d$any_tx <- with(d, ifelse(is.na(liver_tx) & is.na(hsct), NA,
                        as.integer((!is.na(liver_tx) & liver_tx == 1) | (!is.na(hsct) & hsct == 1))))
    d$modality <- factor(with(d, ifelse(is.na(ert_received) | is.na(any_tx), NA,
      ifelse(ert_received == 0 & any_tx == 0, "none",
      ifelse(ert_received == 1 & any_tx == 0, "ERT",
      ifelse(ert_received == 0 & any_tx == 1, "transplant", "both"))))),
      levels = c("none", "ERT", "transplant", "both")); d }
  pctf <- function(a, b) { d <- (b - a) / a * 100; d[!is.finite(d)] <- NA; d }
  db$HDL_pct <- pctf(db$HDL_dx_n, db$HDL_post_n)

  BL <- c(cirrhosis_bl = "Cirrhosis", steatosis_bl = "Hepatic steatosis",
          transaminases_raised_bl = "Raised transaminases", dyslipidemia_bl = "Dyslipidaemia",
          hepatosplenomegaly_bl = "Hepatosplenomegaly", ascvd_bl = "Atherosclerotic CVD")
  brows <- list()
  for (ph in c("CESD", "WD")) for (v in names(BL)) {
    d <- db[db$phenotype %in% ph & !is.na(db[[v]]), ]
    brows[[length(brows) + 1]] <- data.frame(phenotype = ph, feature = BL[[v]], variable = v,
      n_recorded = nrow(d), n_present = sum(d[[v]] == 1),
      pct_present = if (nrow(d)) round(100 * mean(d[[v]] == 1), 1) else NA,
      pct_of_phenotype_recorded = round(100 * nrow(d) / sum(db$phenotype %in% ph), 1))
  }
  BLD <- do.call(rbind, brows)
  write.csv(BLD, lald_out("M10_baseline_burden.csv"), row.names = FALSE)

  OUT <- list(c("ALT_pct", "ALT percent change"), c("AST_pct", "AST percent change"),
              c("LDL_pct", "LDL-C percent change"), c("HDL_pct", "HDL-C percent change"),
              c("TG_pct", "TG percent change"), c("TC_pct", "Total cholesterol percent change"))
  rrows <- list()
  C <- db[db$phenotype %in% "CESD", ]
  for (o in OUT) {
    v <- o[1]; lab <- o[2]
    for (m in c("none", "ERT", "transplant", "both")) {
      s <- C[C$modality %in% m & !is.na(C[[v]]), ]
      rrows[[length(rrows) + 1]] <- data.frame(outcome = lab, variable = v, modality = m,
        n = nrow(s),
        estimate = if (nrow(s)) round(median(s[[v]]), 1) else NA,
        estimator = "median, lo/hi = IQR",
        lo = if (nrow(s) >= 2) round(quantile(s[[v]], .25), 1) else NA,
        hi = if (nrow(s) >= 2) round(quantile(s[[v]], .75), 1) else NA,
        p = NA_real_,
        modelled = m %in% c("none", "ERT"))
    }
    d <- C[C$modality %in% c("none", "ERT") & !is.na(C[[v]]), ]
    d$ert <- as.integer(d$modality == "ERT")
    est <- NA; lo <- NA; hi <- NA; pv <- NA; nn <- nrow(d)
    if (sum(d$ert == 1) >= 3 && sum(d$ert == 0) >= 3) {
      ts <- theil_sen(d$ert, d[[v]])
      bc <- tryCatch(bca_ci(d[, c(v, "ert")], function(dat, i) {
        s <- dat[i, ]; a <- s[[1]][s[[2]] == 1]; b <- s[[1]][s[[2]] == 0]
        if (length(a) < 2 || length(b) < 2) return(NA_real_); median(outer(a, b, "-")) }, R = 2000),
        error = function(e) NULL)
      est <- ts$slope; pv <- ts$p
      if (!is.null(bc)) { lo <- bc$lcl; hi <- bc$ucl }
    }
    rrows[[length(rrows) + 1]] <- data.frame(outcome = lab, variable = v, modality = "ERT vs none",
      n = nn, estimate = round(est, 1),
      estimator = "Theil-Sen shift, lo/hi = BCa 95% CI, p unadjusted (exploratory)",
      lo = round(lo, 1), hi = round(hi, 1), p = pv, modelled = TRUE)
  }
  RES <- do.call(rbind, rrows)
  write.csv(RES, lald_out("M10_response_by_modality.csv"), row.names = FALSE)

  code2 <- function(x) ifelse(is.na(x), NA, as.integer(x == 2))
  CODED <- list(c("steatosis_outcome_coded", "Steatosis improved"),
                c("fibrosis_outcome_coded", "Fibrosis improved"),
                c("hepatomegaly_outcome_coded", "Hepatomegaly improved"),
                c("lipid_outcome_coded", "Lipid response favourable"))
  crows <- list()
  for (cc in CODED) {
    v <- cc[1]; lab <- cc[2]
    for (m in c("none", "ERT", "transplant", "both")) {
      s <- C[C$modality %in% m & !is.na(C[[v]]), ]
      crows[[length(crows) + 1]] <- data.frame(outcome = lab, variable = v, modality = m,
        n = nrow(s), improved = sum(code2(s[[v]]) == 1, na.rm = TRUE),
        pct = if (nrow(s)) round(100 * mean(code2(s[[v]]) == 1, na.rm = TRUE), 1) else NA)
    }
  }
  COD <- do.call(rbind, crows)
  write.csv(COD, lald_out("M10_coded_outcomes_by_modality.csv"), row.names = FALSE)

  INK <- "#17181c"; INK2 <- "#4e5158"; MUTE <- "#82858c"
  COLM <- c(none = "#8a8f98", ERT = "#2f5d50", transplant = "#9a7a1f", both = "#a8452f")
})

local({
  set.seed(42)
  db <- read.csv(lald_in("ipd"), na.strings = c("", "NA"))
  BELOW <- c("abstract", "title", "partial")
  src <- unique(db[, c("source_R_ID", "artefact_level")])
  stopifnot(!anyDuplicated(src$source_R_ID))
  n_src <- nrow(src); n_part <- nrow(db)
  below <- sum(src$artefact_level %in% BELOW)
  part_below <- sum(db$artefact_level %in% BELOW)
  S <- data.frame(
    quantity = c("contributing sources", "sources read in full text", "sources extracted below full text",
                 "participants from sources below full text"),
    n = c(n_src, n_src - below, below, part_below),
    denominator = c(n_src, n_src, n_src, n_part))
  S$pct <- round(100 * S$n / S$denominator, 1)
  write.csv(S, lald_out("M12_availability_summary.csv"), row.names = FALSE)
  OUT <- list("Mortality" = c("mortality", NA), "Liver transplantation" = c("liver_tx", NA),
              "ALT percent change" = c("ALT_dx_n", "ALT_post_n"), "AST percent change" = c("AST_dx_n", "AST_post_n"),
              "LDL-C percent change" = c("LDL_dx_n", "LDL_post_n"), "HDL-C percent change" = c("HDL_dx_n", "HDL_post_n"),
              "TG percent change" = c("TG_dx_n", "TG_post_n"), "Total cholesterol % change" = c("TC_dx_n", "TC_post_n"),
              "Steatosis regression" = c("steatosis_outcome_coded", NA))
  B <- do.call(rbind, lapply(names(OUT), function(lab) {
    v <- OUT[[lab]]
    k <- !is.na(db[[v[1]]]); if (!is.na(v[2])) k <- k & !is.na(db[[v[2]]])
    e <- db[k, ]
    lo <- e[e$artefact_level %in% BELOW, ]
    pct <- round(100 * nrow(lo) / nrow(e), 1)
    data.frame(outcome = lab, n_participants = nrow(e), n_sources = length(unique(e$source_R_ID)),
               median_participants_per_source = round(median(as.vector(table(e$source_R_ID))), 1),
               n_participants_below_fulltext = nrow(lo), pct_participants_below_fulltext = pct,
               n_sources_below_fulltext = length(unique(lo$source_R_ID)),
               exposure = if (pct >= 15) "high" else if (pct >= 5) "moderate" else "low")
  }))
  B <- B[order(-B$pct_participants_below_fulltext, method = "radix"), ]
  write.csv(B, lald_out("M12_availability_by_outcome.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db  <- readRDS(lald_cache("M0_db_derived.rds"))
  ipd <- read.csv(lald_in("ipd"), na.strings = c("", "NA"))
  agg <- read.csv(lald_in("aggregate"), check.names = FALSE, na.strings = c("", "NA"))
  lum <- agg[agg$R_ID == "R166", ]
  row <- function(arm, outcome) { r <- lum[lum$arm == arm & lum$outcome == outcome, ]; stopifnot(nrow(r) == 1); r }

  num <- function(x) suppressWarnings(as.numeric(x))
  db$modality <- with(db, ifelse(num(ert_received) %in% 1 & (num(liver_tx) %in% 1 | num(hsct) %in% 1), "both",
                          ifelse(num(ert_received) %in% 1, "ERT",
                          ifelse(num(liver_tx) %in% 1 | num(hsct) %in% 1, "transplant",
                          ifelse(num(ert_received) %in% 0 & num(liver_tx) %in% 0 & num(hsct) %in% 0, "none", NA)))))
  ipd_or <- function(ph, a_mod, b_mod) {
    d <- db[db$phenotype %in% ph & !is.na(db$mortality) & db$modality %in% c(a_mod, b_mod), ]
    if (!nrow(d)) return(NULL)
    d$x <- as.integer(d$modality == a_mod)
    if (length(unique(d$x)) < 2) return(NULL)
    f <- lald_firth(mortality ~ x, d)
    out <- unname(c(nrow(d), exp(f$coef[2]), exp(f$ci[2, 1]), exp(f$ci[2, 2]), f$p[2]))
    names(out) <- c("n", "or", "lcl", "ucl", "p")
    out
  }
  ext <- list()
  add <- function(src, comparison, meas, est, lcl, ucl, p, n, note)
    ext[[length(ext) + 1]] <<- data.frame(source = src, comparison = comparison, measure = meas,
      estimate = est, lcl = lcl, ucl = ucl, p = p, n = n, note = note)

  w  <- row("HCT with pre-transplant ERT", "survival, counts")
  wo <- row("HCT without pre-transplant ERT", "survival, counts")
  a <- w$n - w$events; b <- w$events; cc <- wo$n - wo$events; d <- wo$events
  or <- (a * d) / (b * cc)
  se <- sqrt(1/a + 1/b + 1/cc + 1/d)
  add("R166 Lum 2023 (EBMT, external)", "HCT + pre-transplant ERT vs HCT alone", "OR for death",
      round(or, 4), round(exp(log(or) - 1.96 * se), 4), round(exp(log(or) + 1.96 * se), 4),
      num(row("HCT without pre-transplant ERT", "overall survival, 5 y")$p_value), w$n + wo$n,
      sprintf("from the verified counts %d/%d and %d/%d survived; the paper's text says 5-y OS 80%% and its Figure 1b says 90%%, the counts say 80%%",
              w$events, w$n, wo$events, wo$n))
  r <- ipd_or("WD", "both", "transplant")
  if (!is.null(r)) add("This review (IPD)", "ERT + transplant vs transplant alone", "OR for death",
                       round(r["or"], 4), round(r["lcl"], 4), round(r["ucl"], 4), r["p"], r["n"],
                       "the 'both' cell was declared unmodellable in advance (n = 2); shown for completeness only")
  r <- ipd_or("WD", "ERT", "transplant")
  if (!is.null(r)) add("This review (IPD)", "ERT vs transplant", "OR for death",
                       round(r["or"], 4), round(r["lcl"], 4), round(r["ucl"], 4), r["p"], r["n"], "primary Q1 contrast")
  gl <- ipd$graft_loss_after_hsct[!is.na(ipd$graft_loss_after_hsct)]
  add("This review (IPD)", "HSCT graft loss or death", "proportion",
      round(100 * sum(gl == 1) / length(gl), 1), NA, NA, NA, length(gl),
      sprintf("%d/%d recipients with the field recorded", sum(gl == 1), length(gl)))
  gf <- row("HCT, whole cohort", "graft failure")
  add("R166 Lum 2023 (EBMT, external)", "HCT graft failure", "proportion",
      round(100 * gf$events / gf$n, 1), NA, NA, NA, gf$n,
      sprintf("%d/%d; all six died despite three second transplants", gf$events, gf$n))
  mo <- row("HCT, whole cohort", "mortality, counts")
  add("R166 Lum 2023 (EBMT, external)", "HCT mortality, whole cohort", "proportion",
      round(100 * mo$events / mo$n, 1), NA, NA, NA, mo$n,
      sprintf("%d/%d; 8 transplant-related, 6 disease-related, 2 unknown", mo$events, mo$n))
  E <- do.call(rbind, ext)
  write.csv(E, lald_out("M15_external_comparison.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db  <- readRDS(lald_cache("M0_db_derived.rds"))
  agg <- read.csv(lald_in("aggregate"), check.names = FALSE)

  fp <- function(a, b, c, d) {
    m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    if (any(is.na(m)) || any(m < 0)) return(NA_real_)
    tryCatch(fisher.test(m)$p.value, error = function(e) NA_real_)
  }

  fragility <- function(a, n1, c, n2, alpha = 0.05) {
    b <- n1 - a; d <- n2 - c
    p0 <- fp(a, b, c, d)
    if (is.na(p0)) return(list(p0 = NA, fi = NA, dir = NA, arm = NA))
    sig <- p0 < alpha
    best <- NA_integer_; barm <- NA_character_; bdir <- NA_character_
    for (arm in 1:2) {
      N <- if (arm == 1) n1 else n2
      for (k in 1:N) {
        for (s in c(+1, -1)) {
          if (arm == 1) { aa <- a + s * k; if (aa < 0 || aa > n1) next
                          p1 <- fp(aa, n1 - aa, c, d) }
          else          { cc <- c + s * k; if (cc < 0 || cc > n2) next
                          p1 <- fp(a, b, cc, n2 - cc) }
          if (is.na(p1)) next
          flipped <- if (sig) (p1 >= alpha) else (p1 < alpha)
          if (flipped && (is.na(best) || k < best)) {
            best <- k; barm <- paste0("arm", arm)
            bdir <- if (s > 0) "add events" else "remove events"
          }
        }
        if (!is.na(best) && best <= k) break
      }
    }
    list(p0 = p0, fi = best, dir = bdir, arm = barm)
  }

  rows <- list()
  add <- function(source, outcome, phenotype, comparison, a, n1, c, n2) {
    f <- fragility(a, n1, c, n2)
    if (is.na(f$p0)) return(invisible())
    sig <- f$p0 < 0.05
    smaller <- min(n1, n2)
    rows[[length(rows) + 1]] <<- data.frame(
      source = source, outcome = outcome, phenotype = phenotype, comparison = comparison,
      events_arm1 = a, n_arm1 = n1, events_arm2 = c, n_arm2 = n2,
      fisher_p = signif(f$p0, 4),
      significant = sig,
      index = if (sig) "fragility index" else "reverse fragility index",
      value = f$fi, direction = f$dir, arm_changed = f$arm,
      fi_over_smaller_arm = if (is.na(f$fi)) NA else round(f$fi / smaller, 3),
      note = "")
  }

  num <- function(x) suppressWarnings(as.numeric(x))
  cnt <- function(d, out, expo, ev = 1) sum(num(d[[out]]) %in% ev & num(d[[expo]]) %in% 1, na.rm = TRUE)

  for (ph in c("CESD", "WD")) {
    d <- db[db$phenotype %in% ph & !is.na(db$mortality) & !is.na(db$ert_received), ]
    if (nrow(d)) add("IPD (this review)", "Mortality", ph, "ERT vs no ERT",
        sum(d$mortality == 1 & d$ert_received == 1), sum(d$ert_received == 1),
        sum(d$mortality == 1 & d$ert_received == 0), sum(d$ert_received == 0))
  }
  db$modality <- with(db, ifelse(num(ert_received) %in% 1 & (num(liver_tx) %in% 1 | num(hsct) %in% 1), "both",
                          ifelse(num(ert_received) %in% 1, "ERT",
                          ifelse(num(liver_tx) %in% 1 | num(hsct) %in% 1, "transplant",
                          ifelse(num(ert_received) %in% 0 & num(liver_tx) %in% 0 & num(hsct) %in% 0, "none", NA)))))
  for (ph in c("CESD", "WD")) for (pair in list(c("ERT", "none"), c("transplant", "none"), c("ERT", "transplant"))) {
    d <- db[db$phenotype %in% ph & !is.na(db$mortality) & db$modality %in% pair, ]
    if (!nrow(d) || length(unique(d$modality)) < 2) next
    a1 <- d$modality == pair[1]; a2 <- d$modality == pair[2]
    add("IPD (this review)", "Mortality", ph, paste(pair[1], "vs", pair[2]),
        sum(d$mortality[a1] == 1), sum(a1), sum(d$mortality[a2] == 1), sum(a2))
  }
  for (o in c("ALT_norm", "AST_norm", "steatosis_outcome_coded")) {
    if (!o %in% names(db)) next
    d <- db[db$phenotype %in% "CESD" & !is.na(db[[o]]) & !is.na(db$ert_received), ]
    if (!nrow(d)) next
    ev <- if (grepl("coded", o)) 2 else 1
    lab <- c(ALT_norm = "ALT normalization", AST_norm = "AST normalization",
             steatosis_outcome_coded = "Steatosis regression")[o]
    add("IPD (this review)", lab, "CESD", "ERT vs no ERT",
        sum(num(d[[o]]) == ev & d$ert_received == 1), sum(d$ert_received == 1),
        sum(num(d[[o]]) == ev & d$ert_received == 0), sum(d$ert_received == 0))
  }

  A <- agg[agg$tier %in% "A" & !is.na(agg$events) & !is.na(agg$n), ]
  pick <- function(rid, pat, arm) {
    s <- A[A$R_ID == rid & grepl(pat, A$outcome, ignore.case = TRUE) &
           grepl(arm, A$arm, ignore.case = TRUE) & !is.na(A$events), ]
    if (!nrow(s)) return(NULL)
    c(ev = as.numeric(s$events[1]), n = as.numeric(s$n[1]))
  }
  e <- pick("R166", "survival, counts", "with pre-transplant ERT")
  u <- pick("R166", "survival, counts", "without pre-transplant ERT")
  if (!is.null(e) && !is.null(u))
    add("R166 Lum 2023 (EBMT)", "Death after HCT", "WD", "HCT + ERT vs HCT alone",
        e["n"] - e["ev"], e["n"], u["n"] - u["ev"], u["n"])
  e <- pick("R024", "12-month survival", "^ERT"); u <- pick("R024", "12-month survival", "Untreated")
  if (!is.null(e) && !is.null(u))
    add("R024 Jones 2017 / VITAL", "12-month survival", "WD", "ERT vs untreated historical control",
        e["ev"], e["n"], u["ev"], u["n"])
  e <- pick("R068", "mortality", "^ERT"); u <- pick("R068", "mortality", "Untreated")
  if (!is.null(e) && !is.null(u))
    add("R068 Lacaille 2026 (French)", "Mortality", "WD", "ERT vs untreated",
        e["ev"], e["n"], u["ev"], u["n"])

  R <- do.call(rbind, rows)
  R$note <- ifelse(!R$significant, "not significant; the value is the reverse fragility index",
             ifelse(R$value <= 2, "one or two patients carry this result",
             ifelse(R$value <= 5, "modest fragility", "robust to several reversals")))
  write.csv(R, lald_out("M16_fragility.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db <- readRDS(lald_cache("M0_db_derived.rds"))

  ev <- function(rr) {
    if (is.na(rr) || rr <= 0) return(NA_real_)
    r <- if (rr < 1) 1 / rr else rr
    r + sqrt(r * (r - 1))
  }
  ev_ci <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return(NA_real_)
    if (lo <= 1 && hi >= 1) return(1)
    ev(if (hi < 1) hi else lo)
  }
  rr_ci <- function(a, n1, c, n2) {
    if (a == 0 || c == 0) { a <- a + 0.5; c <- c + 0.5; n1 <- n1 + 1; n2 <- n2 + 1 }
    r1 <- a / n1; r2 <- c / n2; rr <- r1 / r2
    se <- sqrt(1/a - 1/n1 + 1/c - 1/n2)
    c(rr = rr, lo = exp(log(rr) - 1.96 * se), hi = exp(log(rr) + 1.96 * se))
  }

  rows <- list()
  add <- function(source, outcome, phenotype, comparison, a, n1, c, n2, note = "") {
    if (n1 == 0 || n2 == 0) return(invisible())
    zero <- (a == 0 || c == 0)
    k <- rr_ci(a, n1, c, n2)
    rows[[length(rows) + 1]] <<- data.frame(
      source = source, outcome = outcome, phenotype = phenotype, comparison = comparison,
      events_exposed = a, n_exposed = n1, risk_exposed = round(a / n1, 4),
      events_unexposed = c, n_unexposed = n2, risk_unexposed = round(c / n2, 4),
      risk_ratio = round(k["rr"], 4), rr_lcl = round(k["lo"], 4), rr_ucl = round(k["hi"], 4),
      e_value_point = round(ev(k["rr"]), 2),
      e_value_ci = round(ev_ci(k["lo"], k["hi"]), 2),
      continuity_correction = zero,
      note = note)
  }

  num <- function(x) suppressWarnings(as.numeric(x))
  for (ph in c("CESD", "WD")) {
    d <- db[db$phenotype %in% ph & !is.na(db$mortality) & !is.na(db$ert_received), ]
    if (nrow(d)) add("IPD (this review)", "Mortality", ph, "ERT vs no ERT",
        sum(d$mortality == 1 & d$ert_received == 1), sum(d$ert_received == 1),
        sum(d$mortality == 1 & d$ert_received == 0), sum(d$ert_received == 0),
        if (sum(d$mortality == 1 & d$ert_received == 1) == 0) "zero events in the ERT arm; 0.5 continuity correction applied, so the point E-value is a lower bound" else "")
  }
  db$modality <- with(db, ifelse(num(ert_received) %in% 1 & (num(liver_tx) %in% 1 | num(hsct) %in% 1), "both",
                          ifelse(num(ert_received) %in% 1, "ERT",
                          ifelse(num(liver_tx) %in% 1 | num(hsct) %in% 1, "transplant",
                          ifelse(num(ert_received) %in% 0 & num(liver_tx) %in% 0 & num(hsct) %in% 0, "none", NA)))))
  for (ph in c("CESD", "WD")) for (pr in list(c("ERT","none"), c("transplant","none"), c("ERT","transplant"))) {
    d <- db[db$phenotype %in% ph & !is.na(db$mortality) & db$modality %in% pr, ]
    if (!nrow(d) || length(unique(d$modality)) < 2) next
    i1 <- d$modality == pr[1]; i2 <- d$modality == pr[2]
    add("IPD (this review)", "Mortality", ph, paste(pr[1], "vs", pr[2]),
        sum(d$mortality[i1] == 1), sum(i1), sum(d$mortality[i2] == 1), sum(i2))
  }
  for (o in c("ALT_norm", "AST_norm", "steatosis_outcome_coded")) {
    if (!o %in% names(db)) next
    d <- db[db$phenotype %in% "CESD" & !is.na(db[[o]]) & !is.na(db$ert_received), ]
    if (!nrow(d)) next
    e <- if (grepl("coded", o)) 2 else 1
    lab <- c(ALT_norm = "ALT normalization", AST_norm = "AST normalization",
             steatosis_outcome_coded = "Steatosis regression")[o]
    add("IPD (this review)", lab, "CESD", "ERT vs no ERT",
        sum(num(d[[o]]) == e & d$ert_received == 1), sum(d$ert_received == 1),
        sum(num(d[[o]]) == e & d$ert_received == 0), sum(d$ert_received == 0))
  }
  add("R166 Lum 2023 (EBMT, external)", "Death after HCT", "WD", "HCT + ERT vs HCT alone",
      1, 5, 15, 19, "external corroboration; not pooled with the IPD")

  E <- do.call(rbind, rows)
  write.csv(E, lald_out("M17_evalue.csv"), row.names = FALSE)
})
