if (!exists("lald_out")) source("01_prepare.R")

local({
  set.seed(42)
  db <- readRDS(lald_cache("M0_db_derived.rds"))

  CONT <- c(age_dx_y = "Age at diagnosis, y", onset_age_y = "Age at onset, y",
            fu_years = "Follow-up, y", TC_dx_n = "Total cholesterol, mg/dL",
            LDL_dx_n = "LDL-C, mg/dL", HDL_dx_n = "HDL-C, mg/dL", TG_dx_n = "Triglycerides, mg/dL",
            ALT_dx_n = "ALT, U/L", AST_dx_n = "AST, U/L", severity_index = "Baseline severity index (0-4)")
  BIN  <- c(dyslipidemia_bl = "Dyslipidaemia", hepatosplenomegaly_bl = "Hepatosplenomegaly",
            steatosis_bl = "Hepatic steatosis", cirrhosis_bl = "Cirrhosis",
            growth_failure_bl = "Growth failure", gi_symptoms_bl = "GI symptoms",
            adrenal_calcification = "Adrenal calcification", ascvd_bl = "ASCVD",
            statin = "Statin", LLT_harmonized = "Any lipid-lowering therapy",
            liver_tx = "Liver transplantation", hsct = "HSCT", mortality = "Died")

  fmt_med <- function(x) { x <- x[!is.na(x)]
    if (!length(x)) return("-")
    sprintf("%.1f [%.1f-%.1f] (n=%d)", median(x), quantile(x, .25), quantile(x, .75), length(x)) }
  fmt_bin <- function(x) { x <- x[!is.na(x)]
    if (!length(x)) return("-")
    sprintf("%d/%d (%.0f%%)", sum(x == 1), length(x), 100 * mean(x == 1)) }

  rows <- list()
  PHEN <- c("CESD", "WD", "unassigned")
  for (ph in PHEN) {
    ALL <- if (ph == "unassigned") db[is.na(db$phenotype), ] else db[db$phenotype %in% ph, ]
    d <- ALL[!is.na(ALL$ert_received), ]
    e <- d[d$ert_received == 1, ]; n <- d[d$ert_received == 0, ]
    rows[[length(rows) + 1]] <- data.frame(phenotype = ph, variable = "n", type = "header",
        overall = as.character(nrow(ALL)),
        ert = as.character(nrow(e)), nonert = as.character(nrow(n)),
        test = "", stat = NA, p = NA)
    rows[[length(rows) + 1]] <- data.frame(phenotype = ph, variable = "ERT status not recorded",
        type = "header", overall = as.character(nrow(ALL) - nrow(d)),
        ert = "", nonert = "", test = "", stat = NA, p = NA)
    for (v in names(CONT)) {
      bm <- brunner_munzel(e[[v]], n[[v]])
      rows[[length(rows) + 1]] <- data.frame(phenotype = ph, variable = CONT[[v]], type = "continuous",
          overall = fmt_med(ALL[[v]]), ert = fmt_med(e[[v]]), nonert = fmt_med(n[[v]]),
          test = "Brunner-Munzel", stat = round(bm$phat, 3), p = bm$p)
    }
    for (v in names(BIN)) {
      a <- sum(e[[v]] == 1, na.rm = TRUE); b <- sum(e[[v]] == 0, na.rm = TRUE)
      cc <- sum(n[[v]] == 1, na.rm = TRUE); dd <- sum(n[[v]] == 0, na.rm = TRUE)
      ft <- if (min(a + b, cc + dd) > 0 && (a + cc) > 0)
        fisher.test(matrix(c(a, b, cc, dd), 2))$p.value else NA
      rows[[length(rows) + 1]] <- data.frame(phenotype = ph, variable = BIN[[v]], type = "binary",
          overall = fmt_bin(ALL[[v]]), ert = fmt_bin(e[[v]]), nonert = fmt_bin(n[[v]]),
          test = "Fisher exact", stat = NA, p = ft)
    }
    fmt_sex <- function(x) if (!sum(!is.na(x))) "-" else
      sprintf("%d/%d (%.0f%%)", sum(x == "male", na.rm = TRUE), sum(!is.na(x)),
              100 * mean(x[!is.na(x)] == "male"))
    rows[[length(rows) + 1]] <- data.frame(phenotype = ph, variable = "Male", type = "binary",
        overall = fmt_sex(ALL$sex),
        ert = sprintf("%d/%d (%.0f%%)", sum(e$sex == "male", na.rm = TRUE), sum(!is.na(e$sex)),
                      100 * mean(e$sex[!is.na(e$sex)] == "male")),
        nonert = sprintf("%d/%d (%.0f%%)", sum(n$sex == "male", na.rm = TRUE), sum(!is.na(n$sex)),
                      100 * mean(n$sex[!is.na(n$sex)] == "male")),
        test = "Fisher exact", stat = NA,
        p = tryCatch(fisher.test(matrix(c(sum(e$sex == "male", na.rm = TRUE), sum(e$sex == "female", na.rm = TRUE),
                                          sum(n$sex == "male", na.rm = TRUE), sum(n$sex == "female", na.rm = TRUE)), 2))$p.value,
                     error = function(x) NA))
  }
  T1 <- do.call(rbind, rows)
  T1$p_fmt <- ifelse(is.na(T1$p), "", ifelse(T1$p < 0.001, "<0.001", sprintf("%.3f", T1$p)))
  write.csv(T1, lald_out("M1_table1.csv"), row.names = FALSE)

  s <- T1[T1$phenotype == "CESD" & T1$variable %in%
          c("n","Age at diagnosis, y","Follow-up, y","LDL-C, mg/dL","ALT, U/L","Cirrhosis",
            "Any lipid-lowering therapy","Died","Male"), ]
  s <- T1[T1$phenotype == "WD" & T1$variable %in%
          c("n","Age at diagnosis, y","Follow-up, y","Adrenal calcification","Growth failure","Died","Male"), ]
})

local({
  set.seed(42)
  db <- readRDS(lald_cache("M0_db_derived.rds"))
  COVS <- c("age_dx_y", "severity_index", "LLT_harmonized")

  BIN <- list(c("O1", "Mortality", "mortality", "CESD"), c("O1", "Mortality", "mortality", "WD"),
              c("O2x", "AST normalization", "AST_norm", "CESD"),
              c("O5", "Steatosis regression", "O5_steatosis", "CESD"))
  CON <- list(c("O2b", "ALT % change", "ALT_pct"), c("O2c", "AST % change", "AST_pct"),
              c("O3", "LDL-C % change", "LDL_pct"), c("O3b", "TG % change", "TG_pct"))

  est_bin <- function(d, v) {
    d$.y <- d[[v]]
    d <- d[stats::complete.cases(d[, c(".y", "ert_received")]), ]
    if (nrow(d) < 8 || length(unique(d$.y)) < 2 || length(unique(d$ert_received)) < 2) return(NULL)
    f <- tryCatch(lald_firth(.y ~ ert_received, d), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    j <- which(f$names == "ert_received")
    list(n = nrow(d), est = exp(f$coef[j]), lcl = exp(f$ci[j, 1]), ucl = exp(f$ci[j, 2]), p = f$p[j])
  }
  est_con <- function(d, v) {
    d <- d[!is.na(d[[v]]) & !is.na(d$ert_received), ]
    if (sum(d$ert_received == 1) < 3 || sum(d$ert_received == 0) < 3) return(NULL)
    ts <- theil_sen(d$ert_received, d[[v]])
    list(n = nrow(d), est = ts$slope, lcl = NA, ucl = NA, p = ts$p)
  }
  row_of <- function(scen, id, lab, ph, r, note = "") {
    if (is.null(r)) return(data.frame(scenario = scen, outcome_id = id, outcome = lab, phenotype = ph,
      n = NA, estimate = NA, lcl = NA, ucl = NA, p = NA, note = paste("not estimable", note)))
    data.frame(scenario = scen, outcome_id = id, outcome = lab, phenotype = ph,
               n = r$n, estimate = r$est, lcl = r$lcl, ucl = r$ucl, p = r$p, note = note)
  }
  run_block <- function(scen, data, note = "") {
    out <- list()
    for (b in BIN) out[[length(out) + 1]] <-
      row_of(scen, b[1], b[2], b[4], est_bin(data[data$phenotype %in% b[4], ], b[3]), note)
    for (cc in CON) out[[length(out) + 1]] <-
      row_of(scen, cc[1], cc[2], "CESD", est_con(data[data$phenotype %in% "CESD", ], cc[3]), note)
    do.call(rbind, out)
  }

  R <- list()
  R[[1]] <- run_block("primary (all rows, unadjusted Firth / Theil-Sen)", db)

  d2 <- db
  zfields <- c("ert_received" = "Tx_ERT", "statin" = "ESTATINAS", "liver_tx" = "TRANSPLANTE_HEPATICO",
               "hsct" = "TRANSPLANTE_HEMATOPOYETICO", "dyslipidemia_bl" = "DISLIPIDEMIA_BASAL",
               "hepatosplenomegaly_bl" = "HEPATOESPLENOMEGALIA_BASAL", "steatosis_bl" = "ESTEATOSIS_BASAL",
               "cirrhosis_bl" = "CIRROSIS_BASAL", "ascvd_bl" = "ECVA_BASAL",
               "growth_failure_bl" = "RETRASO_CRECIMIENTO", "gi_symptoms_bl" = "GASTROINTESTINAL",
               "adrenal_calcification" = "CALCIFICACION_ADRENALES")
  for (an in names(zfields)) {
    sf <- zfields[[an]]
    hit <- grepl(paste0("(^|;)", sf, "(;|$)"), ifelse(is.na(d2$rule_coded_zero_fields), "", d2$rule_coded_zero_fields))
    d2[[an]][hit] <- NA
  }
  d2$severity_index <- with(d2, {
    m <- cbind(dyslipidemia_bl, hepatosplenomegaly_bl, steatosis_bl, cirrhosis_bl)
    ifelse(rowSums(!is.na(m)) == 0, NA, rowSums(m, na.rm = TRUE)) })
  R[[2]] <- run_block("S2 rule-coded zeros set missing", d2)

  R[[3]] <- run_block("S3 full-text rows only", db[db$artefact_level == "full_text", ])
  R[[4]] <- run_block("S5 risk of bias Low or Moderate", db[db$rob_overall %in% c("Low", "Moderate"), ])
  R[[5]] <- run_block("S6 published 1990 or later", db[db$era_post1990 %in% 1, ])
  R[[6]] <- run_block("S7 untreated-only comparator",
                      db[db$ert_received %in% 1 | (db$ert_received %in% 0 & db$LLT_harmonized %in% 0), ])
  heavy <- names(sort(table(db$source_R_ID), decreasing = TRUE))[1:5]
  for (s in heavy)
    R[[length(R) + 1]] <- run_block(paste0("S4 leave-one-source-out: ", s),
                                    db[db$source_R_ID != s, ], paste("dropped", s))
  RES <- do.call(rbind, R)

  cl <- list()
  for (cc in CON) {
    d <- db[db$phenotype %in% "CESD" & !is.na(db[[cc[3]]]) & !is.na(db$ert_received), ]
    if (nrow(d) < 12) next
    rownames(d) <- as.character(seq_len(nrow(d)))
    f <- lm(as.formula(paste(cc[3], "~ ert_received")), data = d)
    clv <- setNames(d$source_R_ID, rownames(d))
    V <- tryCatch(vcov_cluster(f, clv), error = function(e) NULL)
    b <- coef(f)["ert_received"]
    se_i <- summary(f)$coefficients["ert_received", 2]
    se_c <- if (!is.null(V)) sqrt(V["ert_received", "ert_received"]) else NA
    cl[[length(cl) + 1]] <- data.frame(
      scenario = "S1 cluster-robust SE on source paper", outcome_id = cc[1], outcome = cc[2],
      phenotype = "CESD", n = nrow(d), n_clusters = length(unique(d$source_R_ID)),
      estimate = b, se_iid = se_i, se_cluster = se_c,
      ratio = se_c / se_i,
      lcl = b - 1.96 * se_c, ucl = b + 1.96 * se_c,
      p = 2 * pnorm(-abs(b / se_c)))
  }
  CL <- do.call(rbind, cl)
  write.csv(RES, lald_out("M5_sensitivity_results.csv"), row.names = FALSE)
  write.csv(CL,  lald_out("M5_cluster_robust.csv"), row.names = FALSE)
})

local({
  set.seed(42)
  db <- readRDS(lald_cache("M8_db_modality.rds"))
  C  <- db[db$phenotype %in% "CESD", ]
  num <- function(x) suppressWarnings(as.numeric(x))

  OUT <- list(c("ALT_dx_n","ALT_post_n","ALT","U/L"),
              c("AST_dx_n","AST_post_n","AST","U/L"),
              c("TC_dx_n","TC_post_n","Total cholesterol","mg/dL"),
              c("LDL_dx_n","LDL_post_n","LDL-C","mg/dL"),
              c("HDL_dx_n","HDL_post_n","HDL-C","mg/dL"),
              c("TG_dx_n","TG_post_n","Triglycerides","mg/dL"))

  rows <- list()
  for (o in OUT) {
    a <- o[1]; b <- o[2]; lab <- o[3]; unit <- o[4]
    d <- C[!is.na(num(C[[a]])) & !is.na(num(C[[b]])) & !is.na(num(C$ert_received)), ]
    if (!nrow(d)) next
    d$pct <- 100 * (num(d[[b]]) - num(d[[a]])) / num(d[[a]])
    e <- d[num(d$ert_received) == 1, ]; u <- d[num(d$ert_received) == 0, ]
    if (nrow(e) < 3 || nrow(u) < 3) next
    ts <- theil_sen(num(d$ert_received), d$pct)
    bc <- tryCatch(bca_ci(data.frame(y = d$pct, x = num(d$ert_received)), function(dat, i) {
            s <- dat[i, ]; y1 <- s$y[s$x == 1]; y0 <- s$y[s$x == 0]
            if (length(y1) < 2 || length(y0) < 2) return(NA_real_)
            median(outer(y1, y0, "-")) }, R = 2000), error = function(e) NULL)
    llt <- num(u$LLT_harmonized); llt <- llt[!is.na(llt)]
    fu  <- num(d$fu_years); fu <- fu[!is.na(fu)]
    rows[[length(rows) + 1]] <- data.frame(
      outcome = lab, unit = unit,
      n_ert = nrow(e), median_pct_ert = round(median(e$pct), 1),
      median_baseline_ert = round(median(num(e[[a]])), 0),
      median_followup_ert = round(median(num(e[[b]])), 0),
      n_untreated = nrow(u), median_pct_untreated = round(median(u$pct), 1),
      median_baseline_untreated = round(median(num(u[[a]])), 0),
      median_followup_untreated = round(median(num(u[[b]])), 0),
      controlled_shift_pp = round(ts$slope, 1),
      shift_lcl = if (is.null(bc)) NA else round(bc$lcl, 1),
      shift_ucl = if (is.null(bc)) NA else round(bc$ucl, 1),
      p_unadjusted = signif(ts$p, 3),
      pct_of_within_arm_change = if (median(e$pct) == 0) NA else
          round(100 * ts$slope / median(e$pct), 0),
      llt_in_untreated = if (length(llt)) sprintf("%d/%d (%.0f%%)", sum(llt), length(llt),
                                                  100 * mean(llt)) else "-",
      median_followup_years = if (length(fu)) round(median(fu), 1) else NA,
      registered = "no - exploratory")
  }
  R <- do.call(rbind, rows)

  S <- C[C$modality %in% c("none", "ERT"), ]
  S$ert_received <- as.integer(S$modality == "ERT")
  srows <- list()
  for (o in OUT) {
    a <- o[1]; b <- o[2]; lab <- o[3]
    d <- S[!is.na(num(S[[a]])) & !is.na(num(S[[b]])), ]
    if (!nrow(d)) next
    d$pct <- 100 * (num(d[[b]]) - num(d[[a]])) / num(d[[a]])
    e <- d[d$ert_received == 1, ]; u <- d[d$ert_received == 0, ]
    if (nrow(e) < 3 || nrow(u) < 3) next
    ts <- theil_sen(d$ert_received, d$pct)
    bc <- tryCatch(bca_ci(data.frame(y = d$pct, x = d$ert_received), function(dat, i) {
            z <- dat[i, ]; y1 <- z$y[z$x == 1]; y0 <- z$y[z$x == 0]
            if (length(y1) < 2 || length(y0) < 2) return(NA_real_)
            median(outer(y1, y0, "-")) }, R = 2000), error = function(e) NULL)
    srows[[length(srows) + 1]] <- data.frame(
      outcome = lab, sens_n_ert = nrow(e), sens_n_untreated = nrow(u),
      sens_pct_ert = round(median(e$pct), 1), sens_pct_untreated = round(median(u$pct), 1),
      sens_shift_pp = round(ts$slope, 1),
      sens_lcl = if (is.null(bc)) NA else round(bc$lcl, 1),
      sens_ucl = if (is.null(bc)) NA else round(bc$ucl, 1),
      sens_p = signif(ts$p, 3),
      sens_share = if (median(e$pct) == 0) NA else round(100 * ts$slope / median(e$pct), 0),
      stringsAsFactors = FALSE)
  }
  R <- merge(R, do.call(rbind, srows), by = "outcome", all.x = TRUE, sort = FALSE)
  R <- R[match(vapply(OUT, `[`, "", 3), R$outcome), ]
  write.csv(R, lald_out("M18_controlled_change.csv"), row.names = FALSE)

  R$q_bh <- signif(p.adjust(R$p_unadjusted, method = "BH"), 3)
  write.csv(R, lald_out("M18_controlled_change.csv"), row.names = FALSE)
})

local({
  set.seed(42)

  db <- read.csv(lald_in("ipd"), na.strings = c("", "NA"))
  U  <- db[db$ert_received %in% 1, ]
  N  <- db[db$ert_received %in% 0, ]

  fmt_p <- function(p) if (is.na(p)) "-" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
  bm_p  <- function(a, b) {
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(NA_real_)
    brunner_munzel(a, b)$p
  }
  med_q <- function(x, type, digits, scale = 1) {
    x <- x[!is.na(x)]
    if (!length(x)) return("-")
    q <- quantile(x, c(0.5, 0.25, 0.75), type = type, names = FALSE)
    if (scale != 1) q <- round(q) / scale
    sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f) [n=%d]"), q[1], q[2], q[3], length(x))
  }
  bin <- function(x, level = 1) {
    x <- x[!is.na(x)]
    if (!length(x)) return("-")
    k <- sum(x == level)
    sprintf("%d/%d (%.1f)", k, length(x), 100 * k / length(x))
  }
  fisher_p <- function(v, level = 1) {
    a <- U[[v]][!is.na(U[[v]])]; b <- N[[v]][!is.na(N[[v]])]
    if (!length(a) || !length(b)) return(NA_real_)
    fisher.test(matrix(c(sum(a == level), sum(a != level), sum(b == level), sum(b != level)), 2))$p.value
  }

  row_cont <- function(label, v) data.frame(variable = label,
    total = med_q(db[[v]], 6, 1), users = med_q(U[[v]], 6, 1), nonusers = med_q(N[[v]], 6, 1),
    p = fmt_p(bm_p(U[[v]], N[[v]])))
  row_bin <- function(label, v, level = 1) data.frame(variable = label,
    total = bin(db[[v]], level), users = bin(U[[v]], level), nonusers = bin(N[[v]], level),
    p = fmt_p(fisher_p(v, level)))

  T1 <- rbind(
    row_cont("Age at diagnosis, years", "age_dx_y"),
    data.frame(variable = "Years of ERT use", total = "-", users = med_q(U$ert_duration_y, 6, 1),
               nonusers = "-", p = "-"),
    row_cont("Follow up, years", "fu_years"),
    row_bin("Sex male (n, %)", "sex", "male"),
    row_bin("ASCVD (n, %)", "ascvd_bl"),
    row_bin("Dyslipidemia (n, %)", "dyslipidemia_bl"),
    row_bin("Hepatosplenomegaly (n, %)", "hepatosplenomegaly_bl"),
    row_bin("Liver steatosis (n, %)", "steatosis_bl"),
    row_bin("Liver cirrhosis (n, %)", "cirrhosis_bl"),
    row_bin("Liver transplant (n, %)", "liver_tx"),
    row_bin("Post liver transplant relapse (n, %)", "recurrence_post_liver_tx"),
    row_bin("LAL Activity measurement (n, %)", "LAL_measured"),
    row_bin("Hematopoietic stem cell transplantation (n, %)", "hsct"),
    row_bin("Stunted growth (n, %)", "growth_failure_bl"),
    row_bin("Adrenal calcification (n, %)", "adrenal_calcification"),
    row_bin("Gastrointestinal symptoms (n, %)", "gi_symptoms_bl"),
    row_bin("Dead (n, %)", "mortality"))
  ert_age <- function(ph) med_q(U$ert_start_age_y[U$phenotype == ph], 6, 1)
  T1 <- rbind(T1[1, ], data.frame(variable = "Age of start of ERT, years", total = "-",
    users = sprintf("CESD %s; WD %s", ert_age("CESD"), ert_age("WD")), nonusers = "-", p = "-"), T1[-1, ])
  write.csv(T1, lald_out("M19_table1.csv"), row.names = FALSE)

  LIP <- list(c("Total cholesterol, mmol/L", "TC", 38.67), c("HDL-C, mmol/L", "HDL", 38.67),
              c("LDL-C, mmol/L", "LDL", 38.67), c("Triglycerides, mmol/L", "TG", 88.57),
              c("Apolipoprotein B, g/L", "ApoB", 100))
  T2 <- rbind(
    data.frame(variable = "Dyslipidemia (n, %)", total = bin(db$dyslipidemia_bl), users = bin(U$dyslipidemia_bl),
               post = "-", nonusers = bin(N$dyslipidemia_bl), p = fmt_p(fisher_p("dyslipidemia_bl"))),
    data.frame(variable = "Lipid-lowering treatment (n, %)", total = bin(db$LLT_harmonized), users = bin(U$LLT_harmonized),
               post = "-", nonusers = bin(N$LLT_harmonized), p = fmt_p(fisher_p("LLT_harmonized"))),
    data.frame(variable = "Statins (n, %)", total = bin(db$statin), users = bin(U$statin),
               post = "-", nonusers = bin(N$statin), p = fmt_p(fisher_p("statin"))))
  for (l in LIP) {
    dx <- paste0(l[2], "_dx_n"); po <- paste0(l[2], "_post_n"); s <- as.numeric(l[3])
    T2 <- rbind(T2, data.frame(variable = l[1],
      total = med_q(db[[dx]], 6, 2, s), users = med_q(U[[dx]], 6, 2, s), post = med_q(U[[po]], 6, 2, s),
      nonusers = med_q(N[[dx]], 6, 2, s), p = fmt_p(bm_p(U[[dx]], N[[dx]]))))
  }
  write.csv(T2, lald_out("M19_table2.csv"), row.names = FALSE)

  SX <- do.call(rbind, lapply(c("CESD", "WD"), function(ph) {
    g <- db[db$phenotype %in% ph, ]
    data.frame(phenotype = ph,
               variable = c("LAL activity measured (n, %)", "Enzyme replacement therapy (n, %)"),
               value = c(bin(g$LAL_measured), bin(g$ert_received)))
  }))
  write.csv(SX, lald_out("M19_tableS5_S7_cells.csv"), row.names = FALSE)
})
