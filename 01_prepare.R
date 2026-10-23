options(stringsAsFactors = FALSE, scipen = 999)

if (!exists("LALD_ROOT")) LALD_ROOT <- getwd()
LALD_FILES <- list(
  ipd       = "data/lald_ipd.csv",
  aggregate = "data/lald_aggregate.csv",
  reports   = "data/lald_reports.csv")
lald_in <- function(key) {
  p <- file.path(LALD_ROOT, LALD_FILES[[key]])
  if (!file.exists(p)) stop("LALD: cannot find ", p)
  p
}
lald_out <- function(...) {
  d <- file.path(LALD_ROOT, "results"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}
lald_cache <- function(...) {
  d <- file.path(LALD_ROOT, "results", "cache"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}
lald_fig <- function(...) {
  d <- file.path(LALD_ROOT, "figures"); dir.create(d, showWarnings = FALSE, recursive = TRUE)
  file.path(d, ...)
}

lald_firth <- function(formula, data) {
  fit <- logistf::logistf(formula, data = data)
  list(coef = coef(fit), ci = cbind(fit$ci.lower, fit$ci.upper),
       p = fit$prob, n = fit$n, names = names(coef(fit)), engine = "logistf")
}

theil_sen <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 3) return(list(slope = NA, intercept = NA, p = NA, n = n))
  ct <- suppressWarnings(cor.test(x, y, method = "kendall"))
  m <- mblm::mblm(y ~ x, dataframe = data.frame(x = x, y = y), repeated = FALSE)
  list(slope = unname(coef(m)[2]), intercept = unname(coef(m)[1]), p = unname(ct$p.value), n = n)
}

bca_ci <- function(data, statistic, R = 2000, level = 0.95, seed = 42) {
  set.seed(seed)
  bo <- boot::boot(data, statistic, R = R)
  ci <- tryCatch(boot::boot.ci(bo, type = "bca", conf = level),
                 error = function(e) NULL)
  if (is.null(ci) || is.null(ci$bca)) {
    q <- quantile(bo$t[, 1], c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE)
    return(list(est = bo$t0[1], lcl = q[[1]], ucl = q[[2]], type = "percentile_fallback", R = R))
  }
  list(est = bo$t0[1], lcl = ci$bca[4], ucl = ci$bca[5], type = "BCa", R = R)
}

vcov_cluster <- function(fit, cluster) {
  mf <- model.frame(fit)
  cl <- cluster[match(rownames(mf), names(cluster))]
  if (anyNA(cl)) stop("cluster vector does not align with the model frame")
  sandwich::vcovCL(fit, cluster = factor(cl), type = "HC1")
}

rma_re <- function(yi, vi, method = "REML") {
  ok <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[ok]; vi <- vi[ok]
  if (!length(yi)) return(NULL)
  r <- metafor::rma(yi = yi, vi = vi, method = method)
  list(b = as.numeric(r$b), se = r$se, ci.lb = r$ci.lb, ci.ub = r$ci.ub,
       pval = r$pval, tau2 = r$tau2, I2 = r$I2, Q = r$QE, Qp = r$QEp,
       k = r$k, method = method)
}
escalc_PLO <- function(xi, ni, add = 0.5) {
  x <- xi; n <- ni
  bad <- (x == 0 | x == n)
  x[bad] <- x[bad] + add; n[bad] <- n[bad] + 2 * add
  list(yi = log(x / (n - x)), vi = 1 / x + 1 / (n - x))
}
escalc_MD <- function(m1, sd1, n1, m2, sd2, n2) {
  list(yi = m1 - m2, vi = sd1^2 / n1 + sd2^2 / n2)
}
escalc_OR <- function(a, b, c, d, add = 0.5) {
  z <- (a == 0 | b == 0 | c == 0 | d == 0)
  a[z] <- a[z] + add; b[z] <- b[z] + add; c[z] <- c[z] + add; d[z] <- d[z] + add
  list(yi = log((a * d) / (b * c)), vi = 1/a + 1/b + 1/c + 1/d)
}

brunner_munzel <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(list(phat = NA, stat = NA, df = NA, p = NA, n1 = n1, n2 = n2))
  r  <- rank(c(x, y)); r1 <- r[1:n1]; r2 <- r[(n1 + 1):(n1 + n2)]
  m1 <- mean(r1); m2 <- mean(r2)
  rx <- rank(x); ry <- rank(y)
  s1 <- sum((r1 - rx - m1 + (n1 + 1) / 2)^2) / (n1 - 1)
  s2 <- sum((r2 - ry - m2 + (n2 + 1) / 2)^2) / (n2 - 1)
  N  <- n1 + n2
  sigma <- N * (s1 / n2 + s2 / n1)
  if (!is.finite(sigma) || sigma <= 0) return(list(phat = NA, stat = NA, df = NA, p = NA, n1 = n1, n2 = n2))
  stat <- (n1 * n2 * (m2 - m1)) / (N * sqrt(n1 * n2 * sigma / N))
  df <- (s1 / n2 + s2 / n1)^2 /
        ((s1 / n2)^2 / (n1 - 1) + (s2 / n1)^2 / (n2 - 1))
  phat <- (m2 - (n2 + 1) / 2) / n1
  list(phat = phat, stat = stat, df = df,
       p = 2 * pt(-abs(stat), max(df, 1)), n1 = n1, n2 = n2)
}

d_from_logOR <- function(logOR) logOR * sqrt(3) / pi
hedges_g <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  sp <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  d  <- (mean(x) - mean(y)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2 - 2)))
  c(g = J * d, lcl = J * d - 1.96 * J * se, ucl = J * d + 1.96 * J * se)
}
bh <- function(p) p.adjust(p, method = "BH")
fmt_ci <- function(e, l, u, d = 2)
  sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), e, l, u)

local({
  set.seed(42)
  db <- read.csv(lald_in("ipd"), check.names = FALSE, na.strings = c("", "NA"))
  pct <- function(dx, post) { d <- (post - dx) / dx * 100; d[!is.finite(d)] <- NA; d }
  db$ALT_pct <- pct(db$ALT_dx_n, db$ALT_post_n)
  db$AST_pct <- pct(db$AST_dx_n, db$AST_post_n)
  db$LDL_pct <- pct(db$LDL_dx_n, db$LDL_post_n)
  db$TG_pct  <- pct(db$TG_dx_n,  db$TG_post_n)
  db$TC_pct  <- pct(db$TC_dx_n,  db$TC_post_n)
  db$ALT_norm <- with(db, ifelse(!is.na(ALT_dx_n) & !is.na(ALT_post_n) & ALT_dx_n > 40,
                                 as.integer(ALT_post_n <= 40), NA))
  db$AST_norm <- with(db, ifelse(!is.na(AST_dx_n) & !is.na(AST_post_n) & AST_dx_n > 40,
                                 as.integer(AST_post_n <= 40), NA))
  code2 <- function(x) ifelse(is.na(x), NA, as.integer(x == 2))
  db$O5_steatosis <- ifelse(db$steatosis_bl %in% 1, code2(db$steatosis_outcome_coded), NA)
  db$O5_steatosis_unrestricted <- code2(db$steatosis_outcome_coded)
  db$O4_lipid     <- code2(db$lipid_outcome_coded)
  saveRDS(db, lald_cache("M0_db_derived.rds"))
})
