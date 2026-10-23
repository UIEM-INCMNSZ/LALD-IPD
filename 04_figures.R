if (!exists("lald_out")) source("01_prepare.R")

local({
  set.seed(42)

  rep_ <- read.csv(lald_in("reports"), colClasses = "character", na.strings = character(0))
  ipd  <- read.csv(lald_in("ipd"), na.strings = c("", "NA"))
  agg  <- read.csv(lald_in("aggregate"), colClasses = "character", na.strings = character(0))

  SEARCH <- c(PubMed = 611, Scopus = 758, SciELO = 10)
  DUPLICATE_RECORDS <- 335

  n_ident   <- sum(SEARCH)
  n_screen  <- n_ident - DUPLICATE_RECORDS
  n_ipd     <- sum(rep_$arm == "ipd" & rep_$R_ID != "INCMNSZ")
  n_agg     <- sum(rep_$arm == "aggregate")
  n_excl    <- sum(rep_$arm == "excluded")
  n_dup     <- sum(rep_$arm == "duplicate")
  n_assess  <- n_ipd + n_agg + n_excl
  n_sought  <- n_assess + n_dup
  n_ex_scr  <- n_screen - n_sought
  n_inst    <- sum(ipd$source_R_ID == "INCMNSZ")
  n_reach   <- sum(rep_$arm == "aggregate" & rep_$reaches_estimate == "1")
  n_groups  <- length(unique(agg$overlap_group[agg$overlap_group != ""]))
  k <- function(x) format(x, big.mark = ",", trim = TRUE)

  PT <- 96 / 72
  INK <- "#1a1d21"; RULE <- "#9ca3af"; SOFT <- "#f4f5f7"; ACC <- "#1f5563"

  rbox <- function(x, y, w, h, r = 1.2, pad = 0.9, ...) {
    x0 <- x - pad; y0 <- y - pad; x1 <- x + w + pad; y1 <- y + h + pad
    a <- seq(0, pi / 2, length.out = 12)
    px <- c(x1 - r + r * cos(a), x0 + r - r * sin(a), x0 + r - r * cos(a), x1 - r + r * sin(a))
    py <- c(y1 - r + r * sin(a), y1 - r + r * cos(a), y0 + r - r * sin(a), y0 + r - r * cos(a))
    polygon(px, py, ...)
  }
  box <- function(x, y, w, h, lines, fs = 10, bold0 = TRUE, face = "white", edge = INK, lw = 1.1) {
    rbox(x, y, w, h, col = face, border = edge, lwd = lw * PT)
    n <- length(lines); step <- h / (n + 0.6)
    for (i in seq_along(lines))
      text(x + w / 2, y + h - step * (i - 1 + 0.8), lines[i], cex = fs / 12, col = INK,
           font = if (i == 1 && bold0) 2 else 1)
  }
  arrow <- function(x1, y1, x2, y2, lw = 1.2, color = INK) {
    dx <- x2 - x1; dy <- y2 - y1; L <- sqrt(dx^2 + dy^2); ux <- dx / L; uy <- dy / L
    hl <- 1.3; hw <- 0.6
    bx <- x2 - ux * hl; by <- y2 - uy * hl
    segments(x1, y1, bx, by, lwd = lw * PT, col = color)
    polygon(c(x2, bx - uy * hw, bx + uy * hw), c(y2, by + ux * hw, by - ux * hw), col = color, border = color)
  }
  stage <- function(y, label) text(1.5, y, label, srt = 90, cex = 9.2 / 12, col = ACC, font = 2)

  draw <- function() {
    par(mar = c(0, 0, 0, 0), family = FAM)
    plot.new(); plot.window(xlim = c(0, 112), ylim = c(-6, 118), xaxs = "i", yaxs = "i")
    MX <- 11; MW <- 53; SX <- 70; SW <- 32
    stage(104, "IDENTIFICATION")
    box(MX, 100, MW, 13, c("Records identified by search", paste("n =", k(n_ident)),
        paste(sprintf("%s %s", names(SEARCH), k(SEARCH)), collapse = "  ·  ")), face = SOFT)
    box(SX, 100, SW, 13, c("Institutional series", "INCMNSZ chart review", sprintf("1 source, %d patients", n_inst)),
        fs = 9.4, edge = RULE, lw = 0.9)
    arrow(MX + MW / 2, 100, MX + MW / 2, 94)
    box(MX, 81, MW, 13, c("Records after duplicates removed", paste("n =", k(n_screen))))
    box(SX, 81, SW, 13, c("Duplicate records", paste("n =", k(DUPLICATE_RECORDS))), fs = 9.4, edge = RULE, lw = 0.9)
    arrow(MX + MW, 87.5, SX, 87.5, lw = 1.0, color = RULE)
    stage(72, "SCREENING")
    arrow(MX + MW / 2, 81, MX + MW / 2, 75)
    box(MX, 62, MW, 13, c("Records screened on title/abstract", paste("n =", k(n_screen))))
    box(SX, 62, SW, 13, c("Excluded, not LAL-D or", "no patient described", paste("n =", k(n_ex_scr))),
        fs = 9.4, edge = RULE, lw = 0.9)
    arrow(MX + MW, 68.5, SX, 68.5, lw = 1.0, color = RULE)
    arrow(MX + MW / 2, 62, MX + MW / 2, 56)
    box(MX, 43, MW, 13, c("Reports sought for retrieval", paste("n =", k(n_sought))))
    box(SX, 43, SW, 13, c("Further duplicates found", "on full record", paste("n =", k(n_dup))),
        fs = 9.4, edge = RULE, lw = 0.9)
    arrow(MX + MW, 49.5, SX, 49.5, lw = 1.0, color = RULE)
    stage(34, "ELIGIBILITY")
    arrow(MX + MW / 2, 43, MX + MW / 2, 37)
    box(MX, 24, MW, 13, c("Reports assessed for eligibility", paste("n =", k(n_assess))))
    box(SX, 24, SW, 13, c("Excluded", "no patient-level and", "no usable group-level data", paste("n =", k(n_excl))),
        fs = 9.4, edge = RULE, lw = 0.9)
    arrow(MX + MW, 30.5, SX, 30.5, lw = 1.0, color = RULE)
    stage(13, "INCLUDED")
    arrow(MX + MW / 2, 24, MX + MW / 2, 18)
    box(MX, 2, MW, 16, c("INCLUDED  ·  individual participant data",
        sprintf("%d studies (%d retrieved + 1 institutional)", n_ipd + 1, n_ipd),
        sprintf("%d participants", nrow(ipd)),
        sprintf("%d CESD  ·  %d Wolman  ·  %d not reported", sum(ipd$phenotype %in% "CESD"),
                sum(ipd$phenotype %in% "WD"), sum(is.na(ipd$phenotype)))), fs = 9.0, face = SOFT)
    box(SX, 2, SW, 16, c("INCLUDED  ·  aggregate data", sprintf("%d studies", n_agg),
        sprintf("%d reach a reported estimate", n_reach), sprintf("%d overlap groups", n_groups)),
        fs = 9.0, face = SOFT, edge = ACC, lw = 1.4)
    arrow(MX + MW, 10, SX, 10, lw = 1.1, color = ACC)
    segments(SX + SW, 106.5, SX + SW + 3.5, 106.5, col = ACC, lwd = 1.1 * PT)
    segments(SX + SW + 3.5, 106.5, SX + SW + 3.5, -0.5, col = ACC, lwd = 1.1 * PT)
    segments(SX + SW + 3.5, -0.5, MX + MW / 2, -0.5, col = ACC, lwd = 1.1 * PT)
    arrow(MX + MW / 2, -0.5, MX + MW / 2, 1.6, lw = 1.1, color = ACC)
  }

  FAM <- "Arial"
  W <- 9.4; H <- 10
  cairo_pdf(lald_fig("Figure1.pdf"), width = W, height = H, family = FAM); draw(); dev.off()
  cairo_ps(lald_fig("Figure1.eps"), width = W, height = H, family = FAM); draw(); dev.off()
  png(lald_fig("Figure1.png"), width = W, height = H, units = "in", res = 600, type = "cairo", family = FAM); draw(); dev.off()
})

local({
  set.seed(42)
  db  <- readRDS(lald_cache("M0_db_derived.rds"))
  agg <- read.csv(lald_in("aggregate"), check.names = FALSE, na.strings = c("", "NA"))
  agg <- agg[agg$row_kind == "result", ]

  CESD <- db[db$phenotype %in% "CESD", ]; WD <- db[db$phenotype %in% "WD", ]
  pctf <- function(a, b) { d <- (b - a) / a * 100; d[!is.finite(d)] <- NA; d }
  CESD$HDL_pct <- pctf(CESD$HDL_dx_n, CESD$HDL_post_n)

  INK <- "#221f26"; INK2 <- "#544e58"; MUTE <- "#87808c"; RULE <- "#ddd8d9"
  ERTC <- "#6b4e71"; CTRLC <- "#8a8f98"; POOL <- "#9c5326"; GOOD <- "#5b7355"; BAD <- "#a34b3a"

  pr <- function(d, v, ert) { x <- d[[v]][!is.na(d[[v]]) & d$ert_received %in% ert]; c(sum(x == 1), length(x)) }
  me <- function(d, v, ert) { x <- d[[v]][!is.na(d[[v]]) & d$ert_received %in% ert]
                              c(median(x), mean(x), sd(x), length(x)) }
  derived_pct <- function(base, eos) round((eos - base) / base * 100, 1)

  altn <- rbind(
    data.frame(study = "IPD (this review)",            grp = "IPD",        x = pr(CESD, "ALT_norm", 1)[1], n = pr(CESD, "ALT_norm", 1)[2], src = "reported"),
    data.frame(study = "ARISE final (Burton 2022)",    grp = "ARISE",      x = 31, n = 66, src = "reported"),
    data.frame(study = "LAL-CL04 5-yr (Malinova 2020)",grp = "LAL-CL01/04", x = 5, n = 5,  src = "reported"),
    data.frame(study = "Registry 1-yr (D'Antiga 2025)",grp = "Registry",   x = 22, n = 53, src = "reported"))
  astn <- rbind(
    data.frame(study = "IPD (this review)",            grp = "IPD",        x = pr(CESD, "AST_norm", 1)[1], n = pr(CESD, "AST_norm", 1)[2], src = "reported"),
    data.frame(study = "ARISE final (Burton 2022)",    grp = "ARISE",      x = 43, n = 65, src = "reported"),
    data.frame(study = "LAL-CL04 5-yr (Malinova 2020)",grp = "LAL-CL01/04", x = 4, n = 4,  src = "reported"),
    data.frame(study = "Registry 1-yr (D'Antiga 2025)",grp = "Registry",   x = 25, n = 51, src = "reported"))
  ipd_ratio <- function(bv, pv) {
    d <- CESD[CESD$ert_received %in% 1 & !is.na(CESD[[bv]]) & !is.na(CESD[[pv]]), ]
    b <- d[[bv]]; a <- d[[pv]]
    c(val = 100 * median(a - b) / median(b), n = nrow(d))
  }
  lip <- function(bv, pv, arise, cl04_d, cl04_b, cl06_d, cl06_b) {
    i <- ipd_ratio(bv, pv)
    rbind(
      data.frame(study = "IPD (this review)",             grp = "IPD",         val = round(i[["val"]], 1), n = i[["n"]], src = "ratio"),
      data.frame(study = "LAL-CL04 5-yr (Malinova 2020)", grp = "LAL-CL01/04", val = round(100 * cl04_d / cl04_b, 1), n = 7,  src = "ratio"),
      data.frame(study = "LAL-CL06 144wk (Burton 2022)",  grp = "LAL-CL06",    val = round(100 * cl06_d / cl06_b, 1), n = 19, src = "ratio"),
      data.frame(study = "ARISE final (Burton 2022)",     grp = "ARISE",       val = arise, n = 66, src = "reported_percent"))
  }
  ldl <- lip("LDL_dx_n", "LDL_post_n", -24, -20.1, 113.1, -52.6, 160.0)
  hdl <- lip("HDL_dx_n", "HDL_post_n",  26,   8.1,  36.9,   9.8,  31.1)
  tg  <- lip("TG_dx_n",  "TG_post_n",   NA,  -1.8, 104.1, -35.4, 159.4)

  pool_prop <- function(d) { es <- escalc_PLO(d$x, d$n); r <- rma_re(es$yi, es$vi, "REML")
    list(est = plogis(r$b) * 100, lo = plogis(r$ci.lb) * 100, hi = plogis(r$ci.ub) * 100, k = r$k, I2 = r$I2, m = r$method) }
  pool_val <- function(d) { d <- d[!is.na(d$val) & d$src == "ratio", ]
    yi <- d$val; vi <- rep(NA, nrow(d))
    w <- d$n / sum(d$n); est <- sum(w * yi)
    se <- sqrt(sum(w^2) * var(yi) * nrow(d) / max(nrow(d) - 1, 1))
    list(est = est, lo = est - 1.96 * se, hi = est + 1.96 * se, k = nrow(d), I2 = NA, m = "n-weighted") }

  cohen_h <- function(x1, n1, x2, n2) {
    p1 <- x1 / n1; p2 <- x2 / n2
    h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2)); se <- sqrt(1 / n1 + 1 / n2)
    c(h, h - 1.96 * se, h + 1.96 * se)
  }
  or_of <- function(x1, n1, x2, n2) {
    a <- x1; b <- n1 - x1; c2 <- x2; d2 <- n2 - x2
    z <- (a == 0 | b == 0 | c2 == 0 | d2 == 0)
    if (z) { a <- a + .5; b <- b + .5; c2 <- c2 + .5; d2 <- d2 + .5 }
    lo <- log(a * d2 / (b * c2)); se <- sqrt(1/a + 1/b + 1/c2 + 1/d2)
    c(exp(lo), exp(lo - 1.96 * se), exp(lo + 1.96 * se))
  }
  B <- list()
  a1 <- pr(CESD, "ALT_norm", 1); a0 <- pr(CESD, "ALT_norm", 0)
  s1 <- pr(CESD, "AST_norm", 1); s0 <- pr(CESD, "AST_norm", 0)
  wsu <- function(d, ert) { x <- d$mortality[!is.na(d$mortality) & d$ert_received %in% ert]; c(sum(x == 0), length(x)) }
  w1 <- wsu(WD, 1); w0 <- wsu(WD, 0); c1 <- wsu(CESD, 1); c0 <- wsu(CESD, 0)
  mkB <- function(lab, est, lo, hi, kind, extra) data.frame(lab = lab, est = est, lo = lo, hi = hi, kind = kind, extra = extra)
  o <- or_of(a1[1], a1[2], a0[1], a0[2]); h <- cohen_h(a1[1], a1[2], a0[1], a0[2])
  B[[1]] <- mkB(sprintf("ALT normalization - IPD  (h=%.2f)", h[1]), o[1], o[2], o[3], "OR", sprintf("%d/%d vs %d/%d", a1[1], a1[2], a0[1], a0[2]))
  o <- or_of(s1[1], s1[2], s0[1], s0[2]); h <- cohen_h(s1[1], s1[2], s0[1], s0[2])
  B[[2]] <- mkB(sprintf("AST normalization - IPD  (h=%.2f)", h[1]), o[1], o[2], o[3], "OR", sprintf("%d/%d vs %d/%d", s1[1], s1[2], s0[1], s0[2]))
  o <- or_of(w1[1], w1[2], w0[1], w0[2]); h <- cohen_h(w1[1], w1[2], w0[1], w0[2])
  B[[3]] <- mkB(sprintf("WD survival - IPD  (h=%.2f)", h[1]), o[1], o[2], o[3], "OR", sprintf("%d/%d vs %d/%d", w1[1], w1[2], w0[1], w0[2]))
  o <- or_of(c1[1], c1[2], c0[1], c0[2]); h <- cohen_h(c1[1], c1[2], c0[1], c0[2])
  B[[4]] <- mkB(sprintf("CESD survival - IPD  (h=%.2f)", h[1]), o[1], o[2], o[3], "OR", sprintf("%d/%d vs %d/%d", c1[1], c1[2], c0[1], c0[2]))
  cont_B <- function(v, lab) {
    e <- CESD[[v]][!is.na(CESD[[v]]) & CESD$ert_received %in% 1]
    n <- CESD[[v]][!is.na(CESD[[v]]) & CESD$ert_received %in% 0]
    g <- hedges_g(e, n); ts <- theil_sen(CESD$ert_received[!is.na(CESD[[v]])], CESD[[v]][!is.na(CESD[[v]])])
    mkB(sprintf("%s - IPD  (g=%+.2f)", lab, g[1]), ts$slope, NA, NA, "pp", sprintf("p=%.3f", ts$p))
  }
  B[[5]] <- cont_B("LDL_pct", "LDL-C % change")
  B[[6]] <- cont_B("HDL_pct", "HDL-C % change")
  B[[7]] <- cont_B("TG_pct",  "TG % change")
  BB <- do.call(rbind, B)

  DX0 <- 0.06; DX1 <- 0.78
  FAM <- "sans"
  open_dev <- function(dev, f, w, h) {
    cairo <- capabilities("cairo")
    FAM <<- if (dev %in% c("cpdf", "eps", "png")) "Arial" else "sans"
    switch(dev,
      pdf  = pdf(f, width = w, height = h),
      svg  = svg(f, width = w, height = h),
      cpdf = if (cairo) cairo_pdf(f, width = w, height = h, family = FAM) else pdf(f, width = w, height = h),
      eps  = if (cairo) cairo_ps(f, width = w, height = h, family = FAM, fallback_resolution = 300)
             else postscript(f, width = w, height = h, horizontal = FALSE, onefile = FALSE,
                             paper = "special", family = "Helvetica"),
      png  = png(f, width = w, height = h, units = "in", res = 300,
                 type = if (cairo) "cairo" else "quartz", family = FAM),
      stop("unknown device ", dev))
  }
  draw_fig1 <- function(dev, path = NULL) {
    f <- if (is.null(path)) lald_fig(paste0("LALD_figure1_v4.", dev)) else path
    open_dev(dev, f, 15.2, 9.6)
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE), heights = c(2.25, 1))
    par(family = FAM, oma = c(0, 0, 0.6, 0))
    A_rows <- list()

    par(mar = c(4.6, 17.5, 2.4, 11.5))
    blocks <- list(list("ALT normalization", altn, "prop"), list("AST normalization", astn, "prop"),
                   list("LDL-C, change from baseline", ldl, "val"),
                   list("HDL-C, change from baseline", hdl, "val"),
                   list("Triglycerides, change from baseline", tg, "val"))
    rowsN <- sum(sapply(blocks, function(b) nrow(b[[2]]) + 2))
    plot(NA, xlim = c(0, 1), ylim = c(0.2, rowsN + 0.4), axes = FALSE, xlab = "", ylab = "")
    title(main = "A)  Pooled ERT response",
          adj = 0, cex.main = 1.06, line = 0.7, col.main = INK)
    LO <- -70; HI <- 105
    mx <- function(v) DX0 + (v - LO) / (HI - LO) * (DX1 - DX0)
    segments(mx(0), 0.2, mx(0), rowsN + 0.2, col = RULE, lty = 3)
    y <- rowsN
    for (b in blocks) {
      nm <- b[[1]]; d <- b[[2]]; kind <- b[[3]]
      text(-0.925, y, nm, adj = 0, font = 2, cex = 0.85, xpd = NA, col = INK); y <- y - 1
      for (i in seq_len(nrow(d))) {
        lab <- d$study[i]
        if (kind == "prop") { v <- 100 * d$x[i] / d$n[i]; txt <- sprintf("%d/%d   %.0f%%", d$x[i], d$n[i], v) }
        else { v <- d$val[i]
               txt <- if (is.na(v)) "no percent change reported" else
                 sprintf("%+.1f%%  (n=%d)%s", v, d$n[i], ifelse(d$src[i] == "reported_percent", "  *", "")) }
        text(-0.025, y, lab, adj = 1, cex = 0.78, col = INK2, xpd = NA)
        ref <- !is.na(d$src[i]) && d$src[i] == "reported_percent"
        if (!is.na(v)) { segments(mx(0), y, mx(v), y, col = ifelse(ref, CTRLC, ERTC),
                                  lwd = ifelse(ref, 2.2, 3.4), lend = 1)
                         points(mx(v), y, pch = ifelse(ref, 5, 18), cex = 1.15,
                                col = ifelse(ref, CTRLC, ERTC)) }
        text(1.02, y, txt, adj = 0, cex = 0.74, col = INK2, xpd = NA)
        A_rows[[length(A_rows) + 1]] <- data.frame(panel = "A", block = nm, study = lab,
          value = ifelse(is.na(v), NA, round(v, 1)), n = d$n[i], source = d$src[i])
        y <- y - 1
      }
      p <- if (kind == "prop") pool_prop(d) else pool_val(d)
      segments(mx(p$lo), y, mx(p$hi), y, col = POOL, lwd = 1.8)
      points(mx(p$est), y, pch = 18, cex = 1.8, col = POOL)
      text(-0.025, y, sprintf("Pooled  (%s, k=%d)", p$m, p$k), adj = 1, cex = 0.78, font = 2, xpd = NA, col = INK)
      text(1.02, y, sprintf("%+.0f%% [%.0f, %.0f]%s", p$est, p$lo, p$hi,
           ifelse(is.na(p$I2), "", sprintf("   I2 %.0f%%", p$I2))), adj = 0, cex = 0.74, font = 2, xpd = NA, col = INK)
      A_rows[[length(A_rows) + 1]] <- data.frame(panel = "A", block = nm, study = "POOLED",
        value = round(p$est, 1), n = p$k, source = p$m)
      y <- y - 1
    }
    tk <- seq(-60, 100, 20)
    axis(1, at = mx(tk), labels = paste0(tk, "%"), cex.axis = 0.72, col = MUTE, col.axis = INK2, line = -0.3)
    mtext("normalization rate, or change from baseline", side = 1, line = 2.0, cex = 0.73, col = INK2, at = mx(20))
    mtext("lipid rows: percent = median paired change / median baseline, computed identically for every pooled cohort",
          side = 1, line = 3.0, cex = 0.64, col = MUTE, at = mx(20))
    mtext("*  ARISE reports a per-patient percent change, a different estimator; shown open, not pooled",
          side = 1, line = 3.9, cex = 0.64, col = MUTE, at = mx(20))

    par(mar = c(4.6, 17.5, 2.4, 11.5))
    ors <- BB[BB$kind == "OR", ]; pps <- BB[BB$kind == "pp", ]
    k <- nrow(ors) + nrow(pps) + 1
    plot(NA, xlim = c(0, 1), ylim = c(0.2, k + 0.4), axes = FALSE, xlab = "", ylab = "")
    title(main = "B)  ERT versus control, within the IPD", adj = 0, cex.main = 1.06, line = 0.7, col.main = INK)
    lo <- 0.02; hi <- 400
    mo <- function(v) DX0 + (log(pmax(pmin(v, hi), lo)) - log(lo)) / (log(hi) - log(lo)) * (DX1 - DX0)
    segments(mo(1), 0.2, mo(1), k + 0.2, col = RULE, lty = 3)
    yy <- k
    for (i in seq_len(nrow(ors))) {
      r <- ors[i, ]
      text(-0.025, yy, r$lab, adj = 1, cex = 0.8, col = INK2, xpd = NA)
      segments(mo(r$lo), yy, mo(r$hi), yy, lwd = 1.6, col = INK2)
      points(mo(r$est), yy, pch = 18, cex = 1.6, col = ifelse(r$lo > 1, GOOD, ERTC))
      text(1.02, yy, sprintf("OR %.2f [%.2f, %.2f]", r$est, r$lo, r$hi), adj = 0, cex = 0.74, col = INK2, xpd = NA)
      text(1.02, yy - 0.34, r$extra, adj = 0, cex = 0.66, col = MUTE, xpd = NA)
      yy <- yy - 1
    }
    yy <- yy - 0.5
    text(-0.925, yy + 0.45, "continuous outcomes, median difference (percentage points)",
         adj = 0, cex = 0.72, font = 3, col = MUTE, xpd = NA)
    for (i in seq_len(nrow(pps))) {
      r <- pps[i, ]
      text(-0.025, yy, r$lab, adj = 1, cex = 0.8, col = INK2, xpd = NA)
      v <- mo(1) + r$est / 100 * (DX1 - DX0) * 0.9
      segments(mo(1), yy, v, yy, col = ERTC, lwd = 3.4, lend = 1); points(v, yy, pch = 18, cex = 1.15, col = ERTC)
      text(1.02, yy, sprintf("%+.1f pp   %s", r$est, r$extra), adj = 0, cex = 0.74, col = INK2, xpd = NA)
      yy <- yy - 1
    }
    tk <- c(0.02, 0.1, 0.5, 1, 5, 20, 100, 400)
    axis(1, at = mo(tk), labels = tk, cex.axis = 0.72, col = MUTE, col.axis = INK2, line = -0.3)
    mtext("odds ratio, log scale        favours control  <--|-->  favours ERT",
          side = 1, line = 2.0, cex = 0.73, col = INK2, at = mo(1))
    mtext("lower rows: Theil-Sen median differences on a linear offset, not odds",
          side = 1, line = 3.1, cex = 0.66, col = MUTE, at = mo(1))

    surv_panel <- function(ttl, rows) {
      par(mar = c(4.6, 17.5, 2.4, 11.5))
      k <- length(rows)
      plot(NA, xlim = c(0, 1), ylim = c(0.3, k + 0.7), axes = FALSE, xlab = "", ylab = "")
      title(main = ttl, adj = 0, cex.main = 1.06, line = 0.7, col.main = INK)
      ms <- function(v) DX0 + v / 100 * (DX1 - DX0)
      for (g in c(0, 25, 50, 75, 100)) segments(ms(g), 0.3, ms(g), k + 0.5, col = RULE, lty = 3)
      yy <- k
      for (r in rows) {
        v <- 100 * as.numeric(r[[2]]) / as.numeric(r[[3]])
        text(-0.025, yy, r[[1]], adj = 1, cex = 0.82, col = INK2, xpd = NA)
        segments(ms(0), yy, ms(max(v, 0.8)), yy, col = r[[4]], lwd = 9, lend = 1)
        text(1.02, yy, sprintf("%s/%s   %.0f%%", r[[2]], r[[3]], v), adj = 0, cex = 0.76, col = INK2, xpd = NA)
        yy <- yy - 1
      }
      axis(1, at = ms(seq(0, 100, 25)), labels = paste0(seq(0, 100, 25), "%"),
           cex.axis = 0.72, col = MUTE, col.axis = INK2, line = -0.3)
      mtext("alive at last reported follow-up", side = 1, line = 2.0, cex = 0.73, col = INK2, at = ms(50))
    }
    surv_panel("C)  CESD survival", list(
      list("IPD, ERT", c1[1], c1[2], ERTC), list("IPD, non-ERT", c0[1], c0[2], CTRLC)))
    surv_panel("D)  Wolman disease survival", list(
      list("IPD, ERT", w1[1], w1[2], ERTC), list("IPD, non-ERT", w0[1], w0[2], CTRLC),
      list("VITAL + CL08 (Vijay 2021), ERT", 15, 19, ERTC),
      list("LAL-1-NH01 untreated (Jones 2016)", 0, 21, CTRLC)))

    dev.off()
    do.call(rbind, A_rows)
  }
  A_data <- draw_fig1("cpdf", lald_fig("Figure2.pdf"))
  invisible(draw_fig1("eps", lald_fig("Figure2.eps")))
  invisible(draw_fig1("png", lald_fig("Figure2.png")))
  write.csv(rbind(cbind(A_data),
                  data.frame(panel = "B", block = BB$kind, study = BB$lab,
                             value = round(BB$est, 3), n = NA, source = BB$extra)),
            lald_out("figure2_data.csv"), row.names = FALSE)

  sc <- function(x, by) x / by
  mk_model <- function(kind, yv, data, terms, labels) {
    d <- data
    for (t in names(terms)) d[[t]] <- terms[[t]](d)
    vars <- c(yv, names(terms))
    d <- d[stats::complete.cases(d[, vars]), vars]
    if (nrow(d) < 10) return(NULL)
    keep <- names(terms)[sapply(names(terms), function(v) length(unique(d[[v]])) > 1)]
    f <- as.formula(paste(yv, "~", paste(keep, collapse = " + ")))
    if (kind == "binary") {
      fit <- tryCatch(lald_firth(f, d), error = function(e) NULL); if (is.null(fit)) return(NULL)
      j <- match(keep, fit$names)
      out <- data.frame(term = keep, label = unname(labels[keep]),
                        est = exp(fit$coef[j]), lo = exp(fit$ci[j, 1]), hi = exp(fit$ci[j, 2]))
      ev <- sum(d[[yv]] == 1)
      list(tab = out, kind = "OR", n = nrow(d), ev = ev,
           fitstat = sprintf("n=%d, events=%d", nrow(d), ev))
    } else {
      fit <- lm(f, data = d); cf <- summary(fit)$coefficients; ci <- confint(fit)
      out <- data.frame(term = keep, label = unname(labels[keep]),
                        est = cf[keep, 1], lo = ci[keep, 1], hi = ci[keep, 2])
      list(tab = out, kind = "pp", n = nrow(d), ev = NA,
           fitstat = sprintf("n=%d, R2=%.2f, adj R2=%.2f", nrow(d),
                             summary(fit)$r.squared, summary(fit)$adj.r.squared))
    }
  }
  LAB <- c(ert_received = "ERT (yes vs no)", age10 = "Age at diagnosis (per 10 y)",
           fu1 = "Follow-up (per year)", alt50 = "Baseline ALT (per 50 U/L)",
           ast50 = "Baseline AST (per 50 U/L)", ldl50 = "Baseline LDL-C (per 1.29 mmol/L, 50 mg/dL)",
           hdl10 = "Baseline HDL-C (per 0.26 mmol/L, 10 mg/dL)", tg100 = "Baseline TG (per 1.13 mmol/L, 100 mg/dL)",
           LLT_harmonized = "Lipid-lowering therapy", onsetm = "Age at onset (per month)", agedxm = "Age at diagnosis (per month)")
  TERMS <- list(ert_received = function(d) d$ert_received,
                age10 = function(d) sc(d$age_dx_y, 10), fu1 = function(d) d$fu_years,
                alt50 = function(d) sc(d$ALT_dx_n, 50), ast50 = function(d) sc(d$AST_dx_n, 50),
                ldl50 = function(d) sc(d$LDL_dx_n, 50), hdl10 = function(d) sc(d$HDL_dx_n, 10),
                tg100 = function(d) sc(d$TG_dx_n, 100), LLT_harmonized = function(d) d$LLT_harmonized,
                onsetm = function(d) d$onset_age_y * 12,
                agedxm = function(d) d$age_dx_y * 12)
  pick_terms <- function(...) TERMS[c(...)]
  WD2 <- WD; WD2$alive <- ifelse(is.na(WD2$mortality), NA, 1L - WD2$mortality)
  MODELS <- list(
    list("ALT normalization",  "binary", "ALT_norm", CESD, pick_terms("ert_received","age10","fu1","alt50")),
    list("AST normalization",  "binary", "AST_norm", CESD, pick_terms("ert_received","age10","fu1","ast50")),
    list("LDL-C (% change)",   "cont",   "LDL_pct",  CESD, pick_terms("ert_received","age10","fu1","ldl50","LLT_harmonized")),
    list("HDL-C (% change)",   "cont",   "HDL_pct",  CESD, pick_terms("ert_received","age10","fu1","hdl10","LLT_harmonized")),
    list("TG (% change)",      "cont",   "TG_pct",   CESD, pick_terms("ert_received","age10","fu1","tg100","LLT_harmonized")),
    list("Wolman survival",    "binary", "alive",    WD2,  pick_terms("ert_received","agedxm","onsetm")))
  FITS <- lapply(MODELS, function(m) { r <- mk_model(m[[2]], m[[3]], m[[4]], m[[5]], LAB)
                                       if (!is.null(r)) r$title <- m[[1]]; r })
  draw_fig2 <- function(dev, path = NULL) {
    f <- if (is.null(path)) lald_fig(paste0("LALD_figure2_v4.", dev)) else path
    open_dev(dev, f, 13.6, 8.8)
    par(mfrow = c(3, 2), family = FAM, oma = c(0, 0, 0.6, 0))
    for (r in FITS) {
      if (is.null(r)) { par(mar = c(4, 12, 3, 6)); plot.new(); next }
      par(mar = c(4.2, 13.5, 3.0, 8.5))
      tb <- r$tab; k <- nrow(tb)
      if (r$kind == "OR") {
        lo <- max(min(tb$lo, na.rm = TRUE) * 0.6, 1e-3); hi <- min(max(tb$hi, na.rm = TRUE) * 1.6, 1e5)
        mp <- function(v) (log(pmax(pmin(v, hi), lo)) - log(lo)) / (log(hi) - log(lo))
        ref <- 1
      } else {
        lo <- min(tb$lo) - 5; hi <- max(tb$hi) + 5
        mp <- function(v) (v - lo) / (hi - lo); ref <- 0
      }
      plot(NA, xlim = c(0, 1), ylim = c(0.4, k + 0.6), axes = FALSE, xlab = "", ylab = "")
      title(main = r$title, adj = 0, cex.main = 1.02, line = 1.5, col.main = INK)
      mtext(r$fitstat, side = 3, line = 0.35, adj = 0, cex = 0.68, col = MUTE)
      segments(mp(ref), 0.4, mp(ref), k + 0.5, col = RULE, lty = 3)
      for (i in seq_len(k)) {
        y <- k - i + 1; isert <- tb$term[i] == "ert_received"
        text(-0.03, y, tb$label[i], adj = 1, cex = 0.78, xpd = NA,
             col = if (isert) INK else INK2, font = if (isert) 2 else 1)
        segments(mp(tb$lo[i]), y, mp(tb$hi[i]), y, lwd = if (isert) 1.9 else 1.3,
                 col = if (isert) ERTC else INK2)
        points(mp(tb$est[i]), y, pch = 18, cex = if (isert) 1.7 else 1.2,
               col = if (isert) ERTC else INK2)
        txt <- if (r$kind == "OR") sprintf("%.2f [%.2f, %.2f]", tb$est[i], tb$lo[i], tb$hi[i])
               else sprintf("%+.1f [%+.1f, %+.1f]", tb$est[i], tb$lo[i], tb$hi[i])
        text(1.03, y, txt, adj = 0, cex = 0.72, xpd = NA, col = if (isert) INK else INK2,
             font = if (isert) 2 else 1)
      }
      if (r$kind == "OR") { tk <- c(0.01, 0.1, 0.5, 1, 2, 10, 100, 1000); tk <- tk[tk > lo & tk < hi]
        axis(1, at = mp(tk), labels = tk, cex.axis = 0.7, col = MUTE, col.axis = INK2, line = -0.4)
        mtext("odds ratio, log scale", side = 1, line = 1.9, cex = 0.7, col = INK2)
      } else { tk <- pretty(c(lo, hi), 5); tk <- tk[tk > lo & tk < hi]
        axis(1, at = mp(tk), labels = tk, cex.axis = 0.7, col = MUTE, col.axis = INK2, line = -0.4)
        mtext("percentage points per unit", side = 1, line = 1.9, cex = 0.7, col = INK2) }
    }
    dev.off()
  }
  draw_fig2("cpdf", lald_fig("Figure3.pdf"))
  draw_fig2("eps", lald_fig("Figure3.eps"))
  draw_fig2("png", lald_fig("Figure3.png"))
  F2 <- do.call(rbind, lapply(FITS, function(r) if (is.null(r)) NULL else
    data.frame(model = r$title, kind = r$kind, n = r$n, r$tab, fitstat = r$fitstat)))
  write.csv(F2, lald_out("figure3_data.csv"), row.names = FALSE)
})
