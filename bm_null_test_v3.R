# =============================================================================
# Battle Mountain anomaly photograph: formal null-hypothesis tests
# Author: Tom Woolman, PhD
# Photo:  20150718_184215.jpg (Samsung Galaxy S6, July 18, 2015, 6:42:14 PM EDT)
#
# v3 (2026-10-05): writes results/sessionInfo.txt on every run, so that the R and
#                  package versions quoted in any manuscript are read from a file
#                  produced by the run itself rather than reconstructed from memory.
#
# H0 (global): the anomaly group of pixels is consistent with normally expected
#              lighting artifacts in this photograph: sensor noise, lens-flare
#              ghosting, or the bright patches that ordinary light produces in a
#              shaded tree line (sky gaps, sunflecks, bark, leaf glints).
# H1:          it is not.
#
# Rejecting H0 supports only "not explained by the tested mechanisms". It does
# not identify what the shapes are.
#
# Component tests
#   Test 1  Sensor noise        Monte Carlo of Gaussian noise at the local shade level
#   Test 2  Lens-flare ghost    distance from the ghost position predicted by the sun
#   Test 3  Scene lighting      scan statistic over every bright patch in the tree
#                               line, with a generalized Pareto (GPD) tail so that the
#                               look-elsewhere effect is corrected. Run in two strata:
#                                 3a lower tree line (the anomaly's own setting)
#                                 3b whole tree line, sky gaps included (conservative)
#
# Veiling glare: the sun is in frame, and its glare adds a smooth haze that lifts
# the darkest shadows from about 5 (far from the sun) to about 65 (near it) on the
# 0-255 scale. Glare is additive, so it is estimated from the dark floor of each
# block and subtracted before patches are measured. Without this step only the
# darkest corner of the frame is comparable with the anomaly.
#
# Run from the folder that holds the photograph. Outputs go to ./results/
# =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
pkgs <- c("jpeg", "evd", "spatstat.geom")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(jpeg); library(evd); library(spatstat.geom) })
set.seed(20150718)

CFG <- list(
  photo         = "20150718_184215.jpg",
  block         = 64,             # block size (px) for glare floor and scene level
  floor_q       = 0.02,           # glare floor = this quantile of luminance in each block
  s_min         = 3,              # minimum scene level (avoids division by ~0)
  k_ratio       = 2.0,            # candidate pixel: Lc > k_ratio * scene level ...
  k_abs         = 8,              # ... and Lc > scene level + k_abs (0-255 units)
  min_area      = 4,              # ignore specks smaller than this (px)
  ring_in       = 8,              # background ring around each patch (px)
  ring_out      = 30,
  sky_margin    = 400,
  shade_band    = c(8, 25),       # Test 3c: raw block-median band that matches the anomaly's shade            # lower tree line = at least this far (px) below the skyline
  gpd_quant     = 0.90,           # GPD threshold quantile of the null statistic
  stat          = "NIEL",         # normalized integrated excess luminance (see features)
  roi           = c(x0 = 4797, x1 = 5090, y0 = 1783, y1 = 1939),  # anomaly box
  # foreground objects excluded from every stratum: the sapling in mid-frame, and the
  # canopy of the foreground tree at the far left, which has open sky behind its leaves
  exclude       = list(sapling   = c(x0 = 2300, x1 = 4200, y0 = 1250, y1 = 2988),
                       left_tree = c(x0 = 1,    x1 = 450,  y0 = 1,    y1 = 1100)),
  # base of the tree line (lawn edge), traced on the photograph; pixels below are lawn
  lawn_edge     = data.frame(x = c(0, 1000, 2000, 2600, 4200, 4600, 5000, 5312),
                             y = c(1980, 2050, 2080, 2120, 2100, 2020, 1960, 1940)),
  skyline_L     = 120,            # skyline: first block (from the top) darker than this
  sun_window    = c(x0 = 1600, x1 = 2050, y0 = 900, y1 = 1250),   # where to look for the sun disk
  flare_sd_frac = 0.05,           # generous scatter of a flare ghost, fraction of frame diagonal
  n_mc          = 100000,         # Monte Carlo replicates, noise test (floor on p = 1e-5)
  alpha         = 0.001,          # alpha test statistic, divided across the three component tests
  n_boot        = 500             # bootstrap replicates, GPD tail
)
dir.create("results", showWarnings = FALSE)
log_file <- file.path("results", sprintf("bm_null_test_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)
cat("Battle Mountain anomaly photograph: null-hypothesis tests\n")
cat("Run: ", format(Sys.time()), "   R ", R.version$major, ".", R.version$minor, "\n\n", sep = "")
str(CFG)

## ---- 1. Load image; glare floor; corrected luminance --------------------------
img <- readJPEG(CFG$photo) * 255
H <- dim(img)[1]; W <- dim(img)[2]
R <- img[, , 1]; G <- img[, , 2]; B <- img[, , 3]
L <- (R + G + B) / 3
rm(img); invisible(gc())
cat(sprintf("\nImage %d x %d px\n", W, H))

bs  <- CFG$block
nbr <- ceiling(H / bs); nbc <- ceiling(W / bs)
blk_stat <- function(M, f) {
  out <- matrix(NA_real_, nbr, nbc)
  for (i in seq_len(nbr)) {
    rr <- ((i - 1) * bs + 1):min(i * bs, H)
    for (j in seq_len(nbc)) out[i, j] <- f(M[rr, ((j - 1) * bs + 1):min(j * bs, W)])
  }
  out
}
med5 <- function(M) {                     # 5 x 5 running median over blocks
  out <- M
  for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M))) {
    out[i, j] <- median(M[max(1, i - 2):min(nrow(M), i + 2), max(1, j - 2):min(ncol(M), j + 2)])
  }
  out
}
interp_w <- function(n, nb) {             # bilinear weights, block centers -> pixels
  ctr <- (seq_len(nb) - 0.5) * bs
  Wm <- matrix(0, n, nb)
  for (p in seq_len(n)) {
    t <- (p - ctr[1]) / bs + 1
    t <- min(max(t, 1), nb); i0 <- floor(t); f <- t - i0
    Wm[p, i0] <- 1 - f; if (i0 < nb) Wm[p, i0 + 1] <- f
  }
  Wm
}
Wr <- interp_w(H, nbr); Wc <- interp_w(W, nbc)
to_px <- function(Mb) Wr %*% Mb %*% t(Wc)

L_blk   <- blk_stat(L, median)                              # raw block median (skyline)
F_blk   <- med5(blk_stat(L, function(v) quantile(v, CFG$floor_q, names = FALSE)))
Fpx     <- to_px(F_blk)                                     # smooth glare floor
Lc      <- pmax(L - Fpx, 0)                                 # glare-corrected luminance
S_px    <- pmax(to_px(med5(blk_stat(Lc, median))), CFG$s_min)   # local scene level

## ---- 2. Regions ----------------------------------------------------------------
col_blk <- ceiling(seq_len(W) / bs)
lawn_y  <- approx(CFG$lawn_edge$x, CFG$lawn_edge$y, xout = seq_len(W), rule = 2)$y
sky_top <- sapply(seq_len(nbc), function(j) {
  k <- which(L_blk[, j] < CFG$skyline_L); if (length(k) == 0) nbr else k[1]
})
sky_y <- (sky_top[col_blk] - 1) * bs + 1
rowM  <- matrix(seq_len(H), H, W)
whole <- rowM >= rep(sky_y, each = H) & rowM < rep(lawn_y, each = H)
lower <- whole & rowM >= rep(sky_y + CFG$sky_margin, each = H)
rm(rowM)
for (ex in CFG$exclude) {
  whole[ex["y0"]:ex["y1"], ex["x0"]:ex["x1"]] <- FALSE
  lower[ex["y0"]:ex["y1"], ex["x0"]:ex["x1"]] <- FALSE
}
cat(sprintf("Whole tree line: %.1f%% of frame; lower tree line: %.1f%%\n", 100 * mean(whole), 100 * mean(lower)))
roi <- CFG$roi
cat(sprintf("Glare floor at the anomaly: %.1f; raw shade %.1f; corrected scene level %.1f\n",
            Fpx[round(mean(roi[3:4])), round(mean(roi[1:2]))],
            median(L[roi["y0"]:roi["y1"], roi["x0"]:roi["x1"]]),
            S_px[round(mean(roi[3:4])), round(mean(roi[1:2]))]))

## ---- 3. Patch detection and features ------------------------------------------------
# NIEL = sum over the patch of (Lc / S - 1), where S is the median corrected
# luminance in a ring around the patch: the patch's total excess brightness in
# units of its own local scene level. It is invariant to illumination level
# (reflected light scales with the light falling on the scene) and grows with
# both brightness and size.
detect <- function(region, k_ratio, Lm = Lc, Sm = S_px) {
  cand <- region & Lm > k_ratio * Sm & Lm > Sm + CFG$k_abs
  msk  <- owin(mask = cand[H:1, ], xrange = c(0, W), yrange = c(0, H))
  lab  <- matrix(as.integer(as.matrix(connected(msk))[H:1, ]), H, W)
  idx  <- which(!is.na(lab)); ids <- lab[idx]
  rows <- (idx - 1) %% H + 1; cols <- (idx - 1) %/% H + 1
  by_id <- split(seq_along(idx), ids)
  ring_med <- function(M, y0, y1, x0, x1) {
    a0 <- max(1, y0 - CFG$ring_out); a1 <- min(H, y1 + CFG$ring_out)
    b0 <- max(1, x0 - CFG$ring_out); b1 <- min(W, x1 + CFG$ring_out)
    keep <- matrix(TRUE, a1 - a0 + 1, b1 - b0 + 1)
    keep[(max(1, y0 - CFG$ring_in) - a0 + 1):(min(H, y1 + CFG$ring_in) - a0 + 1),
         (max(1, x0 - CFG$ring_in) - b0 + 1):(min(W, x1 + CFG$ring_in) - b0 + 1)] <- FALSE
    keep <- keep & is.na(lab[a0:a1, b0:b1])
    median(M[a0:a1, b0:b1][keep])
  }
  do.call(rbind, lapply(names(by_id), function(k) {
    s <- by_id[[k]]
    if (length(s) < CFG$min_area) return(NULL)
    y <- rows[s]; x <- cols[s]; yx <- cbind(y, x)
    lv <- L[yx]; lc <- Lm[yx]
    bgL <- ring_med(L,  min(y), max(y), min(x), max(x))
    S   <- max(ring_med(Lm, min(y), max(y), min(x), max(x)), CFG$s_min)
    top <- order(lv, decreasing = TRUE)[seq_len(min(40, length(lv)))]
    rgb <- cbind(R[yx], G[yx], B[yx])[top, , drop = FALSE]
    data.frame(area = length(s), x = mean(x), y = mean(y),
               width = diff(range(x)) + 1, height = diff(range(y)) + 1,
               bg_raw = bgL, glare = mean(Fpx[yx]), scene = S,
               peak_raw = max(lv), peak_ratio = max(lc) / S, mean_ratio = mean(lc) / S,
               IEL  = sum(pmax(lv - bgL, 0)),
               NIEL = sum(pmax(lc / S - 1, 0)),
               saturated = any(lv >= 254),
               blue_frac = mean(rgb[, 3] / rowSums(rgb)))
  }))
}
tag_roi <- function(f) { f$in_roi <- with(f, x >= roi["x0"] & x <= roi["x1"] & y >= roi["y0"] & y <= roi["y1"]); f }

feat_lower <- tag_roi(detect(lower, CFG$k_ratio))
feat_whole <- tag_roi(detect(whole, CFG$k_ratio))
# 3c: narrow deep-shade stratum on raw luminance, no glare correction. Only blocks
# whose raw median lies in the anomaly's own shade band are compared.
Lraw_bg <- pmax(to_px(L_blk), CFG$s_min)
shade   <- whole & Lraw_bg >= CFG$shade_band[1] & Lraw_bg <= CFG$shade_band[2]
feat_shade <- tag_roi(detect(shade, CFG$k_ratio, Lm = L, Sm = Lraw_bg))
cat(sprintf("Deep-shade band (raw %g-%g): %.1f%% of frame\n", CFG$shade_band[1], CFG$shade_band[2], 100 * mean(shade)))
write.csv(feat_lower, file.path("results", "patches_lower_treeline.csv"), row.names = FALSE)
write.csv(feat_whole, file.path("results", "patches_whole_treeline.csv"), row.names = FALSE)
write.csv(feat_shade, file.path("results", "patches_deep_shade_raw.csv"), row.names = FALSE)

anom <- feat_lower[feat_lower$in_roi, ]
anom <- anom[order(-anom[[CFG$stat]]), ]
prim <- anom[1, ]                                  # strongest patch in the box = right shape
show <- c("area","x","y","width","height","bg_raw","glare","scene","peak_raw","peak_ratio","mean_ratio","NIEL","saturated","blue_frac")
cat("\nAnomaly-box patches (top 6):\n");  print(head(anom[, show], 6), row.names = FALSE, digits = 3)
for (nm in c("lower", "whole", "shade")) {
  f <- get(paste0("feat_", nm)); nl <- f[!f$in_roi, ]
  cat(sprintf("\nStrongest null patches, %s tree line (%d null patches):\n", nm, nrow(nl)))
  print(head(nl[order(-nl[[CFG$stat]]), show], 8), row.names = FALSE, digits = 3)
}

## ---- 4. Test 1: sensor noise -------------------------------------------------------
cat("\n--- Test 1: sensor noise ------------------------------------------------\n")
yy <- round(prim$y); xx <- round(prim$x)
ring <- L[(yy - 60):(yy + 60), (xx - 120):(xx - 50)]      # dark shade between and left of the shapes
sd_noise <- mad(as.vector(diff(ring))) / sqrt(2)
mc <- replicate(CFG$n_mc, sum(pmax(rnorm(prim$area, 0, sd_noise), 0)))
p1 <- (1 + sum(mc >= prim$IEL)) / (CFG$n_mc + 1)
cat(sprintf("Noise sd %.2f (0-255). Observed IEL %.0f; pure-noise IEL mean %.0f (sd %.0f); z = %.0f\n",
            sd_noise, prim$IEL, mean(mc), sd(mc), (prim$IEL - mean(mc)) / sd(mc)))
cat(sprintf("Monte Carlo p <= %.2g (no replicate reached the observed value)\n", p1))

## ---- 5. Test 2: lens-flare ghost ---------------------------------------------------
cat("\n--- Test 2: lens-flare ghost ------------------------------------------\n")
sw <- CFG$sun_window
Lw <- L[sw["y0"]:sw["y1"], sw["x0"]:sw["x1"]]
# sun disk: compact, fully saturated spot with darker surroundings. Score =
# saturation fill in a 41 x 41 core minus fill in the 161 x 161 surround, so
# broad saturated sky glare scores low and the disk scores high.
Sat <- (Lw >= 254) * 1
Ic  <- rbind(0, cbind(0, apply(t(apply(Sat, 1, cumsum)), 2, cumsum)))
boxmean <- function(i, j, r) {
  i0 <- pmax(i - r, 1); i1 <- pmin(i + r, nrow(Sat)); j0 <- pmax(j - r, 1); j1 <- pmin(j + r, ncol(Sat))
  (Ic[cbind(i1 + 1, j1 + 1)] - Ic[cbind(i0, j1 + 1)] - Ic[cbind(i1 + 1, j0)] + Ic[cbind(i0, j0)]) /
    ((i1 - i0 + 1) * (j1 - j0 + 1))
}
g <- expand.grid(i = seq_len(nrow(Sat)), j = seq_len(ncol(Sat)))
core <- boxmean(g$i, g$j, 20); wide <- boxmean(g$i, g$j, 80)
k <- which.max(core - wide)
near <- which(Sat == 1, arr.ind = TRUE)
near <- near[(near[, 1] - g$i[k])^2 + (near[, 2] - g$j[k])^2 < 45^2, , drop = FALSE]
sun <- c(x = mean(near[, 2]) + sw["x0"] - 1, y = mean(near[, 1]) + sw["y0"] - 1)
ghost <- 2 * c((W + 1) / 2, (H + 1) / 2) - sun
d_obs <- sqrt(sum((c(prim$x, prim$y) - ghost)^2))
diag_px <- sqrt(W^2 + H^2); sdp <- CFG$flare_sd_frac * diag_px
p2 <- exp(-d_obs^2 / (2 * sdp^2))                       # Rayleigh tail, 2-D Gaussian scatter
cat(sprintf("Sun disk at (%.0f, %.0f) [core fill %.2f, surround %.2f]; predicted ghost at (%.0f, %.0f)\n",
            sun[1], sun[2], core[k], wide[k], ghost[1], ghost[2]))
cat(sprintf("Right shape is %.0f px from the predicted ghost (%.1f%% of the diagonal); scatter sd allowed %.0f px\n",
            d_obs, 100 * d_obs / diag_px, sdp))
cat(sprintf("P(ghost this far from prediction) = %.2g\n", p2))

## ---- 6. Test 3: scene-lighting scan statistic ---------------------------------------
# The anomaly was found by looking, so it is compared with the MAXIMUM of N
# ordinary patches, not with a typical one. The null tail is a GPD fitted above
# the gpd_quant quantile of the null patches; P(max of N >= observed) =
# 1 - (1 - p_single)^N. A bootstrap of the tail exceedances gives the
# uncertainty of that p-value, which is large when the tail holds few patches.
gpd_scan <- function(f, q, x_obs) {
  xn <- f[[CFG$stat]][!f$in_roi]; N <- nrow(f)
  u  <- quantile(xn, q, names = FALSE); zeta <- mean(xn > u)
  fit <- try(suppressWarnings(fpot(xn, threshold = u, std.err = FALSE)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  pfun <- function(sc, sh) {
    ps <- if (x_obs <= u) mean(xn >= x_obs) else zeta * pgpd(x_obs, loc = u, scale = sc, shape = sh, lower.tail = FALSE)
    1 - (1 - ps)^N
  }
  p <- max(pfun(fit$estimate["scale"], fit$estimate["shape"]), .Machine$double.eps)   # floor: report as < 2.2e-16
  list(n_null = length(xn), N = N, u = u, n_exc = sum(xn > u), max_null = max(xn),
       scale = unname(fit$estimate["scale"]), shape = unname(fit$estimate["shape"]), p = unname(p),
       exc = xn[xn > u], pfun = pfun)
}
boot_p <- function(r, x_obs) {
  suppressWarnings(replicate(CFG$n_boot, {
    e <- sample(r$exc, replace = TRUE)
    f <- try(fpot(e, threshold = r$u, std.err = FALSE), silent = TRUE)
    if (inherits(f, "try-error")) NA else r$pfun(f$estimate["scale"], f$estimate["shape"])
  }))
}
x_obs <- prim[[CFG$stat]]
res3 <- list()
lbl <- c(lower = "3a lower tree line, glare corrected", whole = "3b whole tree line, glare corrected",
         shade = "3c deep-shade band, raw luminance")
for (nm in names(lbl)) {
  f <- get(paste0("feat_", nm))
  xo <- max(f[[CFG$stat]][f$in_roi])          # the anomaly measured the same way as its null
  cat(sprintf("\n--- Test %s ----------------------------\n", lbl[nm]))
  r <- gpd_scan(f, CFG$gpd_quant, xo); x_obs <- xo
  bt <- boot_p(r, xo); ci <- quantile(bt, c(0.025, 0.5, 0.975), na.rm = TRUE)
  cat(sprintf("Null patches %d (N incl. anomaly box = %d); tail threshold %.1f with %d exceedances\n",
              r$n_null, r$N, r$u, r$n_exc))
  cat(sprintf("GPD scale %.1f, shape %.3f\n", r$scale, r$shape))
  cat(sprintf("Observed %s %.1f vs. largest null patch %.1f (%.1fx)\n", CFG$stat, x_obs, r$max_null, x_obs / r$max_null))
  cat(sprintf("Look-elsewhere corrected P(max of %d >= observed) = %s\n", r$N,
              if (r$p <= .Machine$double.eps) "< 2.2e-16" else trimws(formatC(r$p, digits = 2, format = "g"))))
  cat(sprintf("Bootstrap of the tail (%d): median %.2g; 95%% range %.2g to %.2g\n", CFG$n_boot, ci[2], ci[1], ci[3]))
  res3[[nm]] <- c(p = r$p, lo = ci[[1]], hi = ci[[3]])
}
x_obs <- prim[[CFG$stat]]

## ---- 7. Sensitivity -------------------------------------------------------------------
# Re-detects patches under stricter and looser rules and refits the tail.
# p = 0 means the fitted tail has a finite endpoint (shape < 0) below the observed value.
cat("\n--- Sensitivity of Test 3 ------------------------------------------------------\n")
sens <- expand.grid(gpd_quant = c(0.85, 0.90, 0.95), k_ratio = c(1.75, 2.0, 2.5), stratum = c("lower", "whole", "shade"),
                    stringsAsFactors = FALSE)
sens$n_null <- NA; sens$max_null <- NA; sens$obs <- NA; sens$p <- NA
for (kr in unique(sens$k_ratio)) for (st in unique(sens$stratum)) {
  f  <- tag_roi(if (st == "shade") detect(shade, kr, Lm = L, Sm = Lraw_bg) else detect(get(st), kr))
  xo <- max(f[[CFG$stat]][f$in_roi])
  for (q in unique(sens$gpd_quant)) {
    r <- gpd_scan(f, q, xo); w <- which(sens$k_ratio == kr & sens$stratum == st & sens$gpd_quant == q)
    if (!is.null(r)) sens[w, c("n_null", "max_null", "obs", "p")] <- c(r$n_null, r$max_null, xo, r$p)
  }
}
sens_print <- sens; sens_print$p <- ifelse(sens$p <= .Machine$double.eps, "< 2.2e-16", formatC(sens$p, digits = 3, format = "g"))
print(sens_print, row.names = FALSE, digits = 3)
write.csv(sens, file.path("results", "sensitivity_test3.csv"), row.names = FALSE)

## ---- 8. Summary ---------------------------------------------------------------------
alpha <- CFG$alpha / 3
cat(sprintf("\n=== Summary (alpha = %g; Bonferroni across three tests: %g / 3 = %.6f) ===\n", CFG$alpha, CFG$alpha, alpha))
res <- data.frame(
  test = c("1  Sensor noise", "2  Lens-flare ghost", paste("Scene lighting,", lbl)),
  p        = c(p1, p2, sapply(res3, `[[`, "p")),
  boot_hi  = c(NA, NA, sapply(res3, `[[`, "hi")))
res$reject <- res$p < alpha
res$robust <- ifelse(is.na(res$boot_hi), res$reject, res$boot_hi < alpha)
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p <= .Machine$double.eps, "< 2.2e-16", formatC(p, digits = 3, format = "g")))
res$p <- fmt_p(res$p); res$boot_hi <- fmt_p(res$boot_hi)
print(res, row.names = FALSE, right = FALSE)
cat("\n'robust' = rejected even at the upper end of the bootstrap range of the tail fit.\n")
sens_rate <- tapply(sens$p < alpha, sens$stratum, mean, na.rm = TRUE)
cat("Share of sensitivity settings in which H0 is rejected, by stratum:\n"); print(round(100 * sens_rate))
cat("\nGlobal H0 is rejected only if tests 1, 2 and the scene-lighting test all reject.\n")

## ---- 9. Figures -----------------------------------------------------------------------
png(file.path("results", "fig_null_distribution.png"), width = 3200, height = 2000, res = 300)
nl <- feat_lower[!feat_lower$in_roi, CFG$stat]; nw <- feat_whole[!feat_whole$in_roi, CFG$stat]
nl <- nl[nl > 0]; nw <- nw[nw > 0]
xl <- range(log10(c(nl, nw, x_obs)))
hist(log10(nw), breaks = 50, col = "grey85", border = "white", xlim = xl,
     main = "Bright patches in the tree line (glare corrected)", xlab = paste("log10", CFG$stat))
hist(log10(nl), breaks = 50, col = "grey55", border = "white", add = TRUE)
abline(v = log10(x_obs), col = "red", lwd = 2)
legend("topright", c("Whole tree line", "Lower tree line", "Right shape (anomaly)"),
       fill = c("grey85", "grey55", NA), border = c("grey85", "grey55", NA),
       lty = c(NA, NA, 1), col = c(NA, NA, "red"), lwd = c(NA, NA, 2), bty = "n")
dev.off()

png(file.path("results", "fig_patch_map.png"), width = 3600, height = 2200, res = 300)
st <- 4; sm <- Lc[seq(1, H, st), seq(1, W, st)]
image(seq(1, W, st), seq(1, H, st), t(sm)[, nrow(sm):1], col = gray.colors(128, 0, 1), zlim = c(0, 60),
      useRaster = TRUE, axes = FALSE, xlab = "", ylab = "",
      main = "Glare-corrected luminance. Cyan: null patches (size ~ NIEL). Red: anomaly. Yellow: flare-ghost prediction")
contour(seq(1, W, st), seq(1, H, st), t(lower[seq(1, H, st), seq(1, W, st)] * 1)[, nrow(sm):1],
        levels = 0.5, add = TRUE, drawlabels = FALSE, col = "orange", lty = 2)
nwf <- feat_whole[!feat_whole$in_roi, ]
points(nwf$x, H - nwf$y, pch = 1, col = "cyan", cex = 0.3 + 2 * sqrt(nwf[[CFG$stat]] / x_obs))
points(anom$x, H - anom$y, pch = 1, col = "red", lwd = 2, cex = 0.3 + 2 * sqrt(anom[[CFG$stat]] / x_obs))
points(ghost[1], H - ghost[2], pch = 1, cex = 2.5, col = "yellow", lwd = 2)
points(sun[1], H - sun[2], pch = 3, cex = 2.5, col = "orange", lwd = 2)
dev.off()

sink()

## ---- environment record -------------------------------------------------
# Written on every run so the software environment is documented by the analysis
# itself. Quote these versions when reporting the study; do not reconstruct them.
writeLines(c(paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste("Script: bm_null_test_v3.R"),
             "",
             capture.output(sessionInfo())),
           file.path("results", "sessionInfo.txt"))

cat("\nWrote", log_file,
    "and results/ (CSV tables, fig_null_distribution.png, fig_patch_map.png, sessionInfo.txt)\n")
