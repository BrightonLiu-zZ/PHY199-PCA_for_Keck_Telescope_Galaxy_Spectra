library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(minpack.lm)
library(tidyverse)

LIST_PATH   <- "list_file_9.txt"
MATRIX_PATH <- "flux_matrix_9.txt"   # your full matrix file (not the 10-row sample)
SCORES_PATH <- "pc_scores_9.csv"     # has columns: row_index, PC1, PC2, PC3

list_lines <- read_lines(LIST_PATH)

# Parse with a regex that captures filename + numbers
# (robust to spaces after colons/commas)
extract_fields <- function(x) {
  # Pull "Spec1D File", "Redshift", "lambda_min", "lambda_max"
  spec  <- str_match(x, "Spec1D File:\\s*([^,]+)")[,2]
  z     <- str_match(x, "Redshift:\\s*([0-9eE+\\.-]+)")[,2]
  lmin  <- str_match(x, "lambda_min:\\s*([0-9eE+\\.-]+)")[,2]
  lmax  <- str_match(x, "lambda_max:\\s*([0-9eE+\\.-]+)")[,2]
  tibble(fname = spec, z = as.numeric(z),
         lam_min_obs = as.numeric(lmin),
         lam_max_obs = as.numeric(lmax))
}

meta <- bind_rows(lapply(list_lines, extract_fields)) %>%
  filter(!is.na(z), !is.na(lam_min_obs), !is.na(lam_max_obs))

# Compute rest-frame coverage
meta <- meta %>%
  mutate(
    lam_min_rest = lam_min_obs/(1+z),
    lam_max_rest = lam_max_obs/(1+z),
    row_index    = row_number() - 1L     # 0-based to match the CSV
  )

# Each row = one spectrum, columns = flux samples
flux_matrix <- as.matrix(read.table(MATRIX_PATH, header = FALSE))

pc_scores <- read_csv(SCORES_PATH, show_col_types = FALSE) %>%
  select(row_index, PC3)

meta <- meta %>%
  left_join(pc_scores, by = "row_index")

# Sanity check
stopifnot(nrow(meta) == nrow(flux_matrix))
stopifnot(!any(is.na(meta$PC3)))

g_const <- function(l, A, mu, sig, C) {
  C + A * exp(- (l - mu)^2 / (2 * sig^2))
}

oii_doublet <- function(l, A1, A2, mu1, mu2, sig, C) {
  C + A1 * exp(- (l - mu1)^2 / (2 * sig^2)) +
    A2 * exp(- (l - mu2)^2 / (2 * sig^2))
}

OII1  <- 3726.03
OII2  <- 3728.82
OIII  <- 5008.240
sqrt2pi <- sqrt(2*pi)

wavelengths <- seq(3683, 5300, by = 0.3)

# Checks
stopifnot(ncol(flux_matrix) == length(wavelengths))

# Loop over spectra and fit lines
records <- list()

for (i in seq_len(nrow(meta))) {
  row <- meta[i, ]
  
  # Require both lines present in rest-frame span
  if (!(row$lam_min_rest < 3706 && row$lam_max_rest > 5028)) next
  
  spec <- as.numeric(flux_matrix[i, ])
  
  # [O III] single Gaussian (window: ±10 Å)
  msk3 <- wavelengths > (OIII - 10) & wavelengths < (OIII + 10)
  x3   <- wavelengths[msk3]; y3 <- spec[msk3]
  if (length(x3) < 5) next
  
  # Start values (like Python p0)
  p0_OIII <- list(A = max(y3, na.rm = TRUE) - min(y3, na.rm = TRUE),
                  mu = OIII, sig = 1.0, C = median(y3, na.rm = TRUE))
  fit3 <- try(nlsLM(y3 ~ g_const(x3, A, mu, sig, C),
                    start = p0_OIII,
                    control = nls.lm.control(maxiter = 200)),
              silent = TRUE)
  if (inherits(fit3, "try-error")) next
  
  co3 <- coef(fit3)
  F_OIII <- co3["A"] * co3["sig"] * sqrt2pi
  
  # [O II] doublet (window: [3726.03-8, 3728.82+8] Å)
  msk2 <- wavelengths > (OII1 - 8) & wavelengths < (OII2 + 8)
  x2   <- wavelengths[msk2]; y2 <- spec[msk2]
  if (length(x2) < 5) next
  
  # Starts/bounds (match Python)
  amp_guess <- (max(y2, na.rm = TRUE) - min(y2, na.rm = TRUE))/2
  p0_OII <- list(A1 = amp_guess, A2 = amp_guess, mu1 = OII1, mu2 = OII2,
                 sig = 1.0, C = median(y2, na.rm = TRUE))
  
  lower_OII <- c(A1 = 0,   A2 = 0,   mu1 = OII1 - 1, mu2 = OII2 - 1,
                 sig = 0.3, C = -Inf)
  upper_OII <- c(A1 = Inf, A2 = Inf, mu1 = OII1 + 1, mu2 = OII2 + 1,
                 sig = 3.0, C =  Inf)
  
  fit2 <- try(nlsLM(y2 ~ oii_doublet(x2, A1, A2, mu1, mu2, sig, C),
                    start = p0_OII, lower = lower_OII, upper = upper_OII,
                    control = nls.lm.control(maxiter = 400)),
              silent = TRUE)
  if (inherits(fit2, "try-error")) next
  
  co2 <- coef(fit2)
  F_OII <- (co2["A1"] + co2["A2"]) * co2["sig"] * sqrt2pi
  
  # Flux cuts
  if (F_OII <= 50 || F_OIII <= 80) next
  
  records[[length(records) + 1L]] <- tibble::tibble(
    idx = i - 1L,                # 0-based index like Python row_index
    PC3 = row$PC3,
    flux_OII  = as.numeric(F_OII),
    flux_OIII = as.numeric(F_OIII),
    ratio     = as.numeric(F_OII / F_OIII)
  )
}

results <- dplyr::bind_rows(records)
#cat("Spectra kept:", nrow(results), "\n")

# Remove PC3 outliers and fit y ~ x (ratio vs PC3)
LOW_CUT  <- -300
HIGH_CUT <-  500

# Keep only numeric rows we need
clean_base <- results %>%
  dplyr::select(PC3, ratio) %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::filter(is.finite(PC3), is.finite(ratio))

# Remove PC3 outliers
clean <- clean_base %>%
  dplyr::filter(PC3 >= LOW_CUT, PC3 <= HIGH_CUT)

cat(sprintf("Removed %d outliers by PC3; remaining for fit: %d\n",
            nrow(clean_base) - nrow(clean), nrow(clean)))

stopifnot(nrow(clean) >= 3)

# Linear model in original ratio space
fit <- lm(ratio ~ PC3, data = clean)

# Build prediction grid for a smooth line + 95% CI
x_line <- data.frame(PC3 = seq(min(clean$PC3), max(clean$PC3), length.out = 300))
pred   <- predict(fit, newdata = x_line, interval = "confidence", level = 0.95)
x_line$fit <- pred[, "fit"]
x_line$lwr <- pred[, "lwr"]
x_line$upr <- pred[, "upr"]

# Equation text
co   <- coef(fit)
m    <- unname(co["PC3"])
b    <- unname(co["(Intercept)"])
r2   <- summary(fit)$r.squared
n    <- nrow(clean)
eq_s <- sprintf("y = %.3g·x + %.3g\nR² = %.3f\nN = %d", m, b, r2, n)


# Fit in log space (log10 ratio) and plot on log y-axis
df <- clean_base %>%
  filter(PC3 >= LOW_CUT, PC3 <= HIGH_CUT, ratio > 0)

x_rng <- range(df$PC3)
fit_log <- lm(log10(ratio) ~ PC3, data = df)

# Build smooth line in original scale from log fit
x_line2 <- data.frame(PC3 = seq(x_rng[1], x_rng[2], length.out = 300))
ylog    <- predict(fit_log, newdata = x_line2)          # log10
y_line  <- 10^(ylog)

K <- 1

base <- results %>%
  dplyr::select(PC3, ratio) %>%
  dplyr::mutate(across(everything(), as.numeric)) %>%
  dplyr::filter(is.finite(PC3), is.finite(ratio),
                PC3 >= LOW_CUT, PC3 <= HIGH_CUT,
                ratio > 0)

# Robust Z with MAD in log10 space
ylog <- log10(base$ratio)
med  <- median(ylog, na.rm = TRUE)
madv <- mad(ylog, constant = 1.4826, na.rm = TRUE)  # scaled-to-sigma

if (isTRUE(all.equal(madv, 0))) {
  q1  <- quantile(ylog, 0.25, na.rm = TRUE)
  q3  <- quantile(ylog, 0.75, na.rm = TRUE)
  iqr <- max(q3 - q1, 1e-6)
  robust_z <- abs((ylog - median(ylog, na.rm = TRUE)) / (iqr / 1.349))
} else {
  robust_z <- abs((ylog - med) / madv)
}

k <- min(K, max(0, nrow(base) - 3))
drop_idx <- order(robust_z, decreasing = TRUE)[seq_len(k)]
df2 <- if (k > 0) base[-drop_idx, , drop = FALSE] else base

# Fit: log10(ratio) ~ PC3
fit2 <- lm(log10(ratio) ~ PC3, data = df2)
co2  <- coef(fit2); m <- unname(co2["PC3"]); b <- unname(co2["(Intercept)"])
r2   <- summary(fit2)$r.squared

# Build line + 95% CI in log space, then back-transform
x_line3 <- data.frame(PC3 = seq(min(df2$PC3), max(df2$PC3), length.out = 300))
pred2   <- predict(fit2, newdata = x_line3, interval = "confidence", level = 0.95)

y_line      <- 10^(pred2[,"fit"])
ci_lower    <- 10^(pred2[,"lwr"])
ci_upper    <- 10^(pred2[,"upr"])

eq_s2 <- sprintf("log10(y) = %.3g·x + %.3g   R² = %.3f", m, b, r2)

p_A <- ggplot(df2, aes(x = PC3, y = ratio)) + # <--- ASSIGNED TO p_A
  geom_point(size = 1.8, alpha = 0.7) +
  geom_line(data = cbind(x_line3, ratio = y_line),
            aes(x = PC3, y = ratio), linewidth = 1, inherit.aes = FALSE) +
  geom_ribbon(data = cbind(x_line3, lwr = ci_lower, upr = ci_upper),
              aes(x = PC3, ymin = lwr, ymax = upr),
              alpha = 0.2, inherit.aes = FALSE) +
  scale_y_log10() +
  labs(title = "PC 3 vs log([O II]/[O III])",
       x = "PC 3 Amplitude",
       y = "log10([O II] / [O III])") +
  annotate("label", x = min(df2$PC3), y = max(df2$ratio),
           label = eq_s2, hjust = 0, vjust = 1,
           label.size = 0.3, alpha = 0.9) +
  theme_minimal()


LIST_PATH   <- "C:/PHY199/list_file_9.txt"
MATRIX_PATH <- "C:/PHY199/flux_matrix_9.txt"

lines <- readLines(LIST_PATH)
rows  <- lapply(lines[nzchar(trimws(lines))], function(ln) {
  kvs <- strsplit(ln, ",")[[1]]
  kvs <- lapply(kvs, function(kv) {
    parts <- strsplit(kv, ":", fixed = TRUE)[[1]]
    key <- trimws(parts[1]); val <- trimws(paste(parts[-1], collapse=":"))
    setNames(list(val), key)
  })
  as_tibble(do.call(c, kvs))
})
meta <- bind_rows(rows) |>
  rename(fname = `Spec1D File`,
         z = Redshift,
         lam_min_obs = lambda_min,
         lam_max_obs = lambda_max) |>
  mutate(across(c(z, lam_min_obs, lam_max_obs), as.numeric),
         lam_min_rest = lam_min_obs/(1+z),
         lam_max_rest = lam_max_obs/(1+z))

# flux matrix (rows = spectra, columns = wavelength pixels). Space-separated, no header.
flux_matrix <- as.matrix(read.table(MATRIX_PATH, header = FALSE))
# wavelength grid used in the notebook:
wavelengths <- seq(3683, 5300.0, by = 0.3)

# PCA on rows-as-samples (same as sklearn PCA with centering only)
# (prcomp centers by default; we keep scale.=FALSE to match sklearn’s default)
pc <- prcomp(flux_matrix, center = TRUE, scale. = FALSE)
scores <- pc$x[, 1:3, drop = FALSE]  # PC1..PC3

# Attach PC3 (and PC1, PC2) to meta for later plotting
meta <- meta |>
  mutate(PC1 = scores[,1],
         PC2 = scores[,2],
         PC3 = scores[,3])

HBETA   <- 4862.68   # Angstrom, Hβ rest-frame line center
WIN_HB  <- 12        # half-window (Å) around Hβ used for fitting
MIN_F_HB <- 0.0      # minimum integrated Hβ flux to keep

# Single Gaussian + constant continuum
gauss_const <- function(lmbda, A, mu, sigma, C) {
  C + A * exp(-(lmbda - mu)^2 / (2 * sigma^2))
}

# Measure "flux_Hb" in a ±WIN_HB window, keeping rows where Hβ is covered in rest-frame
hb_records <- list()

for (i in seq_len(nrow(meta))) {
  if (!(meta$lam_min_rest[i] < HBETA - 20 && meta$lam_max_rest[i] > HBETA + 20)) next
  
  spec <- flux_matrix[i, ]
  msk  <- wavelengths > (HBETA - WIN_HB) & wavelengths < (HBETA + WIN_HB)
  x <- wavelengths[msk]; y <- spec[msk]
  if (length(x) < 12) next
  
  # crude initials like the notebook
  m0_0 <- median(y)
  A0   <- max(0, max(y) - m0_0)
  p0   <- list(A = A0, mu = HBETA, sigma = 1.2, C = m0_0)
  
  # bounded fit via nlsLM
  fit <- try(nlsLM(y ~ gauss_const(x, A, mu, sigma, C),
                   start = p0,
                   lower = c(A = 0,    mu = HBETA - 3, sigma = 0.4, C = -Inf),
                   upper = c(A = Inf,  mu = HBETA + 3, sigma = 4.0, C =  Inf),
                   control = nls.lm.control(maxiter = 300)),
             silent = TRUE)
  if (inherits(fit, "try-error")) next
  
  cf <- coef(fit)
  # integrated flux = area under Gaussian = A * sigma * sqrt(2*pi)
  F_Hb <- as.numeric(cf["A"] * cf["sigma"] * sqrt(2*pi))
  if (F_Hb <= MIN_F_HB) next
  
  hb_records[[length(hb_records)+1]] <- tibble(
    idx = i, PC1 = meta$PC1[i], PC2 = meta$PC2[i], PC3 = meta$PC3[i],
    flux_Hb = F_Hb,
    A = cf["A"], mu = cf["mu"], sigma = cf["sigma"], C = cf["C"]
  )
}

hb <- bind_rows(hb_records)
#cat("Measured simple Hβ for", nrow(hb), "spectra.\n")

```

```{r}
HBETA   <- 4862.68         # Å
WIN_HB  <- 12              # half-window for fitting
SNR_MIN <- 3.0             # keep same as notebook
MU_TOL  <- 3.0             # |mu_em - HBETA| must be ≤ MU_TOL

# Model in centered coordinates: x = lambda - HBETA
# hbeta(x) = (m0 + m1*x) + A_em*exp(-(x - mu0)^2/(2*sig_em^2)) - A_abs*exp(-x^2/(2*sig_abs^2))
hbeta_centered <- function(x, A_em, mu0, sig_em, A_abs, sig_abs, m0, m1) {
  cont <- m0 + m1 * x
  em   <- A_em * exp(- (x - mu0)^2 / (2 * sig_em^2))
  absb <- A_abs * exp(- (x)^2      / (2 * sig_abs^2))
  cont + em - absb
}

# robust noise estimate (MAD → σ)
mad_sigma <- function(v) 1.4826 * stats::mad(v, center = stats::median(v, na.rm=TRUE), constant = 1, na.rm = TRUE)

hb_rows <- list()

for (i in seq_len(nrow(meta))) {
  # need rest-frame coverage to include Hβ ± 20 Å
  if (!(is.finite(meta$lam_min_rest[i]) && is.finite(meta$lam_max_rest[i]))) next
  if (!(meta$lam_min_rest[i] < HBETA - 20 && meta$lam_max_rest[i] > HBETA + 20)) next
  
  y_all <- flux_matrix[i, ]
  x_all <- wavelengths
  
  # fitting window
  msk <- x_all > (HBETA - WIN_HB) & x_all < (HBETA + WIN_HB)
  x <- x_all[msk]
  y <- y_all[msk]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 20) next
  
  # center x around Hβ for numerical stability
  xc <- x - HBETA
  
  # --- initial continuum via OLS on the whole window
  cont_fit <- try(stats::lm(y ~ xc), silent = TRUE)
  if (inherits(cont_fit, "try-error")) next
  m0_0 <- unname(coef(cont_fit)[1]); m1_0 <- unname(coef(cont_fit)[2])
  
  y_cont   <- m0_0 + m1_0 * xc
  y_resid  <- y - y_cont
  
  # emission/absorption amplitude starts from residuals
  Aem0 <- max(0, max(y_resid, na.rm = TRUE))
  Aab0 <- max(0, max(-y_resid, na.rm = TRUE))
  
  # reasonable sigma starts
  sig_em0 <- 1.5
  sig_ab0 <- 3.5
  
  # multi-start for emission center (in Å, relative to HBETA)
  mu0_starts <- c(0.0, +0.8, -0.8)
  
  # bounds (same as your notebook intent; absorption centered at 0 in xc)
  lower <- c(A_em = 0,   mu0 = -3, sig_em = 0.4, A_abs = 0,   sig_abs = 1.0,  m0 = -Inf, m1 = -Inf)
  upper <- c(A_em = Inf, mu0 =  3, sig_em = 4.0, A_abs = Inf, sig_abs = 8.0,  m0 =  Inf, m1 =  Inf)
  
  best <- NULL
  best_sse <- Inf
  
  for (mu0_0 in mu0_starts) {
    p0 <- list(A_em = Aem0, mu0 = mu0_0, sig_em = sig_em0,
               A_abs = Aab0, sig_abs = sig_ab0, m0 = m0_0, m1 = m1_0)
    fit <- try(minpack.lm::nlsLM(
      y ~ hbeta_centered(xc, A_em, mu0, sig_em, A_abs, sig_abs, m0, m1),
      start = p0, lower = lower, upper = upper,
      control = minpack.lm::nls.lm.control(maxiter = 500)
    ), silent = TRUE)
    if (inherits(fit, "try-error")) next
    
    # compute SSE on in-window data
    yhat <- try(predict(fit), silent = TRUE)
    if (inherits(yhat, "try-error")) next
    if (length(yhat) != length(y) || any(!is.finite(yhat))) next
    sse <- sum((y - yhat)^2, na.rm = TRUE)
    if (is.finite(sse) && sse < best_sse) {
      best <- fit; best_sse <- sse
    }
  }
  
  if (is.null(best)) next
  cf <- coef(best)
  need <- c("A_em","mu0","sig_em","A_abs","sig_abs","m0","m1")
  if (any(!is.finite(cf[need]))) next
  
  # back to absolute mu_em
  mu_em <- as.numeric(HBETA + cf["mu0"])
  
  # integrated net flux (area emission − area absorption)
  area_em <- as.numeric(cf["A_em"] * cf["sig_em"] * sqrt(2*pi))
  area_ab <- as.numeric(cf["A_abs"] * cf["sig_abs"] * sqrt(2*pi))
  F_Hb_net <- area_em - area_ab
  
  # noise from *sidebands* (|xc| in [6, 12] Å), to avoid line core
  sb <- (abs(xc) >= 6) & (abs(xc) <= 12)
  noise <- if (any(sb)) mad_sigma(y[sb] - (cf["m0"] + cf["m1"] * xc[sb])) else mad_sigma(y - (cf["m0"] + cf["m1"] * xc))
  if (!is.finite(noise) || noise <= 0) next
  
  snr <- as.numeric(cf["A_em"] / noise)
  if (!is.finite(F_Hb_net) || !is.finite(snr) || !is.finite(mu_em)) next
  
  # quality cuts (same as notebook semantics)
  if (F_Hb_net <= 0) next
  if (snr < SNR_MIN) next
  if (abs(mu_em - HBETA) > MU_TOL) next
  
  hb_rows[[length(hb_rows) + 1L]] <- tibble::tibble(
    idx = i,
    PC1 = meta$PC1[i], PC2 = meta$PC2[i], PC3 = meta$PC3[i],
    F_Hb_net = F_Hb_net,
    A_em = cf["A_em"], mu_em = mu_em, sig_em = cf["sig_em"],
    A_abs = cf["A_abs"], sig_abs = cf["sig_abs"],
    m0 = cf["m0"], m1 = cf["m1"], snr = snr
  )
}

hb <- dplyr::bind_rows(hb_rows)

library(ggplot2)
library(dplyr)
library(tidyr)

# population-style standardization (ddof = 0, to mirror sklearn)
z_pop <- function(v) {
  mu <- mean(v, na.rm = TRUE)
  sd_pop <- sqrt(mean((v - mu)^2, na.rm = TRUE))
  (v - mu) / sd_pop
}

DF    <- hb
Y_COL <- "F_Hb_net"

df1 <- DF %>%
  select(PC1, PC2, PC3, all_of(Y_COL)) %>%
  drop_na() %>%
  mutate(
    PC1_z = z_pop(PC1),
    PC2_z = z_pop(PC2),
    PC3_z = z_pop(PC3)
  )

# multivariable regression on standardized PCs
fit_multi <- lm(reformulate(c("PC1_z","PC2_z","PC3_z"), Y_COL), data = df1)
coefs <- coef(fit_multi)
a0 <- unname(coefs[1])
a1 <- unname(coefs["PC1_z"])
a2 <- unname(coefs["PC2_z"])
a3 <- unname(coefs["PC3_z"])

# x-axis = linear combination using learned weights
df1 <- df1 %>%
  mutate(x_pred = a1*PC1_z + a2*PC2_z + a3*PC3_z)

# simple y ~ x fit for plotting line + 95% CI
fit_line <- lm(reformulate("x_pred", Y_COL), data = df1)

xnew <- tibble(x_pred = seq(min(df1$x_pred), max(df1$x_pred), length.out = 300))
pred <- predict(fit_line, newdata = xnew, se.fit = TRUE)
tval <- qt(0.975, df = df.residual(fit_line))
xnew <- xnew %>%
  mutate(y = pred$fit,
         y_lo = pred$fit - tval*pred$se.fit,
         y_hi = pred$fit + tval*pred$se.fit)

r2_line <- summary(fit_line)$r.squared

p_B <- ggplot(df1, aes(x = x_pred, y = .data[[Y_COL]])) + # <--- ASSIGNED TO p_B
  geom_point(size = 1.5, alpha = 0.7, stroke = 0.2, shape = 21, fill = "white") +
  geom_line(data = xnew, aes(x = x_pred, y = y), inherit.aes = FALSE) +
  geom_ribbon(data = xnew, aes(x = x_pred, ymin = y_lo, ymax = y_hi),
              alpha = 0.2, inherit.aes = FALSE) +
  labs(
    x = "a·PC1_z + b·PC2_z + c·PC3_z",
    y = "F(Hβ)",
    title = "Linear Combination of standardized PC1, PC2, and PC3 vs F(Hβ)"
  ) +
  annotate(
    "text",
    x = min(df1$x_pred), y = max(df1[[Y_COL]]),
    hjust = 0, vjust = 1,
    label = sprintf("F(Hβ) = %.3g %+ .3g·PC1_z %+ .3g·PC2_z %+ .3g·PC3_z\nR² = %.3f",
                    a0, a1, a2, a3, r2_line),
    size = 3.5
  ) +
  theme_minimal()

DF    <- hb
Y_COL <- "F_Hb_net"

df2 <- DF %>%
  select(PC1, PC2, all_of(Y_COL)) %>%
  drop_na() %>%
  mutate(
    PC1_z = z_pop(PC1),
    PC2_z = z_pop(PC2)
  )

fit_multi_2 <- lm(reformulate(c("PC1_z","PC2_z"), Y_COL), data = df2)
coefs2 <- coef(fit_multi_2)
b0 <- unname(coefs2[1])
b1 <- unname(coefs2["PC1_z"])
b2 <- unname(coefs2["PC2_z"])

df2 <- df2 %>%
  mutate(x_pred = b1*PC1_z + b2*PC2_z)

fit_line_2 <- lm(reformulate("x_pred", Y_COL), data = df2)

xnew2 <- tibble(x_pred = seq(min(df2$x_pred), max(df2$x_pred), length.out = 300))
pred2 <- predict(fit_line_2, newdata = xnew2, se.fit = TRUE)
tval2 <- qt(0.975, df = df.residual(fit_line_2))
xnew2 <- xnew2 %>%
  mutate(y = pred2$fit,
         y_lo = pred2$fit - tval2*pred2$se.fit,
         y_hi = pred2$fit + tval2*pred2$se.fit)

r2_line_2 <- summary(fit_line_2)$r.squared

p_C <- ggplot(df2, aes(x = x_pred, y = .data[[Y_COL]])) + # <--- ASSIGNED TO p_C
  geom_point(size = 1.5, alpha = 0.7, stroke = 0.2, shape = 21, fill = "white") +
  geom_line(data = xnew2, aes(x = x_pred, y = y), inherit.aes = FALSE) +
  geom_ribbon(data = xnew2, aes(x = x_pred, ymin = y_lo, ymax = y_hi),
              alpha = 0.2, inherit.aes = FALSE) +
  labs(
    x = "a·PC1_z + b·PC2_z",
    y = "F(Hβ)",
    title = "Linear Combination of standardized PC1 and PC2 vs F(Hβ)"
  ) +
  annotate(
    "text",
    x = min(df2$x_pred), y = max(df2[[Y_COL]]),
    hjust = 0, vjust = 1,
    label = sprintf("F(Hβ) = %.3g %+ .3g·PC1_z %+ .3g·PC2_z\nR² = %.3f",
                    b0, b1, b2, r2_line_2),
    size = 3.5
  ) +
  theme_minimal()

M   <- flux_matrix
lam <- wavelengths

ROW_KEY <- if ("idx" %in% names(hb)) "idx" else "matrix_row"
hb2 <- hb
if (!("idx" %in% names(hb2))) hb2$matrix_row <- seq_len(nrow(hb2))

targets <- c(25, 400, 850)

chosen_idx <- integer()
for (t in targets) {
  remaining_rows <- setdiff(seq_len(nrow(hb2)), chosen_idx)
  if (!length(remaining_rows)) break
  remaining <- hb2[remaining_rows, , drop = FALSE]
  j <- which.min(abs(remaining[[Y_COL]] - t))
  orig_row <- remaining_rows[j]
  chosen_idx <- c(chosen_idx, orig_row)
}
ex <- hb2[chosen_idx, c("PC1","PC2","PC3", Y_COL, ROW_KEY), drop = FALSE]

SMOOTH_SIGMA_PIX <- 3.0
gauss1d <- function(y, sigma_pix) {
  rad <- ceiling(4.0 * sigma_pix)
  xx  <- seq(-rad, rad)
  ker <- exp(-(xx^2) / (2 * sigma_pix^2)); ker <- ker / sum(ker)
  as.numeric(stats::filter(y, ker, sides = 2, circular = FALSE))
}
specs_sm <- lapply(as.integer(ex[[ROW_KEY]]), function(ridx) gauss1d(M[ridx, ], SMOOTH_SIGMA_PIX))

robust_span <- function(y) diff(stats::quantile(y, c(0.05, 0.95), na.rm = TRUE))
offset_step <- if (length(specs_sm)) 0.8 * stats::median(vapply(specs_sm, robust_span, numeric(1))) else 0
offsets <- seq_along(specs_sm) * offset_step

plot_df <- dplyr::bind_rows(lapply(seq_along(specs_sm), function(i) {
  tibble::tibble(wavelength = lam,
                 flux = specs_sm[[i]] + offsets[i],
                 label = sprintf("F(Hβ)=%.0f, PC1=%.2f, PC2=%.2f, PC3=%.2f",
                                 ex[[Y_COL]][i], ex$PC1[i], ex$PC2[i], ex$PC3[i]))
}))

elem_lines <- c(
  `3727.0921`="OII", `3771`="Hι", `3798`="Hθ", `3835.40`="Hη",
  `3869`="NeIII", `3889.06`="Hζ", `3970.08`="Hε", `4102.89`="Hδ",
  `4341.68`="Hγ", `4862.68`="Hβ", `4932.603`="OIII", `4960.295`="OIII",
  `5008.240`="OIII", `3934.777`="K", `3969.588`="H", `4305.61`="G",
  `5176.7`="Mg", `5269.5`="FeI"
) |> as.list()

p_D <- ggplot(plot_df, aes(x = wavelength, y = flux, color = label)) + # <--- ASSIGNED TO p_D
  geom_line(linewidth = 0.01, show.legend = TRUE) +
  scale_color_discrete(name = NULL) +
  labs(x = "Wavelength (Å)", y = "Flux  (Gaussian-smoothed)",
       title = "Example spectra (stacked, Gaussian-smoothed)\nnear F(Hβ) ≈ 25, 400, 850") +
  theme_minimal()

lam_min <- min(lam); lam_max <- max(lam)
for (wl_str in names(elem_lines)) {
  wl <- as.numeric(wl_str)
  if (!is.na(wl) && wl >= lam_min && wl <= lam_max) {
    p_D <- p_D + geom_vline(xintercept = wl, linetype = "dashed", linewidth = 0.05, alpha = 0.6) + # <--- MODIFIED p_D
      annotate("text", x = wl, y = 250,
               label = paste0(" ", elem_lines[[wl_str]]),
               angle = 90, vjust = 0, hjust = 0, size = 3)
  }
}



lam <- wavelengths
stopifnot(is.numeric(lam), length(lam) == ncol(flux_matrix))

pcs_scores <- readr::read_csv("pc_scores_9.csv", show_col_types = FALSE)

# helper to pick column by common name variants
pick_col <- function(df, patterns) {
  nms <- names(df)
  hit <- nms[Reduce(`|`, lapply(patterns, function(p) grepl(p, nms, ignore.case = TRUE)))]
  if (length(hit) == 0) NA_character_ else hit[1]
}
c_PC1 <- pick_col(pcs_scores, c("^PC1$", "^PC_?1$"))
c_PC2 <- pick_col(pcs_scores, c("^PC2$", "^PC_?2$"))
c_PC3 <- pick_col(pcs_scores, c("^PC3$", "^PC_?3$"))
stopifnot(!is.na(c_PC1), !is.na(c_PC2), !is.na(c_PC3))

PC1_scores_csv <- pcs_scores[[c_PC1]]
PC2_scores_csv <- pcs_scores[[c_PC2]]
PC3_scores_csv <- pcs_scores[[c_PC3]]

# Sanity: these are per-spectrum
stopifnot(length(PC1_scores_csv) == nrow(flux_matrix),
          length(PC2_scores_csv) == nrow(flux_matrix),
          length(PC3_scores_csv) == nrow(flux_matrix))

# prcomp centers variables (columns) by default; scale.=FALSE to match sklearn default
pc <- prcomp(flux_matrix, center = TRUE, scale. = FALSE)


L1 <- pc$rotation[,1]; L2 <- pc$rotation[,2]; L3 <- pc$rotation[,3]
sd1 <- pc$sdev[1];     sd2 <- pc$sdev[2];     sd3 <- pc$sdev[3]

# Sign-align to your CSV scores
s1 <- sign(cor(pc$x[,1], PC1_scores_csv, use = "pairwise.complete.obs")); if (s1 == 0) s1 <- 1
s2 <- sign(cor(pc$x[,2], PC2_scores_csv, use = "pairwise.complete.obs")); if (s2 == 0) s2 <- 1
s3 <- sign(cor(pc$x[,3], PC3_scores_csv, use = "pairwise.complete.obs")); if (s3 == 0) s3 <- 1

# This is the flux contribution of a 1-SD change in each PC score.
L1_flux <- s1 * L1 * sd1
L2_flux <- s2 * L2 * sd2
L3_flux <- s3 * L3 * sd3

plot_df_PC <- dplyr::bind_rows(
  tibble::tibble(wavelength = lam, flux = L1_flux, label = "PC1"),
  tibble::tibble(wavelength = lam, flux = L2_flux, label = "PC2"),
  tibble::tibble(wavelength = lam, flux = L3_flux, label = "PC3")
)

elem_lines <- c(
  `3727.0921`="OII", `3771`="Hι", `3798`="Hθ", `3835.40`="Hη",
  `3869`="NeIII", `3889.06`="Hζ", `3970.08`="Hε", `4102.89`="Hδ",
  `4341.68`="Hγ", `4862.68`="Hβ", `4932.603`="OIII", `4960.295`="OIII",
  `5008.240`="OIII", `3934.777`="K", `3969.588`="H", `4305.61`="G",
  `5176.7`="Mg", `5269.5`="FeI"
) |> as.list()

p_PC_flux <- ggplot(plot_df_PC, aes(x = wavelength, y = flux, color = label)) +
  geom_line(linewidth = 0.01, show.legend = TRUE) +
  scale_color_discrete(name = NULL) +
  labs(
    x = "Wavelength (Å)",
    y = "Flux",
    title = "Stacked PC1, PC2, and PC3"
  ) +
  theme_minimal()

lam_min <- min(lam, na.rm = TRUE); lam_max <- max(lam, na.rm = TRUE)
y_top   <- max(plot_df_PC$flux, na.rm = TRUE)

for (wl_str in names(elem_lines)) {
  wl <- as.numeric(wl_str)
  if (!is.na(wl) && wl >= lam_min && wl <= lam_max) {
    p_PC_flux <- p_PC_flux +
      geom_vline(xintercept = wl, linetype = "dashed", linewidth = 0.01, alpha = 0.6) +
      annotate("text", x = wl, y = 80,
               label = paste0(" ", elem_lines[[wl_str]]),
               angle = 90, vjust = 0, hjust = 0, size = 3)
  }
}

library(patchwork)

# Arrange the four plots: (p_A | p_B) / (p_C | p_D)
# Use plot_layout to give more width to each column
panel <-p_PC_flux / (p_A | p_B) / (p_C | p_D) +
  plot_layout(widths = c(1, 1)) + # Ensures both columns have equal width, adjust if one needs to be wider
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

print(panel)

# Save Overleaf-ready outputs
pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave("PHY199_Figure1.pdf", plot = panel,
       width = 15, height = 10, units = "in", device = pdf_dev) # Increased width
ggsave("PHY199_Figure1.png", plot = panel,
       width = 15, height = 10, units = "in", dpi = 300) # Increased width



