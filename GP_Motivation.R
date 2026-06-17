# ==============================================================================
# GP_Motivation_Variogramme_nach_Residualanalyse.R
# Motivation von Gaußprozessen mit Variogrammen AUFBAUEND auf Residualanalyse.R
# PDF-Ausgabe: NUR PLOTS
# ==============================================================================

# ==============================================================================
# 0) PAKETE UND CHECKS
# ==============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Bitte installiere ggplot2: install.packages('ggplot2')")
}

library(ggplot2)

required_objects <- c(
  "df_val", "df_train", "emulator", "time", "lon", "lat",
  "time_stats", "bias_df", "rmse_df", "sd_df"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]

if (length(missing_objects) > 0) {
  stop(
    "Diese Objekte fehlen: ", paste(missing_objects, collapse = ", "),
    "\nBitte zuerst ausführen:",
    "\n  source('vegemul_gam_df.R')",
    "\n  source('Residualanalyse.R')"
  )
}

required_df_val_cols <- c(
  "tf", "pred", "resid_response", "time_index",
  "lon_index", "lat_index", "lon", "lat"
)

missing_cols <- required_df_val_cols[!required_df_val_cols %in% names(df_val)]

if (length(missing_cols) > 0) {
  stop(
    "In df_val fehlen diese Spalten: ", paste(missing_cols, collapse = ", "),
    "\nBitte prüfen, ob Residualanalyse.R erfolgreich durchgelaufen ist."
  )
}

if (!"resid_response" %in% names(df_train)) {
  stop("In df_train fehlt resid_response. Bitte zuerst Residualanalyse.R ausführen.")
}

if (!exists("theme_simple")) {
  theme_simple <- theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

set.seed(123)

# ==============================================================================
# 1) HILFSFUNKTION: EMPIRISCHES SEMIVARIOGRAMM
# ==============================================================================

empirical_variogram <- function(data,
                                coord_cols,
                                value_col,
                                n_pairs = 200000,
                                n_bins = 18,
                                max_dist_quantile = 0.95,
                                standardize_coords = FALSE) {
  
  dat <- data[, c(coord_cols, value_col), drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) < 10) {
    stop("Zu wenige vollständige Datenpunkte für ein Variogramm.")
  }
  
  coords <- as.matrix(dat[, coord_cols, drop = FALSE])
  values <- dat[[value_col]]
  
  if (standardize_coords) {
    coord_center <- apply(coords, 2, mean, na.rm = TRUE)
    coord_sd <- apply(coords, 2, sd, na.rm = TRUE)
    coord_sd[is.na(coord_sd) | coord_sd == 0] <- 1
    coords <- scale(coords, center = coord_center, scale = coord_sd)
  }
  
  n <- nrow(coords)
  
  i <- sample.int(n, n_pairs, replace = TRUE)
  j <- sample.int(n, n_pairs, replace = TRUE)
  
  keep <- i != j
  i <- i[keep]
  j <- j[keep]
  
  d <- sqrt(rowSums((coords[i, , drop = FALSE] - coords[j, , drop = FALSE])^2))
  gamma <- 0.5 * (values[i] - values[j])^2
  
  max_dist <- as.numeric(quantile(d, probs = max_dist_quantile, na.rm = TRUE))
  
  keep <- is.finite(d) &
    is.finite(gamma) &
    d > 0 &
    d <= max_dist
  
  d <- d[keep]
  gamma <- gamma[keep]
  
  breaks <- seq(0, max_dist, length.out = n_bins + 1)
  bin <- cut(d, breaks = breaks, include.lowest = TRUE)
  
  out <- aggregate(
    cbind(distance = d, semivariance = gamma) ~ bin,
    FUN = mean,
    na.rm = TRUE
  )
  
  n_per_bin <- aggregate(gamma ~ bin, FUN = length)
  names(n_per_bin)[2] <- "n_pairs"
  
  out <- merge(out, n_per_bin, by = "bin")
  out <- out[order(out$distance), ]
  rownames(out) <- NULL
  
  return(out)
}

# ==============================================================================
# 2) DATEN FÜR VARIOGRAMME
# ==============================================================================

# Räumliches Variogramm der mittleren Residuen pro Gitterzelle
spatial_bias_dat <- bias_df
names(spatial_bias_dat)[names(spatial_bias_dat) == "value"] <- "mean_resid"

# Räumliches Variogramm der Rohresiduen auf einer Stichprobe
max_raw_points <- 6000

raw_idx <- sample(
  seq_len(nrow(df_val)),
  min(max_raw_points, nrow(df_val))
)

spatial_raw_dat <- df_val[
  raw_idx,
  c("lon", "lat", "resid_response"),
  drop = FALSE
]

# Zeitliches Variogramm der zeitlich gemittelten Residuen
temporal_dat <- data.frame(
  time_val = time_stats$time,
  mean_resid = time_stats$mean_resid
)

# Klimaraum-Variogramm
possible_covariates <- c("precip", "tcold", "twarm", "co2")
covariates <- possible_covariates[possible_covariates %in% names(df_val)]

climate_vario_available <- length(covariates) >= 2

if (climate_vario_available) {
  
  max_climate_points <- 6000
  
  climate_idx <- sample(
    seq_len(nrow(df_val)),
    min(max_climate_points, nrow(df_val))
  )
  
  climate_dat <- df_val[
    climate_idx,
    c(covariates, "resid_response"),
    drop = FALSE
  ]
}

# ==============================================================================
# 3) VARIOGRAMME BERECHNEN
# ==============================================================================

spatial_bias_variogram <- empirical_variogram(
  data = spatial_bias_dat,
  coord_cols = c("lon", "lat"),
  value_col = "mean_resid",
  n_pairs = 150000,
  n_bins = 18,
  standardize_coords = FALSE
)

spatial_raw_variogram <- empirical_variogram(
  data = spatial_raw_dat,
  coord_cols = c("lon", "lat"),
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = FALSE
)

long_variogram <- empirical_variogram(
  data = spatial_raw_dat,
  coord_cols = c("lon"),
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = FALSE
)

lat_variogram <- empirical_variogram(
  data = spatial_raw_dat,
  coord_cols = c("lat"),
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = FALSE
)

temporal_variogram <- empirical_variogram(
  data = temporal_dat,
  coord_cols = c("time_val"),
  value_col = "mean_resid",
  n_pairs = 100000,
  n_bins = 14,
  standardize_coords = FALSE
)

if (climate_vario_available) {
  
  climate_variogram <- empirical_variogram(
    data = climate_dat,
    coord_cols = covariates,
    value_col = "resid_response",
    n_pairs = 200000,
    n_bins = 18,
    standardize_coords = TRUE
  )
}

# ==============================================================================
# 4) KENNZAHLEN FÜR KONSOLENAUSGABE
# ==============================================================================

r_val <- df_val$resid_response

bias_global <- mean(r_val, na.rm = TRUE)
rmse_global <- sqrt(mean(r_val^2, na.rm = TRUE))
mae_global <- mean(abs(r_val), na.rm = TRUE)
sd_global <- sd(r_val, na.rm = TRUE)

if (!exists("acf_lag1")) {
  acf_lag1 <- acf(time_stats$mean_resid, lag.max = 1, plot = FALSE)$acf[2]
}

estimate_range_heuristic <- function(vario_df) {
  
  max_gamma <- max(vario_df$semivariance, na.rm = TRUE)
  sill_level <- 0.95 * max_gamma
  
  idx <- which(vario_df$semivariance >= sill_level)
  
  if (length(idx) == 0) {
    return(NA_real_)
  }
  
  return(vario_df$distance[min(idx)])
}

spatial_bias_range <- estimate_range_heuristic(spatial_bias_variogram)
spatial_raw_range <- estimate_range_heuristic(spatial_raw_variogram)
long_range <- estimate_range_heuristic(long_variogram)
lat_range <- estimate_range_heuristic(lat_variogram)
temporal_range <- estimate_range_heuristic(temporal_variogram)

if (climate_vario_available) {
  climate_range <- estimate_range_heuristic(climate_variogram)
}

# ==============================================================================
# 5) PDF ERZEUGEN: NUR PLOTS
# ==============================================================================

pdf("Semivariogramme.pdf", width = 12, height = 8)

# ------------------------------------------------------------------------------
# Plot 1: Räumliches Semivariogramm des mittleren Residual-Bias
# ------------------------------------------------------------------------------

p_spatial_bias <- ggplot(
  spatial_bias_variogram,
  aes(x = distance, y = semivariance)
) +
  geom_point(aes(size = n_pairs), alpha = 0.75) +
  geom_line(linewidth = 1) +
  labs(
    title = "Spatial semivariogram of the mean residual bias",
    x = "Spatial distance in lon/lat units",
    y = expression(gamma(h) == frac(1, 2) * mean((r(u) - r(u+h))^2)),
    size = "Number of pairs"
  ) +
  theme_simple

print(p_spatial_bias)

# ------------------------------------------------------------------------------
# Plot 2: Räumliches Semivariogramm der Validierungsresiduen
# ------------------------------------------------------------------------------

p_spatial_raw <- ggplot(
  spatial_raw_variogram,
  aes(x = distance, y = semivariance)
) +
  geom_point(aes(size = n_pairs), alpha = 0.75) +
  geom_line(linewidth = 1) +
  labs(
    title = "Spatial semivariogram of validation residuals",
    x = "Spatial distance in lon/lat units",
    y = expression(gamma(h)),
    size = "Number of pairs"
  ) +
  theme_simple

print(p_spatial_raw)

# ------------------------------------------------------------------------------
# Plot 3: Semivariogramm entlang der Longitude-Richtung
# ------------------------------------------------------------------------------

p_lon <- ggplot(
  long_variogram,
  aes(x = distance, y = semivariance)
) +
  geom_point(aes(size = n_pairs), alpha = 0.75) +
  geom_line(linewidth = 1) +
  labs(
    title = "Semivariogram along the longitude direction",
    x = "Distance in longitude units",
    y = expression(gamma(h)),
    size = "Number of pairs"
  ) +
  theme_simple

print(p_lon)

# ------------------------------------------------------------------------------
# Plot 4: Semivariogramm entlang der Latitude-Richtung
# ------------------------------------------------------------------------------

p_lat <- ggplot(
  lat_variogram,
  aes(x = distance, y = semivariance)
) +
  geom_point(aes(size = n_pairs), alpha = 0.75) +
  geom_line(linewidth = 1) +
  labs(
    title = "Semivariogram along the latitude direction",
    x = "Distance in latitude units",
    y = expression(gamma(h)),
    size = "Number of pairs"
  ) +
  theme_simple

print(p_lat)

# ------------------------------------------------------------------------------
# Plot 5: Zeitliches Semivariogramm
# ------------------------------------------------------------------------------

p_temporal <- ggplot(
  temporal_variogram,
  aes(x = distance, y = semivariance)
) +
  geom_point(aes(size = n_pairs), alpha = 0.75) +
  geom_line(linewidth = 1) +
  labs(
    title = "Temporal semivariogram of averaged residuals",
    x = "Temporal distance",
    y = expression(gamma(h)),
    size = "Number of pairs"
  ) +
  theme_simple

print(p_temporal)

# ------------------------------------------------------------------------------
# Plot 6: Semivariogramm im Klimaraum
# ------------------------------------------------------------------------------

if (climate_vario_available) {
  
  p_climate <- ggplot(
    climate_variogram,
    aes(x = distance, y = semivariance)
  ) +
    geom_point(aes(size = n_pairs), alpha = 0.75) +
    geom_line(linewidth = 1) +
    labs(
      title = "Semivariogram in standardized climate space",
      subtitle = paste("Covariates:", paste(covariates, collapse = ", ")),
      x = "Distance in standardized covariate space",
      y = expression(gamma(h)),
      size = "Number of pairs"
    ) +
    theme_simple
  
  print(p_climate)
}

dev.off()

# ==============================================================================
# 6) KONSOLENAUSGABE
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("  GP-MOTIVATION MIT VARIOGRAMMEN ABGESCHLOSSEN\n")
cat("  -> Ausgabe: GP_Motivation_Variogramme_nach_Residualanalyse.pdf\n")
cat("  -> Die PDF enthält nur Plots.\n")
cat("====================================================================\n\n")

cat("Globale Residual-Kennzahlen:\n")
cat(sprintf("  Globaler Bias:        %+0.5f\n", bias_global))
cat(sprintf("  Globaler RMSE:         %0.5f\n", rmse_global))
cat(sprintf("  Globaler MAE:          %0.5f\n", mae_global))
cat(sprintf("  SD der Residuen:       %0.5f\n", sd_global))
cat(sprintf("  ACF Lag 1:            %+0.5f\n", as.numeric(acf_lag1)))
cat("\n")

cat("Verwendete Kovariaten für Klimaraum-Variogramm:\n")

if (climate_vario_available) {
  cat("  ", paste(covariates, collapse = ", "), "\n\n")
} else {
  cat("  Keine oder zu wenige Kovariaten gefunden.\n\n")
}

cat("Erzeugte Variogramm-Objekte:\n")
cat("  spatial_bias_variogram\n")
cat("  spatial_raw_variogram\n")
cat("  long_variogram\n")
cat("  lat_variogram\n")
cat("  temporal_variogram\n")

if (climate_vario_available) {
  cat("  climate_variogram\n")
}

cat("\n")

cat("Heuristische Reichweiten:\n")
cat(sprintf("  Räumlicher Bias:        %.3f\n", spatial_bias_range))
cat(sprintf("  Räumliche Rohresiduen:  %.3f\n", spatial_raw_range))
cat(sprintf("  Longitude-Richtung:     %.3f\n", long_range))
cat(sprintf("  Latitude-Richtung:      %.3f\n", lat_range))
cat(sprintf("  Zeitlicher Bias:        %.3f\n", temporal_range))

if (climate_vario_available) {
  cat(sprintf("  Klimaraum:              %.3f\n", climate_range))
}

cat("\n")