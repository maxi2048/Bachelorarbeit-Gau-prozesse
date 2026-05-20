# ==============================================================================
# RESIDUALANALYSE – GAM Vegetation Emulator
# Einfache ggplot2-Version für ersten Überblick
# ==============================================================================

library(ggplot2)

# ==============================================================================
# 0) RESIDUEN BERECHNEN
# ==============================================================================

# Validierung
df_val$resid_response <- df_val$tf - df_val$pred

df_val$resid_pearson <- (df_val$tf - df_val$pred) /
  sqrt(pmax(df_val$pred * (1 - df_val$pred), 1e-8))

# Training
df_train$pred_train <- predict(emulator, newdata = df_train, type = "response")

df_train$resid_response <- df_train$tf - df_train$pred_train

df_train$resid_pearson <- (df_train$tf - df_train$pred_train) /
  sqrt(pmax(df_train$pred_train * (1 - df_train$pred_train), 1e-8))

# Kürzel
r_val  <- df_val$resid_response
rp_val <- df_val$resid_pearson
fit    <- df_val$pred

n_lon <- length(lon)
n_lat <- length(lat)

# Einfaches Theme
theme_simple <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# PDF ÖFFNEN
# ==============================================================================

pdf("residualanalyse_GAM_ggplot2_simple.pdf", width = 12, height = 8)

# ==============================================================================
# 1) VERTEILUNG DER RESIDUEN
# ==============================================================================

p1 <- ggplot() +
  geom_histogram(
    data = df_val,
    aes(x = resid_response, y = after_stat(density)),
    bins = 80,
    fill = "lightblue",
    color = "white",
    alpha = 0.65
  ) +
  
  # Validierung zuerst zeichnen
  geom_density(
    data = subset(df_resid_compare, Datensatz == "Validierung"),
    aes(x = resid, color = Datensatz),
    linewidth = 1.1,
    linetype = "solid"
  ) +
  
  # Training danach zeichnen, damit es sichtbar darüber liegt
  geom_density(
    data = subset(df_resid_compare, Datensatz == "Training"),
    aes(x = resid, color = Datensatz),
    linewidth = 1.3,
    linetype = "dashed"
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = summary_text,
    hjust = 1.05,
    vjust = 1.1,
    size = 3.6,
    label.size = 0.3,
    fill = "white",
    color = "black"
  ) +
  
  scale_color_manual(
    values = c(
      "Training" = "darkorange",
      "Validierung" = "steelblue"
    )
  ) +
  
  labs(
    title = expression("Empirical Distribution of " * tf - widehat(tf)),
    subtitle = "Histogramm: Testdata | Density: Train and Test",
    x = "Residuen",
    y = "Dichte",
    color = "Datensatz"
  ) +
  theme_simple

print(p1)

p4 <- ggplot(df_compare, aes(x = Datensatz, y = resid, fill = Datensatz)) +
  geom_boxplot(outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Boxplot: Training vs. Validierung",
    x = "",
    y = "Residuen"
  ) +
  theme_simple +
  theme(legend.position = "none")

print(p4)

# ==============================================================================
# 3) Q-Q-PLOTS
# ==============================================================================

p5 <- ggplot(df_val, aes(sample = resid_response)) +
  stat_qq(alpha = 0.3, size = 0.7) +
  stat_qq_line(color = "red", linewidth = 1) +
  labs(
    title = "Q-Q-Plot der Response-Residuen",
    x = "Theoretische Quantile",
    y = "Empirische Quantile"
  ) +
  theme_simple

print(p5)

# p6 <- ggplot(df_val, aes(sample = resid_pearson)) +
#   stat_qq(alpha = 0.3, size = 0.7) +
#   stat_qq_line(color = "red", linewidth = 1) +
#   labs(
#     title = "Q-Q-Plot der Pearson-Residuen",
#     x = "Theoretische Quantile",
#     y = "Empirische Quantile"
#   ) +
#   theme_simple
# 
# print(p6)

# ==============================================================================
# 4) RESIDUEN VS. FITTED
# ==============================================================================

p7 <- ggplot(df_val, aes(x = pred, y = resid_response)) +
  geom_point(alpha = 0.08, size = 0.5) +
  
  # Rand des Parallelogramms
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = -1),
               inherit.aes = FALSE,
               color = "blue", linetype = "dotted", linewidth = 1) +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0),
               inherit.aes = FALSE,
               color = "blue", linetype = "dotted", linewidth = 1) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 1),
               inherit.aes = FALSE,
               color = "blue", linetype = "dotted", linewidth = 1) +
  geom_segment(aes(x = 1, y = -1, xend = 1, yend = 0),
               inherit.aes = FALSE,
               color = "blue", linetype = "dotted", linewidth = 1) +
  
  # Referenzlinien
  geom_smooth(method = "lm", se = FALSE, color = "orange", linewidth = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  
  labs(
    title = "Residuen vs. vorhergesagte Werte",
    subtitle = expression("Blauer Rand: zulässiger Bereich für " ~ epsilon == tf - widehat(tf)),
    x = expression("Vorhergesagter Wert " ~ widehat(tf)),
    y = expression("Residuum " ~ epsilon)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
  theme_simple

print(p7)


# Heteroskedastizität nachweisen
df_val$pred_bin <- cut(
  df_val$pred,
  breaks = seq(0, 1, by = 0.1),
  include.lowest = TRUE
)

# Empirische Varianz und Standardabweichung pro Intervall
hetero_stats <- aggregate(
  resid_response ~ pred_bin,
  data = df_val,
  FUN = function(x) c(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    var = var(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    mae = mean(abs(x), na.rm = TRUE)
  )
)

# aggregate gibt eine Matrix-Spalte zurück, daher umformen
hetero_stats <- data.frame(
  pred_bin = hetero_stats$pred_bin,
  hetero_stats$resid_response
)

hetero_stats$varcoef <- hetero_stats$sd/hetero_stats$mean

print(hetero_stats)

# summary(lm(resid_response ~ pred, df_val))



# p8 <- ggplot(df_val, aes(x = pred, y = sqrt(abs(resid_response)))) +
#   geom_point(alpha = 0.08, size = 0.5) +
#   geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
#   labs(
#     title = "Scale-Location-Plot",
#     x = "Vorhergesagter Wert",
#     y = "sqrt(|Residuum|)"
#   ) +
#   theme_simple
# 
# print(p8)

p9 <- ggplot(df_val, aes(x = tf, y = pred)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(
    title = "Beobachtet vs. vorhergesagt",
    x = "Beobachtet: tf",
    y = "Vorhergesagt"
  ) +
  theme_simple

print(p9)

# ==============================================================================
# 5) RESIDUEN VS. KOVARIATEN
# ==============================================================================

covariates <- c("precip", "tcold", "twarm")

for (cov_name in covariates) {
  
  # Grenzen des zentralen 99%-Bereichs
  q_99 <- quantile(
    df_val[[cov_name]],
    probs = c(0.025, 0.925),
    na.rm = TRUE
  )
  
  p <- ggplot(df_val, aes(x = .data[[cov_name]], y = resid_response)) +
    geom_point(alpha = 0.08, size = 0.5) +
    geom_smooth(method = "loess", se = FALSE, color = "orange", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    
    # Vertikale Linien für zentralen 99%-Bereich
    geom_vline(xintercept = q_99[1], color = "blue", linetype = "dotted", linewidth = 1) +
    geom_vline(xintercept = q_99[2], color = "blue", linetype = "dotted", linewidth = 1) +
    
    labs(
      title = paste("Residuen vs.", cov_name),
      subtitle = paste0(
        "Blaue Linien: zentraler 95%-Bereich der Daten [",
        round(q_99[1], 2), ", ",
        round(q_99[2], 2), "]"
      ),
      x = cov_name,
      y = "Residuum"
    ) +
    theme_simple
  
  print(p)
}

# ==============================================================================
# 6) RESIDUEN NACH KOVARIATEN-DEZILEN
# ==============================================================================

for (cov_name in covariates) {
  
  q_breaks <- quantile(df_val[[cov_name]], probs = seq(0, 1, 0.1), na.rm = TRUE)
  q_breaks <- unique(q_breaks)
  
  if (length(q_breaks) > 2) {
    
    df_val$cov_decile_temp <- cut(
      df_val[[cov_name]],
      breaks = q_breaks,
      include.lowest = TRUE
    )
    
    p <- ggplot(df_val, aes(x = cov_decile_temp, y = resid_response)) +
      geom_boxplot(fill = "lightblue", outlier.alpha = 0.15) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      labs(
        title = paste("Residuen nach Dezilen von", cov_name),
        x = paste(cov_name, "Dezil"),
        y = "Residuum"
      ) +
      theme_simple +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)
  }
}

df_val$cov_decile_temp <- NULL

# ==============================================================================
# 7) RESIDUEN NACH TREE-COVER-KLASSEN UND BREITENGRAD
# ==============================================================================

df_val$tf_class <- cut(
  df_val$tf,
  breaks = seq(0, 1, by = 0.1),
  include.lowest = TRUE
)

p <- ggplot(df_val, aes(x = tf_class, y = resid_response)) +
  geom_boxplot(fill = "darkseagreen2", outlier.alpha = 0.15) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Residuen nach beobachtetem Tree Cover",
    x = "Tree-Cover-Klasse",
    y = "Residuum"
  ) +
  theme_simple +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)

df_val$lat_band <- cut(
  df_val$lat,
  breaks = seq(min(df_val$lat, na.rm = TRUE),
               max(df_val$lat, na.rm = TRUE),
               length.out = 13),
  include.lowest = TRUE
)

p <- ggplot(df_val, aes(x = lat_band, y = resid_response)) +
  geom_boxplot(fill = "lightblue", outlier.alpha = 0.15) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Residuen nach Breitengrad-Bändern",
    x = "Breitengrad-Band",
    y = "Residuum"
  ) +
  theme_simple +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

print(p)

# ==============================================================================
# 8) ZEITLICHE ANALYSE
# ==============================================================================

val_times <- sort(unique(df_val$time_index))

time_stats <- data.frame(
  time_index = val_times,
  time_val = time[val_times],
  mean_resid = NA,
  rmse = NA,
  sd_resid = NA,
  mae = NA
)

for (k in seq_along(val_times)) {
  
  t <- val_times[k]
  sub <- df_val$resid_response[df_val$time_index == t]
  
  time_stats$mean_resid[k] <- mean(sub, na.rm = TRUE)
  time_stats$rmse[k]       <- sqrt(mean(sub^2, na.rm = TRUE))
  time_stats$sd_resid[k]   <- sd(sub, na.rm = TRUE)
  time_stats$mae[k]        <- mean(abs(sub), na.rm = TRUE)
}

p <- ggplot(time_stats, aes(x = time_val, y = mean_resid)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Zeitlicher Bias",
    x = "Zeit",
    y = "Durchsch. Residuum"
  ) +
  theme_simple

print(p)

p <- ggplot(time_stats, aes(x = time_val, y = rmse)) +
  geom_line(color = "tomato", linewidth = 1) +
  geom_point(color = "tomato", size = 1.5) +
  labs(
    title = "Zeitlicher RMSE",
    x = "Zeit",
    y = "RMSE"
  ) +
  theme_simple

print(p)

p <- ggplot(time_stats, aes(x = time_val, y = sd_resid)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(color = "purple", size = 1.5) +
  labs(
    title = "Zeitliche Residual-Streuung",
    x = "Zeit",
    y = "SD der Residuen"
  ) +
  theme_simple

print(p)

p <- ggplot(time_stats, aes(x = time_val, y = mae)) +
  geom_line(color = "darkorange", linewidth = 1) +
  geom_point(color = "darkorange", size = 1.5) +
  labs(
    title = "Zeitlicher MAE",
    x = "Zeit",
    y = "MAE"
  ) +
  theme_simple

print(p)

# ==============================================================================
# 9) ACF DER ZEITLICH GEMITTELTEN RESIDUEN
# ==============================================================================

acf_obj <- acf(time_stats$mean_resid, plot = FALSE)

acf_df <- data.frame(
  lag = as.numeric(acf_obj$lag),
  acf = as.numeric(acf_obj$acf)
)

ci <- 1.96 / sqrt(nrow(time_stats))

p <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-ci, ci), color = "red", linetype = "dashed") +
  labs(
    title = "ACF der zeitlich gemittelten Residuen",
    x = "Lag",
    y = "ACF"
  ) +
  theme_simple

print(p)

# ==============================================================================
# 10) RÄUMLICHE FELDER: BIAS, RMSE, SD, ANZAHL
# ==============================================================================

cat("Berechne räumliche Residualfelder ...\n")

bias_field <- array(NA, dim = c(n_lon, n_lat))
rmse_field <- array(NA, dim = c(n_lon, n_lat))
sd_field   <- array(NA, dim = c(n_lon, n_lat))
n_field    <- array(NA, dim = c(n_lon, n_lat))

cell_key <- paste0(df_val$lon_index, "_", df_val$lat_index)
cells_split <- split(df_val$resid_response, cell_key)

for (key in names(cells_split)) {
  
  idx <- strsplit(key, "_")[[1]]
  i <- as.integer(idx[1])
  j <- as.integer(idx[2])
  
  vals <- cells_split[[key]]
  
  if (length(vals) > 1) {
    bias_field[i, j] <- mean(vals, na.rm = TRUE)
    rmse_field[i, j] <- sqrt(mean(vals^2, na.rm = TRUE))
    sd_field[i, j]   <- sd(vals, na.rm = TRUE)
    n_field[i, j]    <- length(vals)
  }
}

field_to_df <- function(field) {
  data.frame(
    lon = rep(lon, times = n_lat),
    lat = rep(lat, each = n_lon),
    value = as.vector(field)
  )
}

bias_df <- field_to_df(bias_field)
rmse_df <- field_to_df(rmse_field)
sd_df   <- field_to_df(sd_field)
n_df    <- field_to_df(n_field)

zlim_bias <- max(abs(bias_df$value), na.rm = TRUE)

p <- ggplot(bias_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-zlim_bias, zlim_bias),
    na.value = "transparent"
  ) +
  coord_fixed() +
  labs(
    title = "Räumlicher Bias",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "Bias"
  ) +
  theme_simple

print(p)

p <- ggplot(rmse_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "transparent"
  ) +
  coord_fixed() +
  labs(
    title = "Räumlicher RMSE",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "RMSE"
  ) +
  theme_simple

print(p)

p <- ggplot(sd_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_gradient(
    low = "white",
    high = "darkblue",
    na.value = "transparent"
  ) +
  coord_fixed() +
  labs(
    title = "Räumliche Residual-SD",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "SD"
  ) +
  theme_simple

print(p)

p <- ggplot(n_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    na.value = "transparent"
  ) +
  coord_fixed() +
  labs(
    title = "Anzahl Validierungspunkte pro Gitterzelle",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "n"
  ) +
  theme_simple

print(p)

# ==============================================================================
# 11) RÄUMLICHE RESIDUENKARTEN FÜR EINIGE ZEITPUNKTE
# ==============================================================================

cat("Erstelle räumliche Residuenkarten für ausgewählte Zeitpunkte ...\n")

n_snaps <- min(9, length(val_times))
snap_times <- val_times[round(seq(1, length(val_times), length.out = n_snaps))]

snap_df <- df_val[df_val$time_index %in% snap_times, ]

snap_df$time_label <- paste0("t = ", round(time[snap_df$time_index]))

zlim_snap <- quantile(abs(snap_df$resid_response), 0.99, na.rm = TRUE)

p <- ggplot(snap_df, aes(x = lon, y = lat, fill = resid_response)) +
  geom_raster() +
  facet_wrap(~ time_label, ncol = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-zlim_snap, zlim_snap),
    na.value = "transparent"
  ) +
  coord_fixed() +
  labs(
    title = "Räumliche Residuen für ausgewählte Zeitpunkte",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "Residuum"
  ) +
  theme_simple

print(p)

# ==============================================================================
# 12) HOVMÖLLER-DIAGRAMME: ZEIT × BREITENGRAD UND ZEIT × LÄNGENGRAD
# ==============================================================================

cat("Berechne einfache Hovmöller-Diagramme ...\n")

lat_hov <- aggregate(
  resid_response ~ time_index + lat_index,
  data = df_val,
  FUN = mean,
  na.rm = TRUE
)

lat_hov$time_val <- time[lat_hov$time_index]
lat_hov$lat_val  <- lat[lat_hov$lat_index]

p <- ggplot(lat_hov, aes(x = time_val, y = lat_val, fill = resid_response)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "transparent"
  ) +
  labs(
    title = "Hovmöller-Diagramm: Zeit × Breitengrad",
    x = "Zeit",
    y = "Breitengrad",
    fill = "Residuum"
  ) +
  theme_simple

print(p)

lon_hov <- aggregate(
  resid_response ~ time_index + lon_index,
  data = df_val,
  FUN = mean,
  na.rm = TRUE
)

lon_hov$time_val <- time[lon_hov$time_index]
lon_hov$lon_val  <- lon[lon_hov$lon_index]

p <- ggplot(lon_hov, aes(x = time_val, y = lon_val, fill = resid_response)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "transparent"
  ) +
  labs(
    title = "Hovmöller-Diagramm: Zeit × Längengrad",
    x = "Zeit",
    y = "Längengrad",
    fill = "Residuum"
  ) +
  theme_simple

print(p)

# ==============================================================================
# PDF SCHLIESSEN
# ==============================================================================

dev.off()

# ==============================================================================
# 13) KONSOLEN-AUSGABE
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("  RESIDUALANALYSE ABGESCHLOSSEN\n")
cat("  -> Ausgabe: residualanalyse_GAM_ggplot2_simple.pdf\n")
cat("====================================================================\n\n")

cat("--- Globale Statistiken: Validierung, Response-Residuen ---\n")
cat(sprintf("  Bias / Mean Residuum:  %+.5f\n", mean(r_val, na.rm = TRUE)))
cat(sprintf("  RMSE:                  %.5f\n", sqrt(mean(r_val^2, na.rm = TRUE))))
cat(sprintf("  MAE:                   %.5f\n", mean(abs(r_val), na.rm = TRUE)))
cat(sprintf("  SD Residuen:           %.5f\n", sd(r_val, na.rm = TRUE)))

skewness <- mean((r_val - mean(r_val, na.rm = TRUE))^3, na.rm = TRUE) /
  sd(r_val, na.rm = TRUE)^3

cat(sprintf("  Schiefe:              %+.4f\n", skewness))

expl_var <- 100 * (1 - mean(r_val^2, na.rm = TRUE) / var(df_val$tf, na.rm = TRUE))

cat(sprintf("  Erklärte Varianz:      %.2f%%\n", expl_var))

acf_lag1 <- acf(time_stats$mean_resid, lag.max = 1, plot = FALSE)$acf[2]

cat("\n--- Temporale Autokorrelation ---\n")
cat(sprintf("  ACF Lag 1 der zeitlich gemittelten Residuen: %.4f\n", acf_lag1))

cat("\n====================================================================\n")