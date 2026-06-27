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

r <- df_val$resid_response
r <- r[is.finite(r)]

mu_hat <- mean(r)
sigma2_hat_ml <- mean((r - mu_hat)^2)
sigma_hat_ml <- sqrt(sigma2_hat_ml)

mu_hat
sigma_hat_ml

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



df_val$pred_bin <- cut(
  df_val$pred,
  breaks = seq(0, 1, by = 0.1),
  include.lowest = TRUE
)


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

hetero_stats

hetero_stats <- data.frame(
  pred_bin = hetero_stats$pred_bin,
  hetero_stats$resid_response
)

hetero_stats$varcoef <- hetero_stats$sd/hetero_stats$mean

print(hetero_stats)


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


covariates <- c("precip", "tcold", "twarm", "co2")

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


acf_obj <- acf(time_stats$mean_resid, plot = FALSE)

acf_df <- data.frame(
  lag = as.numeric(acf_obj$lag),
  acf = as.numeric(acf_obj$acf)
)

ci <- 1.96 / sqrt(nrow(time_stats))

p <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0) +
  # geom_hline(yintercept = c(-ci, ci), color = "red", linetype = "dashed") +
  labs(
    title = "ACF der zeitlich gemittelten Residuen",
    x = "Lag",
    y = "ACF"
  ) +
  theme_simple

print(p)
  

# Gemittelt über Zeit

bias_df <- df_val %>%
  group_by(lon_index, lat_index, lon, lat) %>%
  summarise(
    bias = mean(resid_response, na.rm = TRUE),
    rmse = sqrt(mean(resid_response^2, na.rm = TRUE)),
    sd   = sd(resid_response, na.rm = TRUE),
    n    = sum(!is.na(resid_response)),
    .groups = "drop"
  ) %>% 
  mutate(
    lon_plot = ifelse(lon > 180, lon - 360, lon)
  )


# Auch mit Zeit

bias_time_df <- df_val %>%
  group_by(time, lon_index, lat_index, lon, lat) %>%
  summarise(
    bias = mean(resid_response, na.rm = TRUE),
    rmse = sqrt(mean(resid_response^2, na.rm = TRUE)),
    sd   = sd(resid_response, na.rm = TRUE),
    n    = sum(!is.na(resid_response)),
    .groups = "drop"
  ) %>%
  mutate(
    lon_plot = ifelse(lon > 180, lon - 360, lon)
  )


#Maximum der empirischen Semivarianz als geschätzten Plateau und suchst die 
# erste Distanz, bei der 95% davon erreicht sind.

estimate_range_heuristic <- function(vario_df) {
  
  vario_df <- vario_df[order(vario_df$distance), ]
  
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

