library("geosphere")

required_objects <- c(
  "df_val", "df_train", "emulator", "time", "lon", "lat",
  "time_stats", "bias_df", "rmse_df", "sd_df"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
print(missing_objects)

required_df_val_cols <- c(
  "tf", "pred", "resid_response", "time_index",
  "lon_index", "lat_index", "lon", "lat"
)

missing_cols <- required_df_val_cols[!required_df_val_cols %in% names(df_val)]
print(missing_cols)


# Durchschnittlicher Quadratischer Abstand von A und B mit Abstand h 
# Coord Cols sind sowas wie Lon Lat, Lon x Lat, Lon x Lat x Time (Vektor)
# Value Col Abstand über Tree Ffraction, resid etc.
# Unter Umständen hohe Laufzeit, daher n_pairs bestimmt wieviel gesamplet wird 
# n_bins bedeutet Wieviele "Distanzgruppen" es gibt 
# max_quant sollen alle Distanzen berücksichtigt werden oder nur bis zu einem bestimmten Quantil
# standardize_coords, sollen die "Koordinaten" standardisiert werden oder dürfen Äpfel mit Birnen verglichen werdne

empirical_variogram <- function(data,
                                coord_cols,
                                value_col,
                                n_pairs = 200000,
                                n_bins = 18,
                                max_dist_quantile = 0.95,
                                standardize_coords = FALSE,
                                dist_method = "euclidean") {
  
  dat <- data[, c(coord_cols, value_col), drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  if (nrow(dat) < 10) {
    stop("Zu wenige vollständige Datenpunkte für ein Variogramm.")
  }
  
  coords <- as.matrix(dat[, coord_cols, drop = FALSE])
  values <- dat[[value_col]]
  
  if (standardize_coords && dist_method == "haversine") {
    stop("standardize_coords = TRUE ist mit dist_method = 'haversine' nicht sinnvoll.")
  }
  
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
  
  # ----------------------------------------------------------
  # Abstand berechnen
  # ----------------------------------------------------------
  
  if (dist_method == "euclidean") {
    
    d <- sqrt(
      rowSums(
        (coords[i, , drop = FALSE] - coords[j, , drop = FALSE])^2
      )
    )
    
    distance_name <- "distance"
    
  } else if (dist_method == "haversine") {
    
    if (length(coord_cols) != 2) {
      stop("Für dist_method = 'haversine' müssen genau zwei Koordinatenspalten angegeben werden: Longitude und Latitude.")
    }
    
    # Wichtig:
    # coord_cols muss in der Reihenfolge c('lon', 'lat') oder c('lon_plot', 'lat') sein.
    # geosphere::distHaversine erwartet: Longitude, Latitude.
    
    d <- geosphere::distHaversine(
      p1 = coords[i, , drop = FALSE],
      p2 = coords[j, , drop = FALSE]
    ) / 1000
    
    distance_name <- "distance_km"
    
  } else {
    
    stop("Unbekannte dist_method. Erlaubt sind: 'euclidean' oder 'haversine'.")
    
  }
  
  # ----------------------------------------------------------
  # Semivarianz berechnen
  # ----------------------------------------------------------
  
  gamma <- 0.5 * (values[i] - values[j])^2
  
  max_dist <- as.numeric(
    quantile(d, probs = max_dist_quantile, na.rm = TRUE)
  )
  
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
  
  names(out)[names(out) == "distance"] <- distance_name
  
  return(out)
}


spatial_bias_dat <- bias_df

df_val$lon_plot <- ifelse(df_val$lon > 180, df_val$lon - 360, df_val$lon)

spatial_bias_dat_t <- bias_df %>%
  transmute(
    lon = lon_plot,
    lat = lat,
    mean_resid = bias
  ) %>%
  filter(complete.cases(.))

# Gemittelt über Zeit

bias_df <- df_val %>%
  group_by(lon_index, lat_index, lon, lat) %>%
  summarise(
    bias = mean(resid_response, na.rm = TRUE),
    rmse = sqrt(mean(resid_response^2, na.rm = TRUE)),
    sd   = sd(resid_response, na.rm = TRUE),
    n    = sum(!is.na(resid_response)),
    .groups = "drop"
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

spatial_raw_dat <- df_val[raw_idx, , drop = FALSE] %>%
  mutate(
    lon_plot = ifelse(lon > 180, lon - 360, lon)
  ) %>%
  transmute(
    lon = lon_plot,
    lat = lat,
    resid_response = resid_response
  ) %>%
  filter(complete.cases(.))


# Wie verändert sich der Bias über eine größere Sphärische 
spatial_bias_variogram <- empirical_variogram(
  data = spatial_bias_dat_t,
  coord_cols = c("lon", "lat"),
  value_col = "mean_resid",
  n_pairs = 150000,
  n_bins = 18,
  standardize_coords = FALSE,
  dist_method = "haversine"
)


# Keine Kugel weil nur Longitudinal
long_variogram <- empirical_variogram(
  data = spatial_raw_dat,
  coord_cols = c("lon"),
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = FALSE,
  dist_method = "euclidean"
)


lat_variogram <- empirical_variogram(
  data = spatial_raw_dat,
  coord_cols = c("lat"),
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = FALSE,
  dist_method = "euclidean"
)

temporal_variogram <- empirical_variogram(
  data = temporal_dat,
  coord_cols = c("time_val"),
  value_col = "mean_resid",
  n_pairs = 100000,
  n_bins = 14,
  standardize_coords = FALSE,
  dist_method = "euclidean"
)

# Vielleicht auch relevant, aber Interpretation schwierig
climate_variogram <- empirical_variogram(
  data = climate_dat,
  coord_cols = covariates,
  value_col = "resid_response",
  n_pairs = 200000,
  n_bins = 18,
  standardize_coords = TRUE,
  dist_method = "euclidean"
)





p_spatial_bias <- ggplot(
  spatial_bias_variogram,
  aes(x = distance_km, y = semivariance)
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


p_spatial_raw <- ggplot(spatial_raw_variogram,
                        aes(x = distance, y = semivariance)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Räumliches Variogramm der Rohresiduen",
    subtitle = "Stichprobe aus df_val, sphärische Distanz",
    x = "Distanz in km",
    y = "Semivarianz"
  ) +
  theme_minimal()

print(p_spatial_raw)



p_long <- ggplot(long_variogram,
                 aes(x = distance, y = semivariance)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Longitude-Variogramm der Rohresiduen",
    x = "Abstand in Longitude-Graden",
    y = "Semivarianz"
  ) +
  theme_minimal()

print(p_long)


p_lat <- ggplot(lat_variogram,
                aes(x = distance, y = semivariance)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Latitude-Variogramm der Rohresiduen",
    x = "Abstand in Latitude-Graden",
    y = "Semivarianz"
  ) +
  theme_minimal()

print(p_lat)


p_temporal <- ggplot(temporal_variogram,
                     aes(x = distance, y = semivariance)) +
  geom_point() +
  geom_line() +
  labs(
    title = "Zeitliches Variogramm der mittleren Residuen",
    x = "Zeitlicher Abstand",
    y = "Semivarianz"
  ) +
  theme_minimal()

print(p_temporal)



  
  p_climate <- ggplot(climate_variogram,
                      aes(x = distance, y = semivariance)) +
    geom_point() +
    geom_line() +
    labs(
      title = "Klimaraum-Variogramm der Rohresiduen",
      subtitle = paste("Kovariaten:", paste(covariates, collapse = ", ")),
      x = "standardisierte Klimaraum-Distanz",
      y = "Semivarianz"
    ) +
    theme_minimal()
  
  print(p_climate)


