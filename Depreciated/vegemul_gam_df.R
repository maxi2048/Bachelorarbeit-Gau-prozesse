# ============================================================
# Emulator im Dataframe-Format
# ============================================================

# Pakete laden
library(mgcv)
library(ncdf4)
library(ggplot2)
library(maps)


# ============================================================
# 1. Hilfsfunktionen
# ============================================================

# ------------------------------------------------------------
# IPCC-Farbskala von gelb nach grün
# ------------------------------------------------------------

col_magnitude_ipcc_green <- function(n = 21) {
  grDevices::colorRampPalette(
    c(
      rgb(255, 255, 204, maxColorValue = 255),
      rgb(194, 230, 153, maxColorValue = 255),
      rgb(120, 198, 121, maxColorValue = 255),
      rgb(49, 163, 84, maxColorValue = 255),
      rgb(0, 104, 55, maxColorValue = 255)
    ),
    space = "rgb"
  )(n)
}


# ------------------------------------------------------------
# 3D-Array in langes Dataframe umwandeln
#
# Erwartung:
# response hat Dimension:
# time x lon x lat
#
# Jede Zeile im Ergebnis ist eine Kombination aus:
# timestep, lon, lat
# ------------------------------------------------------------

make_emulator_dataframe <- function(response, var_list, lon, lat, time) {
  
  df <- expand.grid(
    timestep = seq_along(time),
    lon_id   = seq_along(lon),
    lat_id   = seq_along(lat)
  )
  
  df$time <- time[df$timestep]
  df$lon  <- lon[df$lon_id]
  df$lat  <- lat[df$lat_id]
  
  # Zielvariable hinzufügen
  df$tf <- as.vector(response)
  
  # Kovariablen hinzufügen
  for (v in names(var_list)) {
    df[[v]] <- as.vector(var_list[[v]])
  }
  
  return(df)
}


# ------------------------------------------------------------
# Trainingsdaten erzeugen
# ------------------------------------------------------------

make_training_data <- function(df, response_name, predictor_names, cal_timesteps) {
  
  train_data <- df[df$timestep %in% cal_timesteps, ]
  
  needed_vars <- c(response_name, predictor_names)
  
  train_data <- train_data[
    complete.cases(train_data[, needed_vars]),
  ]
  
  return(train_data)
}


# ------------------------------------------------------------
# Validierungsdaten erzeugen
# ------------------------------------------------------------

make_validation_data <- function(df, response_name, predictor_names, val_timesteps) {
  
  val_data <- df[df$timestep %in% val_timesteps, ]
  
  needed_vars <- c(response_name, predictor_names)
  
  val_data <- val_data[
    complete.cases(val_data[, needed_vars]),
  ]
  
  return(val_data)
}


# ------------------------------------------------------------
# GAM-Emulator trainieren
# ------------------------------------------------------------

train_emulator_df <- function(train_data, response_name, predictor_names, family) {
  
  form <- as.formula(
    paste0(
      response_name,
      " ~ te(",
      paste(predictor_names, collapse = ", "),
      ")"
    )
  )
  
  emulator <- mgcv::gam(
    formula = form,
    data = train_data,
    family = family
  )
  
  return(emulator)
}


# ------------------------------------------------------------
# Emulator auswerten
#
# Hier werden auch die Residuen berechnet.
# Residuum = beobachteter Wert - vorhergesagter Wert
# ------------------------------------------------------------

evaluate_emulator_df <- function(emulator, val_data, response_name) {
  
  val_data$prediction <- predict(
    emulator,
    newdata = val_data,
    type = "response"
  )
  
  val_data$residual <- val_data[[response_name]] - val_data$prediction
  
  rmse <- sqrt(mean(val_data$residual^2, na.rm = TRUE))
  
  mae <- mean(abs(val_data$residual), na.rm = TRUE)
  
  expl_var <- 1 - mean(val_data$residual^2, na.rm = TRUE) /
    var(val_data[[response_name]], na.rm = TRUE)
  
  results <- data.frame(
    rmse = rmse,
    mae = mae,
    Expl_var = expl_var
  )
  
  return(
    list(
      results = results,
      validation_data = val_data
    )
  )
}


# ------------------------------------------------------------
# RMSE pro Gitterzelle berechnen
# ------------------------------------------------------------

calculate_rmse_field_df <- function(validation_data) {
  
  rmse_field_df <- aggregate(
    residual ~ lon + lat,
    data = validation_data,
    FUN = function(x) sqrt(mean(x^2, na.rm = TRUE))
  )
  
  names(rmse_field_df)[names(rmse_field_df) == "residual"] <- "rmse"
  
  return(rmse_field_df)
}


# ------------------------------------------------------------
# Mittlere Residuen pro Gitterzelle berechnen
# ------------------------------------------------------------

calculate_residual_field_df <- function(validation_data) {
  
  residual_field_df <- aggregate(
    residual ~ lon + lat,
    data = validation_data,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  
  names(residual_field_df)[names(residual_field_df) == "residual"] <- "mean_residual"
  
  return(residual_field_df)
}


# ------------------------------------------------------------
# Plot: Tree Cover für einen bestimmten Zeitschritt
# ------------------------------------------------------------

plot_treecover <- function(df, timestep_to_plot = 1) {
  
  plot_data <- df[df$timestep == timestep_to_plot, ]
  
  world_map <- map_data("world")
  
  ggplot() +
    geom_raster(
      data = plot_data,
      aes(x = lon, y = lat, fill = tf)
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    coord_quickmap() +
    scale_fill_gradientn(
      colours = col_magnitude_ipcc_green()
    ) +
    labs(
      title = paste("Tree Cover, timestep =", timestep_to_plot),
      x = "Longitude",
      y = "Latitude",
      fill = "Tree cover"
    ) +
    theme_minimal()
}


# ------------------------------------------------------------
# Plot: Vorhersage für einen bestimmten Zeitschritt
# ------------------------------------------------------------

plot_prediction <- function(validation_data, timestep_to_plot) {
  
  plot_data <- validation_data[
    validation_data$timestep == timestep_to_plot,
  ]
  
  world_map <- map_data("world")
  
  ggplot() +
    geom_raster(
      data = plot_data,
      aes(x = lon, y = lat, fill = prediction)
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    coord_quickmap() +
    scale_fill_gradientn(
      colours = col_magnitude_ipcc_green()
    ) +
    labs(
      title = paste("Vorhersage, timestep =", timestep_to_plot),
      x = "Longitude",
      y = "Latitude",
      fill = "Prediction"
    ) +
    theme_minimal()
}


# ------------------------------------------------------------
# Plot: Residuen für einen bestimmten Zeitschritt
# ------------------------------------------------------------

plot_residuals <- function(validation_data, timestep_to_plot) {
  
  plot_data <- validation_data[
    validation_data$timestep == timestep_to_plot,
  ]
  
  world_map <- map_data("world")
  
  ggplot() +
    geom_raster(
      data = plot_data,
      aes(x = lon, y = lat, fill = residual)
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    coord_quickmap() +
    labs(
      title = paste("Residuen, timestep =", timestep_to_plot),
      subtitle = "Residuum = beobachtet - vorhergesagt",
      x = "Longitude",
      y = "Latitude",
      fill = "Residuum"
    ) +
    theme_minimal()
}


# ------------------------------------------------------------
# Plot: RMSE pro Gitterzelle
# ------------------------------------------------------------

plot_rmse_field <- function(rmse_field_df) {
  
  world_map <- map_data("world")
  
  ggplot() +
    geom_raster(
      data = rmse_field_df,
      aes(x = lon, y = lat, fill = rmse)
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    coord_quickmap() +
    labs(
      title = "RMSE pro Gitterzelle",
      x = "Longitude",
      y = "Latitude",
      fill = "RMSE"
    ) +
    theme_minimal()
}


# ------------------------------------------------------------
# Plot: mittlere Residuen pro Gitterzelle
# ------------------------------------------------------------

plot_mean_residual_field <- function(residual_field_df) {
  
  world_map <- map_data("world")
  
  ggplot() +
    geom_raster(
      data = residual_field_df,
      aes(x = lon, y = lat, fill = mean_residual)
    ) +
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      fill = NA,
      color = "black",
      linewidth = 0.2
    ) +
    coord_quickmap() +
    labs(
      title = "Mittlere Residuen pro Gitterzelle",
      subtitle = "positive Werte: Modell unterschätzt; negative Werte: Modell überschätzt",
      x = "Longitude",
      y = "Latitude",
      fill = "mittleres Residuum"
    ) +
    theme_minimal()
}


# ============================================================
# 2. Daten laden
# ============================================================

# ------------------------------------------------------------
# Tree Cover laden
# ------------------------------------------------------------

nc <- nc_open("MPI_V3_G_TreeCover_mean_time.nc")

lon <- ncvar_get(nc, "lon")
lat <- rev(ncvar_get(nc, "lat"))
time <- ncvar_get(nc, "time")[1:40]

tf <- aperm(
  ncvar_get(nc, "cover_fract")[, 48:1, 1:40],
  c(3, 1, 2)
)

nc_close(nc)


# ------------------------------------------------------------
# Niederschlag laden
# ------------------------------------------------------------

nc <- nc_open("MPI_V3_G_Pann_mean_time.nc")

time_tmp <- ncvar_get(nc, "time")

precip <- aperm(
  ncvar_get(nc, "precip")[, 48:1, which(time_tmp %in% time)],
  c(3, 1, 2)
)

nc_close(nc)


# ------------------------------------------------------------
# Kalte Temperatur laden
# ------------------------------------------------------------

nc <- nc_open("MPI_V3_G_tcold_mean_time.nc")

time_tmp <- ncvar_get(nc, "time")

tcold <- aperm(
  ncvar_get(nc, "temp2")[, 48:1, which(time_tmp %in% time)],
  c(3, 1, 2)
)

nc_close(nc)


# ------------------------------------------------------------
# Warme Temperatur laden
# ------------------------------------------------------------

nc <- nc_open("MPI_V3_G_twarm_mean_time.nc")

time_tmp <- ncvar_get(nc, "time")

twarm <- aperm(
  ncvar_get(nc, "temp2")[, 48:1, which(time_tmp %in% time)],
  c(3, 1, 2)
)

nc_close(nc)


# ------------------------------------------------------------
# CO2-Daten laden
# ------------------------------------------------------------

co2table <- read.table(
  "CO2_stack_156K_spline_V2.tab",
  header = TRUE,
  sep = "\t",
  skip = 13
)

co2_data <- data.frame(
  age = -co2table$Age..ka.BP. * 1000,
  co2 = co2table$CO2..µmol.mol.
)

co2_timeseries <- sapply(
  seq(-20000, -500, by = 500),
  function(t) {
    mean(
      co2_data$co2[
        co2_data$age >= t &
          co2_data$age < t + 500
      ],
      na.rm = TRUE
    )
  }
)

co2 <- aperm(
  array(
    rep(co2_timeseries, each = length(lon) * length(lat)),
    dim = dim(tf)[c(2, 3, 1)]
  ),
  c(3, 1, 2)
)


# ============================================================
# 3. Variablen zusammenbauen
# ============================================================

response <- tf

var_list <- list(
  precip = precip,
  tcold  = tcold,
  twarm  = twarm,
  co2    = co2
)


# ============================================================
# 4. In Dataframe umwandeln
# ============================================================

df <- make_emulator_dataframe(
  response = response,
  var_list = var_list,
  lon = lon,
  lat = lat,
  time = time
)

# Kontrolle
print(head(df))
print(dim(df))
print(summary(df))


# ============================================================
# 5. Trainings- und Validierungszeitschritte festlegen
# ============================================================

training_data_fraction <- 0.5

cal_timesteps <- seq(
  1,
  length(time),
  by = 1 / training_data_fraction
)

if (training_data_fraction == 1) {
  val_timesteps <- cal_timesteps
} else {
  val_timesteps <- setdiff(seq_along(time), cal_timesteps)
}

print(cal_timesteps)
print(val_timesteps)


# ============================================================
# 6. Trainings- und Validierungsdaten erzeugen
# ============================================================

response_name <- "tf"

predictor_names <- c(
  "precip",
  "tcold",
  "twarm",
  "co2"
)

train_data <- make_training_data(
  df = df,
  response_name = response_name,
  predictor_names = predictor_names,
  cal_timesteps = cal_timesteps
)

val_data <- make_validation_data(
  df = df,
  response_name = response_name,
  predictor_names = predictor_names,
  val_timesteps = val_timesteps
)

print(dim(train_data))
print(dim(val_data))


# ============================================================
# 7. Emulator trainieren
# ============================================================

# Hinweis:
# Falls tf Anteile zwischen 0 und 1 enthält, ist quasibinomial()
# oft robuster als binomial().
#
# Falls du wirklich binomial() willst, ersetze quasibinomial()
# einfach durch binomial().

emulator <- train_emulator_df(
  train_data = train_data,
  response_name = response_name,
  predictor_names = predictor_names,
  family = quasibinomial()
)

summary(emulator)


# ============================================================
# 8. Emulator evaluieren
# ============================================================

eval_out <- evaluate_emulator_df(
  emulator = emulator,
  val_data = val_data,
  response_name = response_name
)

results <- eval_out$results
validation_data <- eval_out$validation_data

print(results)
print(head(validation_data))


# ============================================================
# 9. Residuen ansehen
# ============================================================

# Residuum = beobachtet - vorhergesagt
# positiv: Modell hat zu niedrig vorhergesagt
# negativ: Modell hat zu hoch vorhergesagt

summary(validation_data$residual)

hist(
  validation_data$residual,
  breaks = 50,
  main = "Histogramm der Residuen",
  xlab = "Residuum = beobachtet - vorhergesagt"
)


# ============================================================
# 10. RMSE pro Gitterzelle berechnen
# ============================================================

rmse_field_df <- calculate_rmse_field_df(
  validation_data = validation_data
)

print(head(rmse_field_df))
print(summary(rmse_field_df$rmse))


# ============================================================
# 11. Mittlere Residuen pro Gitterzelle berechnen
# ============================================================

residual_field_df <- calculate_residual_field_df(
  validation_data = validation_data
)

print(head(residual_field_df))
print(summary(residual_field_df$mean_residual))


# ============================================================
# 12. Plots
# ============================================================

# Ursprünglicher Tree Cover für timestep 1
plot_treecover(
  df = df,
  timestep_to_plot = 1
)

# Vorhersage für ersten Validierungszeitpunkt
plot_prediction(
  validation_data = validation_data,
  timestep_to_plot = val_timesteps[1]
)

# Residuen für ersten Validierungszeitpunkt
plot_residuals(
  validation_data = validation_data,
  timestep_to_plot = val_timesteps[1]
)

# RMSE pro Gitterzelle
plot_rmse_field(
  rmse_field_df = rmse_field_df
)

# Mittlere Residuen pro Gitterzelle
plot_mean_residual_field(
  residual_field_df = residual_field_df
)


# ============================================================
# 13. Ergebnisse optional abspeichern
# ============================================================

write.csv(
  df,
  file = "emulator_full_dataframe.csv",
  row.names = FALSE
)

write.csv(
  train_data,
  file = "emulator_training_data.csv",
  row.names = FALSE
)

write.csv(
  validation_data,
  file = "emulator_validation_predictions_residuals.csv",
  row.names = FALSE
)

write.csv(
  rmse_field_df,
  file = "emulator_rmse_field.csv",
  row.names = FALSE
)

write.csv(
  residual_field_df,
  file = "emulator_mean_residual_field.csv",
  row.names = FALSE
)


# ============================================================
# 14. Kurze Kontrolle am Ende
# ============================================================

cat("\nFertig.\n")
cat("Gesamtes Dataframe:\n")
print(dim(df))

cat("\nTrainingsdaten:\n")
print(dim(train_data))

cat("\nValidierungsdaten mit prediction und residual:\n")
print(dim(validation_data))

cat("\nGesamtgüte:\n")
print(results)

cat("\nSpalten in validation_data:\n")
print(names(validation_data))