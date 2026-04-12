library(mgcv)
library(ncdf4)
library(ggplot2)
library(reshape2)

cat("=====================================================\n")
cat(" Start des Skripts: GAM-Emulator mit ausführlicher Ausgabe\n")
cat("=====================================================\n\n")
cat("Ziel dieses Skripts:\n")
cat("1. Klima- und Vegetationsdaten laden\n")
cat("2. Datenstruktur sichtbar machen\n")
cat("3. Ein GAM trainieren\n")
cat("4. Vorhersagen evaluieren\n")
cat("5. Relevante Plots erzeugen\n")
cat("6. Residuen untersuchen\n\n")

# ----------------------------------------------------
# Hilfsfunktion: kurze Zusammenfassung eines Arrays
# ----------------------------------------------------
summarize_array <- function(x, name) {
  cat("---------------------------------------------\n")
  cat("Zusammenfassung für:", name, "\n")
  cat("Dimension:", paste(dim(x), collapse = " x "), "\n")
  cat("Anzahl NA:", sum(is.na(x)), "\n")
  cat("Minimum:", suppressWarnings(min(x, na.rm = TRUE)), "\n")
  cat("Maximum:", suppressWarnings(max(x, na.rm = TRUE)), "\n")
  cat("Mittelwert:", suppressWarnings(mean(x, na.rm = TRUE)), "\n")
  cat("Standardabweichung:", suppressWarnings(sd(as.vector(x), na.rm = TRUE)), "\n")
  cat("---------------------------------------------\n\n")
}

# ----------------------------------------------------
# Hilfsfunktion: Karte eines Felds plotten
# ----------------------------------------------------
plot_field <- function(field2d, lon, lat, title, fill_label, file_name = NULL) {
  df <- melt(field2d)
  names(df) <- c("lon_idx", "lat_idx", "value")
  df$lon <- lon[df$lon_idx]
  df$lat <- lat[df$lat_idx]
  
  p <- ggplot(df, aes(x = lon, y = lat, fill = value)) +
    geom_tile() +
    coord_fixed() +
    labs(title = title, x = "Longitude", y = "Latitude", fill = fill_label) +
    theme_minimal()
  
  print(p)
  
  if (!is.null(file_name)) {
    ggsave(filename = file_name, plot = p, width = 8, height = 4.5)
    cat("Plot gespeichert unter:", file_name, "\n")
  }
}

# ----------------------------------------------------
# Hilfsfunktion: Zeitreihe eines Punkts plotten
# ----------------------------------------------------
plot_timeseries_point <- function(values, time, title, ylab, file_name = NULL) {
  df <- data.frame(time = time, value = values)
  
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line() +
    geom_point(size = 1) +
    labs(title = title, x = "Zeit", y = ylab) +
    theme_minimal()
  
  print(p)
  
  if (!is.null(file_name)) {
    ggsave(filename = file_name, plot = p, width = 8, height = 4.5)
    cat("Plot gespeichert unter:", file_name, "\n")
  }
}

# ----------------------------------------------------
# Trainingsfunktion
# ----------------------------------------------------
train_emulator <- function(response, var_list, family, cal_timesteps, do_plots = TRUE) {
  cat("=====================================================\n")
  cat(" train_emulator() startet\n")
  cat("=====================================================\n\n")
  
  cat("Was passiert jetzt?\n")
  cat("- Die 3D-Arrays werden in einen langen Datensatz umgebaut.\n")
  cat("- Jede Zeile entspricht einer gültigen Beobachtung an einem Ort und Zeitpunkt.\n")
  cat("- Danach wird ein GAM geschätzt.\n\n")
  
  cat("Verwendete Trainings-Zeitschritte:\n")
  print(cal_timesteps)
  cat("\n")
  
  cat("Kovariaten im Modell:\n")
  print(names(var_list))
  cat("\n")
  
  valid_idx <- which(!is.na(response[cal_timesteps,,]))
  cat("Anzahl gültiger Response-Beobachtungen im Training:", length(valid_idx), "\n\n")
  
  covar <- do.call(cbind.data.frame, lapply(var_list, function(x) {
    x[cal_timesteps,,][valid_idx]
  }))
  
  names(covar) <- names(var_list)
  
  formula <- formula(
    paste0("tf ~ te(", paste(names(covar), collapse = ","), ")")
  )
  
  cat("Automatisch erzeugte Modellformel:\n")
  print(formula)
  cat("\n")
  
  data <- cbind.data.frame(
    tf = response[cal_timesteps,,][valid_idx],
    covar
  )
  
  cat("Dimension des Trainingsdatensatzes:\n")
  print(dim(data))
  cat("\nErste 6 Zeilen des Datensatzes:\n")
  print(head(data))
  cat("\nZusammenfassung des Datensatzes:\n")
  print(summary(data))
  cat("\n")
  
  if (do_plots) {
    cat("Erzeuge erste diagnostische Trainings-Plots ...\n")
    
    p1 <- ggplot(data, aes(x = tf)) +
      geom_histogram(bins = 50) +
      labs(title = "Verteilung der Zielvariable tf", x = "tf", y = "Häufigkeit") +
      theme_minimal()
    print(p1)
    
    if ("precip" %in% names(data)) {
      p2 <- ggplot(data, aes(x = precip, y = tf)) +
        geom_point(alpha = 0.2) +
        labs(title = "tf gegen precip", x = "precip", y = "tf") +
        theme_minimal()
      print(p2)
    }
    
    if ("tcold" %in% names(data)) {
      p3 <- ggplot(data, aes(x = tcold, y = tf)) +
        geom_point(alpha = 0.2) +
        labs(title = "tf gegen tcold", x = "tcold", y = "tf") +
        theme_minimal()
      print(p3)
    }
    
    if ("twarm" %in% names(data)) {
      p4 <- ggplot(data, aes(x = twarm, y = tf)) +
        geom_point(alpha = 0.2) +
        labs(title = "tf gegen twarm", x = "twarm", y = "tf") +
        theme_minimal()
      print(p4)
    }
  }
  
  cat("Jetzt wird das GAM mit mgcv::gam() geschätzt ...\n\n")
  emulator <- mgcv::gam(formula, data = data, family = family)
  
  cat("Modell erfolgreich geschätzt.\n\n")
  cat("Modellzusammenfassung:\n")
  print(summary(emulator))
  cat("\n")
  
  cat("AIC des Modells:", AIC(emulator), "\n\n")
  
  if (do_plots) {
    cat("Erzeuge GAM-Diagnostik aus mgcv ...\n")
    plot(emulator, pages = 1)
    gam.check(emulator)
  }
  
  cat("train_emulator() beendet.\n\n")
  emulator
}

# ----------------------------------------------------
# Evaluationsfunktion
# ----------------------------------------------------
evaluate_emulator <- function(emulator, response, var_list, val_timesteps,
                              lon = NULL, lat = NULL, time = NULL, do_plots = TRUE) {
  cat("=====================================================\n")
  cat(" evaluate_emulator() startet\n")
  cat("=====================================================\n\n")
  
  cat("Was passiert jetzt?\n")
  cat("- Es wird ein Validierungsdatensatz aufgebaut.\n")
  cat("- Das Modell macht Vorhersagen.\n")
  cat("- Fehlermaße und Residuen werden berechnet.\n")
  cat("- Zusätzlich werden räumliche Fehlerkarten erzeugt.\n\n")
  
  cat("Verwendete Validierungs-Zeitschritte:\n")
  print(val_timesteps)
  cat("\n")
  
  valid_idx <- which(!is.na(response[val_timesteps,,]))
  cat("Anzahl gültiger Response-Beobachtungen in der Validierung:", length(valid_idx), "\n\n")
  
  covar <- do.call(cbind.data.frame, lapply(var_list, function(x) {
    x[val_timesteps,,][valid_idx]
  }))
  names(covar) <- names(var_list)
  
  data <- cbind.data.frame(
    tf = response[val_timesteps,,][valid_idx],
    covar
  )
  
  cat("Dimension des Validierungsdatensatzes:\n")
  print(dim(data))
  cat("\nErste 6 Zeilen des Validierungsdatensatzes:\n")
  print(head(data))
  cat("\n")
  
  cat("Berechne Vorhersagen auf der Antwortskala ...\n\n")
  predict_data <- predict(emulator, type = "response", newdata = data)
  residuals_data <- data$tf - predict_data
  
  cat("Erste 6 Vorhersagen:\n")
  print(head(predict_data))
  cat("\nErste 6 Residuen:\n")
  print(head(residuals_data))
  cat("\n")
  
  datafield_predict <- response[val_timesteps,,]
  datafield_predict[which(!is.na(datafield_predict))] <- predict_data
  
  residual_field <- response[val_timesteps,,]
  residual_field[which(!is.na(residual_field))] <- residuals_data
  
  rmse_field <- array(NA, dim = dim(datafield_predict)[2:3])
  mae_field <- array(NA, dim = dim(datafield_predict)[2:3])
  
  for (i in 1:dim(rmse_field)[1]) {
    for (j in 1:dim(rmse_field)[2]) {
      obs_ij <- response[val_timesteps, i, j]
      pred_ij <- datafield_predict[, i, j]
      
      rmse_field[i, j] <- suppressWarnings(
        sqrt(mean((obs_ij - pred_ij)^2, na.rm = TRUE))
      )
      
      mae_field[i, j] <- suppressWarnings(
        mean(abs(obs_ij - pred_ij), na.rm = TRUE)
      )
    }
  }
  
  results <- data.frame(
    rmse = sqrt(mean((predict_data - data$tf)^2)),
    mae = mean(abs(predict_data - data$tf)),
    Expl_var = 1 - mean((predict_data - data$tf)^2) / var(data$tf)
  )
  
  cat("Globale Gütemaße:\n")
  print(results)
  cat("\nZusammenfassung der Residuen:\n")
  print(summary(residuals_data))
  cat("\n")
  
  if (do_plots) {
    cat("Erzeuge Evaluations-Plots ...\n")
    
    df_pred <- data.frame(observed = data$tf, predicted = predict_data)
    p_obs_pred <- ggplot(df_pred, aes(x = observed, y = predicted)) +
      geom_point(alpha = 0.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      labs(title = "Beobachtet vs. vorhergesagt", x = "Beobachtet", y = "Vorhergesagt") +
      theme_minimal()
    print(p_obs_pred)
    
    p_res_hist <- ggplot(data.frame(residual = residuals_data), aes(x = residual)) +
      geom_histogram(bins = 50) +
      labs(title = "Histogramm der Residuen", x = "Residual", y = "Häufigkeit") +
      theme_minimal()
    print(p_res_hist)
    
    p_res_pred <- ggplot(data.frame(predicted = predict_data, residual = residuals_data),
                         aes(x = predicted, y = residual)) +
      geom_point(alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Residuen gegen Vorhersage", x = "Vorhergesagt", y = "Residual") +
      theme_minimal()
    print(p_res_pred)
    
    if (!is.null(lon) && !is.null(lat)) {
      plot_field(datafield_predict[1,,], lon, lat,
                 title = "Vorhersagefeld für den ersten Validierungszeitpunkt",
                 fill_label = "Prediction")
      
      plot_field(residual_field[1,,], lon, lat,
                 title = "Residualfeld für den ersten Validierungszeitpunkt",
                 fill_label = "Residual")
      
      plot_field(rmse_field, lon, lat,
                 title = "RMSE-Karte über alle Validierungszeitpunkte",
                 fill_label = "RMSE")
    }
    
    if (!is.null(time)) {
      i_mid <- round(dim(response)[2] / 2)
      j_mid <- round(dim(response)[3] / 2)
      
      plot_timeseries_point(
        values = response[val_timesteps, i_mid, j_mid],
        time = time[val_timesteps],
        title = "Beobachtete tf-Zeitreihe an einem mittleren Gitterpunkt",
        ylab = "tf"
      )
      
      plot_timeseries_point(
        values = datafield_predict[, i_mid, j_mid],
        time = time[val_timesteps],
        title = "Vorhergesagte tf-Zeitreihe an einem mittleren Gitterpunkt",
        ylab = "Vorhersage"
      )
      
      plot_timeseries_point(
        values = residual_field[, i_mid, j_mid],
        time = time[val_timesteps],
        title = "Residual-Zeitreihe an einem mittleren Gitterpunkt",
        ylab = "Residual"
      )
    }
  }
  
  cat("evaluate_emulator() beendet.\n\n")
  
  list(
    results = results,
    rmse_field = rmse_field,
    mae_field = mae_field,
    residual_field = residual_field,
    predictions = datafield_predict,
    pointwise_residuals = residuals_data
  )
}

# ----------------------------------------------------
# Daten laden
# ----------------------------------------------------
cat("=====================================================\n")
cat(" Daten werden geladen\n")
cat("=====================================================\n\n")

nc <- nc_open("MPI_V3_G_TreeCover_mean_time.nc")
lon <- ncvar_get(nc, "lon")
lat <- rev(ncvar_get(nc, "lat"))
time <- ncvar_get(nc, "time")[1:40]
tf <- aperm(ncvar_get(nc, "cover_fract")[,48:1,1:40], c(3,1,2))
nc_close(nc)

cat("Tree Cover geladen.\n")
cat("Länge lon:", length(lon), "\n")
cat("Länge lat:", length(lat), "\n")
cat("Länge time:", length(time), "\n\n")
summarize_array(tf, "tf (Tree Cover)")

nc <- nc_open("MPI_V3_G_Pann_mean_time.nc")
time_tmp <- ncvar_get(nc, "time")
precip <- aperm(ncvar_get(nc, "precip")[,48:1,which(time_tmp %in% time)], c(3,1,2))
nc_close(nc)
cat("Niederschlag geladen.\n\n")
summarize_array(precip, "precip")

nc <- nc_open("MPI_V3_G_tcold_mean_time.nc")
time_tmp <- ncvar_get(nc, "time")
tcold <- aperm(ncvar_get(nc, "temp2")[,48:1,which(time_tmp %in% time)], c(3,1,2))
nc_close(nc)
cat("tcold geladen.\n\n")
summarize_array(tcold, "tcold")

nc <- nc_open("MPI_V3_G_twarm_mean_time.nc")
time_tmp <- ncvar_get(nc, "time")
twarm <- aperm(ncvar_get(nc, "temp2")[,48:1,which(time_tmp %in% time)], c(3,1,2))
nc_close(nc)
cat("twarm geladen.\n\n")
summarize_array(twarm, "twarm")

co2table <- read.table("CO2_stack_156K_spline_V2.tab", header = TRUE, sep = "\t", skip = 13)
co2_timeseries <- data.frame(age = -co2table$Age..ka.BP.*1000,
                             co2 = co2table$CO2..µmol.mol.)
co2_timeseries <- sapply(seq(-20000, -500, by = 500), function(t) {
  mean(co2_timeseries$co2[which(co2_timeseries$age >= t & co2_timeseries$age < t + 500)])
})
co2 <- aperm(array(rep(co2_timeseries, each = length(lon) * length(lat)),
                   dim = dim(tf)[c(2,3,1)]), c(3,1,2))
cat("CO2-Zeitreihe geladen und auf das Raumgitter erweitert.\n\n")
summarize_array(co2, "co2")

# ----------------------------------------------------
# Erste Plots zur Datenexploration
# ----------------------------------------------------
cat("=====================================================\n")
cat(" Erste explorative Plots\n")
cat("=====================================================\n\n")

plot_field(tf[1,,], lon, lat,
           title = "Tree Cover für den ersten Zeitpunkt",
           fill_label = "tf")

plot_field(precip[1,,], lon, lat,
           title = "Niederschlag für den ersten Zeitpunkt",
           fill_label = "precip")

plot_field(tcold[1,,], lon, lat,
           title = "tcold für den ersten Zeitpunkt",
           fill_label = "tcold")

plot_field(twarm[1,,], lon, lat,
           title = "twarm für den ersten Zeitpunkt",
           fill_label = "twarm")

mid_i <- round(dim(tf)[2] / 2)
mid_j <- round(dim(tf)[3] / 2)
cat("Mittlerer Gitterpunkt für Zeitreihenplots: i =", mid_i, ", j =", mid_j, "\n\n")

plot_timeseries_point(tf[, mid_i, mid_j], time,
                      title = "Tree Cover im Zeitverlauf an einem mittleren Gitterpunkt",
                      ylab = "tf")

plot_timeseries_point(precip[, mid_i, mid_j], time,
                      title = "Niederschlag im Zeitverlauf an einem mittleren Gitterpunkt",
                      ylab = "precip")

# ----------------------------------------------------
# Training / Validierung festlegen
# ----------------------------------------------------
cat("=====================================================\n")
cat(" Trainings- und Validierungsaufteilung\n")
cat("=====================================================\n\n")

training_data_fraction <- 0.5
cal_timesteps <- seq(1, length(time), by = 1 / training_data_fraction)
if (training_data_fraction == 1) {
  val_timesteps <- cal_timesteps
} else {
  val_timesteps <- (1:length(time))[!(1:length(time)) %in% cal_timesteps]
}

cat("training_data_fraction =", training_data_fraction, "\n")
cat("Trainings-Zeitschritte:\n")
print(cal_timesteps)
cat("Validierungs-Zeitschritte:\n")
print(val_timesteps)
cat("\n")

# ----------------------------------------------------
# Modellinput festlegen
# ----------------------------------------------------
# var_list <- list(precip = precip, tcold = tcold, twarm = twarm, co2 = co2)
var_list <- list(precip = precip, tcold = tcold, twarm = twarm)
response <- tf

cat("=====================================================\n")
cat(" Modellinput\n")
cat("=====================================================\n\n")
cat("Zielvariable: tf\n")
cat("Verwendete Kovariaten:\n")
print(names(var_list))
cat("\nHinweis: CO2 ist im Skript vorbereitet, aktuell aber nicht im Modell enthalten.\n\n")

# ----------------------------------------------------
# Modell trainieren und evaluieren
# ----------------------------------------------------
emulator <- train_emulator(
  response = response,
  var_list = var_list,
  family = binomial,
  cal_timesteps = cal_timesteps,
  do_plots = TRUE
)

eval_out <- evaluate_emulator(
  emulator = emulator,
  response = response,
  var_list = var_list,
  val_timesteps = val_timesteps,
  lon = lon,
  lat = lat,
  time = time,
  do_plots = TRUE
)

# ----------------------------------------------------
# Abschließende Ausgabe
# ----------------------------------------------------
cat("=====================================================\n")
cat(" Skriptende: wichtigste Ergebnisse\n")
cat("=====================================================\n\n")
cat("Globale Evaluationsmetriken:\n")
print(eval_out$results)
cat("\n")

cat("Kurze Interpretation:\n")
cat("- RMSE misst die mittlere quadratische Abweichung. Kleine Werte sind besser.\n")
cat("- MAE misst die mittlere absolute Abweichung. Kleine Werte sind besser.\n")
cat("- Expl_var nahe 1 bedeutet: Das Modell erklärt einen großen Teil der Variation.\n\n")

cat("Objekte, die jetzt verfügbar sind:\n")
cat("- emulator               -> das geschätzte GAM\n")
cat("- eval_out$results       -> globale Gütemaße\n")
cat("- eval_out$rmse_field    -> räumliche RMSE-Karte\n")
cat("- eval_out$mae_field     -> räumliche MAE-Karte\n")
cat("- eval_out$residual_field-> Residuen im Raum-Zeit-Format\n")
cat("- eval_out$predictions   -> Vorhersagen im Raum-Zeit-Format\n\n")

cat("Nächster methodischer Schritt für deine Bachelorarbeit:\n")
cat("Untersuche die Residuen auf räumliche und zeitliche Struktur.\n")
cat("Wenn dort noch deutliche Muster zu sehen sind, motiviert das einen Gaußprozess auf den Residuen.\n")


residual_field <- eval_out$residual_field

residuals_vec <- eval_out$pointwise_residuals

ggplot(data.frame(res = residuals_vec), aes(res)) +
  geom_histogram(bins = 60) +
  labs(title = "Histogramm der Residuen") +
  theme_minimal()
