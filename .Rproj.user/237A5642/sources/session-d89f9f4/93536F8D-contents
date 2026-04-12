library(mgcv)
library(ncdf4)

cat("=== Skriptstart ===\n")
cat("Dieses Skript lädt Klima- und Vegetationsdaten, trainiert ein GAM und evaluiert die Vorhersagegüte.\n\n")

train_emulator <- function(response, var_list, family, cal_timesteps){
  cat("--- train_emulator() startet ---\n")
  cat("Es wird jetzt der Trainingsdatensatz aus den 3D-Arrays aufgebaut.\n")
  cat("Die erwartete Datenstruktur ist: (Zeit, Lon, Lat).\n")
  cat("Verwendete Trainings-Zeitschritte:", paste(cal_timesteps, collapse = ", "), "\n")
  cat("Anzahl Kovariaten:", length(var_list), "\n")
  cat("Kovariatennamen:", paste(names(var_list), collapse = ", "), "\n\n")
  
  covar <- do.call(cbind.data.frame, lapply(var_list, function(x){
    cat("Extrahiere Kovariate mit Dimension:", paste(dim(x), collapse = " x "), "\n")
    x[cal_timesteps,,][which(!is.na(response[cal_timesteps,,]))]
  }))
  
  names(covar) <- names(var_list)
  cat("Die Kovariaten wurden zu einem Data Frame zusammengefügt.\n")
  cat("Dimension von covar:", paste(dim(covar), collapse = " x "), "\n\n")
  
  formula <- formula(
    paste0("tf ~ te(", paste(names(covar), collapse = ","), ")")
  )
  cat("Die Modellformel wurde erzeugt:\n")
  print(formula)
  cat("\nBedeutung: Die Zielvariable tf wird als glatte gemeinsame Funktion der Kovariaten modelliert.\n\n")
  
  data <- cbind.data.frame(
    tf = response[cal_timesteps,,][which(!is.na(response[cal_timesteps,,]))],
    covar
  )
  cat("Der Trainingsdatensatz wurde aufgebaut.\n")
  cat("Dimension von data:", paste(dim(data), collapse = " x "), "\n")
  cat("Erste Zeilen des Trainingsdatensatzes:\n")
  print(utils::head(data))
  cat("\nJetzt wird das GAM mit mgcv::gam() geschätzt.\n")
  cat("Verwendete family:\n")
  print(family)
  cat("\n")
  
  emulator <- mgcv::gam(formula, data = data, family = family)
  
  cat("Das GAM wurde erfolgreich trainiert.\n")
  cat("Kurze Modellzusammenfassung:\n")
  print(summary(emulator))
  cat("--- train_emulator() beendet ---\n\n")
  
  emulator
}

evaluate_emulator <- function(emulator, response, var_list, val_timesteps) {
  cat("--- evaluate_emulator() startet ---\n")
  cat("Jetzt wird das trainierte Modell auf den Validierungs-Zeitschritten ausgewertet.\n")
  cat("Verwendete Validierungs-Zeitschritte:", paste(val_timesteps, collapse = ", "), "\n\n")
  
  covar <- do.call(cbind.data.frame, lapply(var_list, function(x){
    cat("Extrahiere Validierungswerte für Kovariate mit Dimension:", paste(dim(x), collapse = " x "), "\n")
    x[val_timesteps,,][which(!is.na(response[val_timesteps,,]))]
  }))
  
  names(covar) <- names(var_list)
  cat("Die Kovariaten für die Validierung wurden aufgebaut.\n")
  cat("Dimension von covar:", paste(dim(covar), collapse = " x "), "\n\n")
  
  data <- cbind.data.frame(
    tf = response[val_timesteps,,][which(!is.na(response[val_timesteps,,]))],
    covar
  )
  
  cat("Der Validierungsdatensatz wurde aufgebaut.\n")
  cat("Dimension von data:", paste(dim(data), collapse = " x "), "\n")
  cat("Erste Zeilen des Validierungsdatensatzes:\n")
  print(utils::head(data))
  cat("\nJetzt werden Vorhersagen mit predict(..., type = 'response') erzeugt.\n")
  cat("Das bedeutet: Vorhersagen auf der Antwortskala der Zielvariablen.\n\n")
  
  predict_data <- predict(emulator, type = "response", newdata = data)
  
  cat("Die Vorhersagen wurden berechnet.\n")
  cat("Länge von predict_data:", length(predict_data), "\n")
  cat("Erste Vorhersagen:\n")
  print(utils::head(predict_data))
  cat("\nJetzt werden die Vorhersagen zurück in die ursprüngliche Raum-Zeit-Struktur geschrieben.\n")
  
  datafield_predict <- response[val_timesteps,,]
  datafield_predict[which(!is.na(datafield_predict))] <- predict_data
  
  cat("Dimension von datafield_predict:", paste(dim(datafield_predict), collapse = " x "), "\n")
  cat("\nNun wird für jeden Gitterpunkt ein RMSE über die Validierungszeiten berechnet.\n")
  
  rmse_field <- array(NA, dim = dim(datafield_predict)[2:3])
  for (i in 1:dim(rmse_field)[1]) {
    for (j in 1:dim(rmse_field)[2]) {
      rmse_field[i,j] <- suppressWarnings(
        sqrt(mean((response[val_timesteps,i,j] - datafield_predict[,i,j])^2))
      )
    }
  }
  
  cat("Die räumliche RMSE-Karte wurde berechnet.\n")
  cat("Dimension von rmse_field:", paste(dim(rmse_field), collapse = " x "), "\n\n")
  
  results <- data.frame(
    rmse = sqrt(mean((predict_data - data$tf)^2)),
    mae = mean(abs(predict_data - data$tf)),
    Expl_var = 1 - mean((predict_data - data$tf)^2) / var(data$tf)
  )
  
  cat("Globale Gütemaße wurden berechnet:\n")
  print(results)
  cat("--- evaluate_emulator() beendet ---\n\n")
  
  list(
    results = results,
    rmse_field = rmse_field
  )
}

cat("--- Lade Daten aus NetCDF-Dateien ---\n")

nc <- nc_open("MPI_V3_G_TreeCover_mean_time.nc")
cat("Datei MPI_V3_G_TreeCover_mean_time.nc geöffnet.\n")
lon <- ncvar_get(nc, "lon")
lat <- rev(ncvar_get(nc, "lat"))
time <- ncvar_get(nc, "time")[1:40]
tf <- aperm(ncvar_get(nc, "cover_fract")[,48:1,1:40], c(3,1,2))
cat("Tree cover wurde geladen.\n")
cat("Länge lon:", length(lon), "| Länge lat:", length(lat), "| Länge time:", length(time), "\n")
cat("Dimension von tf (Zeit x Lon x Lat):", paste(dim(tf), collapse = " x "), "\n\n")
nc_close(nc)

nc <- nc_open("MPI_V3_G_Pann_mean_time.nc")
cat("Datei MPI_V3_G_Pann_mean_time.nc geöffnet.\n")
time_tmp <- ncvar_get(nc, "time")
precip <- aperm(ncvar_get(nc, "precip")[,48:1,which(time_tmp %in% time)], c(3,1,2))
cat("Niederschlag wurde geladen.\n")
cat("Dimension von precip:", paste(dim(precip), collapse = " x "), "\n\n")
nc_close(nc)

nc <- nc_open("MPI_V3_G_tcold_mean_time.nc")
cat("Datei MPI_V3_G_tcold_mean_time.nc geöffnet.\n")
time_tmp <- ncvar_get(nc, "time")
tcold <- aperm(ncvar_get(nc, "temp2")[,48:1,which(time_tmp %in% time)], c(3,1,2))
cat("Kalte Temperaturvariable tcold wurde geladen.\n")
cat("Dimension von tcold:", paste(dim(tcold), collapse = " x "), "\n\n")
nc_close(nc)

nc <- nc_open("MPI_V3_G_twarm_mean_time.nc")
cat("Datei MPI_V3_G_twarm_mean_time.nc geöffnet.\n")
time_tmp <- ncvar_get(nc, "time")
twarm <- aperm(ncvar_get(nc, "temp2")[,48:1,which(time_tmp %in% time)], c(3,1,2))
cat("Warme Temperaturvariable twarm wurde geladen.\n")
cat("Dimension von twarm:", paste(dim(twarm), collapse = " x "), "\n\n")
nc_close(nc)

cat("--- Lade und verarbeite CO2-Zeitreihe ---\n")
co2table <- read.table("CO2_stack_156K_spline_V2.tab", header = TRUE, sep = "\t", skip = 13)
co2_timeseries <- data.frame(age = -co2table$Age..ka.BP.*1000, co2 = co2table$CO2..µmol.mol.)
co2_timeseries <- sapply(seq(-20000, -500, by = 500), function(t) {
  mean(co2_timeseries$co2[which(co2_timeseries$age >= t & co2_timeseries$age < t + 500)])
})
co2 <- aperm(array(rep(co2_timeseries, each = length(lon) * length(lat)), dim = dim(tf)[c(2,3,1)]), c(3,1,2))
cat("Die CO2-Zeitreihe wurde eingelesen und auf das Raum-Zeit-Gitter erweitert.\n")
cat("Dimension von co2:", paste(dim(co2), collapse = " x "), "\n\n")

cat("--- Erzeuge Trainings- und Validierungs-Zeitschritte ---\n")
training_data_fraction <- 0.5 # use every second time slice for fitting (50%) and other half for validation
cal_timesteps <- seq(1, length(time), by = 1 / training_data_fraction)
if (training_data_fraction == 1) {
  val_timesteps <- cal_timesteps
} else {
  val_timesteps <- (1:length(time))[!(1:length(time)) %in% cal_timesteps]
}
cat("training_data_fraction =", training_data_fraction, "\n")
cat("Trainings-Zeitschritte:", paste(cal_timesteps, collapse = ", "), "\n")
cat("Validierungs-Zeitschritte:", paste(val_timesteps, collapse = ", "), "\n\n")

# var_list <- list(precip=precip, tcold=tcold, twarm=twarm, co2=co2)
var_list <- list(precip = precip, tcold = tcold, twarm = twarm)
response <- tf

cat("--- Definiere Modellinput ---\n")
cat("Es werden folgende Kovariaten verwendet:", paste(names(var_list), collapse = ", "), "\n")
cat("Die Zielvariable ist tf (Vegetations-/Tree-Cover).\n\n")

emulator <- train_emulator(
  response,
  var_list = var_list,
  family = binomial,
  cal_timesteps = cal_timesteps
)

eval_out <- evaluate_emulator(
  emulator,
  response = response,
  var_list = var_list,
  val_timesteps = val_timesteps
)

cat("=== Skriptende ===\n")
cat("Das Modell wurde trainiert und ausgewertet.\n")
cat("Die globalen Evaluationsmetriken stehen in eval_out$results.\n")
cat("Die räumliche RMSE-Karte steht in eval_out$rmse_field.\n")
