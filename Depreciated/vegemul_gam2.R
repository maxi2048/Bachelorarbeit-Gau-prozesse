library(mgcv)
library(ncdf4)
library(fields)
library(maps)
library(parallel)

# -------- Funktion: Dataframe erstellen --------------
make_emulator_df <- function(response, var_list, time, lon, lat) {
  
  grid <- expand.grid(
    time_index = seq_along(time),
    lon_index  = seq_along(lon),
    lat_index  = seq_along(lat)
  )
  
  grid$time <- time[grid$time_index]
  grid$lon  <- lon[grid$lon_index]
  grid$lat  <- lat[grid$lat_index]
  
  # Zielvariable
  grid$tf <- as.vector(response)
  
  # Kovariaten
  for (v in names(var_list)) {
    grid[[v]] <- as.vector(var_list[[v]])
  }
  
  # NA entfernen
  grid <- grid[complete.cases(grid), ]
  rownames(grid) <- NULL
  
  return(grid)
}

# -------- Load data --------------
nc <- nc_open("MPI_V3_G_TreeCover_mean_time.nc")
lon <- ncvar_get(nc,"lon")
lat <- rev(ncvar_get(nc,"lat"))
time <- ncvar_get(nc,"time")[1:40]
tf <- aperm(ncvar_get(nc,"cover_fract")[,48:1,1:40],c(3,1,2))

nc <- nc_open("MPI_V3_G_Pann_mean_time.nc")
time_tmp <- ncvar_get(nc,"time")
precip <- aperm(ncvar_get(nc,"precip")[,48:1,which(time_tmp %in% time)],c(3,1,2))

nc <- nc_open("MPI_V3_G_tcold_mean_time.nc")
time_tmp <- ncvar_get(nc,"time")
tcold <- aperm(ncvar_get(nc,"temp2")[,48:1,which(time_tmp %in% time)],c(3,1,2))

nc <- nc_open("MPI_V3_G_twarm_mean_time.nc")
time_tmp <- ncvar_get(nc,"time")
twarm <- aperm(ncvar_get(nc,"temp2")[,48:1,which(time_tmp %in% time)],c(3,1,2))

co2table <- read.table("CO2_stack_156K_spline_V2.tab",header=TRUE,sep = "\t",skip = 13)
co2_timeseries <- data.frame(age = -co2table$Age..ka.BP.*1000, co2 = co2table$CO2..µmol.mol.)
co2_timeseries <- sapply(seq(-20000,-500,by=500), function(t) 
  mean(co2_timeseries$co2[which(co2_timeseries$age >= t & co2_timeseries$age < t+500)])
)
co2 <- aperm(array(rep(co2_timeseries,each=length(lon)*length(lat)),
                   dim=dim(tf)[c(2,3,1)]),c(3,1,2))

# -------- Variablenliste --------------
var_list <- list(
  precip = precip,
  tcold  = tcold,
  twarm  = twarm,
  co2 = co2   
)

# -------- Dataframe erstellen --------------
df <- make_emulator_df(
  response = tf,
  var_list = var_list,
  time = time,
  lon = lon,
  lat = lat
)

# Überblick
head(df)
str(df)

# -------- Train / Validation Split --------------
training_data_fraction <- 0.5

cal_timesteps <- seq(1, length(time), by = 1 / training_data_fraction)

if (training_data_fraction == 1) {
  val_timesteps <- cal_timesteps
} else {
  val_timesteps <- (1:length(time))[!(1:length(time) %in% cal_timesteps)]
}

df_train <- df[df$time_index %in% cal_timesteps, ]
df_val   <- df[df$time_index %in% val_timesteps, ]

# -------- Modell (GAM) --------------
n_cores <- max(1, parallel::detectCores() - 1)

cl <- parallel::makeCluster(n_cores)

emulator <- mgcv::gam(
  tf ~ te(precip, tcold, twarm),
  data = df_train,
  family = binomial,
  # method = "REML",
  cluster = cl
)

parallel::stopCluster(cl)

# -------- Vorhersage --------------
df_val$pred <- predict(emulator, newdata = df_val, type = "response")

# -------- Evaluation --------------
rmse <- sqrt(mean((df_val$pred - df_val$tf)^2))
mae  <- mean(abs(df_val$pred - df_val$tf))
expl_var <- 1 - mean((df_val$pred - df_val$tf)^2) / var(df_val$tf)

eval_results <- data.frame(
  rmse = rmse,
  mae = mae,
  Expl_var = expl_var
)

print(eval_results)

# -------- RMSE Map erstellen --------------
rmse_field <- array(NA, dim = c(length(lon), length(lat)))

for (i in seq_along(lon)) {
  for (j in seq_along(lat)) {
    
    subset_data <- df_val[df_val$lon_index == i & df_val$lat_index == j, ]
    
    if (nrow(subset_data) > 0) {
      rmse_field[i,j] <- sqrt(mean((subset_data$tf - subset_data$pred)^2))
    }
  }
}

# -------- Visualisierung --------------

col_magnitude_ipcc_green <- function(n=21) {
  grDevices::colorRampPalette(c(
    rgb(255,255,204,maxColorValue = 255),
    rgb(194,230,153,maxColorValue = 255),
    rgb(120,198,121,maxColorValue = 255),
    rgb(49,163,84,maxColorValue = 255),
    rgb(0,104,55,maxColorValue = 255)
  ))(n)
}

lon_180 <- c(lon[49:96]-360, lon[1:48])

fields::image.plot(
  lon_180,
  lat,
  tf[1,c(49:96,1:48),],
  col = col_magnitude_ipcc_green()
)
maps::map(add = TRUE, interior = FALSE)

# -------- Vorhersage --------------
df_val$pred <- predict(emulator, newdata = df_val, type = "response")

# -------- Objekte für Residualanalyse.R vorbereiten --------------
df_train$pred_train <- predict(emulator, newdata = df_train, type = "response")

df_val$resid_response <- df_val$tf - df_val$pred
df_val$resid_pearson <- (df_val$tf - df_val$pred) /
  sqrt(pmax(df_val$pred * (1 - df_val$pred), 1e-8))

df_train$resid_response <- df_train$tf - df_train$pred_train
df_train$resid_pearson <- (df_train$tf - df_train$pred_train) /
  sqrt(pmax(df_train$pred_train * (1 - df_train$pred_train), 1e-8))

df_resid_compare <- rbind(
  data.frame(
    Datensatz = "Training",
    resid = df_train$resid_response
  ),
  data.frame(
    Datensatz = "Validierung",
    resid = df_val$resid_response
  )
)

df_compare <- df_resid_compare

summary_text <- paste0(
  "Validierung\n",
  "Bias = ", round(mean(df_val$resid_response, na.rm = TRUE), 5), "\n",
  "RMSE = ", round(sqrt(mean(df_val$resid_response^2, na.rm = TRUE)), 5), "\n",
  "MAE = ", round(mean(abs(df_val$resid_response), na.rm = TRUE), 5), "\n",
  "SD = ", round(sd(df_val$resid_response, na.rm = TRUE), 5)
)


# 
# df_val$resid_response <- df_val$tf - df_val$pred
# 
# df_val$resid_pearson <- (df_val$tf - df_val$pred) /
#   sqrt(pmax(df_val$pred * (1 - df_val$pred), 1e-8))
# 
# # Training
# df_train$pred_train <- predict(emulator, newdata = df_train, type = "response")
# 
# df_train$resid_response <- df_train$tf - df_train$pred_train
# 
# df_train$resid_pearson <- (df_train$tf - df_train$pred_train) /
#   sqrt(pmax(df_train$pred_train * (1 - df_train$pred_train), 1e-8))

# Kürzel
r_val  <- df_val$resid_response
rp_val <- df_val$resid_pearson
fit    <- df_val$pred

n_lon <- length(lon)
n_lat <- length(lat)

