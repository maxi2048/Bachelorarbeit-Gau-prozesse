library(mgcv)
library(ncdf4)
library(fields)
library(maps)

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

# Optional: CO2 (falls du es wieder reinnehmen willst)
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
  twarm  = twarm
  # co2 = co2   # optional
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
emulator <- mgcv::gam(
  tf ~ te(precip, tcold, twarm),
  data = df_train,
  family = binomial
)

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

### DOne 

eps <- 0.001
df_train$tf_beta <- pmin(pmax(df_train$tf, eps), 1 - eps)
df_val$tf_beta   <- pmin(pmax(df_val$tf, eps), 1 - eps)

emulator_gp <- mgcv::gam(
  tf_beta ~ 
    te(precip, tcold, twarm, k = c(10, 10, 10)) +
    s(lon, lat, time, bs = "gp", k = 200),
  data = df_train,
  family = betar(),
  method = "REML"
)

df_val$pred_gp <- predict(
  emulator_gp,
  newdata = df_val,
  type = "response"
)

rmse_gp <- sqrt(mean((df_val$pred_gp - df_val$tf)^2))
mae_gp  <- mean(abs(df_val$pred_gp - df_val$tf))
expl_var_gp <- 1 - mean((df_val$pred_gp - df_val$tf)^2) / var(df_val$tf)

data.frame(
  rmse = rmse_gp,
  mae = mae_gp,
  Expl_var = expl_var_gp
)

gam.check(emulator_gp)

summary(emulator_gp)