# Ehemalig Vegemul_Gam 
# Daten in Array laden und in Df Datenstruktur bringen
# GAM rechnen lassen und residuen bestimmen


library(mgcv)
library(ncdf4)
library(fields)
library(maps)
library(ggplot2)
library(dplyr)


theme_simple <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )


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


# Lattitude geht von -87.159095 bis 87.159095 Süden Norden
# Longitude geht von 0 bis 360 Osten Westen
# time 40 mal 500 Jahre


covariates <- list(precip = precip,
                   co2    = co2,
                   tcold  = tcold,
                   twarm  = twarm
                   )


grid <- expand.grid(
  time_index = seq_along(time),
  lon_index  = seq_along(lon),
  lat_index  = seq_along(lat)
)

grid$time <- time[grid$time_index]
grid$lon  <- lon[grid$lon_index]
grid$lat  <- lat[grid$lat_index]

# Indexwerte sind bei time 1 - 40, lon 1 - 96, lat 1 - 48


# Zielvariable
grid$tf <- as.vector(tf)

# Kovariaten
for (v in names(covariates)) {
  grid[[v]] <- as.vector(covariates[[v]])
}

grid$lon_plot <- ifelse(grid$lon > 180, grid$lon - 360, grid$lon)

na_map_tf <- grid %>%
  group_by(lon_plot, lat) %>%
  summarise(
    n_na_tf = sum(is.na(tf)),
    .groups = "drop"
  )

head(na_map_tf)

world <- map_data("world")

ggplot() +
  geom_tile(
    data = na_map_tf,
    aes(x = lon_plot, y = lat, fill = n_na_tf)
  ) +
  geom_polygon(
    data = world,
    aes(x = long, y = lat, group = group),
    fill = NA,
    color = "black",
    linewidth = 0.2
  ) +
  coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90)) +
  labs(
    title = "Anzahl fehlender Werte in TreeCover",
    subtitle = "Gezählt über alle 40 Zeitpunkte",
    x = "Longitude",
    y = "Latitude",
    fill = "Anzahl NAs"
  ) +
  theme_simple


# NA entfernen
grid <- grid[complete.cases(grid), ]
rownames(grid) <- NULL

df <- grid



training_data_fraction <- 0.5

cal_timesteps <- seq(1, length(time), by = 1 / training_data_fraction)

if (training_data_fraction == 1) {
  val_timesteps <- cal_timesteps
} else {
  val_timesteps <- (1:length(time))[!(1:length(time) %in% cal_timesteps)]
}

df_train <- df[df$time_index %in% cal_timesteps, ]
df_val   <- df[df$time_index %in% val_timesteps, ]
# 


# -------- Modell (GAM) --------------
# n_cores <- max(1, parallel::detectCores() - 1)
# 
# cl <- parallel::makeCluster(n_cores)
# 
# emulator <- mgcv::gam(
#   tf ~ te(precip, tcold, twarm, co2),
#   data = df_train,
#   family = quasibinomial,
#   cluster = cl
# )
# 
# parallel::stopCluster(cl)
# 
# 
# emulator

saveRDS(emulator, file = "emulator_gam.rds")

readRDS("emulator_gam.rds")

# Interessant! CO2 hebt R^2 nur um 0.012 an
summary(emulator)

df_val$pred <- predict(emulator, newdata = df_val, type = "response")

rmse <- sqrt(mean((df_val$pred - df_val$tf)^2))
mae  <- mean(abs(df_val$pred - df_val$tf))
expl_var <- 1 - mean((df_val$pred - df_val$tf)^2) / var(df_val$tf)

eval_results <- data.frame(
  rmse = rmse,
  mae = mae,
  Expl_var = expl_var
)

print(eval_results)


df_train$pred_train <- predict(emulator, newdata = df_train, type = "response")

df_val$pred <- predict(emulator, newdata = df_val, type = "response")

df_val$resid_response <- df_val$tf - df_val$pred
df_train$resid_response <- df_train$tf - df_train$pred_train

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

# Kleinster und größter Wert aus Covariate, dann 200 Werte dazwischen generieren
# Die restlichen Kovariaten bleiben gleich (Median)
# Predict Vorhersage durch den emulator mit binomial Linkfunktion


new_precip <- data.frame(
  precip = seq(
    min(df_train$precip, na.rm = TRUE),
    max(df_train$precip, na.rm = TRUE),
    length.out = 200
  ),
  tcold = median(df_train$tcold, na.rm = TRUE),
  twarm = median(df_train$twarm, na.rm = TRUE),
  co2   = median(df_train$co2, na.rm = TRUE)
)

pred_precip <- predict(
  emulator,
  newdata = new_precip,
  type = "link",
  se.fit = TRUE
)

new_precip$fit <- plogis(pred_precip$fit)
new_precip$lwr <- plogis(pred_precip$fit - 2 * pred_precip$se.fit)
new_precip$upr <- plogis(pred_precip$fit + 2 * pred_precip$se.fit)

p_precip <- ggplot(new_precip, aes(x = precip, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(
    title = "Partieller Plot für precip",
    subtitle = "tcold, twarm und co2 auf Median fixiert",
    x = "precip",
    y = "geschätzter TreeCover"
  ) +
  theme_simple

print(p_precip)


new_tcold <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold = seq(
    min(df_train$tcold, na.rm = TRUE),
    max(df_train$tcold, na.rm = TRUE),
    length.out = 200
  ),
  twarm = median(df_train$twarm, na.rm = TRUE),
  co2   = median(df_train$co2, na.rm = TRUE)
)

pred_tcold <- predict(
  emulator,
  newdata = new_tcold,
  type = "link",
  se.fit = TRUE
)

# Intervalle 
new_tcold$fit <- plogis(pred_tcold$fit)
new_tcold$lwr <- plogis(pred_tcold$fit - 2 * pred_tcold$se.fit)
new_tcold$upr <- plogis(pred_tcold$fit + 2 * pred_tcold$se.fit)

p_tcold <- ggplot(new_tcold, aes(x = tcold, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(
    title = "Partieller Plot für tcold",
    subtitle = "precip, twarm und co2 auf Median fixiert",
    x = "tcold",
    y = "geschätzter TreeCover"
  ) + theme_simple

print(p_tcold)


new_twarm <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold  = median(df_train$tcold, na.rm = TRUE),
  twarm = seq(
    min(df_train$twarm, na.rm = TRUE),
    max(df_train$twarm, na.rm = TRUE),
    length.out = 200
  ),
  co2 = median(df_train$co2, na.rm = TRUE)
)

pred_twarm <- predict(
  emulator,
  newdata = new_twarm,
  type = "link",
  se.fit = TRUE
)

new_twarm$fit <- plogis(pred_twarm$fit)
new_twarm$lwr <- plogis(pred_twarm$fit - 2 * pred_twarm$se.fit)
new_twarm$upr <- plogis(pred_twarm$fit + 2 * pred_twarm$se.fit)

p_twarm <- ggplot(new_twarm, aes(x = twarm, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(
    title = "Partieller Plot für twarm",
    subtitle = "precip, tcold und co2 auf Median fixiert",
    x = "twarm",
    y = "geschätzter TreeCover"
  ) +
  theme_simple

print(p_twarm)



new_co2 <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold  = median(df_train$tcold, na.rm = TRUE),
  twarm  = median(df_train$twarm, na.rm = TRUE),
  co2 = seq(
    min(df_train$co2, na.rm = TRUE),
    max(df_train$co2, na.rm = TRUE),
    length.out = 200
  )
)

pred_co2 <- predict(
  emulator,
  newdata = new_co2,
  type = "link",
  se.fit = TRUE
)

new_co2$fit <- plogis(pred_co2$fit)
new_co2$lwr <- plogis(pred_co2$fit - 2 * pred_co2$se.fit)
new_co2$upr <- plogis(pred_co2$fit + 2 * pred_co2$se.fit)

p_co2 <- ggplot(new_co2, aes(x = co2, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(linewidth = 1) +
  labs(
    title = "Partieller Plot für co2",
    subtitle = "precip, tcold und twarm auf Median fixiert",
    x = "co2",
    y = "geschätzter TreeCover"
  ) +
  theme_simple

print(p_co2)



new_precip <- data.frame(
  precip = seq(
    min(df_train$precip, na.rm = TRUE),
    max(df_train$precip, na.rm = TRUE),
    length.out = 200
  ),
  tcold = median(df_train$tcold, na.rm = TRUE),
  twarm = median(df_train$twarm, na.rm = TRUE),
  co2   = median(df_train$co2, na.rm = TRUE)
)

pred_precip <- predict(
  emulator,
  newdata = new_precip,
  type = "link",
  se.fit = TRUE
)

new_precip$fit <- plogis(pred_precip$fit)
new_precip$lwr <- plogis(pred_precip$fit - 2 * pred_precip$se.fit)
new_precip$upr <- plogis(pred_precip$fit + 2 * pred_precip$se.fit)

new_precip$variable <- "precip"
new_precip$x <- new_precip$precip

new_tcold <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold = seq(
    min(df_train$tcold, na.rm = TRUE),
    max(df_train$tcold, na.rm = TRUE),
    length.out = 200
  ),
  twarm = median(df_train$twarm, na.rm = TRUE),
  co2   = median(df_train$co2, na.rm = TRUE)
)

pred_tcold <- predict(
  emulator,
  newdata = new_tcold,
  type = "link",
  se.fit = TRUE
)

new_tcold$fit <- plogis(pred_tcold$fit)
new_tcold$lwr <- plogis(pred_tcold$fit - 2 * pred_tcold$se.fit)
new_tcold$upr <- plogis(pred_tcold$fit + 2 * pred_tcold$se.fit)

new_tcold$variable <- "tcold"
new_tcold$x <- new_tcold$tcold

new_twarm <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold  = median(df_train$tcold, na.rm = TRUE),
  twarm = seq(
    min(df_train$twarm, na.rm = TRUE),
    max(df_train$twarm, na.rm = TRUE),
    length.out = 200
  ),
  co2 = median(df_train$co2, na.rm = TRUE)
)

pred_twarm <- predict(
  emulator,
  newdata = new_twarm,
  type = "link",
  se.fit = TRUE
)

new_twarm$fit <- plogis(pred_twarm$fit)
new_twarm$lwr <- plogis(pred_twarm$fit - 2 * pred_twarm$se.fit)
new_twarm$upr <- plogis(pred_twarm$fit + 2 * pred_twarm$se.fit)

new_twarm$variable <- "twarm"
new_twarm$x <- new_twarm$twarm

new_co2 <- data.frame(
  precip = median(df_train$precip, na.rm = TRUE),
  tcold  = median(df_train$tcold, na.rm = TRUE),
  twarm  = median(df_train$twarm, na.rm = TRUE),
  co2 = seq(
    min(df_train$co2, na.rm = TRUE),
    max(df_train$co2, na.rm = TRUE),
    length.out = 200
  )
)

pred_co2 <- predict(
  emulator,
  newdata = new_co2,
  type = "link",
  se.fit = TRUE
)

new_co2$fit <- plogis(pred_co2$fit)
new_co2$lwr <- plogis(pred_co2$fit - 2 * pred_co2$se.fit)
new_co2$upr <- plogis(pred_co2$fit + 2 * pred_co2$se.fit)

new_co2$variable <- "co2"
new_co2$x <- new_co2$co2


partial_all <- bind_rows(
  new_precip,
  new_tcold,
  new_twarm,
  new_co2
)

partial_all$variable <- factor(
  partial_all$variable,
  levels = c("precip", "tcold", "twarm", "co2")
)


p_partial_all <- ggplot(partial_all, aes(x = x, y = fit)) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.2
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ variable,
    scales = "free_x"
  ) +
  labs(
    title = "Partielle Vorhersageplots des GAM",
    subtitle = "Jeweils eine Kovariate wird variiert, alle anderen werden auf ihrem Median fixiert",
    x = "Wert der Kovariate",
    y = "geschätzter TreeCover"
  ) +
  theme_simple

print(p_partial_all)

make_pdp <- function(model, data, variable, grid_length = 200) {
  
  covars <- c("precip", "tcold", "twarm", "co2")
  
  grid <- seq(
    min(data[[variable]], na.rm = TRUE),
    max(data[[variable]], na.rm = TRUE),
    length.out = grid_length
  )
  
  pdp <- lapply(grid, function(z) {
    
    newdata <- data[, covars]
    newdata[[variable]] <- z
    
    pred <- predict(
      model,
      newdata = newdata,
      type = "response"
    )
    
    data.frame(
      variable = variable,
      x = z,
      fit = mean(pred, na.rm = TRUE)
    )
  }) |>
    bind_rows()
  
  pdp
}



pdp_all <- bind_rows(
  make_pdp(emulator, df_train, "precip"),
  make_pdp(emulator, df_train, "tcold"),
  make_pdp(emulator, df_train, "twarm"),
  make_pdp(emulator, df_train, "co2")
)

pdp_all$variable <- factor(
  pdp_all$variable,
  levels = c("precip", "tcold", "twarm", "co2")
)

p_pdp_all <- ggplot(pdp_all, aes(x = x, y = fit)) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ variable,
    scales = "free_x"
  ) +
  labs(
    title = "Partial Dependence Plots des GAM",
    subtitle = "Jeweils eine Kovariate wird variiert, über die beobachteten Werte der übrigen Kovariaten wird gemittelt",
    x = "Wert der Kovariate",
    y = "mittlere vorhergesagte Response"
  ) +
  theme_simple

print(p_pdp_all)

make_pdp_band <- function(model, data, variable, grid_length = 200) {
  
  covars <- c("precip", "tcold", "twarm", "co2")
  
  grid <- seq(
    min(data[[variable]], na.rm = TRUE),
    max(data[[variable]], na.rm = TRUE),
    length.out = grid_length
  )
  
  pdp <- lapply(grid, function(z) {
    
    newdata <- data[, covars]
    newdata[[variable]] <- z
    
    pred <- predict(
      model,
      newdata = newdata,
      type = "response"
    )
    
    data.frame(
      variable = variable,
      x = z,
      fit = mean(pred, na.rm = TRUE),
      lwr = quantile(pred, 0.05, na.rm = TRUE),
      upr = quantile(pred, 0.95, na.rm = TRUE)
    )
  }) |>
    bind_rows()
  
  pdp
}

pdp_all_band <- bind_rows(
  make_pdp_band(emulator, df_train, "precip"),
  make_pdp_band(emulator, df_train, "tcold"),
  make_pdp_band(emulator, df_train, "twarm"),
  make_pdp_band(emulator, df_train, "co2")
)

pdp_all_band$variable <- factor(
  pdp_all_band$variable,
  levels = c("precip", "tcold", "twarm", "co2")
)

p_pdp_all_band <- ggplot(pdp_all_band, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    title = "Partial Dependence Plots des GAM",
    subtitle = "Linie: mittlere Vorhersage; Band: 5%- bis 95%-Quantile der Vorhersagen",
    x = "Wert der Kovariate",
    y = "mittlere vorhergesagte Response"
  ) +
  theme_simple

print(p_pdp_all_band)

