#0. INSTALL AND LOAD PACKAGES
# =============================
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("rJava")
install.packages("viridis")

# Load libraries
library(dismo)
library(terra)
library(sf)
library(ggplot2)
library(maxnet)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgbif)
library(sdmpredictors)
library(rJava)
library(viridis)

Sys.setenv(JAVA_HOME = "C:/Users/israe/AppData/Local/Programs/Eclipse Adoptium/jdk-11.0.27.6-hotspot")

# =============================
# 1. DOWNLOAD OCCURRENCE DATA FROM GBIF
# =============================
name_back <- name_backbone(name = "Rhincodon typus")
taxon_key <- name_back$usageKey

occs_raw <- occ_search(
  taxonKey = taxon_key,
  limit = 2000,
  hasCoordinate = TRUE,
  hasGeospatialIssue = FALSE
)
occs <- occs_raw$data
occs <- occs[!is.na(occs$decimalLongitude) & !is.na(occs$decimalLatitude), ]
occs <- occs[!duplicated(occs[, c("decimalLongitude", "decimalLatitude")]), ]
occs_sf <- st_as_sf(occs, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# =============================
# 2. LOAD CURRENT ENVIRONMENTAL VARIABLES (Bio-ORACLE)
# =============================
selected_layers <- c(
  "BO_sstmean", "BO_salinity", "BO2_dissoxmean_ss",
  "BO2_chlomean_ss", "BO_ph", "BO_bathymean"
)

env_present <- load_layers(
  selected_layers,
  rasterstack = FALSE,
  equalarea = FALSE,
  datadir = "data/present"
)

zip_files <- list.files("data/present", pattern = "\\.zip$", full.names = TRUE)
lapply(zip_files, unzip, exdir = "data/present")

env_files <- list.files("data/present", pattern = "\\.tif$", full.names = TRUE)
env_present <- rast(env_files)

# =============================
# 3. EXTRACT ENVIRONMENTAL VALUES FOR OCCURRENCES
# =============================
occs_vals <- terra::extract(env_present, vect(occs_sf), ID = FALSE)
occs_vals <- na.omit(occs_vals)

# =============================
# 4. GENERATE BACKGROUND POINTS
# =============================
bathy <- env_present$BO_bathymean_lonlat
ocean_mask <- bathy < -20
set.seed(123)
bg_points <- spatSample(ocean_mask, size = 10000, method = "random", na.rm = TRUE, as.points = TRUE)
bg_vals <- terra::extract(env_present, bg_points, ID = FALSE)
bg_vals <- na.omit(bg_vals)

# =============================
# 5. TRAIN MAXNET MODEL
# =============================
p <- rep(1, nrow(occs_vals))
a <- rep(0, nrow(bg_vals))
all_vals <- rbind(occs_vals, bg_vals)
pa <- c(p, a)

keep_rows <- complete.cases(all_vals)
all_vals <- all_vals[keep_rows, ]
pa <- pa[keep_rows]

names(env_present) <- make.names(names(env_present))
colnames(all_vals) <- make.names(colnames(all_vals))

mod <- maxnet(pa, all_vals, maxnet.formula(pa, all_vals, classes = "default"))
eval <- evaluate(p = occs_vals, a = bg_vals, model = mod)
print(paste("AUC:", round(eval@auc, 3)))

env_present_model <- env_present[[colnames(all_vals)]]
pred_present <- predict(env_present_model, mod, type = "logistic", na.rm = TRUE)
pred_present_masked <- mask(pred_present, ocean_mask)
plot(pred_present_masked, col = viridis(100), main = "Habitat suitability (MaxNet-logistic) - Rhincodon typus")
# =============================
# 6. LOAD AND PREPARE FUTURE ENVIRONMENTAL VARIABLES
# =============================
ruta_2050 <- "C:/Users/israe/OneDrive/Desktop/R/Ocean layers 2050"
chl <- rast(file.path(ruta_2050, "Clorophyll 2050.nc"))
oxy <- rast(file.path(ruta_2050, "Dissolved_Oxygen_2050.nc"))
ph  <- rast(file.path(ruta_2050, "pH 2050.nc"))
sal <- rast(file.path(ruta_2050, "Salinity 2050.nc"))
sst <- rast(file.path(ruta_2050, "Temperature 2050.nc"))

sst_proj <- project(sst, bathy)
sal_proj <- project(sal, bathy)
oxy_proj <- project(oxy, bathy)
chl_proj <- project(chl, bathy)
ph_proj  <- project(ph, bathy)

names(sst_proj) <- "BO_sstmean_lonlat"
names(sal_proj) <- "BO_salinity_lonlat"
names(oxy_proj) <- "Present.Surface.Dissolved.oxygen.Mean"
names(chl_proj) <- "Present.Surface.Chlorophyll.Mean"
names(ph_proj)  <- "BO_ph_lonlat"
names(bathy)    <- "BO_bathymean_lonlat"

env_future <- c(sst_proj, sal_proj, oxy_proj, chl_proj, ph_proj, bathy)
names(env_future) <- make.names(names(env_future))

# =============================
# 7. PREDICT FUTURE SUITABILITY
# =============================
env_future_model <- env_future[[colnames(all_vals)]]
pred_future <- predict(env_future_model, mod, type = "logistic", na.rm = TRUE)
pred_future_masked <- mask(pred_future, bathy < -20)
plot(pred_future_masked, col = viridis(100), main = "Suitability 2050 - Rhincodon typus (SSP5-8.5)")
# =============================
# 8. CALCULATE CHANGE IN SUITABILITY
# =============================
change_map <- pred_future_masked - pred_present_masked
plot(change_map, col = rev(terrain.colors(100)), 
     main = "Change in Habitat Suitability (2050 vs Present)")

# =============================
# 9. CHANGE CLASSIFICATION (BASED ON OBSERVED FREQ)
# =============================
loss_cells   <- 22853
stable_cells <- 5629062
gain_cells   <- 83938
total_cells <- loss_cells + stable_cells + gain_cells

loss_pct   <- round(100 * loss_cells / total_cells, 2)
stable_pct <- round(100 * stable_cells / total_cells, 2)
gain_pct   <- round(100 * gain_cells / total_cells, 2)

cat("\U0001F3DD\uFE0F Habitat change for *Rhincodon typus* (threshold = 0.5)\n")
cat("\U0001F53A Gain: ", gain_pct, "%\n")
cat("\U0001F53B Loss: ", loss_pct, "%\n")
cat("\u23F8\uFE0F Stable: ", stable_pct, "%\n")

# =============================
# 10. PIE CHART
# =============================
change_df <- data.frame(
  class = c("Gain", "Loss", "Stable"),
  count = c(gain_cells, loss_cells, stable_cells)
)

ggplot(change_df, aes(x = "", y = count, fill = class)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual(values = c("Gain" = "green3", "Loss" = "red3", "Stable" = "gray60")) +
  labs(title = "Habitat change for Rhincodon typus (SSP5-8.5 - 2050)")
