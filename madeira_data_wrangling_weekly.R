## Exploration and wrangling of Madeira mosquito ovitrap data 

library(data.table); library(dplyr);library(plyr);library(sf); library(terra); library(stringr); library(lubridate);library(tidyverse)

## -- read in datasets

data = fread("ovitrap_data_2012.csv")
locs = fread("ovitrap_locations_2012.csv")
colnames(locs)[2:3] <- c("LAT_WGS", "LONG_WGS")
traps = fread("Madeira traps info.csv")

# remove row with "withdrawal" in values 
data = data[!data$value == "retirada", ]
data$value = as.numeric(data$value)
# set date to normal date format
data$date <- as.Date(data$date, format = "%m/%d/%Y")
data$month = lubridate::month(data$date)
data$year = lubridate::year(data$date)
data$week = lubridate::week(data$date)

head(data)
head(locs)
## combine data 
data = merge(data, locs, by.x = "Id.Novo", by.y = "Id Novo") ## trap data is missing data for trap 75
data = merge(data, traps[, c(1:8,11:14)], by.x = "Id.Novo", by.y = "ID")
head(data)
str(data)
# set date to normal date format
#data$date <- as.Date(data$date, format = "%m/%d/%Y")

data$presence = as.numeric(data$value > 0)
hist(data$presence)

# have a look at what's going on
wsg_crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs"
library(ggplot2)
d_sf = st_as_sf(data, coords = c("LONG_WGS", "LAT_WGS"), crs =sf::st_crs(wsg_crs))
plot(d_sf)
d2 = data[data$value >0 , ]
ggplot(d2, aes(x = date, fill = Id.Novo)) + geom_histogram(binwidth = 7)

ggplot(data, aes(x = Id.Novo, fill = Id.Novo)) +geom_histogram(binwidth = 1)

nrow(locs)

##map locations
coordinates(locs) = ~LONG_WGS + LAT_WGS 
plot(locs)


# read in ESA 10m landcover data 
## -- !! don't need to do this repeatedly !! -- ##
## -- uncomment for the first run 
#esa1 = rast("ESA_WorldCover_10m_2020_v100_N30W018_Map.tif")
#esa2 = rast("ESA_WorldCover_10m_2020_v100_N33W018_Map.tif")
#plot(esa1)
#esa1
#plot(esa2)

# combine rasters into one
#esa = raster::merge(esa1, esa2)
#esa = do.call(merge, esa)
plot(esa)

# crop to just show main islands to match with trap locations 
lnd<-st_bbox(st_as_sf(locs, coords = c("LAT_WGS", "LONG_WGS"), crs =sf::st_crs(wsg_crs)))
buff<-0.05

lnd[1]<-lnd[1]-buff
lnd[2]<-lnd[2]-buff
lnd[3]<-lnd[3]+buff
lnd[4]<-lnd[4]+buff
esa2 = crop(esa, lnd)
plot(esa2)
points(locs, col = "pink", pch  = 16)
# save so we don't have to do this slow process again
writeRaster(esa2, "Madeira_ESA_WorldCover_10m_2020.tif")

# read in 
mad_esa = rast("Madeira_ESA_WorldCover_10m_2020.tif") 

plot(mad_esa)
plot(d_sf$geometry, add = T, col = "blue")


# need to calculate the percentage cover of land cover types for buffers - get functions from NBMP analysis 
esa_vals = terra::extract(mad_esa, d_sf) # just the values at the sites
table(esa_vals$Madeira_ESA_WorldCover_10m_2020)

d_sf$esa_10m = esa_vals$Madeira_ESA_WorldCover_10m_2020
## Land cover codes:
# 10 - tree cover 
# 20 Shrubland
# 30 Grassland 
# 40 Cropland 
# 50 Built up
# 60 Bare/ Sparse vegetation
# 80 Permenant water bodies (inludes the sea)

# See 'aggregate_madeira_esa.R' script for the creation of these rasters that 
esa_vals_50m = terra::extract(esa_prop_50m, d_sf, fun = mean)
esa_vals_100m = terra::extract(esa_prop_100m, d_sf, fun = mean)
esa_vals_500m = terra::extract(esa_prop_500m, d_sf, fun = mean)

d_sf = cbind(d_sf, esa_vals_50m)
d_sf = cbind(d_sf, esa_vals_100m[, 2:ncol(esa_vals_100m)])
d_sf = cbind(d_sf, esa_vals_500m[, 2:ncol(esa_vals_500m)])

# elevation data from https://land.copernicus.eu/imagery-in-situ/eu-dem/eu-dem-v1.1/view
ele = terra::rast("eu_dem_v11_E10N10/eu_dem_v11_E10N10.TIF")
ele = terra::project(ele, mad_esa) # reproject and crop to same crs and extent as land cover raster
ele_vals = terra::extract(ele, d_sf, fun = mean)
d_sf = cbind(d_sf, ele_vals[,2])

coords = st_coordinates(d_sf)
colnames(coords) = c("longitude", "latitude" )
data2 = st_drop_geometry(d_sf)
data2 = cbind(data2, coords)
write.csv(data2, "madeira_data_esa_weekly.csv")

###############################################################################
## -- weather station weekly averages
###############################################################################
data2 = fread("madeira_data_esa_weekly.csv") %>%
  select(-"V1")
head(data2)
 
#locations of weather stations
#w_locs = fread("madeira_weatherstation_locations.csv")
w_w <- fread("madeira_weekly_weather.csv")%>%
  mutate(week_year = paste(week, Ano, sep = "_")) %>%
  filter(!is.na(mean_temp_1.5m)) 

w_w <- w_w[,c(3:9,13:15,22:32)] %>%
  group_by(nome.da.estação) %>%
  mutate(tmean_1lag = lag(mean_temp_1.5m, n = 1, order_by = week_year),
         tmean_2lag = lag(mean_temp_1.5m, n = 2, order_by = week_year),
         tmean_3lag = lag(mean_temp_1.5m, n = 3, order_by = week_year),
         tmean_4lag = lag(mean_temp_1.5m, n = 4, order_by = week_year),
         tmean_6lag = lag(mean_temp_1.5m, n = 6, order_by = week_year),
         tmax_1lag = lag(max_temp_1.5m, n = 1, order_by = week_year),
         tmax_2lag = lag(max_temp_1.5m, n = 2, order_by = week_year),
         tmax_3lag = lag(max_temp_1.5m, n = 3, order_by = week_year),
         tmax_4lag = lag(max_temp_1.5m, n = 4, order_by = week_year),
         tmax_6lag = lag(max_temp_1.5m, n = 6, order_by = week_year),
         tmin_1lag = lag(min_temp_1.5m, n = 1, order_by = week_year),
         tmin_2lag = lag(min_temp_1.5m, n = 2, order_by = week_year),
         tmin_3lag = lag(min_temp_1.5m, n = 3, order_by = week_year),
         tmin_4lag = lag(min_temp_1.5m, n = 4, order_by = week_year),
         tmin_6lag = lag(min_temp_1.5m, n = 6, order_by = week_year),
         hmean_1lag = lag(mean_humid, n = 1, order_by = week_year),
         hmean_2lag = lag(mean_humid, n = 2, order_by = week_year),
         hmean_3lag = lag(mean_humid, n = 3, order_by = week_year),
         hmean_4lag = lag(mean_humid, n = 4, order_by = week_year),
         hmean_6lag = lag(mean_humid, n = 6, order_by = week_year),
         hmax_1lag = lag(max_humid, n = 1, order_by = week_year),
         hmax_2lag = lag(max_humid, n = 2, order_by = week_year),
         hmax_3lag = lag(max_humid, n = 3, order_by = week_year),
         hmax_4lag = lag(max_humid, n = 4, order_by = week_year),
         hmax_6lag = lag(max_humid, n = 6, order_by = week_year),
         hmin_1lag = lag(min_humid, n = 1, order_by = week_year),
         hmin_2lag = lag(min_humid, n = 2, order_by = week_year),
         hmin_3lag = lag(min_humid, n = 3, order_by = week_year),
         hmin_4lag = lag(min_humid, n = 4, order_by = week_year),
         hmin_6lag = lag(min_humid, n = 6, order_by = week_year),
         pmean_1lag = lag(mean_precip_mm, n = 1, order_by = week_year),
         pmean_2lag = lag(mean_precip_mm, n = 2, order_by = week_year),
         pmean_3lag = lag(mean_precip_mm, n = 3, order_by = week_year),
         pmean_4lag = lag(mean_precip_mm, n = 4, order_by = week_year),
         pmean_6lag = lag(mean_precip_mm, n = 6, order_by = week_year),
         psum_1lag = lag(sum_precip_mm, n = 1, order_by = week_year),
         psum_2lag = lag(sum_precip_mm, n = 2, order_by = week_year),
         psum_3lag = lag(sum_precip_mm, n = 3, order_by = week_year),
         psum_4lag = lag(sum_precip_mm, n = 4, order_by = week_year),
         psum_6lag = lag(sum_precip_mm, n = 6, order_by = week_year),
         pdur_1lag = lag(dur_precip, n = 1, order_by = week_year),
         pdur_2lag = lag(dur_precip, n = 2, order_by = week_year),
         pdur_3lag = lag(dur_precip, n = 3, order_by = week_year),
         pdur_4lag = lag(dur_precip, n = 4, order_by = week_year),
         pdur_6lag = lag(dur_precip, n = 6, order_by = week_year),
         wdir_1lag = lag(wind_direction, n = 1, order_by = week_year),
         wdir_2lag = lag(wind_direction, n = 2, order_by = week_year),
         wdir_3lag = lag(wind_direction, n = 3, order_by = week_year),
         wdir_4lag = lag(wind_direction, n = 4, order_by = week_year),
         wdir_6lag = lag(wind_direction, n = 6, order_by = week_year),
         wmean_1lag = lag(mean_wind_intensity, n = 1, order_by = week_year),
         wmean_2lag = lag(mean_wind_intensity, n = 2, order_by = week_year),
         wmean_3lag = lag(mean_wind_intensity, n = 3, order_by = week_year),
         wmean_4lag = lag(mean_wind_intensity, n = 4, order_by = week_year),
         wmean_6lag = lag(mean_wind_intensity, n = 6, order_by = week_year),
         wmax_1lag = lag(max_wind_intensity, n = 1, order_by = week_year),
         wmax_2lag = lag(max_wind_intensity, n = 2, order_by = week_year),
         wmax_3lag = lag(max_wind_intensity, n = 3, order_by = week_year),
         wmax_4lag = lag(max_wind_intensity, n = 4, order_by = week_year),
         wmax_6lag = lag(max_wind_intensity, n = 6, order_by = week_year)
         ) %>%
  ungroup() %>%
  as.data.table()

w_locs <- w_w %>%
  select(c(nome.da.estação, latitude, longitude, altitude)) %>%
  distinct(.keep_all = T)

# find the closest weather station to each trap - solely based on distance
locs = fread("ovitrap_locations_2012.csv") 
colnames(locs)[2:3] <- c("LAT_WGS", "LONG_WGS")

wsg_crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs"
locs_sf <- st_as_sf(locs, coords = c("LONG_WGS", "LAT_WGS"), crs = wsg_crs)
w_locs_sf <- st_as_sf(w_locs, coords = c("longitude","latitude"), crs = wsg_crs )
w_near <- st_distance(locs_sf, w_locs_sf)
# convert matrix to data frame and set column and row names
dist_matrix <- data.frame(w_near)
names(dist_matrix) <- w_locs$nome.da.estação
rownames(dist_matrix) <- locs$`Id Novo`

# find the nearest station and create new data frame

near <- dist_matrix %>% 
  mutate(ID=rownames(.)) %>% 
  pivot_longer(-ID, names_to = "station", values_to = "dist") %>%
  as.data.table()

# add altitude
near$ID <- as.integer(near$ID)
near = left_join(near, w_locs[,c(1,4)], by = join_by(station == nome.da.estação))

data2$week_year = paste(data2$week, data2$year, sep = "_")

# all trap-week combinations
all_weeks <- sort(unique(data2$week_year))
trap_weeks <- data2[, c(1,45,42,14)]
colnames(trap_weeks)[3] <- "trap_altitude"

trap_weeks$trap_altitude[is.na(trap_weeks$trap_altitude)] <- trap_weeks$Z_m[is.na(trap_weeks$trap_altitude)]

# stations active in each week
week_active <- w_w[, .(w_station = unique(nome.da.estação)), by = week_year]

# candidate trap-week-station combinations
candidates <- merge(trap_weeks, near, by.x = "Id.Novo", by.y = "ID", allow.cartesian = TRUE)
colnames(candidates)[5] <- "w_station"
candidates <- merge(
  week_active, candidates,
  by = c("week_year", "w_station"), all.x = TRUE
)

# keep only active traps
candidates <- candidates[!is.na(Id.Novo)]

# order by distance within trap-week
setorder(candidates, Id.Novo, week_year, dist)

# keep top 3 closest per trap-week, add rank
candidates <- candidates[, .SD[1:min(.N, 3)], by = .(Id.Novo, week_year)]
candidates[, rank := seq_len(.N), by = .(Id.Novo, week_year)]

# compute altitude difference
candidates[, alt_diff := abs(altitude - trap_altitude)]

nearest <- candidates[, .SD[which.min(alt_diff)], by = .(Id.Novo, week_year)]

# mark which station is the "nearest" (best by altitude among top 3)
candidates[, nearest := rank == which.min(alt_diff), by = .(Id.Novo, week_year)]
# ensure all trap-week combinations appear (NA when no station active)
near_stat <- merge(
  trap_weeks, candidates, 
  by = c("Id.Novo", "week_year", "trap_altitude"), all.x = TRUE)


# add to main dataframe 

data3 <- left_join(data2, near_stat[near_stat$nearest == TRUE, ], by = c("Id.Novo", "week_year"))

# add weather data
data3 <- left_join(data3, w_w, by = join_by("w_station" == "nome.da.estação", "week_year"), relationship = "many-to-many")
head(data3)

unique(data3$Id.Novo[is.na(data3$mean_temp_1.5m)])
write.csv(data3, "madeira_data_esa_weather_weekly.csv")

################################################################################
## - Socioeconomic data 
################################################################################
data3 <- fread("madeira_data_esa_weather_weekly.csv")
## - add census data -2011 and 2021

census <- fread("Madeira_parish_census_2011_2021.csv")

  data3 <- data3 %>%
    # join census data by parish
    left_join(census, by = c("Freguesia...Parish" = "Parish")) %>%
    # choose census year column depending on year in data3
    mutate(population_census = if_else(year > 2016, cyr_2021, cyr_2011)) %>%
    # optional: drop the extra columns if you don't need them
    select(-cyr_2011, -cyr_2021) %>%
    str_replace_all(population_census, " ", "")
  
# add income data to use as estimate of poverty
# yearly data are only available from 2015 at the municipality level
income <- fread("Madeira_taxincome_2015_2022.csv")

# split data to pre 2015 and post 2015 for ease of joining 
data_15 <- data3[data3$year <2016, ] %>%
  left_join(income[Year < 2016, c(1,5,7)], by = c("Concelho...Municipalitie" = "Municipality"))
  
data_22 <- data3[data3$year > 2015, ] %>%
  left_join(income[Year > 2015,c(1,2,5,7)], by = c("Concelho...Municipalitie" = "Municipality", "year" = "Year"))
# now recombine and save
data4 <- rbind(data_15, data_22) %>%
  select(-V1)

write.csv(data4, "madeira_data_esa_weather_weekly.csv")

## -- disposable income from WorldPop
disp <- rast("PRT_disp_inc_2015.tif")
plot(disp)
# crop to Madeira and porto santo onluy
disp <- crop(disp, mad_esa)

data4 <- fread("madeira_data_esa_weather_weekly.csv")
d_sf = st_as_sf(data4, coords = c("longitude.x", "latitude.x"), crs =sf::st_crs(wsg_crs))
plot(d_sf)

disp_income = terra::extract(disp, d_sf, fun = max)
data4$disp_income_2015 <- disp_income$PRT_disp_inc_2015

################################################################################
## -- ERA5 monthly anomallies
################################################################################

era5_files <- list.files("ERA5_monthlyannomallies/")
era5 <- rast(paste("ERA5_monthlyannomallies/", era5_files, sep = ""))



################################################################################
## -- Make a polygon outline of Madeira for modelling 
################################################################################

mad_esa = rast("Madeira_ESA_WorldCover_10m_2020.tif") 
# removing water from landcover map to create boundaries for islands
mad_esa_2 = mad_esa
mad_esa_2[mad_esa_2 == 80, ] <- NA
mad_esa_2[mad_esa_2 == 0, ] <- NA # also change 0 values to NA 
plot(mad_esa_2)
# set all non NA values to 1
mad_esa_2[!is.na(mad_esa_2), ] <- 1 
mad_poly = as.polygons(mad_esa_2)
plot(mad_poly)
# simplify outline for mesh creation and to remove holes where there's inland water 
mad_poly <- simplifyGeom(mad_poly, tolerance = 0.0009, 
                         preserveTopology = TRUE)
mad_poly <- fillHoles(mad_poly)
plot(mad_poly)
# simplify outline for mesh creation and to remove holes where there's inland water 
mad_poly <- simplifyGeom(mad_poly, tolerance = 0.0009, 
                         preserveTopology = TRUE)
mad_poly <- fillHoles(mad_poly)
plot(mad_poly)

# simplify outline for mesh creation and to remove holes where there's inland water 
mad_poly <- simplifyGeom(mad_poly, tolerance = 0.0009, 
                         preserveTopology = TRUE)
mad_poly <- fillHoles(mad_poly)
plot(mad_poly)

mad_poly = as(mad_poly, "Spatial")
library(rgdal)
writeOGR(mad_poly, ".", "Maderia_polygon", driver = "ESRI Shapefile" )

