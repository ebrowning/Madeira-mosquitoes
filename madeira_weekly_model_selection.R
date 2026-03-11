library(data.table); library(dplyr); library(INLA);library(terra)
#inla.upgrade(testing=TRUE)

## -- This script tests the weather variables and weekly lags in each variable 


################################################################################

# read in data 
data = fread("madeira_data_esa_weather_weekly.csv")

################################################################################
colnames(data)
# scale all covariates for model
data = data %>%
  dplyr::mutate(year_month = as.Date(paste(year, paste(month, "01", sep = "-"), sep = "-")),
                year_week = paste(year, week.x, sep = "_"),
                X_UTM_rescaled = (X_UTM - min(X_UTM))/1000,
                Y_UTM_rescaled = (Y_UTM - min(Y_UTM))/1000, 
                altitude = scale(Z_m),
                elevation = scale(ele_vals...2.),
                tree_50m = scale(proportion_treecover_50m ), 
                shrub_50m = scale(proportion_shrubland_50m ),
                grass_50m = scale(proportion_grassland_50m ),
                bspveg_50m = scale(proportion_baresparseveg_50m),
                crop_50m = scale(proportion_cropland_50m ), 
                builtup_50m = scale(proportion_builtup_50m ), 
                water_50m = scale(proportion_permwaterbodies_50m ), 
                tree_100m = scale(proportion_treecover_100m ), 
                shrub_100m = scale(proportion_shrubland_100m ),
                grass_100m = scale(proportion_grassland_100m ),
                bspveg_100m = scale(proportion_baresparseveg_100m),
                crop_100m = scale(proportion_cropland_100m ), 
                builtup_100m = scale(proportion_builtup_100m ), 
                water_100m = scale(proportion_permwaterbodies_100m ),
                tree_500m = scale(proportion_treecover_500m ), 
                shrub_500m = scale(proportion_shrubland_500m ),
                grass_500m = scale(proportion_grassland_500m ), 
                bspveg_500m = scale(proportion_baresparseveg_500m),
                crop_500m = scale(proportion_cropland_500m ), 
                builtup_500m = scale(proportion_builtup_500m ), 
                water_500m = scale(proportion_permwaterbodies_500m),
                intercept = 1) %>%
  mutate_at(c(62:134), funs(c(scale(.)))) # scales weather vars

colnames(data)

# some further data sorting 
colnames(data)[c(8:10, 136:137)] <- c("agency_responsible", "location", "in_out", "gross_income", "income_tax")

data$population_census <- scale(data$population_census)
data$gross_income <- gsub(" ", "", data$gross_income)
data$gross_income <- scale(as.numeric(data$gross_income))
data$income_tax <- gsub(" ", "", data$income_tax)
data$income_tax <- as.numeric(data$income_tax)

data$agency_responsible = as.factor(data$agency_responsible)
data$in_out = as.factor(data$in_out)
data$ID.year = as.numeric(data$year) - 2012 +1
data$month = as.factor(data$month)
data$month = factor(data$month, levels = c("1", "2", "3", "4", "5", "6", "7",
                                           "8", "9", "10", "11", "12"), 
                    ordered = TRUE)
data$ID.year_week <- as.numeric(as.factor(data$year_week))
# add column for dengue epidemic period - 2012 -2013
data$epidemic <- 0
data$epidemic[data$year < 2014] <- 1
# time since epidemic 
data$time_since_ep <- 0
data$time_since_ep[data$epidemic == 0] <- data$year[data$epidemic == 0] - 2013

head(data)
table(data$year, data$month) # few records at the beginning of 2012 

ggplot(data, aes(x = year_week, y = value, colour = year)) +geom_point() 

data <- dplyr::filter(data, value < 1000)

################################################################################

## -- correlations 

colnames(data)

# all weather - no lags
cor(data[,c(62:74,143)],  method = "spearman", use = "complete.obs")
# humidity
cor(data[,c(61:63, 91:95)],  method = "spearman", use = "complete.obs")
# temperature
cor(data[,c(62:64, 80:84)],  method = "spearman", use = "complete.obs")
# precipitation
cor(data[,c(68:70, 106:110)],  method = "spearman", use = "complete.obs")
# wind
cor(data[,c(69:71, 130:134)],  method = "spearman", use = "complete.obs")

# land cover and socioeconomic
cor(data[,c(29:35,135:137,143)],method = "spearman", use = "complete.obs")
cor(data[,c(36:42,135:137,143)],method = "spearman", use = "complete.obs")

# all 
cor(data[,c(29:35,135:137, 62:74,143)],  method = "spearman", use = "complete.obs")
cor(data[,c(36:42,135:137, 62:74,143)],  method = "spearman", use = "complete.obs")

################################################################################

### fitMetricsINLA: extract and report on metrics of fit - Gibb et al. 2024 Nat Comms

#' @param mod fitted INLA model
#' @param modname name to give model in resulting dataframe

fitMetricsINLA = function(mod, modname="mod"){
  
  fit = data.frame(modname = modname,
                   dic = mod$dic$dic, 
                   waic = mod$waic$waic,
                   waic_neffp = mod$waic$p.eff,
                   #mae = mean(dx$abs_err, na.rm=TRUE),
                   logscore = -mean(log(mod$cpo$cpo), na.rm=TRUE),
                   cpo_fail = sum(mod$cpo$failure == 1 & !is.na(mod$cpo$failure)))
  return(fit)
}

################################################################################

## -- Set up SPDE and stacks for INLA models 

################################################################################

## -- Setting up spatial components
library(sf)
mad_poly = read_sf("Maderia_polygon.shp")

# fit boundary mesh
max.edge = 0.02
mesh_poly2 = inla.mesh.2d(boundary = mad_poly, 
                          max.edge = c(1, 5)*max.edge,
                          cutoff = max.edge/2, # minimum length of triange edges
                          offset = c(2, 4)*max.edge)

mesh_poly2 = inla.mesh.2d(boundary = mad_poly, 
                          max.edge = max.edge,
                          cutoff = max.edge/2, 
                          offset = 2*max.edge)

plot(mesh_poly2)
points(data$longitude.x, data$latitude.x, col = "pink") #add points to visualise

# Set range prioirs for the spde model
prior.median.sd = 0.01; prior.median.range = 0.05

spde2 = inla.spde2.pcmatern(mesh_poly2, prior.range = c(prior.median.range, 0.5),
                            prior.sigma = c(prior.median.sd, 0.1), constr = T)
indexs <- inla.spde.make.index("s", n.spde = spde2$n.spde, n.group = 11)

# time mesh 
k <- 11
mesh.t <- inla.mesh.1d(seq(0 + 0.5 / k, 1 - 0.5 / k, length = k))

locs = cbind(data$longitude.x, data$latitude.x)
A = inla.spde.make.A(mesh = mesh_poly2, loc = locs, group = data$ID.year, group.mesh = mesh.t)

#check colnames again for covriates to be included in effects list
# add intercept column
data$intercept <- 1

# covs to be used in models 
mod_covs <- c("intercept", "Id.Novo", "month","year", "ID.year", "ID.year_week", "week.x", "agency_responsible",
              "agency_responsible", "location", "in_out", "altitude", "elevation", "population_census", "gross_income",
              "max_humid",  "hmax_1lag",  "hmax_2lag",  "hmax_3lag", "hmax_4lag", "hmax_6lag",
              "sum_precip_mm", "psum_1lag", "psum_2lag", "psum_3lag", "psum_4lag", "psum_6lag",
              "min_temp_1.5m", "tmin_1lag", "tmin_2lag", "tmin_3lag",  "tmin_4lag", "tmin_6lag", 
              "max_temp_1.5m", "tmax_1lag", "tmax_2lag", "tmax_3lag", "tmax_4lag", "tmax_6lag",
              "max_wind_intensity", "wmax_1lag", "wmax_2lag", "wmax_3lag", "wmax_4lag", "wmax_6lag",
               "elevation", "time_since_ep",
              "esa_10m" , "tree_50m", "shrub_50m",             
              "grass_50m", "bspveg_50m", "crop_50m", "builtup_50m",
              "water_50m", 
              "tree_100m", "shrub_100m",         
              "grass_100m", "bspveg_100m", "crop_100m", "builtup_100m",
              "water_100m",
              "tree_500m", "shrub_500m",           
              "grass_500m", "bspveg_500m", "crop_500m", "builtup_500m",
              "water_500m" 
)

# Make stacks for and specify hyperparameters for the random walk
stack_1 <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                      data = list(y = data$value),
                      A =list(A,1),
                      effects=list(s=indexs, # spatial
                                   data.frame(data[ ,..mod_covs])))
my.init = NULL

fx1 <- y ~ -1 + intercept + 
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) 

m1_st <- inla(fx1, family = "Poisson", 
              control.family = list(link = "log"),
              # offset = log(nvisits),
              data = inla.stack.data(stack_1), 
              control.fixed=list(prec = 1),
              control.mode=list(restart=T, theta=my.init),
              control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
              inla.mode = 'experimental')

summary(m1_st)

### -- Add land cover covariates - 500m
fx1_lc <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m + 
  builtup_500m * elevation 

m1_st_lc <- inla(fx1_lc, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_lc)

# test iterative removing each one

fx1_lc1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  crop_500m +
  grass_500m +
  tree_500m 

m1_st_lc1 <- inla(fx1_lc1, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_lc1)

fx1_lc2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  builtup_500m +
  grass_500m +
  tree_500m + 
  builtup_500m * elevation 

m1_st_lc2 <- inla(fx1_lc2, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_lc2)

fx1_lc3 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  builtup_500m +
  crop_500m +
  tree_500m + 
  water_500m +
  builtup_500m * elevation 

m1_st_lc3 <- inla(fx1_lc3, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_lc3)

fx1_lc4 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation 

m1_st_lc4 <- inla(fx1_lc4, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_lc4)


################################################################################
# model fits
fitMetricsINLA(m1_st_lc, modname = "full lcs");fitMetricsINLA(m1_st_lc1, modname = "-built");fitMetricsINLA(m1_st_lc2, modname = "-crop");fitMetricsINLA(m1_st_lc3, modname = "-grass");fitMetricsINLA(m1_st_lc4, modname = "-tree")

################################################################################

### -- Test land cover covariates - 100m
fx1_lc_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  builtup_100m +
  crop_100m +
  grass_100m +
  tree_100m + 
  water_100m +
  builtup_100m * elevation 

m1_st_lc_100 <- inla(fx1_lc_100, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_lc_100)

fx1_lc1_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + 
  elevation +
  crop_100m +
  grass_100m +
  tree_100m +
  water_100m

m1_st_lc1_100 <- inla(fx1_lc1_100, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_lc1_100)

fx1_lc2_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  builtup_100m +
  grass_100m +
  tree_100m + 
  water_100m +
  builtup_100m * elevation 

m1_st_lc2_100 <- inla(fx1_lc2_100, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_lc2_100)

fx1_lc3_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  builtup_100m +
  crop_100m +
  tree_100m + 
  water_100m +
  builtup_100m * elevation 

m1_st_lc3_100 <- inla(fx1_lc3_100, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_lc3_100)

fx1_lc4_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  builtup_100m +
  crop_100m +
  grass_100m +
  water_100m +
  builtup_100m * elevation 

m1_st_lc4_100 <- inla(fx1_lc4_100, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_lc4_100)

fx1_lc5_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  builtup_100m +
  crop_100m +
  grass_100m +
  tree_100m + 
  builtup_100m * elevation 

m1_st_lc5_100 <- inla(fx1_lc5_100, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_lc5_100)

################################################################################
# model fits
fitMetricsINLA(m1_st_lc_100, modname = "full lcs");fitMetricsINLA(m1_st_lc1_100, modname = "-built");fitMetricsINLA(m1_st_lc2_100, modname = "-crop");fitMetricsINLA(m1_st_lc3_100, modname = "-grass");fitMetricsINLA(m1_st_lc4_100, modname = "-tree");fitMetricsINLA(m1_st_lc2_100, modname = "-crop");fitMetricsINLA(m1_st_lc3_100, modname = "-grass");fitMetricsINLA(m1_st_lc5_100, modname = "-water")

################################################################################


### -- Add population census and income covariates
fx1_ip <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) + elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_ip <- inla(fx1_ip, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_ip)


################################################################################

### -- test adding max temperature as non-linear - each lag
fx1_t0 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(max_temp_1.5m), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t0 <- inla(fx1_t0, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_t0)

fx1_t1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_1lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t1 <- inla(fx1_t1, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t1)

fx1_t2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_2lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t2 <- inla(fx1_t2, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t2)

fx1_t4 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_4lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t4 <- inla(fx1_t4, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t4)

fx1_t6 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f((tmax_6lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t6 <- inla(fx1_t6, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t6)

################################################################################
# model fits
fitMetricsINLA(m1_st_t0, modname = "t0");fitMetricsINLA(m1_st_t1,modname = "t1");fitMetricsINLA(m1_st_t2, modname = "t2");fitMetricsINLA(m1_st_t4, modname = "t4");fitMetricsINLA(m1_st_t6, modname = "t6")

################################################################################

### -- test adding max temperature as linear - each lag
fx1_t0nl <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation  +
  max_temp_1.5m

m1_st_t0nl <- inla(fx1_t0nl, family = "nbinomial", 
                   control.family = list(link = "log"),
                   data = inla.stack.data(stack_1), 
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.fixed=list(prec = 1),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                   inla.mode = 'experimental')

summary(m1_st_t0nl)

fx1_t1nl <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation  +
  tmax_1lag

m1_st_t1nl <- inla(fx1_t1nl, family = "nbinomial", 
                   control.family = list(link = "log"),
                   data = inla.stack.data(stack_1), 
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.fixed=list(prec = 1),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                   inla.mode = 'experimental')

summary(m1_st_t1nl)

fx1_t2nl <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  tmax_2lag

m1_st_t2nl <- inla(fx1_t2nl, family = "nbinomial", 
                   control.family = list(link = "log"),
                   data = inla.stack.data(stack_1), 
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.fixed=list(prec = 1),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                   inla.mode = 'experimental')

summary(m1_st_t2nl)

fx1_t4nl <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  tmax_4lag

m1_st_t4nl <- inla(fx1_t4nl, family = "nbinomial", 
                   control.family = list(link = "log"),
                   data = inla.stack.data(stack_1), 
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.fixed=list(prec = 1),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                   inla.mode = 'experimental')

summary(m1_st_t4nl)

fx1_t6_nl <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  tmax_6lag

m1_st_t6nl <- inla(fx1_t6_nl, family = "nbinomial", 
                   control.family = list(link = "log"),
                   data = inla.stack.data(stack_1), 
                   control.inla = list(int.strategy='eb', npoints = 21),
                   control.fixed=list(prec = 1),
                   control.mode=list(restart=T, theta=my.init),
                   control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                   inla.mode = 'experimental')

summary(m1_st_t6nl)
################################################################################
# model fits
fitMetricsINLA(m1_st_t0nl, modname = "t0");fitMetricsINLA(m1_st_t1nl,modname = "t1");fitMetricsINLA(m1_st_t2nl, modname = "t2");fitMetricsINLA(m1_st_t4nl, modname = "t4");fitMetricsINLA(m1_st_t6nl, modname = "t6")

# lowest - t2 

# compare with non-linear best models
fitMetricsINLA(m1_st_t2nl, modname = "t2 nl");fitMetricsINLA(m1_st_t4, modname = "t2");fitMetricsINLA(m1_st_t6, modname = "t6")

################################################################################

### -- test adding min temperature as non-linear - each lag
fx1_t0min <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(min_temp_1.5m, model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t0min <- inla(fx1_t0min, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t0min)

fx1_t1min <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmin_1lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t1min <- inla(fx1_t1min, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t1min)

fx1_t2min <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmin_2lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t2min <- inla(fx1_t2min, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t2min)

fx1_t4min <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmin_4lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t4min <- inla(fx1_t4min, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t4min)

fx1_t6min <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmin_6lag), model = "rw2") +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation 

m1_st_t6min <- inla(fx1_t6min, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_t6min)
################################################################################
# model fits
fitMetricsINLA(m1_st_t0min, modname = "t0min");fitMetricsINLA(m1_st_t1min,modname = "t1min");fitMetricsINLA(m1_st_t2min, modname = "t2min");fitMetricsINLA(m1_st_t4min, modname = "t4min");fitMetricsINLA(m1_st_t6min, modname = "t6min")

################################################################################
### -- test adding total precipitation - each lag
fx1_p0 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  sum_precip_mm

m1_st_p0 <- inla(fx1_p0, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p0)

fx1_p1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_1lag

m1_st_p1 <- inla(fx1_p1, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p1)

fx1_p2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_2lag

m1_st_p2 <- inla(fx1_p2, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p2)

fx1_p3 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_3lag

m1_st_p3 <- inla(fx1_p3, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p3)


fx1_p4 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_4lag

m1_st_p4 <- inla(fx1_p4, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p4)

fx1_p6 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_6lag

m1_st_p6 <- inla(fx1_p6, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_p6)
################################################################################
# model fits
fitMetricsINLA(m1_st_p0, modname = "p0");fitMetricsINLA(m1_st_p1, modname = "p1");fitMetricsINLA(m1_st_p2, modname = "p2");fitMetricsINLA(m1_st_p3, modname = "p3");fitMetricsINLA(m1_st_p4, modname = "p4");fitMetricsINLA(m1_st_p6, modname = "p6")

#lowest - p6 but p4 <4 diff in waic dic and log score higher

################################################################################

################################################################################
################################################################################
### -- test adding max humidity - each lag
fx1_h0 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  max_humid

m1_st_h0 <- inla(fx1_h0, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h0)

fx1_h1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_1lag

m1_st_h1 <- inla(fx1_h1, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h1)

fx1_h2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_2lag

m1_st_h2 <- inla(fx1_h2, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h2)

fx1_h3 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_3lag

m1_st_h3 <- inla(fx1_h3, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h3)

fx1_h4 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_4lag

m1_st_h4 <- inla(fx1_h4, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h4)

fx1_h6 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_6lag

m1_st_h6 <- inla(fx1_h6, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h6)



################################################################################
fitMetricsINLA(m1_st_h0, modname = "h0");fitMetricsINLA(m1_st_h1, modname = "h1");fitMetricsINLA(m1_st_h2, modname = "h2");fitMetricsINLA(m1_st_h3, modname = "h3");fitMetricsINLA(m1_st_h4, modname = "h4");fitMetricsINLA(m1_st_h6, modname = "h6")

# lowest h0

################################################################################

fx1_h0_6 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  max_humid +
  hmax_2lag

m1_st_h0_6 <- inla(fx1_h0_6, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_h0_6)

fitMetricsINLA(m1_st_h2, modname = "h2") ; fitMetricsINLA(m1_st_h0_6, modname = "h0_6")

# not better 
################################################################################
### -- test adding max wind speed - each lag
fx1_w0 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  max_wind_intensity

m1_st_w0 <- inla(fx1_w0, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_w0)

fx1_w1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  wmax_1lag

m1_st_w1 <- inla(fx1_w1, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_w1)

fx1_w2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  wmax_2lag

m1_st_w2 <- inla(fx1_w2, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_w2)

fx1_w4 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  wmax_4lag

m1_st_w4 <- inla(fx1_w4, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_w4)

fx1_w6 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  elevation +
  population_census +
  gross_income + 
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  wmax_6lag

m1_st_w6 <- inla(fx1_w6, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_w6)

################################################################################
fitMetricsINLA(m1_st_w0, modname = "w0");fitMetricsINLA(m1_st_w1, modname = "w1");fitMetricsINLA(m1_st_w2, modname = "w2");fitMetricsINLA(m1_st_w4, modname = "w4");fitMetricsINLA(m1_st_w6, modname = "w6")

#lowest - w0

################################################################################
################################################################################

# model with best weather lag vars - non-linear temp
hyper_rwt <- list(theta = list(prior = "pc.prec", param = c(1, 0.1)))

fx1_bw1 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_6lag), model = "rw2", hyper = hyper_rwt) +
  elevation +
  population_census +
  gross_income +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_6lag +
  psum_6lag +
  max_wind_intensity 

m1_st_bw1 <- inla(fx1_bw1, family = "nbinomial", 
                 control.family = list(link = "log"),
                 data = inla.stack.data(stack_1), 
                 control.inla = list(int.strategy='eb', npoints = 21),
                 control.fixed=list(prec = 1),
                 control.mode=list(restart=T, theta=my.init),
                 control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                 inla.mode = 'experimental')

summary(m1_st_bw1)

## model with non-linear temperature
hyper_rwt <- list(theta = list(prior = "pc.prec", param = c(0.5, 0.1)))

fx1_bw1_nlt <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_6lag), model = "rw2", hyper = hyper_rwt) +
  elevation +
  population_census +
  gross_income +
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  hmax_6lag +
  psum_6lag +
  max_wind_intensity 

m1_st_bw1_nlt <- inla(fx1_bw1_nlt, family = "nbinomial", 
                  control.family = list(link = "log"),
                  data = inla.stack.data(stack_1), 
                  control.inla = list(int.strategy='eb', npoints = 21),
                  control.fixed=list(prec = 1),
                  control.mode=list(restart=T, theta=my.init),
                  control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                  inla.mode = 'experimental')

summary(m1_st_bw1_nlt)

fitMetricsINLA(m1_st_bw1, modname = "best lags lt");fitMetricsINLA(m1_st_bw1_nlt, modname = "best lags nlt");

################################################################################
################################################################################

# model with best weather lag vars and 100m buffer
hyper_rwt <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))

fx1_bw_100 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_6lag), model = "rw2", hyper = hyper_rwt) +
  elevation +
  population_census +
  gross_income +
  builtup_100m +
  crop_100m +
  grass_100m +
  builtup_100m * elevation +
  hmax_6lag +
  psum_6lag +
  max_wind_intensity


m1_st_bw_100 <- inla(fx1_bw_100, family = "nbinomial", 
                  control.family = list(link = "log"),
                  data = inla.stack.data(stack_1), 
                  control.inla = list(int.strategy='eb', npoints = 21),
                  control.fixed=list(prec = 1),
                  control.mode=list(restart=T, theta=my.init),
                  control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                  inla.mode = 'experimental')

summary(m1_st_bw_100)

## -- model fits 

fitMetricsINLA(m1_st, modname = "st base");fitMetricsINLA(m1_st_i, modname = "lcovs"); fitMetricsINLA(m1_st_ip, modname = "lcovs + pop"); 
fitMetricsINLA(m1_st_bw1, modname = "lcs + weather 500m"); fitMetricsINLA(m1_st_bw_100, modname = "lcs + weather 100m")

################################################################################

# model with best weather lag vars without wind intensity

fx1_bw1_2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_6lag), model = "rw2", hyper = hyper_rwt) +
  elevation +
  population_census +
  gross_income +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  hmax_6lag +
  psum_6lag 

m1_st_bw1_2 <- inla(fx1_bw1_2, family = "nbinomial", 
                  control.family = list(link = "log"),
                  data = inla.stack.data(stack_1), 
                  control.inla = list(int.strategy='eb', npoints = 21),
                  control.fixed=list(prec = 1),
                  control.mode=list(restart=T, theta=my.init),
                  control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                  inla.mode = 'experimental')

summary(m1_st_bw1_2)

################################################################################
################################################################################

# model with best weather lag vars and 100m buffer without wind intensity

fx1_bw_100_2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(tmax_6lag), model = "rw2", hyper = hyper_rwt) +
  elevation +
  population_census +
  gross_income +
  builtup_100m +
  crop_100m +
  grass_100m +
  builtup_100m * elevation +
  hmax_6lag +
  psum_6lag 

m1_st_bw_100_2 <- inla(fx1_bw_100_2, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_bw_100_2)

## -- model fits
fitMetricsINLA(m1_st_bw1, modname = "lcs + weather 500m"); fitMetricsINLA(m1_st_bw_100, modname = "lcs + weather 100m");fitMetricsINLA(m1_st_bw1_2, modname = "lcs + weather 500m w/o wind"); fitMetricsINLA(m1_st_bw_100_2, modname = "lcs + weather 100m w/o wind")

# with wind is better

################################################################################
################################################################################
################################################################################

## -- Now test spatiotemporal structure - should estimate changes in space at year or month grouping


################################################################################

## -- Set up SPDE and stacks for INLA models 

################################################################################

## -- Setting up spatial components
library(sf)
mad_poly = read_sf("Maderia_polygon.shp")

# fit boundary mesh
max.edge = 0.03
#mesh_poly2 = inla.mesh.2d(boundary = mad_poly, 
#                        cutoff = max.edge/2, 
#                         max.edge = c(1, 5)*max.edge,
#                       offset = c(2, 4)*max.edge)

mesh_poly2 = fmesher::fm_mesh_2d_inla(boundary = mad_poly, 
                                      max.edge = max.edge,
                                      cutoff = max.edge/2, 
                                      offset = 2*max.edge)

plot(mesh_poly2)
points(data$longitude.x, data$latitude.x, col = "pink") #add points to visualise

# Set range prioirs for the spde model
prior.median.sd = 0.01; prior.median.range = 0.05

spde2 = inla.spde2.pcmatern(mesh_poly2, prior.range = c(prior.median.range, 0.5),
                            prior.sigma = c(prior.median.sd, 0.1), constr = T)
indexs <- inla.spde.make.index("s", n.spde = spde2$n.spde, n.group = 10)

# time mesh 
k <- 10
mesh.t <- fmesher::fm_mesh_1d(seq(0 + 0.5 / k, 1 - 0.5 / k, length = k))

locs = cbind(data$longitude.x, data$latitude.x)
A = inla.spde.make.A(mesh = mesh_poly2, loc = locs, group = data$ID.year, group.mesh = mesh.t)

#check colnames again for covriates to be included in effects list
# add intercept column
data$intercept <- 1

# covs to be used in models 
mod_covs <- c("intercept", "Id.Novo", "week.x", "month","year", "ID.year", "agency_responsible",
              "agency_responsible", "location", "in_out", "altitude", "elevation", "population_census",
              "max_humid",  "hmax_1lag",  "hmax_2lag",  "hmax_3lag", "hmax_4lag", "hmax_6lag",
              "sum_precip_mm", "psum_1lag", "psum_2lag", "psum_3lag", "psum_4lag", "psum_6lag",
              "min_temp_1.5m", "tmin_1lag", "tmin_2lag", "tmin_3lag",  "tmin_4lag", "tmin_6lag", 
              "max_temp_1.5m", "tmax_1lag", "tmax_2lag", "tmax_3lag", "tmax_4lag", "tmax_6lag",
              "max_wind_intensity", "wmax_1lag", "wmax_2lag", "wmax_3lag", "wmax_4lag", "wmax_6lag",
              "elevation", "time_since_ep",
              "esa_10m" , "tree_50m", "shrub_50m",             
              "grass_50m", "bspveg_50m", "crop_50m", "builtup_50m",
              "water_50m", 
              "tree_100m", "shrub_100m",             
              "grass_100m", "bspveg_100m", "crop_100m", "builtup_100m",
              "water_100m",
              "tree_500m", "shrub_500m",             
              "grass_500m", "bspveg_500m", "crop_500m", "builtup_500m",
              "water_500m" 
)

# Make stacks for and specify hyperparameters for the random walk
stack_1 <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                      data = list(y = data$value),
                      A =list(A,1),
                      effects=list(s=indexs, # spatial
                                   data.frame(data[ ,..mod_covs])))
my.init = NULL

fx1 <- y ~ -1 + intercept +
  #f(month, model = "rw1", cyclic = TRUE) +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "rw2")) 

m1_st <- inla(fx1, family = "Poisson", 
              control.family = list(link = "log"),
              #offset = log(nvisits),
              data = inla.stack.data(stack_1), 
              control.fixed=list(prec = 1),
              control.mode=list(restart=T, theta=my.init),
              control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
              inla.mode = 'experimental')

summary(m1_st)

fitMetricsINLA(m1_st, modname = "st")

### -- Add covariates

hyper.rw_temp = list(theta = list(prior="pc.prec", param=c(1, 0.001)))

fx1_i <- y ~ -1 + intercept +
  # f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  f(inla.group(max_temp_1.5m), model = "rw2", hyper = hyper.rw_temp) +
  elevation +
  population_census +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  max_humid + 
  hmax_6lag +
  psum_2lag 

m1_st_i <- inla(fx1_i, family = "nbinomial", 
                control.family = list(link = "log"),
                data = inla.stack.data(stack_1), 
                control.inla = list(int.strategy='eb', npoints = 21),
                control.fixed=list(prec = 1),
                control.mode=list(restart=T, theta=my.init),
                control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                inla.mode = 'experimental')

summary(m1_st_i)

fitMetricsINLA(m1_st_i, modname = "st_i")
################################################################################

## -- now use month_year as the time dimension instead of year

data <- data %>% 
  arrange(year_month) %>%
  mutate(month_year_num = as.numeric(as.factor(year_month))) 

n_months <- length(unique(data$month_year_num))

# Set range prioirs for the spde model
prior.median.sd = 0.01; prior.median.range = 0.05

spde2 = inla.spde2.pcmatern(mesh_poly2, prior.range = c(prior.median.range, 0.5),
                            prior.sigma = c(prior.median.sd, 0.1), constr = T)
indexs <- inla.spde.make.index("s", n.spde = spde2$n.spde, n.group = n_months)

# time mesh using number of months since start of survey up to 2021
k <- n_months
mesh.t <- inla.mesh.1d(seq(0 + 0.5 / k, 1 - 0.5 / k, length = k))

A = inla.spde.make.A(mesh = mesh_poly2, loc = locs, group = data$month_year_num, group.mesh = mesh.t)

# Make stacks for and specify hyperparameters for the random walk
stack_2 <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                      data = list(y = data$value),
                      A =list(A,1),
                      effects=list(s=indexs, # spatial
                                   data.frame(data[ ,..mod_covs])))
my.init = NULL

fx2 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "rw2")) 

m1_st2 <- inla(fx2, family = "nbinomial", 
               control.family = list(link = "log"),
               data = inla.stack.data(stack_2), 
               control.inla = list(int.strategy='eb', npoints = 21),
               control.fixed=list(prec = 1),
               control.mode=list(restart=T, theta=my.init),
               control.predictor=list(A = inla.stack.A(stack_2), link = 1, compute=TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), 
               inla.mode = 'experimental',
               verbose = TRUE)

summary(m1_st2)

#save(m1_st2, file ="weekly_spatiotemporal_m1_st2.RData")

#load("test_spatiotemporal_m0.RData")

################################################################################

## -- full dataset with all covariates

fx2_i <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(s, model = spde2, group = s.group, control.group = list(model = "rw2")) +
  f(inla.group(max_temp_1.5m), model = "rw2") +
  elevation +
  population_census +
  builtup_500m +
  crop_500m +
  grass_500m +
  builtup_500m * elevation +
  max_humid + 
  hmax_6lag +
  psum_2lag 

m2_st <- inla(fx2_i, family = "nbinomial", 
              control.family = list(link = "log"),
              data = inla.stack.data(stack_2), 
              control.inla = list(int.strategy='eb', npoints = 21),
              control.fixed=list(prec = 1),
              control.mode=list(restart=T, theta=my.init),
              control.predictor=list(A = inla.stack.A(stack_2), link = 1, compute=TRUE),
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), 
              inla.mode = 'experimental',
              verbose = TRUE)

summary(m2_st)

