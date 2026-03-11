library(data.table); library(dplyr); library(INLA);library(terra);library(ggplot2)
#inla.upgrade(testing=TRUE)
source('funcs_plotting_inla.R')
################################################################################

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

# head(data)
# table(data$year, data$month) # few records at the beginning of 2012 - will start from March 
# # remove Jan and Feb 2012
# #data <- data[data$year_month > "2012-02-01", ]
# # remove 2022 as only data for Jan-March and this may be affecting tail end of model
data <- data[data$year < 2022, ]
# remove extreme and likely erroneous value

# find the number of years each location has been surveyed
#locs_yrs <- data |>
#  group_by(location) |>
#  summarise(n_distinct(year))
# how many locations have been surveyed fewer than 10 years
#table( locs_yrs$location, locs_yrs$`n_distinct(year)`)
# subset to locs which have been surveyed every year = 135
#locs_10 <- unique(locs_yrs$location[locs_yrs$`n_distinct(year)` > 1]) 
#data <- data[which(data$location %in% locs_10), ]
  
# set week 53 as week 52 as few examples and officially only 52 weeks in a year
data$week.x[data$week.x == 53] <- 52
  
# double check no covariance 
lm500m <- lm(value ~ tree_500m + grass_500m + crop_500m + builtup_500m + 
               max_temp_1.5m + hmax_6lag + psum_6lag + max_wind_intensity + population_census + gross_income + elevation,
             data = data)
car::vif(lm500m)
  
################################################################################

## -- Set up SPDE and stacks for INLA models 

################################################################################

## -- Setting up spatial components
library(sf)
mad_poly = read_sf("Maderia_polygon.shp")

# fit boundary mesh
max.edge = 0.01

# with extension
mesh_poly2 = fmesher::fm_mesh_2d_inla(boundary = mad_poly, 
                          max.edge = c(1, 5)*max.edge,
                          cutoff = max.edge/2, 
                          offset = c(2, 4)*max.edge)

# without extension
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
              "agency_responsible", "location", "in_out", "altitude", "elevation", "population_census", "gross_income",
              "max_humid",  "hmax_6lag",
               "psum_6lag",
              "min_temp_1.5m", 
              "max_temp_1.5m", "tmax_1lag", "tmax_2lag", "tmax_3lag", "tmax_4lag", "tmax_6lag",
              "max_wind_intensity", "wmax_1lag", "wmax_2lag", "wmax_3lag", "wmax_4lag", "wmax_6lag",
              "elevation", "time_since_ep",
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

hyper_rwt <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))
hyper.rw_temp = list(theta = list(prior="pc.prec", param=c(3, 0.0001)))
hyper.rw_hum = list(theta = list(prior="pc.prec", param=c(1, 0.001)))
fx1_500 <- y ~ -1 + intercept +
  f(week.x, model = "rw2", cyclic = TRUE, group = ID.year, control.group = list(model = "rw2")) +
  f(ID.year, model = "rw2") +
  f(s, model = spde2, group = s.group, control.group = list(model = "ar1")) +
  #f(as.factor(location), model = "iid", hyper = hyper.rw_temp) +
  f(max_temp_1.5m, model = "rw2", hyper = hyper.rw_temp) +
  f(hmax_6lag, model = "rw2", hyper = hyper.rw_hum) +
  elevation +
  population_census +
  gross_income +
  builtup_500m +
  crop_500m +
  grass_500m +
  tree_500m +
  builtup_500m * elevation +
  psum_6lag +
  max_wind_intensity 

m1_st_500 <- inla(fx1_500, family = "nbinomial", 
                     control.family = list(link = "log"),
                     data = inla.stack.data(stack_1), 
                     #control.inla = list(int.strategy='eb', npoints = 21),
                     control.fixed=list(prec = 1),
                     control.mode=list(restart=T, theta=my.init),
                     control.predictor=list(A = inla.stack.A(stack_1), link = 1, compute=TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                     inla.mode = 'experimental')

summary(m1_st_500)

###############################################################################

## Model posterior predictive distribution for each model 
#### function to sample the posterior, extract the posterior distribution for a random effect (e.g. week)
post_samp_func <- function(mod, rand_eff, mod_r_eff, s = 1000, n = n, obs_vals){
  samples = inla.posterior.sample(1000, mod)
  samples.m <- sapply(samples, function(x) x$latent[1:n])
  samples.m <- inla.link.invlog(samples.m)
  
  ## https://haakonbakkagit.github.io/btopic112.html#5_Draw_posterior_samples
  
  contents = mod$misc$configs$contents
  effect = rand_eff
  id.effect = which(contents$tag == effect)
  ind.effect = contents$start[id.effect]-1 + (1:contents$length[id.effect])
  samples[[1]]$latent[ind.effect, , drop = FALSE]
  samples.effect = lapply(samples, function(x) x$latent[ind.effect])
  s.eff = matrix(unlist(samples.effect), byrow = TRUE, nrow = length(samples.effect))
  colnames(s.eff) = rownames(samples[[1]]$latent)[ind.effect]
  
  uci_list <- apply(s.eff, 2, quantile,probs=c(0.975), na.rm = TRUE) # 95% credible intervals
  lci_list <- apply(s.eff, 2, quantile,probs=c(0.025))
  
  uci_unlist = matrix(unlist(uci_list), byrow = TRUE, nrow = length(uci_list))
  lci_unlist = matrix(unlist(lci_list), byrow = TRUE, nrow = length(lci_list))
  
  muci_list <- apply(s.eff, 2, quantile,probs=c(0.895), na.rm = TRUE) # 79% credible intervals
  mlci_list <- apply(s.eff, 2, quantile,probs=c(0.105))
  
  post_ref <- cbind(mod_r_eff, colMeans(s.eff),lci_unlist, uci_unlist, mlci_list,  muci_list)
  
  colnames(post_ref) = c("mod_mean", "post_mean", "post_95lci", "post_95uci", "post_79lci", "post_79uci")
  return(as.data.frame(post_ref))
  
  # pred.samples=matrix(NA,nrow=dim(samples.m)[1],ncol=s)
  #for (l in 1:dim(samples.m)[1]){
  # sample from binomial 
  # pred.samples[l,]=dbinom(obs_vals[l], 12, samples.m[l,])
  #}
}
m2_st_post1000 <- post_samp_func(mod = m1_st_500, rand_eff = "week.x", mod_r_eff = m1_st_500$summary.random$week.x$mean, s = 10000, n = nrow(data), obs_vals = data$value)

m2_st_post1000 <- m2_st_post1000 %>%
  mutate(Year = as.factor(c(rep(2012, 52),
                          rep(2013, 52),
                          rep(2014, 52),
                          rep(2015, 52),
                          rep(2016, 52),
                          rep(2017, 52),
                          rep(2018, 52),
                          rep(2019, 52),
                          rep(2020, 52),
                          rep(2021, 52),
                          rep(2022, 52))), 
         Week = rep(1:52, 11))

week_trend_p <- ggplot(m2_st_post1000[1:520, ], aes(x = Week, y = post_mean, colour = Year)) + 
  geom_line() +
  geom_ribbon(aes(ymin = post_95lci, ymax = post_95uci, fill = Year), alpha = 0.3) +
  scale_fill_viridis_d() + xlab("Week") + ylab("Posterior week trend") +
  scale_colour_viridis_d() + theme_linedraw() +
  guides(colour = guide_legend(position = "inside"),
         fill = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0.65, 0.15), 
        legend.direction = "horizontal")

#ggsave(filename = "madeira_weekyr_trends.png", dpi = 400, width = 20, height = 15, units = "cm")


m1_st_post1000_yr <- post_samp_func(mod = m1_st_500, rand_eff = "ID.year", mod_r_eff = m1_st_500$summary.random$ID.year$mean, s = 10000, n = nrow(data), obs_vals = data$value)

m1_st_post1000_yr$Year = (2012:2022)
ggplot(m1_st_post1000_yr[-11,], aes(x = Year, y = post_mean)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = post_95lci, ymax = post_95uci), alpha = 0.3) +
  geom_ribbon(aes(ymin = post_79lci, ymax = post_79uci), alpha = 0.5) +
  scale_fill_viridis_d() + xlab("Year") + ylab("Posterior year trend") +
  scale_colour_viridis_d() + theme_linedraw()
################################################################################
################################################################################

yr_trnd = m1_st_post1000_yr %>%
  mutate(mean_rescaled = (inla.link.invlog(post_mean))/inla.link.invlog(post_mean[2]),
         lower_rescaled = (inla.link.invlog(post_95lci))/inla.link.invlog(post_mean[2]),
         upper_rescaled = (inla.link.invlog(post_95uci))/inla.link.invlog(post_mean[2]),
         lower79_rescaled = (inla.link.invlog(post_79lci))/inla.link.invlog(post_mean[2]),
         upper79_rescaled = (inla.link.invlog(post_79uci))/inla.link.invlog(post_mean[2]))

year_change_p <- ggplot(yr_trnd[-11,], aes(x = Year, y = mean_rescaled)) + geom_point() + geom_line() +
  geom_ribbon(aes(ymin = lower_rescaled, ymax = upper_rescaled), alpha = 0.3) +
  geom_ribbon(aes(ymin = lower79_rescaled, ymax = upper79_rescaled), alpha = 0.4) +
   xlab("Year") + ylab("Relative change in annual abundance (baseline = 2013)") +
  theme_linedraw() + geom_hline( yintercept = 1, linetype = "dotted", colour = "grey20")
  #geom_hline( yintercept = yr_trnd$upper_rescaled[2], linetype = "dashed", colour = "grey20")


week_trnd <- m2_st_post1000 %>%
  group_by(Week) %>%
  mutate(mean_rescaled = (inla.link.invlog(post_mean))/inla.link.invlog(post_mean[2]),
         lower_rescaled = (inla.link.invlog(post_95lci))/inla.link.invlog(post_mean[2]),
         upper_rescaled = (inla.link.invlog(post_95uci))/inla.link.invlog(post_mean[2]))
  
week_change_p <- ggplot(week_trnd[1:520,  ], aes(x = Week, y = log10(mean_rescaled), colour = Year)) + geom_line(aes(colour = Year)) +
  facet_grid(Year~.) +
  xlim(1,52) +
  geom_ribbon(aes(ymin = log10(lower_rescaled), ymax = log10(upper_rescaled), fill = Year), alpha = 0.3) +
  scale_fill_viridis_d() + xlab("Week") + ylab("Log10 relative change in weekly abundance (baseline year = 2013)") +
  scale_colour_viridis_d() + theme_linedraw() + geom_hline( yintercept = 0, linetype = "dotted", colour = "grey20") +
  theme(legend.position = "none")


library(ggpubr)
ggarrange(year_change_p, week_trend_p, week_change_p, nrow = 3)

ggsave(filename = "madeira_temporal_trends.png", dpi = 400, height = 35, width = 20, units = "cm")
