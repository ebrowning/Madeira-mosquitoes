library(data.table); library(dplyr); library(INLA)
# load custom functions for plotting and processing INLA and 
source("00_custom_functions_4_inla_outputs.R")
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
                tree_500m = scale(proportion_treecover_500m ), 
                shrub_500m = scale(proportion_shrubland_500m ),
                grass_500m = scale(proportion_grassland_500m ), 
                bspveg_500m = scale(proportion_baresparseveg_500m),
                crop_500m = scale(proportion_cropland_500m ), 
                builtup_500m = scale(proportion_builtup_500m ), 
                water_500m = scale(proportion_permwaterbodies_500m),
                temp_anomaly = scale(era5_sfc_temperature_2m), 
                prec_anomaly = scale(era5_sfc_precipitation), 
                hum_anomaly = scale(era5_sfc_humidity),
                population_census = scale(population_census), 
                disp_income = scale(disp_income_2015),
                intercept = 1) %>%
  mutate_at(c(62:134), funs(c(scale(.)))) # scales weather vars

colnames(data)

# some further data sorting 
colnames(data)[c(8:10, 136:137)] <- c("agency_responsible", "location", "in_out", "gross_income", "income_tax")

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

data <- dplyr::filter(data, value < 1000) ## one erroneous? 

# # remove 2022 as only data for Jan-March and this may be affecting tail end of model
data <- data[data$year < 2022, ]
# remove extreme and likely erroneous value

# set week 53 as week 52 as few examples and officially only 52 weeks in a year
data$week.x[data$week.x == 53] <- 52




# covs to be used in models 
cv_covs <- c("value", "intercept", "Freguesia...Parish", "Id.Novo", "week.x", "month","year", "ID.year", 
              "location", "longitude.x", "latitude.x",
             "elevation", 
              "hmax_6lag", "psum_6lag",
             "max_temp_1.5m","max_wind_intensity",
              "tree_500m",             
              "grass_500m", "crop_500m", "builtup_500m",
             "temp_anomaly", "prec_anomaly", "disp_income"
)

ddf <- data[ ,..cv_covs]


# ================= initialises 5fold exclusion =================

# n.b. the below only runs once, the first time this script is run for a specified save_dir
# creates candidate models list, tracker and completed files to track and store model runs

# parameters for setup

n_reps = 12
k_folds = 5 

save_dir <- ("/Users/ellabrowning/Dropbox/MEWAR/Madiera/OOS/") ## change this (directory must already exist)
# create initialise files have not all created already
ll = list.files(save_dir)
 
  # unique IDs vec
  kfold_ids = c()
  
  # ----- create n_reps k-fold sets ------
  
  for(nn in 1:n_reps){
    
    unique_id_n = sample(1:10^6, 1)
    mod_names_n = paste(unique_id_n, "stempOOS", sep="_")
    kfold_ids = c(kfold_ids, unique_id_n)
    
    # partition dataset - parish and week-year combination for spatiotemporal

      folds = ddf %>% dplyr::select(Freguesia...Parish, year, week.x) %>% distinct()
      folds$kfold = kfold_func(folds, k = k_folds)
      write.csv(folds, paste(save_dir, unique_id_n, "_folds.csv", sep=""), row.names=FALSE)
      
      ddf = left_join(ddf, folds)
    
  } # end partitioning loop

## ------------------------------------------------------------------------- ###
  
  # ----- create models dataframe -----
  
  form_base = paste(c("y ~ -1 + intercept",
    "f(week.x, model = 'rw2', cyclic = TRUE, group = ID.year, control.group = list(model = 'rw2'))",
    "f(ID.year, model = 'rw2', hyper = hyper.rw_yr)",
    "f(s, model = spde2, group = s.group, control.group = list(model = 'ar1'))"),
  collapse = " + "
  )
  

# covariate names
  effect_names_500 = 
    c("f(max_temp_1.5m, model = 'rw2', hyper = hyper_rw_temp, scale.model=TRUE, constr=TRUE)",
      "f(hmax_6lag, model = 'rw2', hyper = hyper.rw_hum)",
      "f(prec_anomaly, model = 'ar1', hyper = hyper.ar_precan)",
      "elevation",
      "estat_pop_1km",
      "builtup_500m", 
      "crop_500m",
      "grass_500m",
      "tree_500m",
      "builtup_500m * elevation",
      "psum_6lag",
      "max_wind_intensity")
  
  fx = vector("list", length=length(effect_names_500)+1)
  
  # models
  m1 = paste(effect_names_500, collapse=" + ")
  
  name = "full"

  fx[1] = list(m1)
  
  # 500 m
  for(i in 1:length(effect_names_500)){
    fx[[ i + 1 ]] = paste(effect_names_500[ -i ], collapse=" + ")
  }
  # unlist and add into df
  fx = unlist(fx)
  # create data frame including formulae
  fx = data.frame(modid = 1:length(fx),
                  fx = fx,
                  candidate = c(name, "tmax", "hmax_6lag", "prec_anom", "elevation", "estat_pop_1km", "builtup",
                                "crop_500", "grass", "tree", "builtup*elevation", "psum_6lag", "max_wind_intensity"),
                  formula = paste(form_base, fx, sep=" + "))
  
  bs = data.frame(modid = "baseline", fx = "baseline", candidate="baseline", formula=form_base)
  fx = rbind(fx, bs)
  
  # repeat for each n_reps ID to create full set of models to fit, cross referenced to unique holdout set
  fx_mods = data.frame()
  for(ii in kfold_ids){
    fx_mods = rbind(fx_mods, 
                    fx %>% dplyr::mutate(unique_id = ii,
                                         model_filename = paste(ii, "stempOOS", "nb_model", modid, sep="_")))
  }
  fx_mods$model_identifier = 1:nrow(fx_mods)  
  
  # ------- save initialisation objects ---------
  
  # list of models to fit
  write.csv(fx_mods, paste(save_dir, "models_list_all.csv", sep=""), row.names=FALSE)
  
  # tracker (which models are currently running)
  tracker = data.frame(unique_id = "dummy", model_filename = "dummy", model_identifier = "dummy")
  write.csv(tracker, paste(save_dir, "models_tracker.csv", sep=""), row.names=FALSE)
  
  # completed
  completed = data.frame(unique_id = "dummy", model_filename = "dummy", model_identifier = "dummy", candidate = "dummy",
                         n_models_fitted = "dummy", mae_oos = "dummy", rmse_oos = "dummy")
  write.csv(completed, paste(save_dir, "models_completed.csv", sep=""), row.names=FALSE)
  
  
### ------------------------------------------------------------------------- ###

  # setup SPDE for models 
  my.init = NULL
  library(sf)
  mad_poly = read_sf("Maderia_polygon.shp")
  
  # fit boundary mesh
  max.edge = 0.01

  # without extension
  mesh_poly2 = fmesher::fm_mesh_2d_inla(boundary = mad_poly, 
                                        max.edge = max.edge,
                                        cutoff = max.edge/2, 
                                        offset = 2*max.edge)
  
  plot(mesh_poly2)
  points(ddf$longitude.x, ddf$latitude.x, col = "pink") #add points to visualise
  
  # Set range prioirs for the spde model
  prior.median.sd = 0.01; prior.median.range = 0.05
  
  spde2 = inla.spde2.pcmatern(mesh_poly2, prior.range = c(prior.median.range, 0.5),
                              prior.sigma = c(prior.median.sd, 0.1), constr = T)
  indexs <- inla.spde.make.index("s", n.spde = spde2$n.spde, n.group = 10)
  
  # time mesh 
  k <- 10
  mesh.t <- fmesher::fm_mesh_1d(seq(0 + 0.5 / k, 1 - 0.5 / k, length = k))
  
  locs = cbind(ddf$longitude.x, ddf$latitude.x)
  A = inla.spde.make.A(mesh = mesh_poly2, loc = locs, group = ddf$ID.year, group.mesh = mesh.t)
  
  hyper_rw_temp = list(theta = list(prior="pc.prec", param=c(3, 0.0001)))
  hhyper.rw_yr = list(theta = list(prior="pc.prec", param=c(20, 0.0000001)))
  hyper.ar_precan = list(rho = list(prior="pc.prec", param=c(0.5, 0.000000000001)))

### ------------------------------------------------------------------------- ###

# ---------------- chooses and fits model under k-fold CV -------------------- #

# check currently running and completed models
tracker = read.csv(paste(save_dir, "models_tracker.csv", sep="")) %>%
  dplyr::select(-c(X))
#tracker <- tracker[,-1]
completed = read.csv(paste(save_dir, "models_completed.csv", sep=""))

# list of models to fit
fx = read.csv(paste(save_dir, "models_list_all.csv", sep="")) %>%
  dplyr::filter(!model_identifier %in% c(tracker$model_identifier, completed$model_identifier)) 

  
  # select model and fit if > 0 models left in list
  if(nrow(fx) > 0){
    
    # select 1 model
    fx = fx %>% sample_n(1)
    
    # append selected model to tracker (removed after model completed)
    to_append = fx[ , c("unique_id", "model_filename", "model_identifier")]
    tracker = rbind(tracker, to_append)
    write.csv(tracker, file=paste(save_dir, "models_tracker.csv", sep=""), 
                col.names=FALSE, row.names=FALSE, sep=",")
    
  # add correct folds information to dd
folds = fread(paste(save_dir, fx$unique_id, "_folds.csv", sep=""))
ddf = left_join(ddf, folds)

  # fit each of the models in dataframe fx (by default 1)
for(i in 1:nrow(fx)){
  
  # formula
  fx_i = fx[ i, ]
  form_i = formula(as.vector(fx_i$formula))
  
  # fit INLA model nested in tryCatch
  e = simpleError("error fitting")
  
  # the folds to fit
  folds_seq = unique(ddf$kfold)[ order(unique(ddf$kfold)) ]
  
  # oos_results data frame
  oos_results = data.frame()
  
  # for each fold
  for(k in folds_seq){
    
    # data and set group k to NA
    dd_i <- ddf %>% dplyr::mutate(y = value) %>%
      dplyr::mutate(y = replace(y, kfold == k, NA))
    
    # Make stacks for and specify hyperparameters for the random walk
    stack_i <- inla.stack(tag = 'est', # name tag of the stack (e.g. here est = estimating)
                          data = list(y = dd_i$y),
                          A =list(A,1),
                          effects=list(s=indexs, # spatial
                                       data.frame(dd_i[ ,..cv_covs])))
    
    print("")
    print("==========================================")
    print(Sys.time())
    print("==========================================")
    print("")
    
    
    mod_i <- tryCatch(
      inla(form_i, family = "nbinomial", 
                      control.family = list(link = "log"),
                      data = inla.stack.data(stack_i), 
                      control.inla = list(int.strategy='eb', npoints = 21),
                      control.fixed=list(prec = 1),
                      control.mode=list(restart=T, theta=my.init),
                      control.predictor=list(A = inla.stack.A(stack_i), link = 1, compute=TRUE),
                      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),# set config to TRUE if want to do INLA::inla.posterior.sample()
                      inla.mode = 'experimental'),
      error = function(e) e
    )
    
    # try again if failed
    if(inherits(mod_i, "error")){
      mod_i = tryCatch(
        inla(form_i, family = "nbinomial", 
             control.family = list(link = "log"),
             data = inla.stack.data(stack_i), 
             control.inla = list(int.strategy='eb', npoints = 21),
             control.fixed=list(prec = 1),
             control.mode=list(restart=T, theta=my.init),
             control.predictor=list(A = inla.stack.A(stack_i), link = 1, compute=TRUE),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),# set config to TRUE if want to do INLA::inla.posterior.sample()
             inla.mode = 'experimental'),
        error = function(e) e
      )
    }
    
    # if failed again write timeout to result
    if(inherits(mod_i, "error")){
      
      ex = fx_i; ex$result = "error in fitting"
      ex_file_name = paste(k, "_", type, "_nb_err_", fx_i$model_identifier, ".csv", sep="")
      write.csv(ex, paste(save_dir, "errors/", ex_file_name, sep=""), row.names=FALSE)
      
      # otherwise calculate and save fit metrics and model
    } else{
      
      fm = fitMetricsINLA(mod_i, data=dd_i, modname=fx_i$fx, inla.mode="experimental")
      res_i = cbind(fx_i, fm)
      fm_file_name = paste(k, "_", "stempOOS", "_nb_fitmetrics_", fx_i$model_identifier, ".csv", sep="")
      write.csv(res_i, paste(save_dir, "fitmetrics/", fm_file_name, sep=""), row.names=FALSE)
      
      # save model
      #save(mod_i, file=paste(save_dir, "models/", paste(k, fx_i$model_filename, sep="_"), sep=""))
      
      # extract observed and fitted and save
      dd_o = ddf %>%
        dplyr::select(Id.Novo, year, week.x, value, kfold) %>%
        dplyr::mutate(oos = ifelse(kfold == k, TRUE, FALSE),
                      model = fx_i$candidate,
                      holdout_id = fx_i$unique_id,
                      model_identifier = fx_i$model_identifier,
                      holdout_type = "stempOOS",
                      mean = mod_i$summary.linear.predictor$mean[1:nrow(ddf)],
                      lower = mod_i$summary.linear.predictor$`0.025quant`[1:nrow(ddf)],
                      upper = mod_i$summary.linear.predictor$`0.975quant`[1:nrow(ddf)]) %>%
        dplyr::filter(oos == TRUE)
      write.csv(dd_o, paste(save_dir, "models/", paste(k, "output", fx_i$model_identifier, fx_i$model_filename, ".csv", sep="_"), sep=""), row.names=FALSE)
      
      # add to growing OOS results dataframe
      oos_results = rbind(oos_results, dd_o, fill = TRUE)
      
    } # end of model failed if-else 
      
  } # end of kfold loop
  
  
  # ============ final operations ==============
  
  print("Saving results")
  
  # create dataframe to add to "completed" csv
  completed_i = data.frame(unique_id = fx_i$unique_id, 
                           model_filename = fx_i$model_filename, 
                           model_identifier = fx_i$model_identifier, 
                           candidate = fx_i$candidate,
                           n_models_fitted = n_distinct(oos_results$kfold), 
                           mae_oos = NA, 
                           rmse_oos = NA)
  
  # calculate predicted and residual error (run through link function)
  oos_results$predicted = exp(oos_results$value + oos_results$mean)
  oos_results$resid = oos_results$predicted - oos_results$value
  
  # calculate summary stats
  completed_i$mae_oos = mean(abs(oos_results$resid), na.rm=TRUE)
  completed_i$rmse_oos = sqrt(mean(oos_results$resid^2, na.rm=TRUE))
  
  # remove from 
  #if(completed_i$n_models_fitted == k_folds){
  
  # append to completed
  write.table(completed_i, file=paste(save_dir, "models_completed.csv", sep=""), 
              append=TRUE, col.names=FALSE, row.names=FALSE, sep=",")
  
  # remove model from tracker
  read.csv(paste(save_dir, "models_tracker.csv", sep="")) %>%
    dplyr::filter(model_identifier != fx_i$model_identifier) %>%
    write.csv(paste(save_dir, "models_tracker.csv", sep=""))
  
  } # end of model fitting loop
  
} # end of if statement

