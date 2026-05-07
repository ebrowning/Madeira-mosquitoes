# Functions and packages for plotting INLA model outputs

# Visualising the covariates
cov_plot <- function(mod_sum_fixed,  ylims, y_text_pos, 
                     cov_labels = c(rownames(mod_sum_fixed)), plot_title, hline = 0){
  
  plot((mod_sum_fixed$mean),ylim = ylims, ylab="Linear predictor", xlab = ' ',
       pch=19, col="black", main = plot_title, cex.lab=1, cex = 1, axes= FALSE)
  axis(side = 1, at = 1:length(rownames(mod_sum_fixed)),labels = cov_labels, cex = 0.5, las = 2, pos = y_text_pos)
  axis(side = 2, at = ylims, las = 0)
  # add points for the credible intervals
  points((mod_sum_fixed$`0.025quant`),col="blue",pch=18, cex = 1)
  points((mod_sum_fixed$`0.975`),col="blue",pch=18, cex = 1)
  # add lables for each covariate
 # text(x = c(1:nrow(mod_sum_fixed)),y = 
  #       y_text_pos , labels = cov_labels,  font = 2, cex = 0.7, srt = 90)
  # add dotted lines to join CIs
  for(i in 1:nrow(mod_sum_fixed)){
    segments(i, (mod_sum_fixed$`0.025quant`[i]),i,(mod_sum_fixed$`0.975`[i]),lwd=2, lty = 3, col = "blue")
    
  }
  abline(h=hline,lty=2,lwd=1,col="gray50")
}

# Plotting the temporal trend - scaled so that the second year = 
## Plot the yearly trend
yr_trend_plot <- function(mod_yr_sum, ylims, plot_title, ylab = 'Probability of observation (baseline =1999)', xlab = 'Year'){
yr_trnd = mod_yr_sum
yr_trnd$mean_rescaled = (inla.link.invlog(yr_trnd$mean))/inla.link.invlog(yr_trnd$mean[2])*100
yr_trnd$lower_rescaled = (inla.link.invlog(yr_trnd$`0.025quant`))/inla.link.invlog(yr_trnd$mean[2])*100
yr_trnd$upper_rescaled = (inla.link.invlog(yr_trnd$`0.975quant`))/inla.link.invlog(yr_trnd$mean[2])*100
plot(yr_trnd$ID, yr_trnd$mean_rescaled, type="l", ylim = ylims, xlab = xlab, 
     ylab = ylab, main = plot_title)
lines(yr_trnd$ID, yr_trnd$lower_rescaled, lty=2)
lines(yr_trnd$ID, yr_trnd$upper_rescaled, lty=2)
abline(h=100, lty=3)
}

# Plotting the posterior spatial field 
# The spatial field 
local.plot.field = function(field, mesh,  xlim = c(80, 670), ylim = c(10, 1100), zlim = c(-6, 6), colours = magma(16),...){
  stopifnot(length(field) == mesh$n)
  # - error when using the wrong mesh
  proj = inla.mesh.projector(mesh, xlim = xlim, 
                             ylim = ylim, dims=c(300, 300))
  # - Can project from the mesh onto a 300x300 plotting grid 
  field.proj = inla.mesh.project(proj, field)
  # - Do the projection
  image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
             xlim = xlim, ylim = ylim, zlim = zlim, col = colours, ...)  
}


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


# ------------------- 2. INLA model fitting function --------------------

### fitINLAModel: fit and return INLA model with specified formula and family
### From Gibb et al. 2024 Nature Communications

#' @param formx inla formula object; i.e. created using formula(y ~ x + f())
#' @param data model data frame
#' @param family likelihood
#' @param config boolean (default FALSE): set config in compute to TRUE for inla.posterior.sample()
#' @param verbose verbose reporting on or off? default FALSE
#' @param return.marginals specify whether model should save and return marginals

fitINLAModel = function(formx, data, family, config=FALSE, verbose=FALSE, return.marginals=FALSE){
  
  return(
    
    INLA::inla(
      
      # formula, data and model family
      formx,
      data = data,
      family = family,
      
      # fixed effects calibration
      control.fixed = list(mean.intercept=0, 
                           prec.intercept=1, # precision 1
                           mean=0, 
                           prec=1), # weakly regularising on fixed effects (sd of 1)
      
      # save predicted values on response scale
      control.predictor = list(compute=FALSE, 
                               link=1),
      
      # items to compute
      control.compute = list(cpo=TRUE, 
                             waic=TRUE, 
                             dic=TRUE, 
                             config=config, # set config to TRUE if want to do INLA::inla.posterior.sample()
                             return.marginals=return.marginals), # do not return marginals unless specified (saves memory)
      
      # configure inla approx
      control.inla = list(strategy='adaptive', # adaptive gaussian
                          cmin=0), # fixing Q factorisation issue https://groups.google.com/g/r-inla-discussion-group/c/hDboQsJ1Mls)
      
      # verbose?
      inla.mode = "experimental", # new version of INLA algorithm (requires R 4.1 and INLA testing version)
      verbose = verbose
    )
  )
  
  # # quick viz of prior sd distribution with precision
  # precx = 1
  # hist(rnorm(100000, 0, sd=sqrt(1/precx)), 100)
}


###################################################################################################

# kfold function from dismo package

kfold_func = function (x, k = 5, by = NULL) 
{
  singlefold <- function(obs, k) {
    if (k == 1) {
      return(rep(1, obs))
    }
    else {
      i <- obs/k
      if (i < 1) {
        stop("insufficient records:", obs, ", with k=", 
             k)
      }
      i <- round(c(0, i * 1:(k - 1), obs))
      times = i[-1] - i[-length(i)]
      group <- c()
      for (j in 1:(length(times))) {
        group <- c(group, rep(j, times = times[j]))
      }
      r <- order(runif(obs))
      return(group[r])
    }
  }
  if (is.vector(x)) {
    if (length(x) == 1) {
      if (x > 1) {
        x <- 1:x
      }
    }
    obs <- length(x)
  }
  else if (inherits(x, "Spatial")) {
    if (inherits(x, "SpatialPoints")) {
      obs <- nrow(coordinates(x))
    }
    else {
      obs <- nrow(x@data)
    }
  }
  else {
    obs <- nrow(x)
  }
  if (is.null(by)) {
    return(singlefold(obs, k))
  }
  by = as.vector(as.matrix(by))
  if (length(by) != obs) {
    stop("by should be a vector with the same number of records as x")
  }
  un <- unique(by)
  group <- vector(length = obs)
  for (u in un) {
    i = which(by == u)
    kk = min(length(i), k)
    if (kk < k) 
      warning("lowered k for by group: ", u, "  because the number of observations was  ", 
              length(i))
    group[i] <- singlefold(length(i), kk)
  }
  return(group)
}


### extractRandomINLA: extract random effect and rename columns

# if effect is grouped/replicated by a factor, automatically assign each subgroup to its grouping factor (labelled 1:n) 
# if BYM model, further partition into u and v components

#' @param summary_random points to model$summary.random$effect_of_interest
#' @param effect_name name to assign to fitted effect in dataframe (can be anything)
#' @param model_is_bym boolean; to specify if model is joint Besag-York-Mollie
#' @param transform specify whether to exponentiate coefficients (i.e. back transform to relative risk)
#' 
extractRandomINLA = function(summary_random, effect_name, model_is_bym=FALSE, transform=FALSE){
  
  # extract model effect
  rf = summary_random %>%
    dplyr::rename("value"=1, "lower"=4, "median"=5, "upper"=6)
  
  # label by grouping factor (if not replicated, group is 1 for all observations)
  rf$group = rep(1:as.vector(table(rf$value)[1]), each=n_distinct(rf$value))
  
  # partition BYM into u and v components
  if(model_is_bym){
    rf$component = rep(c("uv_joint", "u_besag"), each=n_distinct(rf$value)/2)
    rf$value = rep(1:(n_distinct(rf$value)/2), n_distinct(rf$group)*2)
  }
  
  # back transform if specified
  if(transform == TRUE){
    rf[ , 2:7 ] = exp(rf[ , 2:7])
  }
  
  # name and return
  rf$effect = effect_name
  return(rf)
}


### extractFixedINLA: extract fixed effects and rename columns

extractFixedINLA = function(model, model_name="mod", transform=FALSE){
  ff = model$summary.fixed
  ff$param = row.names(ff)
  ff$param[ ff$param == "(Intercept)" ] = "Intercept"
  names(ff)[3:5] = c("lower", "median", "upper")
  if(transform == TRUE){
    ff[ 1:5 ] = exp(ff[ 1:5 ])
  }
  ff
}


