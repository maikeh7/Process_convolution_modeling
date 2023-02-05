library(mgcv)
library(qmap)
library(hetGP)
library(dplyr)
library(MASS)


# NOTE: first three functions have been modified so that post mean x (e.g. random coefficient estimated values) is 
# just extracted from LME model fit rather than calculating by hand (which may cause issues)
# Also much faster this way! 
#########################################################################################################
# Function to adjust x's of PC model for WRF based on PC model run on observed data in training set
# value: modified x's of WRF PC model
# wrfmodel_test = LME model fit for ALL WRF data in historical period! (test years are filtered out later)
# obsmodel_train = LME model fit for observed data in training set ONLY
# post mean of x of WRF is adjusted using sigmas from PC model fit on observed training data
########################################################################################################
adjust_x_crossval = function(wrfmodel_test, obsmodel_train, c_star = 0.95){
  Model_listW = get_sigmas_K(wrfmodel_test) # get sigmas for PC model for all WRF data in historical period
  
  basis_dims = Model_listW$Nvec
  basis_names = Model_listW$kmats
  WRF_sigmas = Model_listW$sigmas
  
  Model_listO = get_sigmas_K(obsmodel_train) # same for observed
  
  Obs_sigmas = Model_listO$sigmas
  
  Sigma_ratios = Obs_sigmas / WRF_sigmas
  
  R_resid = Model_listO$sige / Model_listW$sige
  
  #post cond dist of x for future wrf model fit
  PostXWorig = post_cond_x(wrfmodel_test) 
  
  #now modify post mean of x and residuals
  PostXMod = PostXWorig
  corrX = PostXMod$meanx
  
  # adjust WRF x's by multiplying by sigma ratios. Use adjustment of c_star for 3.5 (which should be the first basis mat!)
  newX_list = list()
  for (i in 1:length(Sigma_ratios)){
    if (i == 1){
      my_new_x = corrX[1:basis_dims[1] ] * ((Sigma_ratios[1]) * c_star)
      newX_list[[i]] = my_new_x
    }else{
      my_new_x = corrX[(1 +  cumsum(basis_dims)[(i-1)]) : cumsum(basis_dims)[i] ] * Sigma_ratios[i]
      newX_list[[i]] = my_new_x
    }
  }
  
  x_adjust = unlist(newX_list)
  
  PostXMod$meanx = x_adjust # replace post mean of X with modified post mean of x
  
  NewWRFPreds = get_preds(PostXMod, Model_listW) # get new fitted values
  
  OldPreds = get_preds(PostXWorig, Model_listW) # get original fitted values
  
  resids = (NewWRFPreds$TMAX - OldPreds$mypreds) * R_resid # resids using non-adjusted post 
  
  #add resids to the preds made using adjusted post mean of x
  CorrWRF = resids + NewWRFPreds$mypreds 
  return(CorrWRF)
}

######################################################
# Get post mean of x via the 'easy' way
######################################################
post_cond_x = function(LME_model){
  
  meanx = as.numeric(unlist(LME_model$coefficients$random))
  
  return(list(meanx=meanx))
}


#######################################################################################
# get predicted values using post mean of x (PostX). Also need to add mean back in!
# Note that getting predicted values this way, the mean must be added back in
# if you use the predict() method on the LME model, you SHOULD NOT add mean back in.
#######################################################################################
get_preds = function(PostX, Model_list){
  Kpred = Model_list$K # basis matrix 
  meanx = PostX$meanx # post mean of x
  mu = Model_list$mu # intercept (mean)--must be added to Kx
  TMAX = Model_list$y # this is just the response variable in the model
  
  mypreds = as.vector(Kpred %*% meanx + mu) #these are the fitted values --> Kx + mu
  return(list(mypreds = mypreds, TMAX = TMAX))
  
}

get_sigmas_K = function(LMEMODEL){
  require(nlme)
  sigmas <- VarCorr(LMEMODEL)
  sigmas <- unique(sigmas[,2])
  sigmas <- as.numeric(sigmas[-1])
  sige = sigmas[length(sigmas)]
  sigmas = sigmas[-length(sigmas)]  
  
  #parse call to get the correct order of basis matrices...less prone to error than hard coding!
  test = as.character(LMEMODEL$call)[4]
  test1 = unlist(strsplit(test, c(",", "=")))
  test2 = sub('.*~', '', test1)
  kmats = sub(" .*", "", test2)
  
  
  Nvec = vector(length=length(kmats))
  for(i in 1:length(kmats)){
    temp = as.matrix(LMEMODEL$data[[kmats[i]]])
    if (ncol(temp) > 1){
      Nvec[i] = ncol(LMEMODEL$data[[kmats[i]]])
    }else{
      Nvec[i] = 1
    }
  }
  
  modeldat = LMEMODEL$data
  y = LMEMODEL$data[,1]
  modeldat = modeldat[ , kmats]
  K = as.matrix(modeldat)
  
  return(list(sigmas = sigmas, sige = sige, Nvec = Nvec, K= K, modeldat = modeldat, 
              mu = LMEMODEL$coefficients$fixed, y= y, kmats = kmats))
}

#######################################################################################
# wrfmodel_test = LME model for WRF test set only
# obsmodel_train = LME model for OBS train set only
# wrfdat = entire WRF dataset. 
# testyears = years in test set for ith fold
# value: all WRF data w/ adjusted TMAX (tmax_norm_VAR1_adjust) and 
# the same except only WRF in test set
#######################################################################################
process_WRF_adjustX = function(wrfdat, testyears, wrfmodel_test = wrfLMEmodel, 
                       obsmodel_train = obsLMEmodel){
  tmax_norm_VAR1_adjust = adjust_x_crossval(wrfmodel_test, obsmodel_train)
  
  wrfdat$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust
  wrftest = filter(wrfdat, YEAR %in% testyears)
  return(list(wrfdat = wrfdat, wrftest = wrftest))
}


#######################################################################################
#######################################################################################
# Function for detrending WRF:
# extract long term trend w/ spline using 5 knots (more knots and it will be too wiggly),
# extract seasonal mean and sd using HetGP model w/ Matern(5,2) covariance function, 
# get constant scalar using median absolute deviation. 
# wrfdat = entire WRF dataset. 
# testyears = years in test set for ith fold
# value: all WRF data w/ adjusted TMAX (tmax_norm_VAR1_adjust) and 
# the same except only WRF in test set
#######################################################################################
# Use HetGP model to extract seasonal means and SDs (we use HetGP in lieu of PC b/c it is just simpler *and*
# **if the Matern(5,2) covariance function is used**, you get pretty much exactly the same seasonal mean as from PC model
# first, a spline is used to extract long term trend. Then, HetGP model w/ Matern(5,2) cov function is fit to detrended data
# Next, WRF is detrended --> subtract, divide by SD, divide by constant scaling factor to ensure sd approx 1 for whole series
# value: processed WRF and processed WRF in test set only
#######################################################################################
process_WRF = function(wrfdat, testyears){
  
  wrfdat$dailyTimeStep = 1:nrow(wrfdat)
  
  wrftrend = gam(meanWRF ~ s(dailyTimeStep, k = 5), data = wrfdat)
  
  Long_term_trend_fut = wrftrend$fitted.values
  wrfdat$Long_term_trend_fut = Long_term_trend_fut
  wrfdat$detrended = wrfdat$meanWRF - Long_term_trend_fut
  
  preWRF = filter(wrfdat, timestep %in% 316:365)$detrended
  postWRF = filter(wrfdat, timestep %in% 1:50)$detrended
  
  num_years = length(unique(wrfdat$YEAR))
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  allWRF = c(preWRF, rep(wrfdat$detrended, 3), postWRF)
  
  # fit het GP using defaults
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Matern5_2")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf = rep(WRFpreds$mean[366:730], num_years)
  Var_scaling_factors_wrf = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], num_years)
  
  wrfdat$Seasonal_trend_wrf = seasonal_fits
  
  wrfdat$Var_scaling_factors_wrf = Var_scaling_factors_wrf
  
  wrfdat$resids = wrfdat$detrended - wrfdat$Seasonal_trend_wrf
  
  wrfdat$tmax_norm = wrfdat$resids / wrfdat$Var_scaling_factors_wrf
  Constant_scalar_wrf = mean(movMAD(wrfdat$tmax_norm, window = 100)) 
  wrfdat$tmax_norm_VAR1 = wrfdat$tmax_norm / Constant_scalar_wrf
  wrfdat$Constant_scalar_wrf = Constant_scalar_wrf
  
  wrftest = filter(wrfdat, YEAR %in% testyears)
  return(list(wrfdat = wrfdat, wrftest = wrftest))
  
}


#######################################################################################
# Function for detrending observed data:
# wrfdat: all WRF data for historical period
# obsdat: all observed data
# trainyears/testyears: years in training and test sets for ith fold
# Detrending: subtract WRF long term trend from observed, 
# extract seasonal mean and sd using HetGP model w/ Matern(5,2) covariance function, 
# get constant scalar using median absolute deviation. 
# value: observed data in train set, observed data in test set, and
# seasonal means, sd's, and constant scalar in a separate dataset
#######################################################################################

process_OBS = function(obsdat, wrfdat, trainyears, testyears, second_detrend = FALSE){
  if (second_detrend == TRUE){
    print("processing Obs--second_detrend is TRUE")
   #take WRF trend out of observed data (this the 'post processed' long term trend--hence the postfix 'new')
  obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut_new
  }else{
    print("processing Obs--second_detrend is FALSE")
    # if we are not adjusting x or do not want second detrending, we don't do further post processing--can use long term trend from original WRF
    obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut
  }
  # filter out observed test data
  obstest = filter(obsdat, YEAR %in% testyears)
  obstrain = filter(obsdat, YEAR %in% trainyears)
  
  # seasonal trend
  # set up data such that seasonal trend can be extracted
  num_years = length(unique(obstrain$YEAR))
  
  preOBS = filter(obstrain, timestep %in% 316:365)$detrended
  postOBS = filter(obstrain, timestep %in% 1:50)$detrended
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime, midtime, aftertime, postdays)
  
  allOBS = c(preOBS, rep(obstrain$detrended, 3), postOBS)
  
  # fit het GP using defaults
  het_Obs <- mleHetGP(alltime, allOBS, covtype = "Matern5_2")
  het_Obs_rebuild = rebuild(het_Obs, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  Obspreds = predict(x = pred_times, object = het_Obs_rebuild)
  
  Seasonal_trend_obs = rep(Obspreds$mean[366:730], num_years)
  Var_scaling_factors_obs = rep(sqrt(Obspreds$sd2 + Obspreds$nugs)[366:730], num_years)
  #plot(Seasonal_trend_obs[1:730])
  #plot(Var_scaling_factors_obs[1:730], col = 'red')
  
  obstrain$Seasonal_trend_obs = Seasonal_trend_obs
  
  obstrain$Var_scaling_factors_obs = Var_scaling_factors_obs
  
  obstrain$resids = obstrain$detrended - obstrain$Seasonal_trend_obs
  
  obstrain$tmax_norm = obstrain$resids / obstrain$Var_scaling_factors_obs
  Constant_scalar_obs = mean(movMAD(obstrain$tmax_norm, window = 100)) 
  obstrain$tmax_norm_VAR1 = obstrain$tmax_norm / Constant_scalar_obs
  obstrain$Constant_scalar_obs = Constant_scalar_obs
  
  obs_params = data.frame(Var_scaling_factors_obs = obstrain$Var_scaling_factors_obs[1:365],
                          Constant_scalar_obs = Constant_scalar_obs,
                          Seasonal_trend_obs = Seasonal_trend_obs[1:365],
                          timestep = 1:365)
  return(list(obstrain = obstrain, obstest = obstest, obs_params = obs_params))
  
}


# This is only used for EQM WITH TEMPORAL ADJUSTMENT **OR** if a second detrending is required
#use adjust_x = TRUE if you want to put the twice detrending back in
process_OBS2 = function(obsdat, wrfdat, trainyears, testyears, second_detrend = TRUE){
  if (second_detrend == TRUE){
    print("processing Obs--second_detrend is TRUE")
  # take WRF trend out of observed data (this the 'post processed' long term trend--hence the postfix 'new')
   obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut_new
  }else{
    print("processing Obs--second_detrend is FALSE")
    # if we are not adjusting x, we don't do further post processing--can use long term trend from original WRF
    obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut
  }
  # filter out observed test data
  obstest = filter(obsdat, YEAR %in% testyears)
  obstrain = filter(obsdat, YEAR %in% trainyears)
  
  # seasonal trend
  # set up data such that seasonal trend can be extracted
  num_years = length(unique(obstrain$YEAR))
  
  preOBS = filter(obstrain, timestep %in% 316:365)$detrended
  postOBS = filter(obstrain, timestep %in% 1:50)$detrended
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime, midtime, aftertime, postdays)
  
  allOBS = c(preOBS, rep(obstrain$detrended, 3), postOBS)
  
  # fit het GP using defaults
  het_Obs <- mleHetGP(alltime, allOBS, covtype = "Matern5_2")
  het_Obs_rebuild = rebuild(het_Obs, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  Obspreds = predict(x = pred_times, object = het_Obs_rebuild)
  
  Seasonal_trend_obs = rep(Obspreds$mean[366:730], num_years)
  Var_scaling_factors_obs = rep(sqrt(Obspreds$sd2 + Obspreds$nugs)[366:730], num_years)
  #plot(Seasonal_trend_obs[1:730])
  #plot(Var_scaling_factors_obs[1:730], col = 'red')
  
  obstrain$Seasonal_trend_obs = Seasonal_trend_obs
  
  obstrain$Var_scaling_factors_obs = Var_scaling_factors_obs
  
  obstrain$resids = obstrain$detrended - obstrain$Seasonal_trend_obs
  
  obstrain$tmax_norm = obstrain$resids / obstrain$Var_scaling_factors_obs
  Constant_scalar_obs = mean(movMAD(obstrain$tmax_norm, window = 100)) 
  obstrain$tmax_norm_VAR1 = obstrain$tmax_norm / Constant_scalar_obs
  obstrain$Constant_scalar_obs = Constant_scalar_obs
  
  obs_params = data.frame(Var_scaling_factors_obs = obstrain$Var_scaling_factors_obs[1:365],
                          Constant_scalar_obs = Constant_scalar_obs,
                          Seasonal_trend_obs = Seasonal_trend_obs[1:365],
                          timestep = 1:365)
  return(list(obstrain = obstrain, obstest = obstest, obs_params = obs_params))
  
}

movMAD = function(vec, window){
  goodlength = length(vec) - window
  resvec = vector(length = goodlength)
  for (i in 1:goodlength){
    resvec[i] = mad(vec[i:(i+(window-1))])
  }
  return(resvec)
}


#######################################################################################
# Optional function for a second detrending of data that may be necessary if adjustment
# of x's results in too much increase in variance
# wrfdat: all WRF data for historical period INCLUDING adjusted WRF: tmax_norm_VAR1
# value: dataset including new variables. Variables ending in 'new' denote
# new estimated values for 2nd detrending; tmax_norm_VAR1_adjust_prime is like tmax_norm_VAR1_adjust, 
# but it has been detrended again.
#######################################################################################
add_adjust_x_postProcess = function(wrfdat, testyears, wrfmodel_test = wrfLMEmodel, 
                        obsmodel_train = obsLMEmodel){
  wrfdat_new = wrfdat
  # adjust x's
  tmax_norm_VAR1_adjust = adjust_x_crossval(wrfmodel_test, obsmodel_train)
  wrfdat_new$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust
  
  # upscale to TMAX scale and then remove long term, seasonal mean and SD trends AGAIN
 
  up1 = (wrfdat_new$tmax_norm_VAR1_adjust * wrfdat_new$Constant_scalar_wrf[1]) *wrfdat_new$Var_scaling_factors_wrf
  up2 = up1 + wrfdat_new$Seasonal_trend_wrf
  up3 = up2 + wrfdat_new$Long_term_trend_fut
  
  # add newly adjusted WRF to wrfdat_new
  wrfdat_new$wrf_adjust=up3
  
  # add long term trend back in and take it out
  tempdf = data.frame(WRF_TMAX = up3, dailyTimeStep = wrfdat_new$dailyTimeStep)
  
  wrftrend = gam(WRF_TMAX ~ s(dailyTimeStep, k = 5), data = tempdf)

  Long_term_trend_fut_new = wrftrend$fitted.values
  wrfdat_new$Long_term_trend_fut_new = Long_term_trend_fut_new
  
  # subtract long term trend again (it is nearly identical to old long term trend, but we do it again)
  wrfdat_new$detrended_adjust = up3 - Long_term_trend_fut_new
  
  num_years= length(unique(wrfdat_new$YEAR))
  
  preWRF = filter(wrfdat_new, timestep %in% 316:365)$detrended_adjust
  postWRF = filter(wrfdat_new, timestep %in% 1:50)$detrended_adjust
  
  num_years = length(unique(wrfdat_new$YEAR))
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  allWRF = c(preWRF, rep(wrfdat_new$detrended_adjust, 3), postWRF)
  
  # fit het GP using defaults
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Matern5_2")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  # make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf_new = rep(WRFpreds$mean[366:730], num_years)
  
  Var_scaling_factors_wrf_new = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], num_years)
  
  wrfdat_new$Seasonal_trend_wrf_new = Seasonal_trend_wrf_new
  
  wrfdat_new$Var_scaling_factors_wrf_new = Var_scaling_factors_wrf_new
  
  wrfdat_new$resids = wrfdat_new$detrended_adjust - wrfdat_new$Seasonal_trend_wrf_new
  
  # now do the detrending process again
  wrfdat_new$tmax_norm_prime = wrfdat_new$resids / wrfdat_new$Var_scaling_factors_wrf_new
  Constant_scalar_wrf_new = mean(movMAD(wrfdat_new$tmax_norm_prime, window = 100)) 
  wrfdat_new$tmax_norm_VAR1_adjust_prime = wrfdat_new$tmax_norm_prime / Constant_scalar_wrf_new
  wrfdat_new$Constant_scalar_wrf_new = Constant_scalar_wrf_new
  

  return(wrfdat_new)
}

process_WRFtest = function(wrfdat){
  
  wrfdat$dailyTimeStep = 1:nrow(wrfdat)
  
  wrftrend = gam(meanWRF ~ s(dailyTimeStep, k = 5), data = wrfdat)
  
  Long_term_trend_fut = wrftrend$fitted.values
  wrfdat$Long_term_trend_fut = Long_term_trend_fut
  wrfdat$detrended = wrfdat$meanWRF - Long_term_trend_fut
  
  # seasonal trendN
  # set up data such that seasonal trend can be extracted
  num_years = length(unique(wrfdat$YEAR))
  
  preWRF = filter(wrfdat, timestep %in% 316:365)$detrended
  postWRF = filter(wrfdat, timestep %in% 1:50)$detrended
  
  num_years = length(unique(wrfdat$YEAR))
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  allWRF =  c(preWRF, rep(wrfdat$detrended, 3), postWRF)
  
  # fit het GP using defaults
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Matern5_2")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf = rep(WRFpreds$mean[366:730], num_years)
  Var_scaling_factors_wrf = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], num_years)
  
  wrfdat$Seasonal_trend_wrf = Seasonal_trend_wrf
  
  wrfdat$Var_scaling_factors_wrf = Var_scaling_factors_wrf
  
  wrfdat$resids = wrfdat$detrended - wrfdat$Seasonal_trend_wrf
  
  wrfdat$tmax_norm = wrfdat$resids / wrfdat$Var_scaling_factors_wrf
  Constant_scalar_wrf = mean(movMAD(wrfdat$tmax_norm, window = 100)) 
  wrfdat$tmax_norm_VAR1 = wrfdat$tmax_norm / Constant_scalar_wrf
  wrfdat$Constant_scalar_wrf = Constant_scalar_wrf
  
  #wrftest = filter(wrfdat, YEAR %in% testyears)
  
  
  return(wrfdat)
  
}
