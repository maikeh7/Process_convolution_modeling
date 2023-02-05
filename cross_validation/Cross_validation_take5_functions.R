library(mgcv)
library(qmap)
library(hetGP)
library(dplyr)
library(MASS)


process_WRFtest = function(wrfdat){
 
  wrfdat$dailyTimeStep = 1:nrow(wrfdat)
  
  wrftrend = gam(meanWRF ~ s(dailyTimeStep, k = 5), data = wrfdat)
  
  Long_term_trend_fut = wrftrend$fitted.values
  wrfdat$Long_term_trend_fut = Long_term_trend_fut
  wrfdat$detrended = wrfdat$meanWRF - Long_term_trend_fut
  
  # seasonal trendN
  # set up data such that seasonal trend can be extracted
  Num_years = length(unique(wrfdat$YEAR))
  
  aftertime = rep(731:1095, Num_years)
  midtime = rep(366:730, Num_years)
  beforetime = rep(1:365, Num_years)
  alltime = c(beforetime,midtime, aftertime)
  
  allWRF = rep(wrfdat$detrended, 3)
  
  # fit het GP using defaults
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Gaussian")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf = rep(WRFpreds$mean[366:730], Num_years)
  Var_scaling_factors_wrf = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], Num_years)
 
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
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Gaussian")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf = rep(WRFpreds$mean[366:730], Num_years)
  Var_scaling_factors_wrf = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], Num_years)
  
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

process_OBS = function(obsdat, wrfdat, trainyears, testyears, adjust_x=FALSE){
  if (adjust_x == TRUE){
    print("processing Obs--adjust_x is TRUE")
  # take WRF trend out of observed data (this the 'post processed' long term trend--hence the postfix 'new')
  obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut_new
  }else{
    print("processing Obs--adjust_x is FALSE")
    # if we are not adjusting x, we don't do further post processing--can use long term trend from original WRF
    obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut
  }
  # filter out observed test data
  obstest = filter(obsdat, YEAR %in% testyears)
  obstrain = filter(obsdat, YEAR %in% trainyears)
  
  
  
  # seasonal trend
  # set up data such that seasonal trend can be extracted
  Num_years = length(unique(obstrain$YEAR))
  
  preOBS = filter(obstrain, timestep %in% 316:365)$detrended
  postOBS = filter(obstrain, timestep %in% 1:50)$detrended
  
  mydays = -rev(0:49)
  predays = rep(mydays, Num_years)
  postdays = rep(1096:1145, Num_years)
  aftertime = rep(731:1095, Num_years)
  midtime = rep(366:730, Num_years)
  beforetime = rep(1:365, Num_years)
  
  alltime = c(predays, beforetime, midtime, aftertime, postdays)
  
  allOBS = c(preOBS, rep(obstrain$detrended, 3), postOBS)
  
  # fit het GP using defaults
  het_Obs <- mleHetGP(alltime, allOBS, covtype = "Gaussian")
  het_Obs_rebuild = rebuild(het_Obs, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  Obspreds = predict(x = pred_times, object = het_Obs_rebuild)
  
  Seasonal_trend_obs = rep(Obspreds$mean[366:730], Num_years)
  Var_scaling_factors_obs = rep(sqrt(Obspreds$sd2 + Obspreds$nugs)[366:730], Num_years)
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

# this function adjusts

add_adjust_x_postProcess = function(wrfdat, testyears, wrfmodel_test = wrfLMEmodel, 
                        obsmodel_train = obsLMEmodel, adjust_x=TRUE){
  wrfdat_new = wrfdat
  # adjust x's
  tmax_norm_VAR1_adjust = adjust_x_crossval_alt(wrfmodel_test, obsmodel_train)
  wrfdat_new$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust
  
  # upscale to TMAX scale and then remove long term, seasonal mean and SD trends again
 
  up1 = (wrfdat_new$tmax_norm_VAR1_adjust * wrfdat_new$Constant_scalar_wrf[1]) *wrfdat_new$Var_scaling_factors_wrf
  up2 = up1 + wrfdat_new$Seasonal_trend_wrf
  up3 = up2 + wrfdat_new$Long_term_trend_fut
  
  # add long term trend back in and take it out
  tempdf = data.frame(WRF_TMAX = up3, dailyTimeStep = wrfdat_new$dailyTimeStep)
  
  wrftrend = gam(WRF_TMAX ~ s(dailyTimeStep, k = 5), data = tempdf)

  Long_term_trend_fut_new = wrftrend$fitted.values
  wrfdat_new$Long_term_trend_fut_new = Long_term_trend_fut_new
  
  # subtract long term trend again (it is nearly identical to old long term trend, but we do it again)
  wrfdat_new$detrended_adjust = up3 - Long_term_trend_fut_new
  
  Num_years= length(unique(wrfdat_new$YEAR))
  
  aftertime = rep(731:1095, Num_years)
  midtime = rep(366:730, Num_years)
  beforetime = rep(1:365, Num_years)
  alltime = c(beforetime,midtime, aftertime)
  
  allWRF = rep(wrfdat_new$detrended_adjust, 3)
  
  # fit het GP using defaults
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Gaussian")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  
  
  Seasonal_trend_wrf_new = rep(WRFpreds$mean[366:730], Num_years)
  
  Var_scaling_factors_wrf_new = rep(sqrt(WRFpreds$sd2 + WRFpreds$nugs)[366:730], Num_years)
  
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


adjust_x_crossval = function(wrfmodel_test, obsmodel_train){
  
  Model_listW = get_sigmas_K(wrfmodel_test)
  Model_listO = get_sigmas_K(obsmodel_train)
  basis_dims = Model_listW$Nvec
  
  R_3pt5 =  Model_listO$sigmas[1] / Model_listW$sigmas[1] *.96
  R_180 = Model_listO$sigmas[2] / Model_listW$sigmas[2] 
  
  R_resid = Model_listO$sige / Model_listW$sige  
  
  #post cond dist of x for future wrf model fit
  PostXWorig = post_cond_x(Model_listW) 
  
  #now modify post mean of x and residuals
  PostXMod = PostXWorig
  corrX = PostXMod$meanx
  
  basis_dims = Model_listW$Nvec
  
  kmat3pt5_len = basis_dims[1]
  kmat180_len = basis_dims[2]
  
  #adjust 3.5 and 180 x's
  newX1 = corrX[1:kmat3pt5_len] * R_3pt5 
  newX2 = corrX[(kmat3pt5_len + 1):(kmat3pt5_len + kmat180_len)] * R_180
  
  newX = c(newX1, newX2)
  
  PostXMod$meanx = newX # replace post mean of X with modified post mean of x
  
  
  NewWRFPreds = get_preds(PostXMod, Model_listW) # get new fitted values
  
  OldPreds = get_preds(PostXWorig, Model_listW) # get original fitted values
  
  resids = (NewWRFPreds$TMAX - OldPreds$mypreds) * R_resid # resids using non-adjusted post 
  #resids = (NewWRFPreds$TMAX - NewWRFPreds$mypreds) * R_resid
  #add resids to the preds made using adjusted post mean of x
  CorrWRF = resids + NewWRFPreds$mypreds 
  return(CorrWRF)
}


adjust_x_crossval_alt = function(wrfmodel_test, obsmodel_train){
  
  Model_listW = get_sigmas_K(wrfmodel_test)
  Model_listO = get_sigmas_K(obsmodel_train)
  
  #post cond dist of x for future wrf model fit
  PostXWorig = post_cond_x(Model_listW) 
  
  OldPreds = get_preds(PostXWorig, Model_listW) #original fitted values
  
  PostXModadj = post_cond_x_alt(wrfmodel_test, obsmodel_train) #adjusted post x
  
  NewWRFpreds = get_preds(PostXModadj, Model_listW) # adjusted fitted values 
  
  resids = (NewWRFpreds$TMAX - OldPreds$mypreds) # Need to or not multiply be r_resid? I think not....?? 
  
  #add resids to the preds made using adjusted post mean of x
  CorrWRF = resids + NewWRFpreds$mypreds 
  return(CorrWRF)
}

post_cond_x_alt = function(wrfmodel_test, obsmodel_train){
  Model_listW = get_sigmas_K(wrfmodel_test)
  Model_listO = get_sigmas_K(obsmodel_train)
  modeldat = Model_listW$modeldat
  mu = Model_listW$mu
  sigmas = Model_listO$sigmas
  sige= Model_listO$sige
  Nvec = Model_listW$Nvec
  K = Model_listW$K
  sigx =  unlist(mapply(rep, sigmas, Nvec))
  y = Model_listW$y
  y = y - mu #must subtract intercept
  # compute the posterior precision for x
  Precx <- 1/sige^2*(t(K)%*%K) + diag(sigx^(-2))
  invPrec = ginv(Precx)
  RHS = as.vector(1/sige^2*(t(K)%*%y))
  meanx = invPrec %*% RHS
  #meanx <- solve(Precx, as.vector(1/sige^2*(t(K)%*%y))) #get the inverse here of precision matrix
  return(list(Precx = Precx, meanx = meanx))
}

post_cond_x = function(Model_list){
  modeldat = Model_list$modeldat
  mu = Model_list$mu
  sigmas = Model_list$sigmas
  sige= Model_list$sige
  Nvec = Model_list$Nvec
  K = Model_list$K
  sigx =  unlist(mapply(rep, sigmas, Nvec))
  y = Model_list$y
  y = y - mu #must subtract intercept
  # compute the posterior precision for x
  Precx <- 1/sige^2*(t(K)%*%K) + diag(sigx^(-2))
  
  meanx <- solve(Precx, as.vector(1/sige^2*(t(K)%*%y))) #get the inverse here of precision matrix
  return(list(Precx = Precx, meanx = meanx))
}

get_preds = function(PostX, Model_list){
  Kpred = Model_list$K
  meanx = PostX$meanx
  invPrecx = PostX$invPrecx
  mu = Model_list$mu #intercept
  TMAX = Model_list$y #this is just the response variable in the model
  sige = Model_list$sige
  mypreds = as.vector(Kpred %*% meanx + mu) #these are the fitted values
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
    Nvec[i] = ncol(LMEMODEL$data[[kmats[i]]])
  }
  
  modeldat = LMEMODEL$data
  y = LMEMODEL$data[,1]
  modeldat = modeldat[ , kmats]
  K = as.matrix(modeldat)
  
  return(list(sigmas = sigmas, sige = sige, Nvec = Nvec, K= K, modeldat = modeldat, 
              mu = LMEMODEL$coefficients$fixed, y= y, kmats = kmats))
}
