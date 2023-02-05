###########################################################################
# instead of using analytical solution, just grab x's from the model
# When using weighted lme models, you end up with singular matrix issues
# also works fine for non-weighted models and is a lot faster!
# when calculating mean x
##########################################################################

post_cond_x = function(LME_model){

  meanx = as.numeric(unlist(LME_model$coefficients$random))

  return(list(meanx=meanx))
}
head(meanx)

########################################################################
#function to get fitted values, given posterior mean object 
# note that preds are calculated using 'post mean x' from the model fit
# (random effects coefficients)
#and return value from get_sigmas_K()
########################################################################

get_preds = function(PostX, Model_list){
  Kpred = Model_list$K
  meanx = PostX$meanx
  mu = Model_list$mu #intercept--must be added to Kx
  TMAX = Model_list$y #this is just the response variable in the model
  
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
modelwrf=mod1
modelobs = obsmodel
R_resid=1
correct_wrf_x = function(modelwrf, modelobs, c_star = 0.9){
  Model_listW = get_sigmas_K(modelwrf) # get sigmas of WRF model fit on future period
  
  basis_dims = Model_listW$Nvec
  basis_names = Model_listW$kmats
  WRF_sigmas = Model_listW$sigmas
  
  Model_listO = get_sigmas_K(modelobs) # same for observed
  
  Obs_sigmas = Model_listO$sigmas
  
  Sigma_ratios = Obs_sigmas / WRF_sigmas
  
  R_resid = Model_listO$sige / Model_listW$sige
  
  #post cond dist of x for future wrf model fit
  PostXWorig = post_cond_x(modelwrf) 
  
  #now modify post mean of x and residuals
  PostXMod = PostXWorig
  corrX = PostXMod$meanx
  
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
par(mfrow=c(1,1))
sta = 100
sto=150
wrf_all[150, ]
oldwrf = mod1$data$Dtmax
plot(sta:sto, oldwrf[sta:sto], col = "blue", type = "l", xlab = "Day", ylab = "Detrended TMAX (°C)")
lines(sta:sto,CorrWRF[sta:sto], col = "red")
lines(sta:sto,oldwrf[sta:sto], col = 'blue')
legend("topright", legend = c(expression(tilde(T*"'")[Mod]),
                              expression(T*"'"[Mod])), lty = 1, col =c("red", "blue") )


detrend_again = function(wrfdat){
  wrfdat_new=wrfdat
  up1 = (wrfdat_new$tmax_norm_VAR1_adjust * wrfdat_new$Constant_scalar_wrf[1]) *wrfdat_new$Var_scaling_factors_wrf
  up2 = up1 + wrfdat_new$Seasonal_trend_wrf
  up3 = up2 + wrfdat_new$Long_term_trend_fut
  
  wrfdat_new$WRF_corr_adjust = up3
  # add long term trend back in and take it out
  tempdf = data.frame(WRF_TMAX = up3, dailyTimeStep = wrfdat_new$dailyTimeStep)
  
  wrftrend = gam(WRF_TMAX ~ s(dailyTimeStep, k = 50), data = tempdf)
  
  Long_term_trend_fut_new = wrftrend$fitted.values
  wrfdat_new$Long_term_trend_fut_new = Long_term_trend_fut_new
  
  # subtract long term trend again (it is nearly identical to old long term trend, but we do it again)
  wrfdat_new$detrended_adjust = up3 - Long_term_trend_fut_new
  
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
  
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Matern5_2")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
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

# plot(wrfdat_new$Seasonal_trend_wrf[1:365], type  ='l')
# lines(wrfdat_new$Seasonal_trend_wrf_new[1:365], col = "red")
# plot(wrfdat_new$Var_scaling_factors_wrf[1:365], type = "l")
# lines(wrfdat_new$Var_scaling_factors_wrf_new[1:365], col = "red")
detrend_again = function(wrfdat){
  wrfdat_new=wrfdat
  up1 = (wrfdat_new$tmax_norm_VAR1_adjust * wrfdat_new$Constant_scalar_wrf[1]) *wrfdat_new$Var_scaling_factors_wrf
  up2 = up1 + wrfdat_new$Seasonal_trend_wrf
  up3 = up2 + wrfdat_new$Long_term_trend_fut
  
  # add long term trend back in and take it out
  tempdf = data.frame(WRF_TMAX = up3, dailyTimeStep = wrfdat_new$dailyTimeStep)
  
  wrftrend = gam(WRF_TMAX ~ s(dailyTimeStep, k = 50), data = tempdf)
  
  Long_term_trend_fut_new = wrftrend$fitted.values
  wrfdat_new$Long_term_trend_fut_new = Long_term_trend_fut_new
  
  # subtract long term trend again (it is nearly identical to old long term trend, but we do it again)
  wrfdat_new$detrended_adjust = up3 - Long_term_trend_fut_new
  
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
  
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Gaussian")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  #make predictions
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
