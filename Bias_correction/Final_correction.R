# get observed data--get spatial average
library(dplyr)
traindf = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/TMAXRCP85_avesV2.csv")
ymd = traindf %>% group_by(MONTH, DAY) %>% summarise(meanT = mean(meanStationTMAX))
ymd$timestep = 1:365
ymd$meanT = NULL
TMAX = dplyr::select(traindf, YEAR, MONTH, DAY, meanStationTMAX)
TMAX = right_join(TMAX, ymd, by = c("MONTH", "DAY"))
TMAX = TMAX %>% arrange(YEAR, MONTH, DAY)

# get WRF for all years, get spatial average; also get WRF at station level
WRF = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF = WRF %>% group_by(YEAR, MONTH, DAY) %>% summarise(meanWRF = mean(test.pred.mu))
WRF = right_join(WRF, ymd, by = c("MONTH", "DAY"))
WRF = WRF %>% arrange(YEAR, MONTH, DAY)
WRF_stationlevel = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF = filter(WRF, YEAR < 2006)


####################################################################################
saveRDS(wrf_all, "wrf_all.Rds")
#WRF=filter(WRF, YEAR < 2006)

# process: detrend, take out seasonal mean and sd , constant scalar
wrf_data = process_WRF(WRF, kt=50)

# historical WRF
wrfhist = wrf_data$wrfhist

# historical and future WRF
wrfdat = wrf_data$wrfdat

# process obs: detrend using WRF trend, take out seasonal means and sd, constant scalar
obs_stuff = process_OBS(TMAX, wrfhist)

obstrain = obs_stuff$obsdat

obs_params = obs_stuff$obs_params

# join wrf data with observed trends
wrf_all = right_join(wrfdat, obs_params, by = "timestep")

# perform the correction at spatially-averaged scale
corr1 = ( wrf_all$tmax_norm_VAR1 * wrf_all$Constant_scalar_obs[1] ) 

corr2 = corr1 * wrf_all$Var_scaling_factors_obs

corr3 = corr2 + wrf_all$Long_term_trend_fut + wrf_all$Seasonal_trend_obs

wrf_all$WRF_corr = corr3

wrf_all$dailyTimeStep = 1:nrow(wrf_all)

# join spatially averaged , corrected wrf dataframe w/ station-level wrf
WRF_stationlevel = right_join(WRF_stationlevel, wrf_all, by = c("YEAR", "MONTH", "DAY"))

#saveRDS(wrf_all, "wrf_all.Rds")
# perform linear correction at station - level
corr_wrf = list()
#params1 = data.frame()
for (m in 1:nrow(wrf_all)){
  if (m %% 10000 == 0){
    print(m)
  }
  mytimestep = m
  tempdf = filter(WRF_stationlevel, dailyTimeStep == mytimestep)
  rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF for day m
  corrMean = tempdf$WRF_corr[1] # corrected spatially averaged WRF for day m
  
  sdObs = tempdf$Var_scaling_factors_obs[1] # scaling factors observed data 
  sdMod = tempdf$Var_scaling_factors_wrf[1]
  #mo = tempdf$MONTH[1]

  a = corrMean - rawMean*(sdObs/sdMod)
  b = sdObs / sdMod

  corrWRF = b*tempdf$test.pred.mu + a
  
  tempdf$TMAX_CORR_TEST = corrWRF
  corr_wrf[[m]] = tempdf
  # yr=tempdf$YEAR[1]
  # mo = tempdf$MONTH[1]
  # day=tempdf$DAY[1]
  # myparams = c(ObsMean, corrMean, a, b, yr, mo, day)
  # params1=rbind(params1, myparams)
}
corr_orig = do.call("rbind", corr_wrf)
saveRDS(corr_orig, "corr_orig.Rds")
names(params1) = c("ObsMean", "corrMean", "a", "b", "YEAR", "MONTH", "DAY")
head(params1)
params1$WRFMean=NULL
head(corr_orig)
corr_hist = filter(corr_orig, YEAR < 2006)
hist(corr_hist$TMAX, breaks = 60, col = 'blue', freq = FALSE)
hist(corr_hist$TMAX_CORR_TEST, breaks = 60, col = "red", freq=FALSE, add = T)


##################################################
# Now do the same thing but with adjusted WRF  
###################################################
obsdat = obs_stuff$obsdat
saveRDS(obsdat, "obsdat.Rds")

library(nlme)
# Future_1976_2005
mod1 = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/WRF_1976_2005.Rds")
# Future_2006_2036
mod2= readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/WRF_2006_2036.Rds")
# Future_2037_2067
mod3 = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/WRF_2037_2067.Rds")
# Future_2068_2099
mod4= readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/WRF_2069_2099.Rds")
# observed historical period
obsmodel = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/OBS_1976_2005.Rds")

# take a quick look at sigmas
 test = get_sigmas_K(mod1)
 test$sige
 test$sigmas
 test2 = get_sigmas_K(mod2)
 test2$sigmas
 test2$sige
 test3 = get_sigmas_K(mod3)
 test3$sigmas
 test3$sige
 test4 = get_sigmas_K(mod4)
 test4$sigmas
 test4$sige
 obtest=get_sigmas_K(obsmodel)
 obtest$sigmas
 obtest$sige

wrf1 = correct_wrf_x(modelwrf = mod1, modelobs = obsmodel)
wrf2 = correct_wrf_x(modelwrf = mod2, modelobs = obsmodel)
wrf3 = correct_wrf_x(modelwrf = mod3, modelobs = obsmodel)
wrf4 = correct_wrf_x(modelwrf = mod4, modelobs = obsmodel)

tmax_norm_VAR1_adjust = c(wrf1, wrf2, wrf3, wrf4)

# process: detrend, take out seasonal mean and sd , constant scalar
wrf_data = process_WRF(WRF)
# historical WRF
wrfhist = wrf_data$wrfhist

# historical and future WRF
wrfdat = wrf_data$wrfdat

# process obs: detrend using WRF trend, take out seasonal means and sd, constant scalar
obs_stuff = process_OBS(TMAX, wrfhist)

obs_params = obs_stuff$obs_params

wrfdat$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust

# here is where you need to detrend again! (b/c seasonal variance increases slightly due to adj of temp dependence)
# IF YOU WANT SECOND DETRENDING USING final_correction_second_detrend.R please!
#wrfdat_new  = detrend_again(wrfdat)

# join wrf data with observed trends
#wrf_all = right_join(wrfdat_new, obs_params, by = "timestep")
wrf_all = right_join(wrfdat, obs_params, by = "timestep")

# perform the correction at spatially-averaged scale
corr1 = ( wrf_all$tmax_norm_VAR1_adjust * wrf_all$Constant_scalar_obs[1] ) 

# multiply by observed seasonal sd
corr2 = corr1 * wrf_all$Var_scaling_factors_obs

corr3 = corr2 + wrf_all$Long_term_trend_fut + wrf_all$Seasonal_trend_obs

wrf_all$WRF_corr_adjust = corr3

wrf_all$dailyTimeStep = 1:nrow(wrf_all)
saveRDS(wrf_all, "wrf_all_adjust.Rds")

# join spatially averaged , corrected wrf dataframe w/ station-level wrf
WRF_stationlevel = right_join(WRF_stationlevel, wrf_all, by = c("YEAR", "MONTH", "DAY"))

# perform linear correction at station - level
corr_wrf = list()

for (m in 1:nrow(wrf_all)){
  if (m %% 10000 == 0){
    print(m)
  }
  
  mytimestep = m
  tempdf = filter(WRF_stationlevel, dailyTimeStep == mytimestep)
  rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF in test period
  corrMean = tempdf$WRF_corr_adjust[1] # corrected spatially averaged WRF in test period
  
  sdObs = tempdf$Var_scaling_factors_obs[1] # Seasonal SD of observed data
  sdMod = tempdf$Var_scaling_factors_wrf[1] # Seasonal SD of WRF data (from HetGP model)
  
  a = corrMean - rawMean*(sdObs/sdMod)
  b = sdObs / sdMod
  
  corrWRF = b*tempdf$test.pred.mu + a
  
  tempdf$TMAX_CORR_TEST = corrWRF
  corr_wrf[[m]] = tempdf
}
corr_df_wrf_adjust = do.call("rbind", corr_wrf)

corr_adjust_hist = filter(corr_df_wrf_adjust, YEAR < 2006)

hist(corr_adjust_hist$TMAX, breaks=60, col = "blue",freq = FALSE)
hist(corr_adjust_hist$TMAX_CORR_TEST, col = "red", breaks = 60, freq = FALSE, add = TRUE)
hist(corr_hist$TMAX_CORR_TEST, col = "green", breaks = 60, freq = FALSE, add = TRUE)
saveRDS(corr_df_wrf_adjust, "corr_df_wrf_adjust.Rds")
saveRDS(wrf_all, "wrf_all_adjust.Rds")

library(mgcv)
movMAD = function(vec, window){
  goodlength = length(vec) - window
  resvec = vector(length = goodlength)
  for (i in 1:goodlength){
    resvec[i] = mad(vec[i:(i+(window-1))])
  }
  return(resvec)
}


process_WRF = function(WRF, histyears = 1976:2005, MyCovariance = "Matern5_2", kt = 50){
  wrfdat=WRF
  #kt = 50
  wrfdat$dailyTimeStep = 1:nrow(WRF)
  
  wrftrend = gam(meanWRF ~ s(dailyTimeStep, k = kt), data = wrfdat)
  
  Long_term_trend_fut = wrftrend$fitted.values
  wrfdat$Long_term_trend_fut = Long_term_trend_fut
  wrfdat$WRF_detrended = wrfdat$meanWRF - Long_term_trend_fut
  preWRF = filter(wrfdat, timestep %in% 316:365)$WRF_detrended
  postWRF = filter(wrfdat, timestep %in% 1:50)$WRF_detrended
  
  num_years = length(unique(wrfdat$YEAR))
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  allWRF = c(preWRF, rep(wrfdat$WRF_detrended, 3), postWRF)
  
  het_WRF <- mleHetGP(alltime, allWRF, covtype = MyCovariance)
 # het_WRF <- mleHetGP(alltime, allWRF, covtype = "Gaussian")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  Seasonal_trend_wrf = rep(WRFpreds$mean[366:730], num_years)
  Var_scaling_factors_wrf = rep(sqrt(WRFpreds$sd2[366:730]+ WRFpreds$nugs[366:730]), num_years)
  
  wrfdat$Seasonal_trend_wrf = Seasonal_trend_wrf
  wrfdat$Var_scaling_factors_wrf = Var_scaling_factors_wrf

  wrfdat$resids = wrfdat$WRF_detrended - wrfdat$Seasonal_trend_wrf
  wrfdat$tmax_norm = wrfdat$resids / wrfdat$Var_scaling_factors_wrf
  
  Constant_scalar_wrf = mean(movMAD(wrfdat$tmax_norm, window = 100)) 
  wrfdat$tmax_norm_VAR1 = wrfdat$tmax_norm / Constant_scalar_wrf
  wrfdat$Constant_scalar_wrf = Constant_scalar_wrf
  
  wrfhist = filter(wrfdat, YEAR %in% histyears)
  #return(wrfdat)
  return(list(wrfdat = wrfdat, wrfhist = wrfhist))
  
}

process_OBS = function(obsdat, wrfdat){
  # take WRF trend out of observed data
  obsdat$detrended = obsdat$meanStationTMAX - wrfdat$Long_term_trend_fut

  preOBS = filter(obsdat, timestep %in% 316:365)$detrended
  postOBS = filter(obsdat, timestep %in% 1:50)$detrended
  
  num_years = length(unique(obsdat$YEAR))
  
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  allOBS = c(preOBS, rep(obsdat$detrended, 3), postOBS)
  
  het_OBS <- mleHetGP(alltime, allOBS, covtype = "Matern5_2")
  
  het_OBS_rebuild = rebuild(het_OBS, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  OBSpreds = predict(x = pred_times, object = het_OBS_rebuild)
  Seasonal_trend_obs = rep(OBSpreds$mean[366:730], num_years)
  
  Var_scaling_factors_obs = rep(sqrt(OBSpreds$sd2[366:730]+ OBSpreds$nugs[366:730]), num_years)
  
  obsdat$Seasonal_trend_obs = Seasonal_trend_obs
  obsdat$Var_scaling_factors_obs = Var_scaling_factors_obs
  
  obsdat$resids = obsdat$detrended - obsdat$Seasonal_trend_obs
  obsdat$tmax_norm = obsdat$resids / obsdat$Var_scaling_factors_obs
  
  Constant_scalar_obs = mean(movMAD(obsdat$tmax_norm, window = 100)) 
  obsdat$tmax_norm_VAR1 = obsdat$tmax_norm / Constant_scalar_obs
  obsdat$Constant_scalar = Constant_scalar_obs
  
  obs_params = data.frame(Var_scaling_factors_obs = obsdat$Var_scaling_factors_obs[1:365],
                          Constant_scalar_obs = Constant_scalar_obs,
                          Seasonal_trend_obs = Seasonal_trend_obs[1:365],
                          timestep = 1:365)
  return(list(obs_params = obs_params, obsdat = obsdat))
  
}
####################################################################################
## IMPORTANT!!
## see Final_correction_modified_functions.R for slightly modified versions of
## post_cond_x, get_preds, and correct_wrf_x--made to deal with weighted LME models
## but also work fine for non weighted lme models and it's faster!
####################################################################################
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

######################################################################
#function to get posterior conditional mean of x
#input is the return value object of get_sigmas_K(LMEMODEL)
# using analytical solution (slower and MAY FAIL)
######################################################################

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
  
  # get precision matrix - leave off for now
  # invPrecx <- solve(Precx)
  return(list(Precx = Precx, meanx = meanx))
}


########################################################################
#function to get fitted values, given posterior mean object 
#and return value from get_sigmas_K()
########################################################################

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

correct_wrf_x = function(modelwrf, modelobs, c_star = 0.97){
  Model_listW = get_sigmas_K(modelwrf) # get sigmas of WRF model fit on future period
  
  basis_dims = Model_listW$Nvec
  basis_names = Model_listW$kmats
  WRF_sigmas = Model_listW$sigmas
  
  Model_listO = get_sigmas_K(modelobs) # same for observed
  
  Obs_sigmas = Model_listO$sigmas
  
  Sigma_ratios = Obs_sigmas / WRF_sigmas
  
  R_resid = Model_listO$sige / Model_listW$sige
  
  #post cond dist of x for future wrf model fit
  PostXWorig = post_cond_x(Model_listW) 
  
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


#IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT
#IMPORTANT--get seasonal means and update datasets
#these are lme fits where only periodic basis matrices are used (14, 40, 60, 90, 180)
#predict to get seasonal mean, then also add the intercept to get the actual seasonal mean
# just updating the datasets here!
# NOTE: in OBS_ALL_Future and WRF_ALL_1976_2005 in ...Final_datasets/Matern I just detrended using the 
# seasonal mean from HetGP which is nearly identical to what you get from PC model
# However, I included the PC model estimates for the seasonal mean in the datasets as well w/ postfix 'PC'
# so this stuff below is not really important
obsdat$Long_term_trend_fut_wrf = wrfdat$Long_term_trend_fut[1:nrow(obsdat)]
wrfseas = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/WRF_MAT_ALL_seasonalOnly.Rds")
obseas = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/OBS_MAT_seasonal.Rds")
wrfsigs = get_sigmas_K(wrfseas)
obsigs = get_sigmas_K(obseas)
wrfmu = wrfsigs$mu
obsmu = obsigs$mu
wrfpreds = predict(wrfseas)
obspreds = predict(obseas)
plot(wrfpreds[1:365], type = "l", ylim = c(-18, 17))
lines(obspreds[1:365], col = "red")

obsdat$Seasonal_trend_obs_PC = obspreds
wrfdat$Seasonal_trend_wrf_PC = wrfpreds


write.csv(obsdat, "C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/OBS_ALL_Future.csv")
write.csv(wrfdat, "C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/WRF_ALL_1976_2099.csv")

testdat = filter(wrfdat, YEAR < 2006)
testdat$Mean_seasonal_trend_PC = predsw
testdat$WRF_seas_detrended = testdat$meanWRF - testdat$Long_term_trend_fut - testdat$Mean_seasonal_trend_PC
write.csv(testdat, "C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/WRF_vartest_data.csv")
write.csv(Otestdat, "C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/OBS_vartest_data.csv")

colnames(obsdat)
Otestdat = obsdat
Obs_seas_detrended = Otestdat$meanStationTMAX - Otestdat$Long_term_trend_fut_wrf - as.numeric(Otestdat$Seasonal_trend_obs_PC)
Otestdat$Obs_seas_detrended = Obs_seas_detrended
pp=Otestdat$meanStationTMAX - Otestdat$Long_term_trend_fut_wrf - Otestdat$Seasonal_trend_obs_PC

head(Otestdat$Obs_seas_detrended)
plot(Otestdat$Var_scaling_factors_obs)

# junk/maybe not necessary
# upscale to TMAX scale and then remove long term, seasonal mean and SD trends again
# may not actually be necessary
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

######################################################################
## Dave's alternative method of correction that does not work so well
######################################################################
#try alternative way of correction....
corr_wrf = list()
params = data.frame()
for (m in 1:10950){
  if (m %% 10000 == 0){
    print(m)
  }
  mytimestep = m
  tempdf = filter(WRF_stationlevel, dailyTimeStep == mytimestep)
  rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF in test period
  WRFMean = tempdf$Seasonal_trend_wrf[1] # corrected spatially averaged WRF in test period
  ObsMean = tempdf$Seasonal_trend_obs[1]
  deltaMean = ObsMean-WRFMean
  corrMean = rawMean + deltaMean
  
  sdObs = tempdf$Var_scaling_factors_obs[1] # scaling factors observed data 
  sdMod = tempdf$Var_scaling_factors_wrf[1]
  
  a = corrMean - rawMean*(sdObs/sdMod)
  b = sdObs / sdMod
  
  corrWRF = b*tempdf$test.pred.mu + a
  
  tempdf$TMAX_CORR_TEST = corrWRF
  corr_wrf[[m]] = tempdf
  yr=tempdf$YEAR[1]
  mo = tempdf$MONTH[1]
  day=tempdf$DAY[1]
  myparams = c(WRFMean, ObsMean, deltaMean, corrMean, a, b, yr, mo, day)
  params=rbind(params, myparams)
}
corr_df_wrf_het_alt  = do.call("rbind", corr_wrf)
names(params) = c("WRFMean", "ObsMean", "deltaMean", "corrMean", "a", "b", "YEAR", "MONTH", "DAY")
head(params)
aves1 = corr_orig %>% group_by(YEAR, MONTH, DAY) %>% summarise(meancorr = mean(TMAX_CORR_TEST),
                                                               meanObs = mean(TMAX), 
                                                               meanMod = mean(test.pred.mu),
                                                               sdcorr = sd(TMAX_CORR_TEST),
                                                               sdObs = sd(TMAX), 
                                                               sdMod = sd(test.pred.mu)
)

aves2 = corr_df_wrf_het_alt %>% group_by(YEAR, MONTH, DAY) %>% 
  summarise(meancorr = mean(TMAX_CORR_TEST), sdcorr= sd(TMAX_CORR_TEST))

aves1_b = aves1 %>% group_by(MONTH) %>% summarise(sdcorr= sd(meancorr),
                                                  sdObs = sd(meanObs),
                                                  sdMod = sd(meanMod))

aves2_b= aves2 %>% group_by(MONTH) %>% summarise(sdcorr= sd(meancorr))

# sd after averaging by ymd
plot(1:nrow(aves1_b), aves1_b$sdObs, type = "l", ylim=c(2.5, 8))
lines(1:nrow(aves1_b), aves1_b$sdMod, col = "red")
lines(1:nrow(aves1_b), aves1_b$sdcorr, col = "blue")
lines(1:nrow(aves1_b), aves2_b$sdcorr, col= "green")
legend("bottomleft", legend = c("obs", "mod", "corr_orig", "corr_alt"), 
       col = c("black", "red","blue", "green"), lwd=1)

# sd before averaging by ymd
plot(1:nrow(aves1), aves1$sdObs, type = "l")
lines(1:nrow(aves1), aves1$sdMod, col = "red")
lines(1:nrow(aves1), aves1$sdcorr, col = "blue")
lines(1:nrow(aves1), aves2$sdcorr, col= "green")


paves1 = params1 %>% group_by(YEAR, MONTH, DAY) %>% summarise(Cmean = mean(corrMean), mean_a = mean(a))
paves = params %>% group_by(YEAR, MONTH, DAY) %>% summarise(Cmean = mean(corrMean), mean_a = mean(a))

pstuff1= paves1 %>% group_by(MONTH, DAY) %>% summarise(corrMeanSD = sd(Cmean), sd_a = sd(mean_a))
pstuff= paves %>% group_by(MONTH, DAY) %>% summarise(corrMeanSD = sd(Cmean), sd_a = sd(mean_a))


plot(1:nrow(pstuff1), pstuff1$corrMeanSD, type = "l", col = "red")
lines(1:nrow(pstuff), pstuff$corrMeanSD, col = "green")
legend("bottomleft", legend = c("SD of orig corrMean aves", "SD of alt corrMean aves"), 
       col = c("red","green"), lwd=1)

plot(pstuff1$sd_a, pstuff$sd_a)
abline(0,1)
summary(paves1$mean_a)
summary(paves$mean_a)


plot(1:nrow(pstuff), pstuff$sd_a, col = "green", type = "l")
lines(1:nrow(pstuff1), pstuff1$sd_a, type = "l", col = "red")
legend(200, .8, legend = c("SD of a aves (orig)", "SD of a aves (alt)"), 
       col = c("red", "green"), lwd=1)



plot(1:nrow(aves1), aves1$sdObs, col = "blue", type = "l")
lines(1:nrow(aves1), aves1$sdMod, col = "green", type="l")
lines(1:nrow(aves1), aves1$sdcorr, col = "darkgreen")
lines(1:nrow(aves1), aves2$sdcorr, col = "purple")
aves1
aves1$sdObs/aves1$sdcorr
head(params1)
head(params$b)

######################################################################
## if you want to detrend after adjusting WRF -- may not be necessary
######################################################################
##############################
# process: detrend, take out seasonal mean and sd , constant scalar
wrf_data = process_WRF(WRF)
# historical WRF
wrfhist = wrf_data$wrfhist

# historical and future WRF
wrfdat = wrf_data$wrfdat

# process obs: detrend using WRF trend, take out seasonal means and sd, constant scalar
obs_stuff = process_OBS(TMAX, wrfhist)

obs_params = obs_stuff$obs_params

wrfdat$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust

# here is where you need to detrend again! (b/c seasonal variance increases slightly due to adj of temp dependence)
wrfdat_new  = detrend_again(wrfdat)

# join wrf data with observed trends
wrf_all = right_join(wrfdat_new, obs_params, by = "timestep")
wrf_all = right_join(wrfdat, obs_params, by = "timestep")

# perform the correction at spatially-averaged scale
# tmax_norm_VAR1_adjust_prime is mean0 var1 WRF with extra detrending from function detrend_again()
# all of the 'new' trends resulting from the second detrending have suffix 'new'
corr1 = ( wrf_all$tmax_norm_VAR1_adjust_prime * wrf_all$Constant_scalar_obs[1] ) 
#corr1 = ( wrf_all$tmax_norm_VAR1_adjust * wrf_all$Constant_scalar_obs[1] ) 

corr2 = corr1 * wrf_all$Var_scaling_factors_obs
#corr2 = corr1 * wrf_all$Var_scaling_factors_obs

corr3 = corr2 + wrf_all$Long_term_trend_fut_new + wrf_all$Seasonal_trend_obs
#corr3 = corr2 + wrf_all$Long_term_trend_fut + wrf_all$Seasonal_trend_obs

wrf_all$WRF_corr_adjust_prime = corr3
#wrf_all$WRF_corr_adjust = corr3

wrf_all$dailyTimeStep = 1:nrow(wrf_all)

#saveRDS(wrf_all, "wrf_all_adjust_prime.Rds")
# join spatially averaged , corrected wrf dataframe w/ station-level wrf
WRF_stationlevel = right_join(WRF_stationlevel, wrf_all, by = c("YEAR", "MONTH", "DAY"))

orig=WRF_stationlevel
WRF_stationlevel=filter(WRF_stationlevel, MONTH %in% c(4,5))

# perform linear correction at station - level
corr_wrf = list()

for (m in 1:nrow(wrf_all)){
  if (m %% 10000 == 0){
    print(m)
  }
  
  mytimestep = m
  tempdf = filter(WRF_stationlevel, dailyTimeStep == mytimestep)
  rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF in test period
  corrMean = tempdf$WRF_corr_adjust_prime[1] # corrected spatially averaged WRF in test period (including second detrending step)
  
  sdObs = tempdf$Var_scaling_factors_obs[1] # scaling factors observed data 
  sdMod = tempdf$Var_scaling_factors_wrf_new[1] # use the 'new' WRF sd obtained from function detrend_again()
  
  a = corrMean - rawMean*(sdObs/sdMod)
  b = sdObs / sdMod
  
  corrWRF = b*tempdf$test.pred.mu + a
  
  tempdf$TMAX_CORR_TEST = corrWRF
  corr_wrf[[m]] = tempdf
}
corr_df_wrf_adjust_het_prime = do.call("rbind", corr_wrf)

tail(corr_df_wrf_adjust_het_prime)
saveRDS(corr_df_wrf_adjust_het_prime, "corr_df_wrf_adjust_het_prime.Rds")