# IMPORTANT IMPORTANT IMPORTANT
# USE THIS IF YOU WANT TO APPLY THE SECOND DETRENDING AFTER ADJUSTMENT OF X
# Otherwise, use Cross_validation_take6
# Note that take6 functions are used 
library(mgcv)
library(qmap)
library(hetGP)
library(dplyr)

source("cross_validation_take6_functions.R")
#use this---we are doing 10 fold cross val
cvdf = read.csv("cvdfConsecYears10.csv")
#cvdf = read.csv("cvdfConsecYears.csv")

traindf = read.csv("TMAXRCP85_avesV2.csv")

ymd = traindf %>% group_by(MONTH, DAY) %>% summarise(meanT = mean(meanStationTMAX))
ymd$timestep = 1:365
ymd$meanT = NULL
TMAX = dplyr::select(traindf, YEAR, MONTH, DAY, meanStationTMAX)
WRF = dplyr::select(traindf, YEAR, MONTH, DAY, meanWRFTMAX)
WRF$meanWRF = WRF$meanWRFTMAX
WRF$meanWRFTMAX = NULL

TMAX = right_join(TMAX, ymd, by = c("MONTH", "DAY"))
WRF = right_join(WRF, ymd, by = c("MONTH", "DAY"))

# station level data
newdf =  readRDS("newDF_with_white_REVISED_3pt5.Rds")

newdf2 = dplyr::select(newdf, YEAR ,MONTH,DAY, TMAX, ID, LAT,LON, ELEV, test.pred.mu)
# use adjust_x = TRUE if you want to do temporal adjustment, otherwise use adjust_x = FALSE
Do_crossval5 = function(stationlevel_dat = newdf2, TMAXdat = TMAX, WRFdat = WRF, cross_val_df = cvdf, adjust_x = TRUE,
                        second_detrend = TRUE){
  
  TMAX = TMAXdat
  WRF = WRFdat
  newdf2 = stationlevel_dat
  cvdf=cross_val_df
  
  # get seasonal means/sds/process WRF. Only need to do this once on entire WRF dataset
  wrfall = process_WRFtest(WRF)
  
  # list to hold MAEs/metrics
  bias_list = list()
  Err_list = list()
  corrwrf_list = list()
  Err_listDOY = list()
  Kfold = unique(cvdf$K)
  for (i in Kfold){
    print(i)
    trainyears = filter(cvdf,  K != i)$YEAR
    
    testyears = filter(cvdf,  K == i)$YEAR
    
    newdf_train = filter(newdf2, YEAR %in% trainyears)
    newdf_test = filter(newdf2, YEAR %in% testyears)
    
    newdf_test_qmap = list()
    for (j in 1:12){
      monthtrain = filter(newdf_train, MONTH == j)
      monthtest = filter(newdf_test, MONTH == j)
      obs = monthtrain$TMAX
      mod = monthtrain$test.pred.mu
      mod_test = monthtest$test.pred.mu
      
      mytf = fitQmapQUANT(obs, mod, wet.day = F, qstep = 0.0001)
      
      test_corr = doQmapQUANT(mod_test, mytf)
      monthtest$EQM_CORR = test_corr
      newdf_test_qmap[[j]] = monthtest
    }
    
    newdf_test_qmap = do.call("rbind", newdf_test_qmap)
    
    if (adjust_x == TRUE){
      # all WRF is processed in exactly the same way. This PC was run on processed WRF data so can 
      # be used for all folds. However, sigmas from PC model run on observed training data is used to adjust x's
      wrfLMEmodel = readRDS("WRF_ALL_forCrossVal.Rds")
      obsLMEmodel = readRDS(paste("OBS_TRAIN_", i, 
                                  ".Rds", sep = ""))
      
      # 
      # adjust WRF x's using sigmas from observed training data in fold k
      # add the column to entire WRF dataset (wrfall from above)
      # adjusted mean 0 var 1 WRF is called tmax_norm_VAR1_adjust
      # **apply second detrending**
      wrf_data_adjust_x = add_adjust_x_postProcess(wrfdat = wrfall, testyears = testyears,
                                                   wrfmodel_test = wrfLMEmodel, 
                                   obsmodel_train = obsLMEmodel)
      
      wrftest = filter(wrf_data_adjust_x, YEAR %in% testyears)
      
      # write.csv(wrftest, paste("WRF_test_sets/wrftest", i, ".csv", sep = ""))
      
      # process observed data--take out WRF long term trend (from second detrending), seasonal mean/sd, make mean 0 var 1
      # if you use second_detrend = TRUE, then the new long term term is subtracting from observed data;
      # otherwise, the long term trend from original WRF is used
      obs_data = process_OBS(obsdat = TMAX, wrfdat = wrf_data_adjust_x, trainyears = trainyears, 
                             testyears = testyears, second_detrend = TRUE)
      
      obstest = obs_data$obstest
      
      obstrain = obs_data$obstrain
      
      obs_params = obs_data$obs_params
      
      # upscale to TMAX scale, apply correction
      wrftest = right_join(wrftest, obs_params, by = "timestep")
    
      corr1 = ( wrftest$tmax_norm_VAR1_adjust_prime * wrftest$Constant_scalar_obs[1] )
      
      corr2 = corr1 * wrftest$Var_scaling_factors_obs
      
      #note that the 'new' long term trend is added back in!
      corr3 = corr2 + wrftest$Long_term_trend_fut_new + wrftest$Seasonal_trend_obs
    }else{
      print("not adjusting x's !")
      ###############################################################################
      
      wrftest = filter(wrfall, YEAR %in% testyears)
      
      # NOTE: only observed data in training set is used to get seasonal trends
      # we use second_detrend = FALSE, b/c there is no adjustment of x's...so second detrending is not relevant
      obs_data = process_OBS(obsdat = TMAX, wrfdat = wrfall, trainyears = trainyears, 
                             testyears = testyears, second_detrend = FALSE)
      
      obstest = obs_data$obstest
      
      obstrain = obs_data$obstrain
      
      obs_params = obs_data$obs_params
      
      wrftest = right_join(wrftest, obs_params, by = "timestep")
      
      # upscale to TMAX scale, apply correction
      corr1 = ( wrftest$tmax_norm_VAR1 * wrftest$Constant_scalar_obs[1] ) 
    
      corr2 = corr1 * wrftest$Var_scaling_factors_obs
    
      corr3 = corr2 + wrftest$Long_term_trend_fut + wrftest$Seasonal_trend_obs
    }
    
    wrftest$WRF_corr = corr3
    
   
    wrftest$dailyTimeStep = 1:nrow(wrftest)
    
    newdf_wrftest = right_join(newdf_test_qmap, wrftest, by = c("YEAR", "MONTH", "DAY"))
    
    corr_wrf = list()
    
    for (m in 1:nrow(wrftest)){
      mytimestep = m
      tempdf = filter(newdf_wrftest, dailyTimeStep == mytimestep)
      rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF in test period
      corrMean = tempdf$WRF_corr[1] # corrected spatially averaged WRF in test period
      
      sdObs = tempdf$Var_scaling_factors_obs[1] # scaling factors observed data (from train period)
      if (second_detrend == TRUE){
        sdMod = tempdf$Var_scaling_factors_wrf_new[1] # note 'new' means this came from 2nd detrending after adjusting x's
      }else{
        sdMod = tempdf$Var_scaling_factors_wrf[1]
      }
      # sdMod = sd(tempdf$test.pred.mu)
      #apply linear correction
      a = corrMean - rawMean*(sdObs/sdMod)
      b = sdObs / sdMod
      
      corrWRF = b*tempdf$test.pred.mu + a
      
      tempdf$TMAX_CORR_TEST = corrWRF
      corr_wrf[[m]] = tempdf
    }
    
    corr_df_wrf  = do.call("rbind", corr_wrf)
    
    
    errdf = data.frame()
    
    for (month in 1:12){
      monthdat = filter(corr_df_wrf, MONTH == month)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      sortMYCORR = quantile(monthdat$TMAX_CORR_TEST, seq(0,1, length.out = 1000))
      ModObs_err = mean(abs(sortObs - sortMod))
      EQM_err = mean(abs(sortObs - sortEQM))
      MY_err = mean(abs(sortObs - sortMYCORR))
      errVec = c(ModObs_err, EQM_err, MY_err)
      errdf = rbind(errdf, errVec)
    }
    names(errdf) = c("Raw_MAE", "EQM_MAE", "NEW_MAE")
    errdf$MONTH  =1:12
    errdf$KFOLD = i
    
    biasdf = data.frame()
    
    for (month in 1:12){
      monthdat = filter(corr_df_wrf, MONTH == month)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      sortMYCORR = quantile(monthdat$TMAX_CORR_TEST, seq(0,1, length.out = 1000))
      ModObs_err = mean(sortObs - sortMod)
      EQM_err = mean(sortObs - sortEQM)
      MY_err = mean(sortObs - sortMYCORR)
      errVec = c(ModObs_err, EQM_err, MY_err)
      biasdf = rbind(biasdf, errVec)
    }
    names(biasdf) = c("Raw_bias", "EQM_bias", "NEW_bias")
    biasdf$MONTH  =1:12
    biasdf$KFOLD = i
    
    
    ymd = corr_df_wrf %>% group_by(MONTH, DAY) %>% summarise(meanTMAX = mean(TMAX))
    ymd$meanTMAX = NULL
    ymd$DOY = 1:365
    corr_df_wrf = right_join(corr_df_wrf, ymd, by = c("MONTH", "DAY"))
    errdfD = data.frame()
    for (d in 1:365){
      monthdat = filter(corr_df_wrf, DOY == d)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      sortMYCORR = quantile(monthdat$TMAX_CORR_TEST, seq(0,1, length.out = 1000))
      ModObs_err = mean(abs(sortObs - sortMod))
      EQM_err = mean(abs(sortObs - sortEQM))
      MY_err = mean(abs(sortObs - sortMYCORR))
      errVec = c(ModObs_err, EQM_err, MY_err)
      errdfD = rbind(errdfD, errVec)
    }
    
    names(errdfD) = c("Raw_MAE", "EQM_MAE", "NEW_MAE")
    errdfD$DOY  = 1:365
    errdfD$KFOLD = i
    
    bias_list[[i]] = biasdf
    
    Err_list[[i]] = errdf
    
    Err_listDOY[[i]] = errdfD
    
    corrwrf_list[[i]] = corr_df_wrf
    
  }
  return(list(Err_list = Err_list, Err_listDOY = Err_listDOY, bias_list = bias_list, 
              corrwrf_list = corrwrf_list))
}
starttime = Sys.time() 
newres_adjust = Do_crossval5(second_detrend = TRUE)
stoptime = Sys.time()

saveRDS(newres_adjust, "CV_newMethod_secondDetrendTRUE.Rds")
print(stoptime-starttime)


