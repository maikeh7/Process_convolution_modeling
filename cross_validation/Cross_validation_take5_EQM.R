###################################################################################################################
# This code is ONLY for cross-validation using EQM WITH TEMPORAL ADJUSTMENT
# use cross_validation_take6 to perform cross-validation via EQM without temporal adjustment (my method can be done with and without temporal adjustment in that code)

library(mgcv)
library(qmap)
library(hetGP)
library(dplyr)
#cvdf = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/cvdf_consecYears.csv")

source("cross_validation_take6_functions.R")
#use this---we are doing 10 fold cross val
cvdf = read.csv("cvdfConsecYears10.csv")

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

Do_crossval5_EQM = function(stationlevel_dat = newdf2, TMAXdat = TMAX, WRFdat = WRF, cross_val_df = cvdf){
  
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
    
      # all WRF is processed in exactly the same way. This PC was run on processed WRF data so can 
      # be used for all folds. However, sigmas from PC model run on observed training data is used to adjust x's
      wrfLMEmodel = readRDS("WRF_ALL_forCrossVal.Rds")
     
      obsLMEmodel = readRDS(paste("OBS_TRAIN_", i,".Rds", sep = ""))
       
      # adjust WRF x's using sigmas from observed training data in fold k
      # add the column to entire WRF dataset (wrfall from above)
      # adjusted mean 0 var 1 WRF is called tmax_norm_VAR1_adjust
      # This is NOT necessary for EQM w/ temporal adjustment b/c we WILL want the 
      # second detrending which occurs in add_adjust_x_postProcess() below
      # wrf_data = process_WRF_adjustX(wrfdat = wrfall, testyears = testyears, 
      #                       wrfmodel_test = wrfLMEmodel, 
      #                       obsmodel_train = obsLMEmodel)
                             
      # Apply second detrending after x's have been adjusted second detrending
      wrf_data_adjust_x = add_adjust_x_postProcess(wrfdat = wrfall, testyears = testyears,
                                                   wrfmodel_test = wrfLMEmodel, 
                                                   obsmodel_train = obsLMEmodel)
     #plot(wrf_data_adjust_x$Var_scaling_factors_wrf[1:365], col = "blue", ylim=c(2.9, 8.5))
     #lines(wrf_data_adjust_x$Var_scaling_factors_wrf_new[1:365], col = "red")
      # filter out test set
      wrftest = filter(wrf_data_adjust_x, YEAR %in% testyears)
      
      #write.csv(wrftest, paste("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Crossval/WRF_test_sets/wrftest", i, ".csv", sep = ""))
      
      # I don't think this is necessary...
      # process observed data--take out WRF long term trend, seasonal mean/sd, make mean 0 var 1
      # obs_data = process_OBS2(obsdat = TMAX, wrfdat = wrf_data_adjust_x, trainyears = trainyears, 
      #                      testyears = testyears, adjust_x = TRUE)
            
      #obstest = obs_data$obstest
      
      #obstrain = obs_data$obstrain
      
      #obs_params = obs_data$obs_params
      
      #wrftest = right_join(wrftest, obs_params, by = "timestep")
      
           # upscale to TMAX scale (no correction here--just use WRF values. We only want to incorporate the modified temporal dependence)
      wrf1 = ( wrf_data_adjust_x$tmax_norm_VAR1_adjust_prime * wrf_data_adjust_x$Constant_scalar_wrf_new[1] )
      
      wrf2 = wrf1 * wrf_data_adjust_x$Var_scaling_factors_wrf_new
      
      #note that the 'new' long term trend is added back in!
      wrf3 = wrf2 + wrf_data_adjust_x$Long_term_trend_fut_new + wrf_data_adjust_x$Seasonal_trend_wrf_new
      
      #ok so wrf3 is adjusted wrf...but not yet corrected
      wrf_data_adjust_x$wrf_adjust = wrf3
      
      # include important variables. Note that the var's ending w/ 'new' are the seasonal means/sd's after adjusting x.
     
      wrf_data_adjust_x = dplyr::select(wrf_data_adjust_x, YEAR, MONTH, DAY, wrf_adjust, 
                                        Seasonal_trend_wrf_new,
                                        Seasonal_trend_wrf, Var_scaling_factors_wrf, 
                                        Var_scaling_factors_wrf_new, dailyTimeStep)
      # newdf2 = station-level WRF data
      newdf_wrftest = right_join(newdf2, wrf_data_adjust_x, by = c("YEAR", "MONTH", "DAY"))
     
      eqm_list = list()
     
      # apply adjustment of x's of WRF to raw WRF at station level
      for (m in 1:nrow(wrf_data_adjust_x)){
        mytimestep = m
        tempdf = filter(newdf_wrftest, dailyTimeStep == mytimestep)
        rawMean = mean(tempdf$test.pred.mu) # spatial mean of raw WRF 
        corrMean = tempdf$wrf_adjust[1] # adjusted spatially averaged WRF 
        
        #sdNew = tempdf$Var_scaling_factors_wrf_new[1] # seasonal sd after second detrending
        #sdOld = tempdf$Var_scaling_factors_wrf[1] #seasonal sd of original WRF before adjustment of x's
        sdNew = 1
        sdOld = 1
        
        # sdMod = sd(tempdf$test.pred.mu)
        #apply linear correction
        a = corrMean - rawMean*(sdNew/sdOld)
        b = sdNew / sdOld
        
        WRF_eqm = b*tempdf$test.pred.mu + a
        
        tempdf$test.pred.mu_adjust = WRF_eqm
        eqm_list[[m]] = tempdf
      }
      
    #plot(wrf_data_adjust_x$Seasonal_trend_wrf[1:365], type = "l")
    #lines(wrf_data_adjust_x$Seasonal_trend_wrf_new[1:365], col = 'red')
    eqm_df  = do.call("rbind", eqm_list)
    
    #hist(eqm_df$test.pred.mu, freq = FALSE, breaks = 60, ylim=c(0, .05), add = TRUE)
    #ist(eqm_df$test.pred.mu_adjust, freq= FALSE, breaks = 60, col = 'red', add = TRUE)
    # Finally, split into train/test sets and carry out EQM as usual on the adjusted WRF data
    newdf_train = filter(eqm_df, YEAR %in% trainyears)
    newdf_test = filter(eqm_df, YEAR %in% testyears)
    
    newdf_test_qmap = list()
    for (j in 1:12){
      monthtrain = filter(newdf_train, MONTH == j)
      monthtest = filter(newdf_test, MONTH == j)
      obs = monthtrain$TMAX
      mod = monthtrain$test.pred.mu_adjust
      mod_test = monthtest$test.pred.mu_adjust
      
      mytf = fitQmapQUANT(obs, mod, wet.day = F, qstep = 0.0001)
      
      test_corr = doQmapQUANT(mod_test, mytf)
      monthtest$EQM_CORR = test_corr
      newdf_test_qmap[[j]] = monthtest
    }
    
    newdf_test_qmap = do.call("rbind", newdf_test_qmap)
    
    
    errdf = data.frame()
    
    for (month in 1:12){
      monthdat = filter(newdf_test_qmap, MONTH == month)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      ModObs_err = mean(abs(sortObs - sortMod))
      EQM_err = mean(abs(sortObs - sortEQM))
      errVec = c(ModObs_err, EQM_err)
      errdf = rbind(errdf, errVec)
    }
    
    names(errdf) = c("Raw_MAE", "EQM_MAE")
    errdf$MONTH  =1:12
    errdf$KFOLD = i

    biasdf = data.frame()
    
    for (month in 1:12){
      monthdat = filter(newdf_test_qmap, MONTH == month)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      ModObs_err = mean(sortObs - sortMod)
      EQM_err = mean(sortObs - sortEQM)
      errVec = c(ModObs_err, EQM_err)
      biasdf = rbind(biasdf, errVec)
    }
    names(biasdf) = c("Raw_bias", "EQM_bias")
    biasdf$MONTH  =1:12
    biasdf$KFOLD = i
    
    
    ymd = newdf_test_qmap %>% group_by(MONTH, DAY) %>% summarise(meanTMAX = mean(TMAX))
    ymd$meanTMAX = NULL
    ymd$DOY = 1:365
    newdf_test_qmap = right_join(newdf_test_qmap, ymd, by = c("MONTH", "DAY"))
    errdfD = data.frame()
    for (d in 1:365){
      monthdat = filter(newdf_test_qmap, DOY == d)
      sortMod = quantile(monthdat$test.pred.mu, seq(0,1, length.out = 1000))
      sortObs = quantile(monthdat$TMAX, seq(0,1, length.out = 1000))
      sortEQM = quantile(monthdat$EQM_CORR, seq(0,1, length.out = 1000))
      ModObs_err = mean(abs(sortObs - sortMod))
      EQM_err = mean(abs(sortObs - sortEQM))
      errVec = c(ModObs_err, EQM_err)
      errdfD = rbind(errdfD, errVec)
    }
    
    names(errdfD) = c("Raw_MAE", "EQM_MAE")
    errdfD$DOY  = 1:365
    errdfD$KFOLD = i
    
    bias_list[[i]] = biasdf
    
    Err_list[[i]] = errdf
    
    Err_listDOY[[i]] = errdfD
    
    
    
  }
  return(list(Err_list = Err_list, Err_listDOY = Err_listDOY, bias_list = bias_list))
}

newres_adjust_EQM = Do_crossval5_EQM()

saveRDS(newres_adjust_EQM, "CV_results_EQM_adjustxTRUE_SDs1.Rds")
print("done")