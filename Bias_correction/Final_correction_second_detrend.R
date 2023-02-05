##################################################
# Do adjustment but include second detrending
###################################################
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

# 
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
hist(wrfhist$tmax_norm_VAR1, breaks = 60, freq = FALSE, xlim = c(-5, 5))
hist(tmax_norm_VAR1_adjust, breaks = 60, freq = FALSE, add = TRUE, col = "blue")

summary(tmax_norm_VAR1_adjust)
summary(wrfhist$tmax_norm_VAR1)
wrfmeans = rollmean(wrfhist$tmax_norm_VAR1, k = 3)
wrfadjust_means =  rollmean(tmax_norm_VAR1_adjust, k =3)
obsmeans = rollmean(obsdat$tmax_norm_VAR1, k = 3)

png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Histograms.png", width = 1500, height=1000, res=250)
par(mfrow=c(3,1), mar=c(1,0,0,0), oma = c(4,4,1,1))
par(mfrow=c(1,1))
hist(obsmeans, breaks = 60, freq = FALSE, col = "gray",xaxt="n", xlim = c(-5, 5),ylim=c(0, .7), xlab = "",
     main="")
text(-4.5, .4, "a", cex=2)
hist(wrfmeans, breaks = 60, freq = FALSE, xlim = c(-5, 5),ylim=c(0, .5), xaxt="n", xlab="", col = 'blue',
     main= "", add = TRUE)
text(-4.5, .4, "b", cex=2)
hist(wrfadjust_means, breaks = 60, freq = FALSE, col = "red", xlim = c(-5, 5) ,ylim=c(0, .5),
     xlab = "", main = "", add = TRUE)
text(-4.5, .4, "c", cex=2)
dev.off()

summary(obsmeans)
summary(wrfmeans)
summary(wrfadjust_means)
dat = rbind(summary(obsmeans), summary(wrfmeans), summary(wrfadjust_means))
colnames(dat) = names(summary(obsmeans))
rownames(dat)  =c("Obs", "WRF", "WRF_adjust")
dat
# process: detrend, take out seasonal mean and sd , constant scalar
nrow(WRF)
wrf_data = process_WRF(WRF)
# historical WRF
wrfhist = wrf_data$wrfhist

# historical and future WRF
wrfdat = wrf_data$wrfdat

# process obs: detrend using WRF trend, take out seasonal means and sd, constant scalar
obs_stuff = process_OBS(TMAX, wrfhist)
obsdat = obs_stuff$obsdat
obs_params = obs_stuff$obs_params

wrfdat$tmax_norm_VAR1_adjust = tmax_norm_VAR1_adjust

# here is where you need to detrend again! 
# (b/c seasonal variance increases slightly due to adjustment of temporal dependence)
wrfdat_new  = detrend_again(wrfdat)

# join wrf data with observed trends
wrf_all = right_join(wrfdat_new, obs_params, by = "timestep")

# wrf is already 'upscaled' to TMAX --> WRF_corr_adjust as part of detrend_again()!

wrf_all$dailyTimeStep = 1:nrow(wrf_all)

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
  sdMod = tempdf$Var_scaling_factors_wrf_new[1] # Seasonal SD of WRF data (AFTER second detrending!)
  
  a = corrMean - rawMean*(sdObs/sdMod)
  b = sdObs / sdMod
  
  corrWRF = b*tempdf$test.pred.mu + a
  
  tempdf$TMAX_CORR_TEST = corrWRF
  corr_wrf[[m]] = tempdf
}
corr_df_wrf_adjust_prime = do.call("rbind", corr_wrf)
saveRDS(corr_df_wrf_adjust_prime, "corr_df_wrf_adjust_prime.Rds")
