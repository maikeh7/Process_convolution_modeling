
library(dplyr)
library(qmap)

mydata = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")

calib_data = filter(mydata, YEAR %in% 1976:2005)
head(obs)
head(mod)
ymd$meanT = NULL
i=12
qmap_corr_df = data.frame()
qmap_corr_df_hist = data.frame()
par(mfrow=c(3,4), mar=c(0,0,0,0), oma = c(4,4,1,1))
for (i in 1:12){
  mymonth = i 
  monthdat = filter(calib_data, MONTH == mymonth)
  obs = monthdat$TMAX
  mod = monthdat$test.pred.mu
  monthtest = filter(mydata, MONTH == mymonth)
  mod_test = monthtest$test.pred.mu
  monthTF = fitQmapQUANT(obs, mod, wet.day = FALSE, qstep = 0.01)
 # mymax = max(c(obs, mod))
  #mymin = min(c(obs, mod))
 # png(filename = paste("MONTH", i , ".png", sep = ""), height = 1000, width = 1000,  res = 200)
  # plot(monthTF$par$modq, monthTF$par$fitq, pch=19, 
  #      ylim=c(-35, 35),
  #      xlim=c(-35, 35), xlab ="", ylab = "", axes=F)
  # polygon(x=c(8,8,12,12), y= c(8,12,12,8), col = "lightblue")
  # points(monthTF$par$modq, monthTF$par$fitq, pch=19)
  # box()
  # if (mymonth >=9){
  #   axis(1)
  # }
  # abline(0,1, col = "red", lwd=2)
  #abline(v=12)
  #abline(v=8)
#  dev.off()
   TMAX_qqmap = doQmapQUANT(mod_test, monthTF)
   TMAX_qqmap_historical = doQmapQUANT(monthdat$test.pred.mu, monthTF)
   monthdat$TMAX_qqmap = TMAX_qqmap_historical
   monthtest$TMAX_qqmap = TMAX_qqmap
   qmap_corr_df = rbind(qmap_corr_df, monthtest)
   qmap_corr_df_hist = rbind(qmap_corr_df_hist, monthdat)
}

?polygon
plot(c(1, 9), 1:2, type = "n")
polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
        col = c("red", "blue"),
        border = c("green", "yellow"),
        lwd = 3, lty = c("dashed", "solid"))
par(mfrow = c(1,1))
# EQM TF plot Month 12
qq_quants = data.frame(Obsquants = monthTF$par$fitq, modquants = monthTF$par$modq)

newX = seq(17.91493, 35, by = .2)
extrapX = newX + 14.88507
doQmapQUANT(c(5, 18, 22), monthTF)
png("Month12EQM.png", width = 1200, height=1200, res=200)
plot(monthTF$par$modq, monthTF$par$fitq, main = paste("MONTH: ", i), xlim = c(-37,50),
     ylim = c(-37, 50), type = "l", lwd = 2, ylab = "Observed quantiles (°C)", xlab = "Model quantiles (°C)")

lines(newX, extrapX, col = "blue", lty=2, lwd=2)
abline(0,1, col = "red")
legend("topleft", legend = c("EQM TF", "Linear extrapolation", "1:1 line"), lwd = 2, lty=c(1,2,1),
       col = c("black", "blue", "red"), cex=.75)
dev.off()

png("Month12EQM.png", width = 1200, height=1200, res=200)
plot(monthTF$par$modq, monthTF$par$fitq, main = paste("MONTH: ", i), xlim = c(10,30),
     ylim = c(10, 40), type = "l", lwd = 2, ylab = "Observed quantiles (°C)", xlab = "Model quantiles (°C)")

lines(newX, extrapX, col = "blue", lty=2, lwd=2)
abline(0,1, col = "red")
legend("topleft", legend = c("EQM TF", "Linear extrapolation", "1:1 line"), lwd = 2, lty=c(1,2,1),
       col = c("black", "blue", "red"), cex=.75)
dev.off()


# EQM TF plot Month 4
newX = seq(28.91308, 50, by = .2)
extrapX = newX + 5.486922
doQmapQUANT(c(15, 35), monthTF)
png("Month4EQM.png", width = 1200, height=1200, res=200)
plot(monthTF$par$modq, monthTF$par$fitq, main = paste("MONTH: ", i), xlim = c(-20,50),
     ylim = c(-20, 50), type = "l", lwd = 2, ylab = "Observed quantiles (°C)", xlab = "Model quantiles (°C)")

lines(newX, extrapX, col = "blue", lty=2, lwd=2)
abline(0,1, col = "red")
legend("topleft", legend = c("EQM TF", "Linear extrapolation", "1:1 line"), lwd = 2, lty=c(1,2,1),
       col = c("black", "blue", "red"), cex=.75)
dev.off()

# These are processed using Matern covariance function
wrf_all = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/wrf_all.Rds")
corr_df_wrf_het = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/corr_orig.Rds")
corr_hist = filter(corr_df_wrf_het, YEAR < 2006)
corr_df_wrf_fut = filter(corr_df_wrf_het, YEAR > 2005)
corr_df = corr_df_wrf_het %>% group_by(YEAR, MONTH, DAY) %>% summarise(meancorr = mean(TMAX_CORR_TEST))
corr_df = right_join(corr_df, ymd, by = c("MONTH", "DAY"))
corr_df = corr_df %>% arrange(YEAR, MONTH, DAY)
head(corr_df)
wrf_all_fut  =filter(wrf_all, YEAR > 2005)
wrf_all_fut = dplyr::select(wrf_all_fut, YEAR, MONTH, DAY, Long_term_trend_fut, timestep)
wrf_all_fut = arrange(wrf_all_fut,YEAR, MONTH, DAY)

wrf_all_fut$corr_WRF = corr_df$meancorr
wrf_all_fut$detrended = wrf_all_fut$corr_WRF - wrf_all_fut$Long_term_trend_fut

#include everything...all years....
corr_df_wrf_het = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/corr_orig.Rds")

wrf_all_allyears = dplyr::select(wrf_all, YEAR, MONTH, DAY, Long_term_trend_fut, timestep)
wrf_all_allyears = arrange(wrf_all_allyears, YEAR, MONTH, DAY)

corr_df = corr_df_wrf_het %>% group_by(YEAR, MONTH, DAY) %>% summarise(meancorr = mean(TMAX_CORR_TEST))

corr_df = right_join(corr_df, ymd, by = c("MONTH", "DAY"))
corr_df = arrange(corr_df) %>% arrange(YEAR, MONTH, DAY)

wrf_all_allyears$corr_WRF = corr_df$meancorr
wrf_all_allyears$detrended = wrf_all_allyears$corr_WRF - wrf_all_allyears$Long_term_trend_fut

wrf_temp = allcorr_hist
head(wrf_temp)
length(alltime)
length(allWRF)
test=rep(wrf_temp$detrended, 3)

get_trends = function(my_wrf_dat){
  wrf_temp = my_wrf_dat
  num_years = length(unique(wrf_temp$YEAR))
  
  #aftertime = rep(731:1095, 30)
  #midtime = rep(366:730, 30)
  #beforetime = rep(1:365, 30)
  mydays = -rev(0:49)
  predays = rep(mydays, num_years)
  postdays = rep(1096:1145, num_years)
  aftertime = rep(731:1095, num_years)
  midtime = rep(366:730, num_years)
  beforetime = rep(1:365, num_years)
  
  alltime = c(predays, beforetime,midtime, aftertime, postdays)
  
  preWRF = filter(wrf_temp, timestep %in% 316:365)$detrended
  postWRF = filter(wrf_temp, timestep %in% 1:50)$detrended
  
  allWRF = c(preWRF, rep(wrf_temp$detrended, 3), postWRF)
  
  het_WRF <- mleHetGP(alltime, allWRF, covtype = "Matern5_2")
  het_WRF_rebuild = rebuild(het_WRF, robust = TRUE)
  
  pred_times = matrix(1:1095, ncol = 1)
  
  WRFpreds = predict(x = pred_times, object = het_WRF_rebuild)
  Seasonal_trend_wrf_corr = rep(WRFpreds$mean[366:730], num_years)
  Var_scaling_factors_wrf_corr = rep(sqrt(WRFpreds$sd2[366:730]+ WRFpreds$nugs[366:730]), num_years)
  
  wrf_temp$Seasonal_trend_wrf_corr = Seasonal_trend_wrf_corr

  wrf_temp$Var_scaling_factors_wrf_corr = Var_scaling_factors_wrf_corr
 
  return(wrf_temp)
}

plot(histres_corr$Seasonal_trend_wrf_corr[1:365])
histres_corr = get_trends(my_wrf_dat = filter(wrf_all_allyears, YEAR < 2006))
earlyres_corr = get_trends(my_wrf_dat = filter(wrf_all_allyears, YEAR < 2036 & YEAR > 2005) )
middleres_corr = get_trends(my_wrf_dat = filter(wrf_all_allyears, YEAR < 2067 & YEAR >= 2036))
lateres_corr = get_trends(my_wrf_dat = filter(wrf_all_allyears, YEAR >= 2067))
alldat = get_trends(my_wrf_dat = wrf_all_allyears)

# seasonal scaling
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/SDCorr.png", width = 1500, height=1200,
    res=250)
plot(histres_corr$Var_scaling_factors_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(2.5, 8),
     ylab = "Seasonal SD (°C)", xlab = "DOY")
#lines(alldat$Var_scaling_factors_wrf_corr[1:365], col = "purple", type = "l", lwd=2)
lines(earlyres_corr$Var_scaling_factors_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_corr$Var_scaling_factors_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_corr$Var_scaling_factors_wrf_corr[1:365], col = "darkgreen", lwd = 2)
lines(wrf_all$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)
legend("bottomleft", legend = c("DMTA 1976-2005", "DMTA 2006-2036", "DMTA 2037-2067", 
                                "DMTA 2068-2099", "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()

png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/SDAllCorr.png", width = 1500, height=1200,
    res=250)
plot(wrf_all$Var_scaling_factors_obs[1:365], type= 'l', lwd = 2, lty=2,
     ylim = c(3, 7), ylab = "Seasonal SD (°C)",
     xlab = "DOY")
lines(alldat$Var_scaling_factors_wrf_corr[1:365], col = "purple", type = "l", lwd=2)
legend("bottomleft", legend = c("DMTA 1976-2099","DMTA 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()

# seasonal means
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/MeansCorr.png",
    width = 1500, height=1200,
    res=250)
plot(histres_corr$Seasonal_trend_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(-18, 17),
     ylab = "Seasonal mean (°C)", xlab = "DOY")
#lines(alldat$Seasonal_trend_wrf_corr[1:365], col = "purple", type = "l", lwd=2)
lines(earlyres_corr$Seasonal_trend_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_corr$Seasonal_trend_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_corr$Seasonal_trend_wrf_corr[1:365], col = "darkgreen", lwd = 2)
lines(wrf_all$Seasonal_trend_obs[1:365], lty= 2, lwd = 2)
legend("topleft", legend = c("DMTA 1976-2005", "DMTA 2006-2036", "DMTA 2037-2067", "DMTA 2068-2099", 
                             "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()

png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/AllMeansCorr.png", 
    width = 1500, height=1200,
    res=250)
plot(alldat$Seasonal_trend_wrf_corr[1:365], col = "purple", type = "l", lwd=4, ylab = "Seasonal mean (°C)",
     xlab = "DOY")
lines(wrf_all$Seasonal_trend_obs[1:365], lty= 2, lwd = 2)
legend("topleft", legend = c("DMTA 1976-2099","Obs 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()
###############################################
#now make the same plot for uncorrected data
###############################################
wrf_all_WRF= dplyr::select(wrf_all, YEAR, MONTH, DAY, meanWRF, timestep, Long_term_trend_fut, Seasonal_trend_wrf,
                           Seasonal_trend_obs, Var_scaling_factors_wrf, Var_scaling_factors_obs )

wrf_all_WRF$detrended = wrf_all_WRF$meanWRF - wrf_all_WRF$Long_term_trend_fut
histres_corrW = get_trends(my_wrf_dat = filter(wrf_all_WRF, YEAR < 2006))
earlyres_corrW = get_trends(my_wrf_dat = filter(wrf_all_WRF, YEAR < 2036 & YEAR > 2005) )
middleres_corrW = get_trends(my_wrf_dat = filter(wrf_all_WRF, YEAR < 2067 & YEAR >= 2036))
lateres_corrW = get_trends(my_wrf_dat = filter(wrf_all_WRF, YEAR >= 2067))
alldatW = get_trends(my_wrf_dat = wrf_all_WRF)

# seasonal scaling -uncorr
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/SDUncorr.png", width = 1500, height=1200,
    res=250)
plot(histres_corrW$Var_scaling_factors_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(2.5, 8),
     ylab = "Seasonal SD (°C)", xlab = "DOY")
lines(earlyres_corrW$Var_scaling_factors_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_corrW$Var_scaling_factors_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_corrW$Var_scaling_factors_wrf_corr[1:365], col = "darkgreen", lwd = 2)
lines(wrf_all_WRF$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)
legend("bottomleft", legend = c("Mod 1976-2005", "Mod 2006-2036", "Mod 2037-2067", "Mod 2068-2099",
                                "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()

png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/AllSDUncorr.png",
    width = 1500, height=1200,
    res=250)
plot(alldatW$Var_scaling_factors_wrf_corr[1:365], col = "purple", type = "l", lwd=2,
      ylab = "Seasonal SD (°C)", xlab = "DOY")
lines(wrf_all_WRF$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)
legend("bottomleft", legend = c("Mod 1976-2099","Obs 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()

# seasonal means
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/MeansUncorr.png",
    width = 1500, height=1200,
    res=250)
plot(histres_corrW$Seasonal_trend_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(-18, 17),
     ylab = "Seasonal mean (°C)", xlab = "DOY")
#lines(alldatW$Seasonal_trend_wrf_corr[1:365], col = "purple", type = "l", lwd=2)
lines(earlyres_corrW$Seasonal_trend_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_corrW$Seasonal_trend_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_corrW$Seasonal_trend_wrf_corr[1:365], col = "darkgreen", lwd = 2)
#lines(wrf_all_WRF$Seasonal_trend_wrf_corr[1:365], lty= 2, lwd = 2)
lines(wrf_all_WRF$Seasonal_trend_obs[1:365], lty=2, lwd = 2)
legend("topleft", legend = c("Mod 1976-2005", "Mod 2006-2036", "Mod 2037-2067", "Mod 2068-2099", "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()

png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/MeansAllUncorr.png",
    width = 1500, height=1200,
    res=250)
plot(alldatW$Seasonal_trend_wrf_corr[1:365], col = "purple", type = "l", lwd=2, ylab = "Seasonal mean (°C)",
     xlab = "DOY")
lines(wrf_all_WRF$Seasonal_trend_obs[1:365], lty=2, lwd = 2)
legend("topleft", legend = c("Mod 1976-2099","Obs 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()

#################
#### EQM ########
#################
wrf_all = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/wrf_all.Rds")

wrf_all = arrange(wrf_all, YEAR, MONTH, DAY)
#filter(wrf_all, dailyTimeStep == 40000)
#plot(wrf_all$Long_term_trend_fut, ylab = "TMAX (°C)", xlab = "Day", type = "l", lwd = 2, xaxt = "n")
#axis(1, at = c(0, 10000, 20000, 30000, 40000), labels = c(1976, 2003, 2030 , 2058, 2085))

eqm_all = qmap_corr_df %>% group_by(YEAR,MONTH,DAY) %>% summarise(meanEQM = mean(TMAX_qqmap))
eqm_all = right_join(eqm_all, ymd, by = c("MONTH", "DAY"))
eqm_all = arrange(eqm_all, YEAR, MONTH, DAY)
eqm_all$Long_term_trend_fut = wrf_all$Long_term_trend_fut
eqm_all$detrended = eqm_all$meanEQM - eqm_all$Long_term_trend_fut


histres_eqm = get_trends(my_wrf_dat = filter(eqm_all, YEAR < 2006))
earlyres_eqm = get_trends(my_wrf_dat = filter(eqm_all, YEAR < 2036 & YEAR > 2005) )
middleres_eqm = get_trends(my_wrf_dat = filter(eqm_all, YEAR < 2067 & YEAR >= 2036))
lateres_eqm = get_trends(my_wrf_dat = filter(eqm_all, YEAR >= 2067))
alldat_eqm = get_trends(my_wrf_dat = eqm_all)

# seasonal scaling EQM
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/SDEQM.png",
    width = 1500, height=1200,
    res=250)
plot(histres_eqm$Var_scaling_factors_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(2.5, 8),
     ylab = "Seasonal SD (°C)", xlab = "DOY")
#lines(alldat_eqm$Var_scaling_factors_wrf_corr[1:365], col = "purple", type = "l", lwd=2)
lines(earlyres_eqm$Var_scaling_factors_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_eqm$Var_scaling_factors_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_eqm$Var_scaling_factors_wrf_corr[1:365], col = "darkgreen", lwd = 2)
lines(wrf_all$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)
legend("bottomleft", legend = c("EQM 1976-2005", "EQM 2006-2037", "EQM 2037-2067", "EQM 2068-2099",
                                "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()

# scaling EQM all
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/SDAllEQM.png",
    width = 1500, height=1200,
    res=250)
plot(alldat_eqm$Var_scaling_factors_wrf_corr[1:365], col = "purple", type = "l", lwd=2, ylab = 
       "Seasonal SD (°C)", xlab = "DOY" , ylim = c(2.5, 7))
lines(wrf_all$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)

legend("bottomleft", legend = c("EQM 1976-2099","Obs 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()

# seasonal means -EQM
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/MeansEQM.png",
    width = 1500, height=1200,
    res=250)
plot(histres_eqm$Seasonal_trend_wrf_corr[1:365], lwd = 2, type = "l", ylim = c(-18, 17),
     ylab = "Seasonal mean  (°C)", xlab = "DOY")

lines(earlyres_eqm$Seasonal_trend_wrf_corr[1:365], col = "red", lwd = 2)
lines(middleres_eqm$Seasonal_trend_wrf_corr[1:365], col = "blue", lwd= 2)
lines(lateres_eqm$Seasonal_trend_wrf_corr[1:365], col = "darkgreen", lwd = 2)
#lines(wrf_all_WRF$Seasonal_trend_wrf_corr[1:365], lty= 2, lwd = 2)
lines(wrf_all$Seasonal_trend_obs[1:365], lty=2, lwd = 2)
legend("topleft", legend = c("EQM 1976-2005", "EQM 2006-2037", "EQM 2037-2067", "EQM 2068-2099",
                             "Obs 1976-2005"),
       col = c("black", "red", "blue", "darkgreen", "black"), lty=c(1,1,1,1,2), cex = .8, lwd=2)
dev.off()


png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/MeansAllEQM.png",
    width = 1500, height=1200,
    res=250)
plot(alldat$Seasonal_trend_wrf_corr[1:365], col = "purple", type = "l", lwd=2, ylab = "Seasonal mean (°C)",
     xlab = "DOY")
lines(wrf_all$Seasonal_trend_obs[1:365], lty=2, lwd = 2)
legend("topleft", legend = c("EQM 1976-2099","Obs 1976-2005"),
       col = c("purple", "black"), lty = c(1,2), cex = .8, lwd = 2)
dev.off()

library(reshape2)
hquants = quantile(historical$test.pred.mu, probs = seq(0,1, length.out = 10000))
fquants = quantile(midfut$test.pred.mu, probs = seq(0,1, length.out = 10000))
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/Trends/HistVSFutqqplot.png",
    width = 1500, height=1200,
    res=250)
plot(hquants, fquants, pch=19, cex = .8, ylim = c(-40, 40), xlab = "Historical model quantiles (°C) (1976-2005)",
     ylab = "Future model quantiles (°C) (2085-2099)")
abline(0,1,  lty=2, col = "blue")
dev.off()
historical = filter(corr_df_wrf, YEAR < 2006)
midfut = filter(corr_df_wrf, YEAR >2029)
# Boxplots 
corr_df_wrf = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/HetGP/corr_df_wrf_het.Rds")

futqmap= filter(qmap_corr_df, YEAR > 2084 )
futNew = filter(corr_df_wrf, YEAR > 2084 )
futqmap1 = dplyr::select(futqmap, YEAR, MONTH, DAY, test.pred.mu, TMAX_qqmap, ID)
#futqmap1 = melt(futqmap1, id.vars = c("YEAR", "MONTH", "DAY"))
futNew = dplyr::select(futNew, YEAR, MONTH, DAY, TMAX_CORR_TEST, ID)
corrdat = right_join(futqmap1, futNew, by = c("YEAR", "MONTH", "DAY", "ID"))
corrdat$ID = NULL
corrmelt = melt(corrdat, id.vars = c("YEAR", "MONTH", "DAY"))

levels(corrmelt$variable) = c("Mod", "EQM", "DMTA")
head(corrmelt)
for (i in 1:12){
  tdat = filter(corrmelt, MONTH == i)
  ggplot(data = tdat, aes(x = factor(DAY), y = value, fill = variable)) +
    geom_boxplot() + 
    ggtitle(paste("MONTH: ", i, "2084-2096")) +
    ylab("TMAX (°C)")
  #geom_point(alpha = .3) 
  ggsave(paste("boxplot", i, ".png", sep = ""))
}

head(histNEW)
histqmap = filter(qmap_corr_df_hist, YEAR < 2006)
histqmap1 = dplyr::select(histqmap, YEAR, MONTH, DAY, test.pred.mu, TMAX_qqmap, ID)
histNEW  =filter(corr_df_wrf, YEAR < 2006)
histNEW = dplyr::select(histNEW, YEAR, MONTH, DAY, TMAX_CORR_TEST, ID)

hist_dat = right_join(histNEW, histqmap1, by = c("YEAR", "MONTH", "DAY", "ID"))

hist_dat$ID = NULL
hist_melt = melt(hist_dat, id.vars = c("YEAR", "MONTH", "DAY"))

levels(hist_melt$variable) = c("DMTA", "Mod", "EQM")
corrmelt$variable_char = as.character(corrmelt$variable)
corrmelt$Type = "Future"
hist_melt$variable_char = as.character(hist_melt$variable)
hist_melt$Type = "Historical"
alldat= rbind(hist_melt, corrmelt)
unique(alldat$variable_char)
head(corrmelt)
alldat$variable = factor(alldat$variable, levels = c("Mod", "EQM", "DMTA"))
colnames(alldat)[4] = "Data"
library(ggplot2)
head(alldat)
unique(tdat$variable)
i=1
for (i in 1:12){
  tdat = filter(alldat, MONTH == i)
  tdat = filter(tdat, Type == "Future")
  ggplot(data = tdat, aes(x = factor(DAY), y = value)) +
    geom_boxplot(aes( fill = Data, col = Data), outlier.colour = NULL) + 
    geom_boxplot(aes( fill = Data), outlier.colour = NA) + 
    scale_fill_manual(values = c("#0072B2", "#CC79A7", "#009E73")) +
    scale_colour_manual(values = c("#0072B2", "#CC79A7", "#009E73")) + 
    ggtitle(paste("MONTH ", i)) +
    ylab("TMAX (°C)") +
    xlab("DAY") +
    theme_minimal()
    ggsave(paste("boxplotFuture", i, ".png", sep = ""), width = 40, height=15, units="cm")
   
  tdat = filter(alldat, MONTH == i)
  tdat = filter(tdat, Type == "Historical")
  ggplot(data = tdat, aes(x = factor(DAY), y = value)) +
    geom_boxplot(aes( fill = Data, col = Data), outlier.colour = NULL) + 
    geom_boxplot(aes( fill = Data), outlier.colour = NA) + 
    scale_fill_manual(values = c("#0072B2", "#CC79A7", "#009E73")) +
    scale_colour_manual(values = c("#0072B2", "#CC79A7", "#009E73")) + 
    ggtitle(paste("MONTH ", i)) +
    ylab("TMAX (°C)") +
    xlab("DAY") +
    theme_minimal()
  
  #facet_wrap(~Type)
  #geom_point(alpha = .3) 
  ggsave(paste("boxplotHistorical", i, ".png", sep = ""), width = 40, height=15, units="cm")
}

#add observed data in case that might be nice
histqmap = filter(qmap_corr_df_hist, YEAR < 2006)
histqmap1 = dplyr::select(histqmap, YEAR, MONTH, DAY, test.pred.mu, TMAX_qqmap, ID)
histNEW  =filter(corr_df_wrf, YEAR < 2006)
histNEW = dplyr::select(histNEW, YEAR, MONTH, DAY, TMAX , TMAX_CORR_TEST, ID)

hist_dat = right_join(histNEW, histqmap1, by = c("YEAR", "MONTH", "DAY", "ID"))
hist_dat$ID = NULL
hist_melt = melt(hist_dat, id.vars = c("YEAR", "MONTH", "DAY"))
head(hist_melt)
levels(hist_melt$variable) = c("Obs", "DMTA", "Mod", "EQM")
hist_melt$variable = factor(hist_melt$variable, levels = c("Obs", "Mod", "EQM", "DMTA"))
head(hist_melt)
colnames(hist_melt)[4] = "Data"
for (i in 1:12){
  tdat = filter(hist_melt, MONTH == i)
 
  ggplot(data = tdat, aes(x = factor(DAY), y = value)) +
    geom_boxplot(aes( fill = Data, col = Data), outlier.colour = NULL) + 
    geom_boxplot(aes( fill = Data), outlier.colour = NA) + 
    scale_fill_manual(values = c("gray", "#0072B2", "#CC79A7", "#009E73")) +
    scale_colour_manual(values = c("gray", "#0072B2", "#CC79A7", "#009E73")) + 
    ggtitle(paste("MONTH ", i)) +
    ylab("TMAX (°C)") +
    xlab("DAY") +
    theme_minimal()
  
  #facet_wrap(~Type)
  #geom_point(alpha = .3) 
  ggsave(paste("boxplotHistoricalWithObs", i, ".png", sep = ""), width = 40, height=15, units="cm")
}


head(wrf_all)

hw = filter(wrf_all, YEAR > 2089)
mean(hw$meanWRF)







#######################
## Just plots of means / sds
#######################
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pie(rep(1, 8), col = colorBlindBlack8)

# get WRF for all years, get spatial average; also get WRF at station level
WRF = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF = WRF %>% filter(YEAR < 2006) %>% group_by(MONTH, DAY) %>% summarise(meanWRF = mean(test.pred.mu), meanTMAX = mean(TMAX))
head(WRF)
WRF$DOY = 1:365
library(reshape2)
meltdat = melt(WRF, id.vars = c("MONTH", "DAY", "DOY"))
colnames(meltdat)[4:5] = c("Data", "Mean")
levels(meltdat$Data) = c("Mod", "Obs")
str(meltdat)
levels(meltdat$Data)
ggplot(data = meltdat, aes(x = DOY, y = Mean, col = Data)) + 
  geom_line(size = .9) +
  scale_color_manual(values = c("#0072B2","#E69F00")) +
  theme_minimal() + 
  ylab("Mean TMAX  (°C)")

ggsave("WRFandObsMeansoverDOY.png")

WRF = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF = WRF %>% filter(YEAR < 2006) %>% group_by(MONTH, DAY) %>% summarise(sdWRF = sd(test.pred.mu), sdTMAX = sd(TMAX))
WRF$DOY = 1:365
meltdat = melt(WRF, id.vars = c("MONTH", "DAY", "DOY"))
colnames(meltdat)[4:5] = c("Data", "SD")
levels(meltdat$Data) = c("Mod", "Obs")
str(meltdat)
levels(meltdat$Data)
ggplot(data = meltdat, aes(x = DOY, y = SD, col = Data)) + 
  geom_line(size = .9) +
  scale_color_manual(values = c("#0072B2","#E69F00")) +
  theme_minimal() + 
  ylab("SD (°C)")

ggsave("WRFandObsSDoverDOY.png")



plot(wrf_all$Var_scaling_factors_wrf)
histdat = filter(wrf_all, YEAR < 2006)
early = filter(wrf_all, YEAR < 2036 & YEAR > 2005) 
middle = filter(wrf_all, YEAR < 2067 & YEAR >= 2036) 
late = filter(wrf_all, YEAR >= 2067)


histres = seasonal_stuff(histdat)
earlyres = seasonal_stuff(early)
middleres = seasonal_stuff(middle)
lateres = seasonal_stuff(late)

histres = get_seasonal_trend(histdat)
earlyres = get_seasonal_trend(early)
middleres = get_seasonal_trend(middle)
lateres = get_seasonal_trend(late)


plot(wrf_all$Var_scaling_factors_wrf[1:365], type = "l" , col = "purple", ylim = c(2, 6), lwd = 2, ylab = 
       "Seasonal Scaling", main = "Raw WRF")
lines(histres$Var_scaling_factors_wrf[1:365], lwd = 2)
lines(wrf_all$Var_scaling_factors_obs[1:365], lty= 2, lwd = 2)
lines(earlyres$Var_scaling_factors_wrf[1:365], col = "red", lwd = 2)
lines(middleres$Var_scaling_factors_wrf[1:365], col = "blue", lwd= 2)
lines(lateres$Var_scaling_factors_wrf[1:365], col = "green", lwd = 2)
plot(wrfdat2$Var_scaling_factors_wrf[1:365], col = "orange", lty =3) #see lines 209+ for wrfdat2
lines(wrf_all$Var_scaling_factors_obs[1:365], col = 'red')
legend("bottomleft", legend = c("WRF < 2006", "WRF > 2005 and <2036", "WRF >= 2036 and < 2076", "WRF > 2067",
                                "ALL Years"),
       col = c("black", "red", "blue", "green", "purple"), lty = c(1,1,1,1, 1), cex = .8)

plot(wrf_all$Seasonal_trend_wrf[1:365], type = "l" , col = "purple",  lwd = 2, ylab = 
       "Seasonal Means")
lines(histres$Seasonal_trend_wrf[1:365], lwd = 2)
lines(wrf_all$Seasonal_trend_obs[1:365], lty= 2, lwd = 2)
lines(earlyres$Seasonal_trend_wrf[1:365], col = "red", lwd = 2)
lines(middleres$Seasonal_trend_wrf[1:365], col = "blue", lwd= 2)
lines(lateres$Seasonal_trend_wrf[1:365], col = "green", lwd = 2)
legend("topleft", legend = c("WRF < 2006", "WRF > 2005 and <2036", "WRF >= 2036 and < 2076", "WRF > 2067",
                             "Obs", "WRF all years"), col = c("black", "red", "blue", "green", "black",
                                                              "purple"), lty = c(1,1,1,1,2,1), cex = .8)
library(viridis)
library(reshape2)
viridiscols = viridis(20)
mydat = data.frame(Historical_model = histres$Seasonal_trend_wrf[1:365],
                   Observed = wrf_all$Seasonal_trend_obs[1:365],
                   Model_2006_2036 = earlyres$Seasonal_trend_wrf[1:365],
                   Model_2037_2066 = middleres$Seasonal_trend_wrf[1:365],
                   Model_2067_2099 = lateres$Seasonal_trend_wrf[1:365],
                   DOY  = 1:365)
meltdat = melt(mydat, id.vars = "DOY")
names(meltdat) = c("DOY", "Data", "Seasonal_trend")
ggplot(data = meltdat, aes(x = DOY, y = Seasonal_trend)) +
  geom_line(aes(col = Data), size=1) +
  scale_color_manual(values = c("#440154FF" ,"#35B779FF","#88CCEE", "#CC6677", "#888888"))


mydat = data.frame(Historical_model = histres$Var_scaling_factors_wrf[1:365],
                   Observed = wrf_all$Var_scaling_factors_obs[1:365],
                   Model_2006_2036 = earlyres$Var_scaling_factors_wrf[1:365],
                   Model_2037_2066 = middleres$Var_scaling_factors_wrf[1:365],
                   Model_2067_2099 = lateres$Var_scaling_factors_wrf[1:365], 
                   Corrected_model = wrfdat2$Var_scaling_factors_wrf[1:365],
                   DOY  = 1:365)
meltdat = melt(mydat, id.vars = "DOY")
names(meltdat) = c("DOY", "Data", "Seasonal_SD")
ggplot(data = meltdat, aes(x = DOY, y = Seasonal_SD)) +
  geom_line(aes(col = Data), size=1) +
  scale_color_manual(values = c("#440154FF" ,"#35B779FF","#88CCEE", "#CC6677", "#888888", "#B8DE29FF"))




head(mydat)

start_year = 1976
end_year = 2005
trend_list = list()
i=1
while (end_year < 2100){
  tempdat = filter(wrf_all, YEAR <= end_year & YEAR >= start_year)
  
  newdat = seasonal_stuff(tempdat)
  yeardat = data.frame(YEAR = start_year, Var_scaling_factors_wrf = newdat$Var_scaling_factors_wrf[1:365])
  trend_list[[i]] = yeardat
  i=i+1
  start_year = start_year + 1
  end_year = end_year + 1
}
all_trend_list = do.call("rbind", trend_list)
range(all_trend_list$YEAR)
plot(all_trend_list$Var_scaling_factors_wrf[1:365], ylim = c(2, 6), type = "l",col = mycols[1])
mycols = rainbow(130)
lcol = 1

for (i in 1976:2070){
  lcol = lcol +1
  tdat = filter(all_trend_list, YEAR == i)
  lines(tdat$Var_scaling_factors_wrf, col = mycols[lcol])
}
legend("bottomleft", legend = c('early', "late"), col = c("red", "purple"), lwd=2)

get_seasonal_trend = function(wrfdat){
  
  # mean_season = data.frame(Seasonal_trend_wrf = wrf_all$Seasonal_trend_wrf[1:365], timestep = 1:365)
  ks=20
  wrfdat=dplyr::select(wrfdat, YEAR, MONTH, DAY, meanWRF, timestep, Long_term_trend_fut,
                       detrended)
  # seasonal trend
  mylm = gam(detrended ~  s(timestep, bs = 'cc', k =  ks), data = wrfdat)
  
  resids = residuals(mylm)
  
  seasonal_fits = mylm$fitted.values 
  
  wrfdat$Seasonal_trend_wrf = seasonal_fits
 
  return(wrfdat)
}




seasonal_stuff = function(wrfdat){
 # mean_season = data.frame(Seasonal_trend_wrf = wrf_all$Seasonal_trend_wrf[1:365], timestep = 1:365)
  ks=20
  wrfdat=dplyr::select(wrfdat, YEAR, MONTH, DAY, meanWRF, timestep, Long_term_trend_fut,
                       detrended,Seasonal_trend_wrf)
  
  resids = wrfdat$detrended - wrfdat$Seasonal_trend_wrf
  
  wrfdat$resids = resids
  
  #model seasonal "variance"
  dftmax = data.frame(absResids = abs(wrfdat$resids), timestep = wrfdat$timestep)
  
  fit2resids = gam(absResids ~ s(timestep, bs = "cc", k = 6), data = dftmax)
  
  Var_scaling_factors_wrf = fit2resids$fitted.values
  
  wrfdat$Var_scaling_factors_wrf = Var_scaling_factors_wrf
  return(wrfdat)
}

wrfdat2 = dplyr::select(wrf_all, YEAR, MONTH, DAY, meanWRF, timestep, Long_term_trend_fut, Seasonal_trend_obs,
                        WRF_corr)

wrfdat2$detrended = wrfdat2$WRF_corr - wrfdat2$Long_term_trend_fut
resids = wrfdat2$detrended - wrfdat2$Seasonal_trend_obs
wrfdat2$resids = resids

#model seasonal "variance"
dftmax = data.frame(absResids = abs(wrfdat2$resids), timestep = wrfdat2$timestep)

fit2resids = gam(absResids ~ s(timestep, bs = "cc", k = 6), data = dftmax)

Var_scaling_factors_wrf = fit2resids$fitted.values

wrfdat2$Var_scaling_factors_wrf = Var_scaling_factors_wrf




i=2
quant_df = data.frame()
for (i in 1:12){
  mymonth = i
  NewMonth = filter(futNew, MONTH == i)
  NewEQM = filter(futqmap, MONTH == i)
  historical_new = filter(histNEW, MONTH == i)
  futWRF_Hiquants = quantile(NewEQM$test.pred.mu, probs = seq(0.9, 1, length.out = 20))
  futqmap_Hiquants = quantile(NewEQM$TMAX_qqmap, probs = seq(.9, 1, length.out = 20))
  futNEW_Hiquants = quantile(NewMonth$TMAX_CORR_TEST, probs = seq(.9, 1, length.out = 20))
 # hist_Hiquants = quantile(historical_new$TMAX_CORR_TEST, probs = seq(.9, 1, length.out = 20))
#  Hist_ModHiquants = quantile(historical_new$test.pred.mu, probs = seq(.9, 1, length.out = 20))
  newdat = data.frame(Mod = futWRF_Hiquants, EQM = futqmap_Hiquants, NEW = futNEW_Hiquants, 
                     MONTH = mymonth)
  quant_df = rbind(quant_df, newdat)
}
eqmPchange = ((quant_df$EQM - quant_df$Mod )/ quant_df$Mod )* 100

NewPchange = ((quant_df$NEW - quant_df$Mod )/ quant_df$Mod )* 100

mydat = data.frame(MONTH = quant_df$MONTH, EQM_change = eqmPchange, NEW_change = NewPchange)
meltdat = melt(mydat, id.vars = "MONTH")
head(meltdat)
ggplot(data = meltdat, aes(x = factor(MONTH), y = value, fill=variable)) +
  geom_boxplot()
EQMdiff = quant_df$EQM - quant_df$Mod
NEWdiff = quant_df$NEW - quant_df$Mod

head(melt_quant)
melt_quant = melt(quant_df, id.vars = "MONTH")
ggplot(data = melt_quant, aes(x = factor(MONTH), y = value, fill=variable)) +
  geom_boxplot()

diff_df = data.frame(EQMdiff = EQMdiff, NEWdiff = NEWdiff, MONTH = quant_df$MONTH)
melt_diff = melt(diff_df, id.vars = "MONTH")
ggplot(data = melt_diff, aes(x = factor(MONTH), y = value, fill=variable)) +
  geom_boxplot() +
  ylab("Difference") +
  xlab("MONTH")


WRFhist = filter(mydata, YEAR < 2006)
sddat = WRFhist %>% group_by(MONTH, DAY) %>% summarise(SDobs = sd(TMAX),
                                                       sdMod = sd(test.pred.mu))

library(ggplot2)
library(reshape2)
viridiscols = viridis(10)
#"gray", "#0072B2", "#CC79A7", "#009E73"
viridiscols
colnames(sddat) = c("MONTH", "DAY", "Obs", "Mod")
meltdat = melt(sddat, id.vars = c("MONTH", "DAY"))
meltdat$DOY = rep(1:365, 2)
colnames(meltdat)[3] = "Data"
colnames(meltdat)[4] = "SD"
meltdat$Data = as.character(meltdat$Data)
sdplot=ggplot(data = meltdat, aes(x = DOY, y=SD,)) + 
  geom_line(aes(col = Data), size = .8) +
  scale_color_manual(values = c("#0072B2" ,"darkgray")) +
  ylab("SD (°C)") +
  theme_minimal()

meandat = WRFhist %>% group_by(MONTH, DAY) %>% summarise(meanobs = mean(TMAX),
                                                       meanMod = mean(test.pred.mu))

colnames(meandat) = c("MONTH", "DAY", "Obs", "Mod")
meltdat = melt(meandat, id.vars = c("MONTH", "DAY"))
meltdat$DOY = rep(1:365, 2)

colnames(meltdat)[3] = "Data"
colnames(meltdat)[4] = "Mean"
meanplot=ggplot(data = meltdat, aes(x = DOY, y=Mean,)) + 
  geom_line(aes(col = Data), size = .8) +
  scale_color_manual(values = c("#0072B2" ,"darkgray"))+
  ylab("TMAX (°C)") +
  theme_minimal()


library(ggpubr)
library(reshape2)
ggarrange(meanplot, sdplot, 
          common.legend = TRUE,
          ncol = 2, legend = "bottom")
ggsave("WRFMeanSDBias.png")
plot(Seasonal_trend_wrf[1:365], type = "l" , xlab = "DOY", ylab = "Mean (°C)", lwd=2, main = "WRF model seasonal mean")
plot(Var_scaling_factors_wrf[1:365], type = "l" , xlab = "DOY", ylab = "SD (°C)", lwd=2, main = "WRF model seasonal SD")
plot(Seasonal_trend_obs[1:365], type = "l" , xlab = "DOY", ylab = "Mean (°C)", lwd=2, main = "Observed seasonal mean")

plot(Var_scaling_factors_obs[1:365], type = "l" , xlab = "DOY", ylab = "SD (°C)", lwd=2, main = "Observed seasonal SD")


meandat = data.frame(Mean_WRF = Seasonal_trend_wrf[1:365], Mean_Obs = Seasonal_trend_obs[1:365],
                     DOY = 1:365)
meanmelt  = melt(meandat, id.vars = "DOY")
colnames(meanmelt)[2] = "Data"
colnames(meanmelt)[3] = "Mean"
levels(meanmelt$Data) = c("Mod", "Obs")
meanmelt$Data <- factor(meanmelt$Data, levels = c("Obs", "Mod"))

SDdat = data.frame(SD_WRF = Var_scaling_factors_wrf[1:365], SD_Obs = Var_scaling_factors_obs[1:365],
                   DOY = 1:365)
sdmelt  = melt(SDdat, id.vars = "DOY")
colnames(sdmelt)[2] = "Data"
colnames(sdmelt)[3] = "SD"
levels(sdmelt$Data) = c("Mod", "Obs")
sdmelt$Data <- factor(sdmelt$Data, levels = c("Obs", "Mod"))



sdplot=ggplot(data = sdmelt, aes(x = DOY, y=SD,)) + 
  geom_line(aes(col = Data), size = 1.2) +
  scale_color_manual(values = c("#0072B2" ,"darkgray")) +
  ylab("SD (°C)") +
  theme_minimal()

meanplot=ggplot(data = meanmelt, aes(x = DOY, y=Mean,)) + 
  geom_line(aes(col = Data), size = 1) +
  scale_color_manual(values = c("#0072B2" ,"darkgray"))+
  ylab("TMAX (°C)") +
  theme_minimal()

ggarrange(meanplot, sdplot, 
          common.legend = TRUE,
          ncol = 2, legend = "bottom")
ggsave("HetGPMeanSD.png")

plot(qmapmeansFut$abs_change, type = "l")
lines(qmapmeansFut$orig_abs_change, col = "red")


corr_hist = filter(corr_df_wrf, YEAR < 2006)
corr_fut = filter(corr_df_wrf, YEAR > 2005)
corrmeans = corr_hist %>% group_by(MONTH, DAY) %>% summarise(meanTMAX = mean(TMAX), meanMod = mean(test.pred.mu),
                                                             meanNEW = mean(TMAX_CORR_TEST))
corrmeans$abs_change = abs(corrmeans$meanNEW - corrmeans$meanMod)

corrmeansFut = corr_fut %>% group_by(MONTH, DAY) %>% summarise(meanMod = mean(test.pred.mu),
                                                               meanNEW = mean(TMAX_CORR_TEST))

corrmeansFut$abs_change = abs(corrmeansFut$meanNEW - corrmeansFut$meanMod)

plot(corrmeansFut$abs_change, type = "l")
lines(corrmeans$abs_change, type = "l", col = "red")
par(mfrow = c(1,1))





futqmap[which.max(futqmap$test.pred.mu), ]
par(mfrow = c(1,1))
range(filter(histqmap, MONTH == 6)$TMAX)
range(filter(histqmap, MONTH == 6)$test.pred.mu)
range(filter(futqmap, MONTH == 6)$test.pred.mu)
range(filter(futqmap, MONTH == 6)$TMAX_qqmap)
head(histqmap)
df1 = histqmap %>% group_by(YEAR, MONTH, DAY) %>% summarise(meanEQM = mean(TMAX_qqmap))
df1 = right_join(df1, ymd, by = c("MONTH", "DAY"))


mytrends = filter(wrf_all, YEAR < 2006)$Long_term_trend_fut
df1$Long_term_trend_fut = mytrends

ks = 18
kv = 6
head(df1)
# take WRF trend out of observed data
df1$detrended = df1$meanEQM - df1$Long_term_trend_fut

# seasonal trend
mylm = gam(detrended ~  s(timestep, bs = 'cc', k = 18), data = df1)

resids = residuals(mylm)

Seasonal_trend = mylm$fitted.values 

df1$resids = resids

df1$Seasonal_trend = Seasonal_trend

#model seasonal "variance"
dftmax = data.frame(absResids = abs(df1$resids), timestep = df1$timestep)

fit2resids = gam(absResids ~ s(timestep, bs = "cc", k = 6), data = dftmax)

Var_scaling_factors = fit2resids$fitted.values

df1$Var_scaling_factors = Var_scaling_factors

df1$tmax_norm = df1$resids / df1$Var_scaling_factors
Constant_scalar = mean(movMAD(df1$tmax_norm, window = 100)) 
df1$tmax_norm_VAR1 = df1$tmax_norm / Constant_scalar
df1$Constant_scalar = Constant_scalar
plot(df1$tmax_norm_VAR1)
write.csv(df1, "EQM_ALL.csv")

stations = filter(WRF_stationlevel, YEAR < 2006)
monthdat = filter()