library(ggmap)
library(dplyr)
setwd("/data/mholthui/climate_analysis/Chapter2/TMAX")

WRF_stationlevel = readRDS("RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF_stationlevel = filter(WRF_stationlevel, YEAR < 2006)

register_google(key = "AIzaSyCjdQJewMV9u-MyVQGlsuul69KLvd9Wvjo")
mymap = get_map(c(lon = -73.14653, lat = 44.94444), zoom = 7, maptype = "toner", source = "stamen")

basemap = ggmap(mymap)


WRF_stationlevel$mydate = paste(WRF_stationlevel$YEAR, WRF_stationlevel$MONTH, 
                                WRF_stationlevel$DAY, sep = "-")
head(WRF_stationlevel)
# the first day of each months for random years
# only end up using months 1-6
#dates = c("1976-1-1", "1980-2-1", "1982-3-1" , "1984-4-1", "1985-5-1", "1987-6-1", "1989-7-1", 
#          "1990-8-1", "1992-9-1", "1994-10-1", "1996-11-1", "2000-12-1")

dates = c("2000-1-1", "2000-2-1", "2000-3-1" , "2000-4-1", "2000-5-1", "2000-6-1", "2000-7-1", 
          "2000-8-1", "2000-9-1", "2000-10-1", "2000-11-1", "2000-12-1")
test1 = as.POSIXct(dates)
dates_posix = strftime(test1, format='%Y-%b-%d')
dates_posix
#dates_1 = dates[1:6]
dates_1 = dates[c(1, 3, 5, 7, 9, 11)]
#dates_2 = dates[c(2, 4, 6, 8, 10, 12)]
tdat = filter(WRF_stationlevel, mydate %in% dates_1)
#tdat = filter(WRF_stationlevel, dailyTimeStep %in% days)
#days = unique(tdat$dailyTimeStep)
date_df = data.frame(mydate = dates_1, dates_posix = dates_posix[c(1, 3, 5, 7, 9, 11)])
#date_df = data.frame(mydate = dates_1[c(2, 4, 6, 8, 10, 12)], dates_posix = dates_posix[c(2, 4, 6, 8, 10, 12)])
date_df$mydate = as.character(date_df$mydate)
date_df
mymin = -19
mymax = 33


library(grid)
library(gridExtra)
tdat = filter(WRF_stationlevel, mydate %in% dates_1)

tdat_sub = tdat %>% dplyr::select(mydate, TMAX, test.pred.mu, LON, LAT)
head(tdat_sub)
colnames(tdat_sub)[3] = "Model"
colnames(tdat_sub)[2] = "Observed"

tdat_melt =melt(tdat_sub, id.vars = c("mydate", "LON", "LAT"))
head(tdat_melt)
range(tdat_melt$value)
mymin = -15
mymax = 33

testmap=basemap + 
  geom_point(data = tdat_melt, mapping = aes(x = LON, y = LAT, col = value), size = 1) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax), breaks = seq(mymin, mymax, by = 2)) +
  #theme(legend.position= "none",
  #      axis.title.x=element_blank(),
  #      axis.text.x=element_blank(),
  #      axis.ticks.x=element_blank()) +
  coord_equal() +
  theme(legend.key.width = unit(2, "cm"), legend.position = "bottom") +
  facet_grid(variable~mydate)

#use this for the plot (stations only)
#unique(tdat_melt$mydate)
Fdat=filter(tdat_melt,variable == "Observed")
Fdat = right_join(Fdat, date_df, by = "mydate")
Fdat$better_date = gsub("2000-", "", Fdat$dates_posix)
Fdat$better_date_factor <- factor(Fdat$better_date, levels = c("Jan-01", "Mar-01", "May-01", "Jul-01", "Sep-01", "Nov-01"))
levels(Fdat$better_date_factor)
head(Fdat)
str(Fdat)

testmap2=basemap + 
  geom_point(data = Fdat,
             mapping = aes(x = LON, y = LAT, col = value), size = 1) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  ylab("LAT") +
  scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax),breaks = seq(mymin, mymax, by = 2)) +
  theme(legend.position= "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), plot.margin=unit(c(1,-1,-2,1), "cm")) +
  coord_equal() +
  facet_grid(~ better_date_factor)

#get the legend from this plot
# testmap=basemap + 
#   geom_point(data = Fdat,
#              mapping = aes(x = LON, y = LAT, col = value), size = 1) +
#   ylim(c(43, 46)) + 
#   xlim(c(-75,-71)) +
#   scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax),breaks = seq(mymin, mymax, by = 2)) +
#   theme(legend.key.width = unit(2, "cm"), legend.position = "bottom", plot.margin=unit(c(1,1,-0.5,1), "cm")) +
#   ylab("LAT") +
#   coord_equal() +
#   facet_grid(~mydate) 

# for legend (rhs)
testmap=basemap + 
  geom_point(data = Fdat,
             mapping = aes(x = LON, y = LAT, col = value), size = 1) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax),breaks = seq(mymin, mymax, by = 4)) +
  theme(legend.key.height = unit(1, "cm"), legend.position = "right", plot.margin=unit(c(1,1,-0.5,1), "cm")) +
  labs(color='TMAX (°C)') +
  ylab("LAT") +
  coord_equal() +
  facet_grid(~ better_date_factor)

leg = get_legend(testmap)



# wrf grid plot
unique(allwrf$mydate)
allwrf = right_join(allwrf, date_df, by = "mydate")
allwrf$better_date = gsub("2000-", "", allwrf$dates_posix)
allwrf$better_date_factor <- factor(allwrf$better_date, levels = c("Jan-01", "Mar-01", "May-01", "Jul-01", "Sep-01", "Nov-01"))
levels(allwrf$better_date_factor)

wrfgrids = basemap +
  geom_tile(data=allwrf, aes(x=LON, y= LAT, fill=TMAX), alpha = .8) +
  scale_fill_gradientn(colors=topo.colors(40), limits = c(mymin, mymax), breaks = seq(mymin, mymax, by = 2)) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  ylab("LAT") +
  xlab("LON") +
  theme(legend.position= "none",  
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin=unit(c(-0.5,-1,1,1), "cm")) +
  coord_equal() +
  facet_wrap(~better_date_factor, nrow = 1, ncol= 6) 
# make the plot and save it  
#use this!
hlay <- rbind(c(1,1,1,1,3),
              c(1,1,1,1,3),
              c(2,2,2,2,3),
              c(2,2,2,2,3))
png("Fig2Dates2000BetterDates.png", height = 800, width = 1800, res = 200)
grid.arrange(testmap2, wrfgrids, leg, layout_matrix = hlay) 
dev.off()

library(marmap)
#now need to get gridded WRF for the same days.
make_sp_df = function(wrfdat){
  mydate = wrfdat$mydat[1]
  dfday = dplyr::select(wrfdat, LON, LAT, TMAX)
  wrfgridded = griddify(dfday, nlon = 63, nlat = 69)
  test_spdf <- as(wrfgridded, "SpatialPixelsDataFrame")
  test_df <- as.data.frame(test_spdf)
  names(test_df) = c("TMAX", "LON", "LAT")
  test_df$mydate = mydate
  return(test_df)
  
}

dates_1
wrf1 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf1 = filter(wrf1, MONTH == 1, DAY ==1)
wrf1$mydate = paste(wrf1$YEAR, wrf1$MONTH, wrf1$DAY, sep = "-")
wrf1sp = make_sp_df(wrf1)
head(wrf1sp)

wrf2 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf2 = filter(wrf2, MONTH == 3, DAY ==1)
wrf2$mydate = paste(wrf2$YEAR, wrf2$MONTH, wrf2$DAY, sep = "-")
wrf2sp = make_sp_df(wrf2)
head(wrf2sp)

wrf3 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf3 = filter(wrf3, MONTH == 5, DAY ==1)

wrf3$mydate = paste(wrf3$YEAR, wrf3$MONTH, wrf3$DAY, sep = "-")
wrf3sp = make_sp_df(wrf3)
head(wrf3sp)


wrf4 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf4 = filter(wrf4, MONTH == 7, DAY ==1)
wrf4$mydate = paste(wrf4$YEAR, wrf4$MONTH, wrf4$DAY, sep = "-")
wrf4sp = make_sp_df(wrf4)
head(wrf4sp)

wrf5 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf5 = filter(wrf5, MONTH == 9, DAY ==1)
wrf5$mydate = paste(wrf5$YEAR, wrf5$MONTH, wrf5$DAY, sep = "-")
wrf5sp = make_sp_df(wrf5)
head(wrf5sp)

wrf6 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf6 = filter(wrf6, MONTH == 11, DAY ==1)
wrf6$mydate = paste(wrf6$YEAR, wrf6$MONTH, wrf6$DAY, sep = "-")
wrf6sp = make_sp_df(wrf6)

allwrf = rbind(wrf1sp, wrf2sp, wrf3sp, wrf4sp, wrf5sp, wrf6sp)
head(allwrf)
#################### other dates
wrf1 = get_wrf_data(1989, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf1 = filter(wrf1, MONTH == 2, DAY ==1)
wrf1$mydate = paste(wrf1$YEAR, wrf1$MONTH, wrf1$DAY, sep = "-")
wrf1sp = make_sp_df(wrf1)
head(wrf1sp)

wrf2 = get_wrf_data(1990, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf2 = filter(wrf2, MONTH == 4, DAY ==1)
wrf2$mydate = paste(wrf2$YEAR, wrf2$MONTH, wrf2$DAY, sep = "-")
wrf2sp = make_sp_df(wrf2)
head(wrf2sp)

wrf3 = get_wrf_data(1992, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf3 = filter(wrf3, MONTH == 6, DAY ==1)

wrf3$mydate = paste(wrf3$YEAR, wrf3$MONTH, wrf3$DAY, sep = "-")
wrf3sp = make_sp_df(wrf3)
head(wrf3sp)


wrf4 = get_wrf_data(1994, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf4 = filter(wrf4, MONTH == 8, DAY ==1)
wrf4$mydate = paste(wrf4$YEAR, wrf4$MONTH, wrf4$DAY, sep = "-")
wrf4sp = make_sp_df(wrf4)
head(wrf4sp)

wrf5 = get_wrf_data(1996, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf5 = filter(wrf5, MONTH == 10, DAY ==1)
wrf5$mydate = paste(wrf5$YEAR, wrf5$MONTH, wrf5$DAY, sep = "-")
wrf5sp = make_sp_df(wrf5)
head(wrf5sp)

wrf6 = get_wrf_data(2000, VAR = "tmax", historical = FALSE, rcp="RCP85")
wrf6 = filter(wrf6, MONTH == 12, DAY ==1)
wrf6$mydate = paste(wrf6$YEAR, wrf6$MONTH, wrf6$DAY, sep = "-")
wrf6sp = make_sp_df(wrf6)




library(mgcv)
#PLOT OF WRF OUT TO FUTURE WITH OBS AND TRENDS
WRF_stationlevel = readRDS("RCP85_TMAX_WRFInterp2GHCND1976_2099.Rds")
WRF_stationlevel = filter(WRF_stationlevel, YEAR < 2006)
wrf_all = readRDS("wrf_all.Rds")
head(WRF_stationlevel)

wrfmeans = WRF_stationlevel %>% group_by(YEAR, MONTH, DAY) %>% summarise(meanWRF = mean(test.pred.mu))
wrfmeans$Time = 1:nrow(wrfmeans)

obsmeans = WRF_stationlevel %>% filter(YEAR < 2006) %>% group_by(YEAR, MONTH, DAY) %>% summarise(meanTMAX = mean(TMAX))
obsmeans$Time = 1:nrow(obsmeans)

obstrend = gam(meanTMAX ~ s(Time, bs = "cr", k = 10), data = obsmeans)
obsfit = obstrend$fitted.values
obsfit = data.frame(value = obsfit, dailyTimeStep = 1:length(obstrend$fitted.values), type = "Observations")

wrf_fit = gam(meanWRF ~ s(Time, k = 50), data = wrfmeans)
wrf_fit = data.frame(value = wrf_fit$fitted.values, dailyTimeStep = 1:length(wrf_fit$fitted.values), type = "Model")
all_fit = rbind(obsfit, wrf_fit)
colnames(all_fit)[2] = "Time"
colnames(all_fit)[3] = "Data"
head(all_fit)

head(wrfmeans)
colnames(wrfmeans)[4] = "Model"
colnames(obsmeans)[4] = "Observations"
meltwrf= melt(wrfmeans, id.vars = c("YEAR", "MONTH", "DAY", "Time"))
meltobs = melt(obsmeans, id.vars = c("YEAR", "MONTH", "DAY",  "Time"))
head(meltobs)
head(meltwrf)
alldat = rbind(meltobs, meltwrf)
colnames(alldat)[5] = "Data"


topocols = topo.colors(40)
col1 = topocols[1]
col2 = topocols[15]
mylabs = c(1976, 2006, 2035, 2065, 2095)
alldat$Data = factor(alldat$Data, levels = c("Observations", "Model"))
obsdat = filter(alldat, Data == "Observations")
head(obsdat)
plot1 = ggplot(data = alldat, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .5) +
  geom_line(data = all_fit, aes(x = Time, y = value, col = Data)) +
  geom_line(data = obsdat, aes(x = Time, y = value), col = col2, alpha = .5) +
  #geom_smooth(method = "loess", span = .9, se= FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 10950, 21900, 32850, 43800), labels= mylabs) +
  theme(legend.key.width = unit(2, "cm"), 
        legend.key.height = unit(0.55, "cm"),
        legend.position = c(0.87, 0.3),
        legend.background = element_rect(fill = "white", color = "black")) +
  ylab("TMAX (°C)")
#ggsave("WRFObsTimeSeries.png")

plot2 = ggplot(data = alldat, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .7) +
  geom_line(data = all_fit, aes(x = Time, y = value, col = Data)) +
  geom_line(data = obsdat, aes(x = Time, y = value), col = col2, alpha = .6) +
  #geom_smooth(method = "loess", span = .9, se= FALSE) + 
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 10950, 21900, 32850, 43800), labels= mylabs, limits = c(0, 11000)) +
  theme(legend.position = "none")+
  ylab("TMAX (°C)")
#ggsave("WRFObsTimeSeriesCALIB.png")

hlay <- rbind(c(1,1,1,1,1),
              c(1,1,1,1,1),
              c(2,2,2,2,2),
              c(2,2,2,2,2))
png("TEST.png", height = 800, width = 1800, res = 200)
grid.arrange(plot1, plot2, layout_matrix = hlay) 
dev.off()
########################################################################
#PLOT OF 5 YEARS OF OBS AND SIM AND SEASONAL EFFECT (SHOWING BIAS)
########################################################################
topocols = topo.colors(40)
col1 = topocols[1]
col2 = topocols[15]

mylabs = c(1976, 2006, 2035, 2065, 2095)
wrf_2 = dplyr::select(wrf_all, YEAR, MONTH, DAY, Var_scaling_factors_obs, Seasonal_trend_obs,
                      Var_scaling_factors_wrf, Seasonal_trend_wrf,
                      Long_term_trend_fut, WRF_detrended, meanWRF)
wrf_2 = filter(wrf_2, YEAR < 2006)
tail(wrf_2)
obsmeans = WRF_stationlevel %>% group_by(YEAR, MONTH, DAY) %>% summarise(meanTMAX = mean(TMAX)) %>%
  filter(YEAR < 2006)
obsmeans$OBS_detrended = obsmeans$meanTMAX - wrf_2$Long_term_trend_fut
range(wrf_2$WRF_detrended)

alldat = data.frame(Observations = obsmeans$OBS_detrended, Model = wrf_2$WRF_detrended, Time = 1:nrow(obsmeans), 
                    YEAR = wrf_2$YEAR, MONTH = wrf_2$MONTH, DAY = wrf_2$DAY)
Malldat = melt(alldat, id.vars = c("YEAR", "MONTH", "DAY", "Time"))
head(Malldat)
colnames(Malldat)[5] = "Data"

#PLOT OF 5 YEARS OF OBS AND SIM AND SEASONAL EFFECT (SHOWING BIAS)
wrf_all = readRDS("wrf_all.Rds")
newdat = data.frame(DOY = 1:365, Model = wrf_all$Seasonal_trend_wrf[1:365],
                    Observations = wrf_all$Seasonal_trend_obs[1:365])
new_melt = melt(newdat, id.vars = "DOY")
colnames(new_melt)[2] = "Data"

newdat_sd = data.frame(DOY = 1:365, Model = wrf_all$Var_scaling_factors_wrf[1:365],
                       Observations = wrf_all$Var_scaling_factors_obs[1:365])
new_melt_sd = melt(newdat_sd, id.vars = "DOY")
colnames(new_melt_sd)[2] = "Data"
head(new_melt_sd)

alldat_sub = filter(Malldat, YEAR < 1982)
new_melt$Time = 1:365
mylabs = c(1976, 1977, 1978, 1979, 1980, 1981)
yearplot_null = ggplot(data = alldat_sub, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .7) +
  geom_line(data = new_melt, aes(x = Time, y = value, col = Data), size=1) +
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460, 1825), labels= mylabs) +
  theme(legend.key.width = unit(3, "cm"), legend.key.height = unit(.2, "cm"), legend.position = c(.85, .22),
        legend.background = element_rect(fill = "white", color = "black")) +
  ylab("TMAX (°C)")

yearplot = ggplot(data = alldat_sub, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .7) +
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460, 1825), labels= mylabs) +
  theme(legend.position = "none") +
  ylab("TMAX (°C)")

seasPlot = ggplot(data = new_melt, aes(x = DOY, y = value, col = Data)) +
  geom_line(alpha = .9, size=1) +
  theme_minimal() +
  scale_color_manual(values = c(col1, col2)) +
  theme(legend.position = "none") +
  ylab("TMAX (°C)")

SD_plot = ggplot(data = new_melt_sd, aes(x = DOY, y = value, col = Data)) +
  geom_line(alpha = .9, size=1) +
  theme_minimal() +
  scale_color_manual(values = c(col1, col2)) +
  theme(legend.position = "none") +
  ylab("SD TMAX (°C)")

leg2 = get_legend(yearplot_null)

png("test2.png", height = 800, width = 1500, res = 200)

grid.arrange(yearplot_null, seasPlot, SD_plot, layout_matrix = hlay) 
dev.off()
hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))

####################################################################################################
png("obsmean.png", height = 800, width = 800, res = 200)
plot(wrf_2$Seasonal_trend_obs[1:365], type = "l", col = "blue", ylab = "Mean", xlab = "DOY", lwd=2, cex.lab = .8)
dev.off()
png("obsSD.png", height = 800, width = 800, res = 200)
plot(wrf_2$Var_scaling_factors_obs[1:365], type = "l", col = "blue", ylab = "SD", xlab = "DOY", lwd=2, cex.lab=.8)
dev.off()

png("WRFmean.png", height = 800, width = 800, res = 200)
plot(wrf_2$Seasonal_trend_wrf[1:365], type = "l", ylab = "Mean", xlab = "DOY", lwd=2, cex.lab = .8)
dev.off()
png("WRFSD.png", height = 800, width = 800, res = 200)
plot(wrf_2$Var_scaling_factors_wrf[1:365], type = "l", ylab = "SD", xlab = "DOY", lwd=2, cex.lab = .8)
dev.off()

#############################################################################
# junk
#############################################################################
wrf_all = readRDS("wrf_all.Rds")
newdat = data.frame(DOY = 1:365, Model = wrf_all$Seasonal_trend_wrf[1:365],
                    Observations = wrf_all$Seasonal_trend_obs[1:365])
new_melt = melt(newdat, id.vars = "DOY")
colnames(new_melt)[2] = "Data"

newdat_sd = data.frame(DOY = 1:365, Model = wrf_all$Var_scaling_factors_wrf[1:365],
                       Observations = wrf_all$Var_scaling_factors_obs[1:365])
new_melt_sd = melt(newdat_sd, id.vars = "DOY")
colnames(new_melt_sd)[2] = "Data"
head(alldat_sub)

head(alldat)
alldat_sub = filter(alldat, YEAR < 1982)
mylabs = c(1976, 1977, 1978, 1979, 1980, 1981)
yearplot_null = ggplot(data = alldat_sub, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .7) +
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460, 1825), labels= mylabs) +
  theme(legend.key.width = unit(2, "cm"), legend.position = c(.87, .1),
        legend.background = element_rect(fill = "white", color = "black")) +
  ylab("TMAX (°C)")

yearplot = ggplot(data = alldat_sub, aes(x = Time, y = value, col = Data)) +
  geom_line(alpha = .7) +
  theme_minimal() +
  scale_color_manual(values = c(col2, col1)) +
  scale_x_continuous(breaks = c(0, 365, 730, 1095, 1460, 1825), labels= mylabs) +
  theme(legend.position = "none") +
  ylab("TMAX (°C)")

seasPlot = ggplot(data = new_melt, aes(x = DOY, y = value, col = Data)) +
  geom_line(alpha = .9, size=1.5) +
  theme_minimal() +
  scale_color_manual(values = c(col1, col2)) +
  theme(legend.position = "none") +
  ylab("TMAX (°C)")

SD_plot = ggplot(data = new_melt_sd, aes(x = DOY, y = value, col = Data)) +
  geom_line(alpha = .9, size=1.5) +
  theme_minimal() +
  scale_color_manual(values = c(col1, col2)) +
  theme(legend.position = "none") +
  ylab("SD (°C)")

leg2 = get_legend(yearplot_null)

png("test2.png", height = 800, width = 1500, res = 200)

grid.arrange(yearplot_null, seasPlot, SD_plot, layout_matrix = hlay) 
dev.off()
hlay <- rbind(c(1,1,1,1),
              c(1,1,1,1),
              c(2,2,3,3),
              c(2,2,3,3))
png("test2.png", height = 800, width = 2000, res = 200)

grid.arrange(yearplot_null, seasPlot, seasPlot, SD_plot, layout_matrix = hlay) 
dev.off()
hlay <- rbind(c(1,1,1,1,2,2),
              c(1,1,1,1,2,2),
              c(3,3,3,4,4,4),
              c(3,3,3,4,4,4))
png("test2.png", height = 800, width = 2000, res = 200)

grid.arrange(yearplot, seasPlot, leg2,
             ncol=2, nrow = 2, heights=unit(c(8,1), c("cm", "cm")),
             widths=unit(c(12,9), c("cm", "cm")))
dev.off()
mymax = ceiling(max(tdat$test.pred.mu, tdat$TMAX))
mymin = floor(min(tdat$test.pred.mu, tdat$TMAX))
mymin=-19
#days = c(2, 24, 66, 91, 90, 126)
m=1
maplist = list()
mylabels = vector()
for (m in 1:6){
  myday = dates_1[m]
  tempdat= filter(WRF_stationlevel, mydate == myday)
  wrfdat = dplyr::select(tempdat, YEAR, MONTH, DAY, LAT, LON, test.pred.mu)
  obsdat = dplyr::select(tempdat, YEAR, MONTH, DAY, LAT, LON, TMAX)
  yr = obsdat$YEAR[1]
  mo = obsdat$MONTH[1]
  dy = obsdat$DAY[1]
  
  wrfmap = basemap + 
    geom_point(data = wrfdat, mapping = aes(x = LON, y = LAT, col = test.pred.mu), size = 2) +
    ylim(c(43, 46)) + 
    xlim(c(-75,-71)) +
    scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax)) +
    coord_equal() +
    theme(legend.position= "none") 
  #ggtitle(paste("Model ", yr,"/", mo, "/", dy))
  maplist[[m]] = wrfmap
  lab=paste("Model ", yr,"/", mo, "/", dy)
  mylabels[m] = lab
  
  obsmap = basemap + 
    geom_point(data = obsdat, mapping = aes(x = LON, y = LAT, col = TMAX), size = 2) +
    ylim(c(43, 46)) + 
    xlim(c(-75,-71)) +
    scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax)) +
    coord_equal() +
    theme(legend.position= "none") 
  #ggtitle(paste("Observed ", yr,"/", mo, "/", dy))
  maplist[[m+6]] = obsmap
  lab = paste("Observed ", yr,"/", mo, "/", dy)
  mylabels[m+6] = lab
}

obsmap = basemap + 
  geom_point(data = obsdat, mapping = aes(x = LON, y = LAT, col = TMAX), size = 2) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  scale_color_gradientn(colors = topo.colors(40), breaks = seq(mymin, mymax, by = 2),
                        limits=c(mymin, mymax)) +
  coord_equal() +
  theme(legend.key.width = unit(2, "cm"), legend.position = "bottom") 
#ggtitle(paste("Observed ", yr,"/", mo, "/", dy))

leg = get_legend(obsmap)

# junk
mymax = ceiling(max(tdat$test.pred.mu, tdat$TMAX))
mymin = floor(min(tdat$test.pred.mu, tdat$TMAX))
mymin=-19
#days = c(2, 24, 66, 91, 90, 126)
m=1
maplist = list()
mylabels = vector()
for (m in 1:6){
  myday = dates_1[m]
  tempdat= filter(WRF_stationlevel, mydate == myday)
  wrfdat = dplyr::select(tempdat, YEAR, MONTH, DAY, LAT, LON, test.pred.mu)
  obsdat = dplyr::select(tempdat, YEAR, MONTH, DAY, LAT, LON, TMAX)
  yr = obsdat$YEAR[1]
  mo = obsdat$MONTH[1]
  dy = obsdat$DAY[1]
  
  wrfmap = basemap + 
    geom_point(data = wrfdat, mapping = aes(x = LON, y = LAT, col = test.pred.mu), size = 2) +
    ylim(c(43, 46)) + 
    xlim(c(-75,-71)) +
    scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax)) +
    coord_equal() +
    theme(legend.position= "none") 
  #ggtitle(paste("Model ", yr,"/", mo, "/", dy))
  maplist[[m]] = wrfmap
  lab=paste("Model ", yr,"/", mo, "/", dy)
  mylabels[m] = lab
  
  obsmap = basemap + 
    geom_point(data = obsdat, mapping = aes(x = LON, y = LAT, col = TMAX), size = 2) +
    ylim(c(43, 46)) + 
    xlim(c(-75,-71)) +
    scale_color_gradientn(colors = topo.colors(40), limits=c(mymin, mymax)) +
    coord_equal() +
    theme(legend.position= "none") 
  #ggtitle(paste("Observed ", yr,"/", mo, "/", dy))
  maplist[[m+6]] = obsmap
  lab = paste("Observed ", yr,"/", mo, "/", dy)
  mylabels[m+6] = lab
}

obsmap = basemap + 
  geom_point(data = obsdat, mapping = aes(x = LON, y = LAT, col = TMAX), size = 2) +
  ylim(c(43, 46)) + 
  xlim(c(-75,-71)) +
  scale_color_gradientn(colors = topo.colors(40), breaks = seq(mymin, mymax, by = 2),
                        limits=c(mymin, mymax)) +
  coord_equal() +
  theme(legend.key.width = unit(2, "cm"), legend.position = "bottom") 
#ggtitle(paste("Observed ", yr,"/", mo, "/", dy))

leg = get_legend(obsmap)