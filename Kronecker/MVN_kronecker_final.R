
# this function efficiently calculates quadratic term of MVN w/ kronecker'd cov matrix
Fast_kronecker_quadratic=function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance){
  p=nrow(C)
  k=nrow(R)
  
  Zmat = matrix(z, nrow = k, ncol = p)
  
  Myvec = as.vector(Ur_t %*% Zmat %*% Uc) 

  eigenvals = marginal_variance * as.vector((outer(Sr, Sc))) + nugget #ok this works for sure
  
  
  InvEigenMat = diag((1 / eigenvals), nrow = p*k, ncol = p*k)
  
  ssqKronecker = (t(Myvec) %*% InvEigenMat %*% Myvec) 
  LogDet = sum(log(eigenvals)) 
  
  return(list(quadratic = ssqKronecker, LogDet = LogDet))
}


LogLik_MVN <- function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance, Temporal_nugget) {
  p=nrow(C) #temporal
  k=nrow(R) #spatial
  n = p*k #total length of z vec
  
  #calculates the quadratic term
  kronResult = Fast_kronecker_quadratic(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance) 
  ssqKronecker = kronResult$quadratic 
  
  LogDet = kronResult$LogDet
  #return log likelihood
  loglike = - .5*LogDet -n/2*log(2 *pi) - .5*ssqKronecker
  
  #SigSquared = 1000000 #penalty term to keep nuggets within correct range. But did not seem to work well!
  # loglike = (-n/2) *log(2 *pi) - .5*LogDet - .5*ssqKronecker - 
  #    .5*SigSquared*(.285^2 - marginal_variance - nugget - Temporal_nugget)^2 -
  #    .5*SigSquared*(0.0807932 - nugget - Temporal_nugget)^2
  
  
  return(loglike)
}

Find_scale_params = function(parameters, distC, distR, z){
  # period for covariance function
  P = 365
  
  #ell controls temporal scale, smaller = less dependence
  ell = parameters[1]
  
  #spatial scale
  scaleR = parameters[2]
  #scaleR = 120
  nugget = parameters[3] 
#  marginal_variance = parameters[3] #try specifying marg var first
  Temporal_nugget = parameters[4]
  marginal_variance = 0.15
  
  # this is the covariance function for the temporal part
  C = exp(-2*ell^2 * (sin(pi * distC/P))^2 ) + diag(Temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
  
  # eigendecomposition of cov matrix
  svdC = svd(C)
  Uc = svdC$u
  Sc = svdC$d
  Uc_t = t(svdC$v)
  
  # spatial covariance matrix (Gaussian) and eigendecomposition
  R = exp(-(distR / scaleR)^2) 
  svdR = svd(R)
  Ur = svdR$u
  Sr = svdR$d
  Ur_t = t(svdR$v)
  
  MylogLik = LogLik_MVN(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance, Temporal_nugget)
  # print(MylogLik)
  return(-MylogLik) #return negative log likelihood
}

# Calculate spatial distance matrix
# will give distance in km
library(geosphere)
Make_SigmaNN = function(stationcoords){
  mat = distm(stationcoords[, c("LON", "LAT")])/1000
  return(mat)
}

# calculate cross-distance matrix of 1km grid coords to station coords
Make_SigmaGrid_N = function(gridcoords, stationcoords){
  mat = distm(gridcoords[, c("LON", "LAT")], coords[, c("LON", "LAT")]) / 1000
  return(mat)
}

#######################################################################################################

# Here is how you get whitened residuals
# you first need to get station level corrected data using 'Final_correction.R'
corr_df_wrf = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/corr_df_wrf_adjust.Rds")
corr_df_wrf_adjust =readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/corr_df_wrf.Rds")
corr_df_hist = filter(corr_df_wrf_adjust, YEAR < 2006)
corr_df_wrf_raw = filter(corr_df_wrf, YEAR < 2006)
head(corr_df_hist)
TMAX_white = ((corr_df_hist$TMAX- corr_df_hist$Long_term_trend_fut - corr_df_hist$Seasonal_trend_obs) / corr_df_hist$Var_scaling_factors_obs ) / corr_df_hist$Constant_scalar_obs[1]
corr_df_hist$TMAX_white = TMAX_white


TMAX_DMTA_white = ((corr_df_hist$TMAX_CORR_TEST -  corr_df_hist$Long_term_trend_fut - corr_df_hist$Seasonal_trend_obs) / corr_df_hist$Var_scaling_factors_obs ) / corr_df_hist$Constant_scalar_obs[1]
corr_df_hist$TMAX_DMTA_white = TMAX_DMTA_white
head(corr_df_hist)
# get whitened residuals
corr_df_hist$Delta_TMAX_white = corr_df_hist$TMAX_DMTA_white - corr_df_hist$TMAX_white

## for non-PC-adjusted WRF
# so this WRF has been corrected but wrf x's have NOT been adjusted
TMAX_white = ((corr_df_wrf_raw$TMAX - corr_df_wrf_raw$Long_term_trend_fut - corr_df_wrf_raw$Seasonal_trend_obs) / corr_df_wrf_raw$Var_scaling_factors_obs ) / corr_df_wrf_raw$Constant_scalar_obs[1]
corr_df_wrf_raw$TMAX_white = TMAX_white

TMAX_WRF_white_noadjust = ((corr_df_wrf_raw$TMAX_CORR_TEST -  corr_df_wrf_raw$Long_term_trend_fut - corr_df_wrf_raw$Seasonal_trend_obs) / corr_df_wrf_raw$Var_scaling_factors_obs ) / corr_df_wrf_raw$Constant_scalar_obs[1]
corr_df_wrf_raw$TMAX_WRF_white_noadjust = TMAX_WRF_white_noadjust

# get whitened residuals
corr_df_wrf_raw$Delta_TMAX_white_noadjust = corr_df_wrf_raw$TMAX_WRF_white_noadjust - corr_df_wrf_raw$TMAX_white


# Set up data
# Please see covariance functions / other important functions in simulateCGA.R
dfall = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/newDF_with_white_REVISED_3pt5.Rds")
# Delta_TMAX_white is: Whitened corrected WRF - whitened station data
monthAves = dfall %>% group_by(MONTH,DAY, ID) %>% summarise(MeanDeltaWhite = mean(Delta_TMAX_white)) %>% 
  arrange(MONTH, DAY, ID)

# get temporal distance matrix
sigmaT = make_temporal_distmat(daysInYear = 365)

#make distance matrix for stations
library(geosphere)
newdfaves = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/newdfaves.csv")
newdfaves = newdfaves %>% arrange(MONTH, DAY, ID)
stations = newdfaves[!(duplicated(newdfaves$ID)), ] %>% arrange(ID)
coords = stations[, c("LON", "LAT")]
SigmaNN = Make_SigmaNN(coords)

# make cross-covariance matrix of 1km grid points to stations
gridcoords = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/latlonGrid.csv")
SigmaGridN = Make_SigmaGrid_N(gridcoords, coords)




# need this for optimization
library(dfoptim)

# you will need this as input to optimization
distC = sigmaT
distR = SigmaNN
z = monthAves$MeanDeltaWhite


#ell controls temporal scale, smaller = less dependence
ell = parameters[1]

#spatial scale
scaleR = parameters[2]
#scaleR = 120
nugget = parameters[3] 
#  marginal_variance = parameters[3] #try specifying marg var first
Temporal_nugget = parameters[4]



# ell = parameters[1]
# scaleR = parameters[2]
# nugget = parameters[3]
# marginal_variance = parameters[4] #but let's set this to the estimate = (.285)^2 = 0.081
# Temporal_nugget = parameters[5]

#0.672604743 47.258957106  0.007379727  0.999999993  0.013868505
testmkb = nmkb(c(2, 100, .01, .001), Find_scale_params,distC = distC, distR=distR, z=z,
               lower = c(.2, 10, 0.00001, .00001), upper = c(100, 1000, .08, .08))


