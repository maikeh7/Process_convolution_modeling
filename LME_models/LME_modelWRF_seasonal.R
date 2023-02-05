
#install.packages("nlme")
library(nlme)
library(dplyr)
#function to construct distance matrix
doyDist <- function(t1,t2=tvals,daysInYear=365){
  # computes the distance in days between the days in the vector t1
  # and the days in the vector t2, using only the day of the year for this.
  # So Dec 30 and Jan 1 are only 2 days apart, even if they're in different
  # years.
  
  halfyear = round(daysInYear/2)
  t1c <- t1%%daysInYear
  t2c <- t2%%daysInYear
  dist1 <- outer(t1c,t2c,'-')
  #
  t1cs <- (t1+halfyear)%%daysInYear
  t2cs <- (t2+halfyear)%%daysInYear
  dist2 <- outer(t1cs,t2cs,'-')
  #
  dist3 <- ifelse(abs(dist1)<abs(dist2),dist1,dist2)
}


#################################################################################
### functions to construct seasonal, random monthly, and random yearly matrices
################################################################################
make_kmats = function(num_years=30, kernel_SD){
  day=seq(1,1825, by = 1)
  if (kernel_SD >= 90){
    print("90")
    seq1=seq(73,1825, by = 73)
    knots1 = length(seq(73, 365, by = 73)) 
  }
  else if (kernel_SD >= 60){
    print("60")
    seq1 = seq(36.5, 1825, by = 36.5)
    knots1 = length(seq(36.5, 365, by = 36.5))
  }
  else if (kernel_SD >= 30){
    print("30")
    seq1 = seq(18.25, 1825, by = 18.25)
    knots1 = length(seq(18.25, 365, by = 18.25))
  }
  else if (kernel_SD >= 14){
    print("14")
    seq1 = seq(9.125, 1825, by = 9.125)
    knots1 = length(seq(9.125, 365, by = 9.125))
  }
  else if (kernel_SD >= 7){
    print("7")
    seq1=seq(4.5625, 1825, by = 4.5625)
    knots1 = length(seq(4.5625, 365, by = 4.5625))
  }
  else if (kernel_SD >= 3.5){
    print("3.5")
    seq1=seq(2.28125, 1825, by = 2.28125)
    knots1 = length(seq(2.28125, 365, by = 2.28125))
  }
  dmat = doyDist(day, seq1,
                 daysInYear=99999)
  kmat = dnorm(dmat , sd = kernel_SD) / dnorm(0, sd = kernel_SD)
  
  testmat = matrix(NA, ncol = knots1, nrow=365)
  
  counter=1
  
  for(i in 731:1095){
    
    all1 = kmat[i, 1:(knots1)] + kmat[i, (knots1+1):(2*knots1)] +  kmat[i, (2*knots1+1):(3*knots1)] + 
      kmat[i, (3*knots1+1):(4*knots1)] + kmat[i, (4*knots1+1):(5*knots1)]
    
    testmat[counter, ] = all1
    counter=counter+1
  }
  
  int2k = max(testmat %*% t(testmat))
  newkmat = testmat / sqrt(int2k)
  
  bigmat = rep(1,num_years) %x% newkmat
  
  return(bigmat)
}


make_longterm_Kmat = function(day, kernel_SD){
  print("making long term kmat")
  boundary_extension = 2.25*kernel_SD
  dmat = doyDist(day,seq(1 - boundary_extension, length(day) + boundary_extension, by = kernel_SD),
                 daysInYear=9999999)
  kmat = dnorm(dmat , sd = kernel_SD) / dnorm(0, sd = kernel_SD)
  int2k = max(kmat%*%t(kmat))
  newkmat = kmat / sqrt(int2k)
  print(dim(newkmat)[2])
  return(newkmat)
  
}
##############################################
#2006:2036
#2037:2067
#2068:2099
#############################################
#also try this on qqmapped data
#traindf = read.csv("WRF_DMTAadjusted.csv")
#traindf = read.csv("EQM_ALL_1976_2099.csv")
traindf = read.csv("WRF_ALL_1976_2099.csv")
traindf = filter(traindf, YEAR %in% 1976:2005)

#specify days to train on

numdays = nrow(traindf)
day = 1:nrow(traindf)

delta = (traindf$meanWRF - traindf$Long_term_trend_fut) / traindf$Var_scaling_factors_wrf
years = length(unique(traindf$YEAR))


kmats180 = make_kmats(num_years = years, kernel_SD = 180)
kmats90 = make_kmats(num_years = years, kernel_SD =  90)
kmats60 = make_kmats(num_years = years, kernel_SD = 60)
kmats30 = make_kmats(num_years = years, kernel_SD = 30)
kmats14 = make_kmats(num_years = years, kernel_SD = 14)


#test on n years 
df2train <- data.frame(Dtmax=delta[1:nrow(traindf)])

#seasonal
df2train$seasonX180=kmats180[1:numdays,]
df2train$seasonX90=kmats90[1:numdays,]
df2train$seasonX60=kmats60[1:numdays,]
df2train$seasonX30=kmats30[1:numdays,]
df2train$seasonX14=kmats14[1:numdays,]
#df2train$seasonX7=kmats7[1:numdays,]
#df2train$seasonX3pt5 =kmats3pt5[1:numdays,]
df2train$sub=1
#df2train$ftimestep = as.factor(traindf$ftimestep)

start_time <- Sys.time()

#run random effects model
vf = traindf$Var_scaling_factors_wrf^2
fit2train <- lme(fixed= Dtmax ~ 1, 
                 random= list(sub = pdIdent(~seasonX14-1),
                 sub = pdIdent(~seasonX30-1),
                 sub = pdIdent(~seasonX60-1),
                 sub = pdIdent(~seasonX90-1),
                 sub = pdIdent(~seasonX180-1)),
                 #weights = ~vf, 
                 data=df2train, 
                 na.action=na.omit)

end_time <- Sys.time()
end_time - start_time
saveRDS(fit2train, "WRF_MAT_HIST_seasonalOnly_NOWEIGHTS.Rds")
#saveRDS(fit2train, "Future_1976_2005_EQM.Rds")
