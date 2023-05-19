
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
make_kmats = function(day, kernel_SD){
  num_knots = ceiling(365 / kernel_SD)
  knot_spacing = 365 / num_knots
  dmats = doyDist(day,seq(1,364, by = knot_spacing))
  kmats <- dnorm(dmats, sd = kernel_SD) / dnorm(0, sd = kernel_SD)
  int2k = max(kmats%*%t(kmats))
  newkmat = kmats / sqrt(int2k)
  print(dim(newkmat)[2])
  return(newkmat)
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
traindf = read.csv("WRF_ALL_crossval.csv")
#traindf = filter(traindf, YEAR %in% 1976:2005)
#specify days to train on

#20 years
#day = 1:7000
numdays = nrow(traindf)
day = 1:nrow(traindf)

delta = traindf$tmax_norm_VAR1


#kmats3pt5 = make_kmats(day = day, kernel_SD = 3.5)
#kmats7 = make_kmats(day = day, kernel_SD = 7)
#kmats14 = make_kmats(day = day, kernel_SD = 14)
#kmats21 = make_kmats(day = day, kernel_SD = 21)
#kmats28 = make_kmats(day = day, kernel_SD = 28)
#kmats35 = make_kmats(day = day, kernel_SD = 35)
#kmats42 = make_kmats(day = day, kernel_SD = 42)
kmatt3pt5 = make_longterm_Kmat(day, kernel_SD = 3.5)
#kmatt7 = make_longterm_Kmat(day, kernel_SD = 7)
#kmatt14 = make_longterm_Kmat(day, kernel_SD = 14)
#kmatt21 = make_longterm_Kmat(day, kernel_SD = 21)
#kmatt28 = make_longterm_Kmat(day, kernel_SD = 28)
kmatt180 = make_longterm_Kmat(day, kernel_SD = 180)
#kmatt365 = make_longterm_Kmat(day, kernel_SD = 365)
#kmatt730 = make_longterm_Kmat(day, kernel_SD = 730)
#test on n years 
df2train <- data.frame(Dtmax=delta[1:nrow(traindf)])
#seasonal K
#df2train$seasonX3pt5=kmats3pt5[1:numdays,]
#df2train$seasonX7=kmats7[1:numdays, ]
#df2train$seasonX14=kmats14[1:numdays, ]
#df2train$seasonX28=kmats28[1:numdays,]
#df2train$seasonX35=kmats35[1:numdays,]
#df2train$seasonX21=kmats21[1:numdays,]
#df2train$seasonX42=kmats42[1:numdays,]
#df2train$seasonX56=kmats56[1:numdays,]
#long term K
df2train$trendt3pt5=kmatt3pt5[1:numdays,]
#df2train$trendt7=kmatt7[1:numdays,]
#df2train$trendt14=kmatt14[1:numdays,]
#df2train$trendt21=kmatt21[1:numdays,]
#df2train$trendt28=kmatt28[1:numdays,]
df2train$trendt180=kmatt180[1:numdays,]
#df2train$trendt365 = kmatt365[1:numdays,]
#df2train$trendt730=kmatt730[1:numdays,]
df2train$sub=1

start_time <- Sys.time()

#run random effects model
fit2train <- lme(fixed= Dtmax ~ 1, 
                 random= list(sub = pdIdent(~trendt3pt5-1),
                 sub = pdIdent(~trendt180-1)),
                  #sub = pdIdent(~trendt365-1),
                  #sub = pdIdent(~trendt730-1)),
                    data=df2train, na.action=na.omit)

end_time <- Sys.time()
end_time - start_time
saveRDS(fit2train, "WRF_ALL_forCrossVal.Rds")
print("DONE")