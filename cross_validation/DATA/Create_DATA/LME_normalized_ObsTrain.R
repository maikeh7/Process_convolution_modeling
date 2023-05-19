#!/usr/bin/Rscript
# use this script to construct training PC models of observed data
# just run the .sh file
args<-commandArgs(trailingOnly=TRUE)

library(nlme)

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
#############################################
run_LME = function(filenumbers = args[1]){
  for(j in 1:length(filenumbers)) {
    print("filenumbers is:")
    print(filenumbers)
    MY_FILE = as.integer(filenumbers[j])
    
    
    file_to_read = paste("OBS_TRAIN_", MY_FILE, ".csv", sep = "")
    print(file_to_read)
    traindf = read.csv(file_to_read)

    numdays = nrow(traindf)

    day = 1:nrow(traindf)

    delta = traindf$tmax_norm_VAR1


    kmatt3pt5 = make_longterm_Kmat(day, kernel_SD = 3.5)

    kmatt180 = make_longterm_Kmat(day, kernel_SD = 180)
    
    #kmatt365 = make_longterm_Kmat(day, kernel_SD = 365)
    
    #kmatt730 = make_longterm_Kmat(day, kernel_SD = 730)

    #test on n years 
    df2train <- data.frame(Dtmax=delta[1:nrow(traindf)])

    #long term K
    df2train$trendt3pt5=kmatt3pt5[1:numdays,]

    df2train$trendt180=kmatt180[1:numdays,]
    
    df2train$trendt365=kmatt365[1:numdays,]
    
    df2train$trendt730=kmatt730[1:numdays,]

    df2train$sub=1


    #run random effects model
    fit2train <- lme(fixed= Dtmax ~ 1, 
                     random= list(sub = pdIdent(~trendt3pt5-1),
                     sub = pdIdent(~trendt180-1),
                      #sub = pdIdent(~trendt365-1),
                      # sub = pdIdent(~trendt730-1)),
                     # sub = pdIdent(~trendt1825-1)),
                        data=df2train, na.action=na.omit)
    file_to_save = paste("OBS_TRAIN_", MY_FILE, ".Rds", sep = "")
    print(file_to_save)
    saveRDS(fit2train, file_to_save)
    print("done")
    }
}
run_LME()
    #saveRDS(fit2train, "/for_testing/LME6basesALLdays.Rds")