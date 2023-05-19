#apply to actual data
#C --> covariance matrix for days (365x365)
#R --> covariance matrix for space (78x78)
#G --> cross covariance matrix for grid points to station locations (69552x78)
#Here are parameter estimates for temporal range, spatial range, nugget, and marginal variance
#nugget = 0.0118931, marginal variance = 0.058755, temporal range 20.5, spatial range = 186

# 1.118438489 10.000007315  0.048818415  0.027361500  0.006105197
library(dplyr)
library(geosphere)

# use this to calculate distance matrix for points
Make_SigmaNN = function(stationcoords){
  mat = distm(stationcoords[, c("LON", "LAT")])/1000
  return(mat)
}

#gridcoords will be the fine scale grid of coords
Make_SigmaGrid_N = function(gridcoords, stationcoords){
  mat = distm(gridcoords[, c("LON", "LAT")], coords[, c("LON", "LAT")]) / 1000
  return(mat)
}

#fine - scale grid -- you need this for post cond mean calculation
gridcoords = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/latlonGrid.csv")

dfall = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/newDF_with_white_REVISED_3pt5.Rds")

# Delta_TMAX_white is: Whitened corrected WRF - whitened station data
monthAves = dfall %>% group_by(MONTH,DAY, ID) %>% summarise(MeanDeltaWhite = mean(Delta_TMAX_white)) %>% 
  arrange(MONTH, DAY, ID)


# make distance matrix for stations
newdfaves = read.csv("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/GP_datasets/Data/newdfaves.csv")
newdfaves = newdfaves %>% arrange(MONTH, DAY, ID)
monthAves = newdfaves %>% group_by(MONTH,DAY, ID) %>% summarise(MeanResid = mean(meanDelta_white)) %>% 
  arrange(MONTH, DAY, ID)
stations = newdfaves[!(duplicated(newdfaves$ID)), ] %>% arrange(ID)
coords = stations[, c("LON", "LAT")]

########################
# set up data
#########################
sigmaT = make_temporal_distmat(daysInYear = 365) # function in simulateCGA.R
distC = sigmaT

SigmaNN = Make_SigmaNN(coords)

SigmaGridN = Make_SigmaGrid_N(gridcoords, coords)

z = monthAves$MeanDeltaWhite
#sd ratio
# newdfaves = newdf %>% group_by(MONTH, DAY, ID) %>% 
#   summarise(SD_Ratio_white = sd(TMAX_CORR_white) / sd(TMAX_white)) %>% 
#   arrange(MONTH, DAY, ID)
#   newdfaves$LogSD_ratio_white = log(newdfaves$SD_Ratio_white)
# z = newdfaves$LogSD_ratio_white


# final (??) params for mean
ell=10
scaleR = 134
nugget = 0.03368421
marginal_variance =   0.070
temporal_nugget = 0.05921053

# can also run log of sd ratio
#log of sd ratio
#ell = 1.5
#scaleR = 92
#nugget = 0.080792
#marginal_variance =  0.080792 / 2
#temporal_nugget = 0

#spatial and temporal covariance matrices
P = 365
# construct temporal cov matrix, C
C = exp(-2*ell^2 * (sin(pi * distC/P))^2 ) + diag(temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
######################################################
#not sure this is necessary!
#Cprime = exp(-2*ell^2 * (sin(pi * distC/P))^2 )
#CminP = matrix(min(Cprime), nrow = nrow(Cprime), ncol = ncol(Cprime))
#C2p = Cprime - CminP
#C3P = C2p / max(C2p)
#Cprime = C3P

# spatial covariance matrix
R = exp(-(distR / scaleR)^2) #+ diag(0.01, nrow = nrow(distR), ncol = ncol(distR))


# cross covariance matrix
G= SigmaGridN
G = exp(-(G / scaleR)^2)

# do eigen decompositions of C and R
svdC = svd(C)
svdR = svd(R)

Uc = svdC$u
Sc = svdC$d
Uc_t = t(svdC$v)

Ur = svdR$u
Sr = svdR$d
Ur_t = t(svdR$v)

#inverse of kronecker of singular values of C and R
InvSSmat = diag(1/(marginal_variance*as.vector(outer(Sr, Sc)) + nugget), nrow = 28470, ncol = 28470)

#reshape z to be 78x365 (no. stations x days in year)
zreshape = matrix(z, nrow=78, ncol=365)

# these steps allow for efficient calculation of conditional posterior mean
res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% solve(t(Uc)))

resmat = matrix(res1, nrow = 78, ncol = 365)

PostMeanTest = marginal_variance * as.vector(G%*%solve(Ur_t) %*% resmat %*% t(Cprime%*%solve(Uc_t)))

#plot the posterior mean on a grid for each day of the year
for (i in 1:365){
  mystart = start
  end = mystart + (69552-1) # there are 69552 grid cells on fine scale grid
  day1 = PostMeanTest[mystart:end]
  mymean = mean(day1)

  gridcoordsTest = gridcoords # coordinates of fine scale grid
  gridcoordsTest$postMean = day1
  coordinates(gridcoordsTest) = ~LON + LAT
   
  ggplot(gridcoordsTest, aes(x = LON, y = LAT, col = postMean)) +
   geom_point() +
   scale_colour_gradientn(colours = bpy.colors(15)) 
   scale_colour_gradientn(colours = bpy.colors(15), limits=c(-.7, .8))
   ggsave(paste("plot", i, ".png", sep = ""))
  
  
  
  #png(paste("postMean", i, ".png", sep = ""), width = 800, height = 1000, res = 200)
  # png("TEST.png")
  #print(spplot(gridcoordsTest, zcol = "postMean", key.space = "right",
  #             zlim = c(-5,5), cuts = mycuts, main = paste("Day ", i)))
  #dev.off()
  #dev.copy(png, paste("TEST", i, ".png", sep = ""))
  #dev.off()
}


####################
##JUNK JUNK JUNK####
####################
plotPostMean(scaleC, scaleR, nugget, marginal_variance)
hist(testmat)
plotPostMean = function(scaleC, scaleR, nugget, marginal_variance){
  #spatial and temporal covariance matrices
  P = 365
  #use this for making temporal cov mat
  C = exp(-2*ell^2 * (sin(pi * distC/P))^2 ) + diag(temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
  Cmin = matrix(min(C), nrow = nrow(C), ncol = ncol(C))
  C2 = C - Cmin
  C3 = C2 / max(C2)
  C=C3 + diag(temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
  range(Cprime)
  range(C)
  image(C)
  plot(C[,180])
  is.positive.definite(Cprime)
  ######################################################
  Cprime = exp(-2*ell^2 * (sin(pi * distC/P))^2 )
  CminP = matrix(min(Cprime), nrow = nrow(Cprime), ncol = ncol(Cprime))
  C2p = Cprime - CminP
  C3P = C2p / max(C2p)
  Cprime = C3P
  
  R = exp(-(distR / scaleR)^2) #+ diag(0.01, nrow = nrow(distR), ncol = ncol(distR))
  
  
  #Rtest = nearPD(R, keepDiag = T)
  #R = as.matrix(Rtest$mat)
  #cross covariance--Sigma_grid_n
  G= SigmaGridN
  G = exp(-(G / scaleR)^2)
  
  #svd C and R
  
  svdC = svd(C)
  svdR = svd(R)
  
  Uc = svdC$u
  Sc = svdC$d
  Uc_t = t(svdC$v)
  
  Ur = svdR$u
  Sr = svdR$d
  Ur_t = t(svdR$v)
  
  #inverse of kronecker of singular values of C and R
  #length(z)
  InvSSmat = diag(1/(marginal_variance*as.vector(outer(Sr, Sc)) + nugget), nrow = 28470, ncol = 28470)
  range(InvSSmat)
  
  #reshape z to be 78x365
  zreshape = matrix(z, nrow=78, ncol=365)
  
  
  res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% solve(t(Uc)))
  
  resmat = matrix(res1, nrow = 78, ncol = 365)
  
  #PostMeanTest = as.vector(G%*%solve(Ur_t) %*% resmat %*% t(C%*%solve(Uc_t)))
  
  PostMeanTest = marginal_variance * as.vector(G%*%solve(Ur_t) %*% resmat %*% t(Cprime%*%solve(Uc_t)))
  head(PostMeanTest)
  range(PostMeanTest)
  range(z)
  PostMeanTest = exp(PostMeanTest)
  hist(PostMeanTest)
  sum(PostMeanTest[PostMeanTest>1000])
  #range(PostMeanTest)
  #overallMean = vector(length = 365)
  pvec1 = vector(length = 365)
  pvec2 = vector(length = 365)
  pvec3 = vector(length = 365)
  pvec4 = vector(length = 365)
  pvec5 = vector(length = 365)
  pvec6 = vector(length = 365)
  pvec7 = vector(length = 365)
  pvec8 = vector(length = 365)
  pvec9 = vector(length = 365)
  start=1
  
  for (i in 1:365){
    mystart = start
    end = mystart + (69552-1)
    day1 = PostMeanTest[mystart:end]
    #print(mean(day1))
    start = end + 1
    pixvalue1 = day1[100]
    pixvalue2 = day1[10000]
    pixvalue3 = day1[20000]
    pixvalue4 = day1[55000]
    pixvalue5 = day1[44030]
    pixvalue6 = day1[5000]
    pixvalue7 = day1[15000]
    pixvalue8 = day1[30500]
    pixvalue9 = day1[60400]
    mymean = mean(day1)
    
    pvec1[i] = pixvalue1
    pvec2[i] = pixvalue2
    pvec3[i] = pixvalue3
    pvec4[i] = pixvalue4
    pvec5[i] = pixvalue5
    pvec6[i] = pixvalue6
    pvec7[i] = pixvalue7
    pvec8[i] = pixvalue8
    pvec9[i] = pixvalue9
    
    #hist(day1)
    gridcoordsTest = gridcoords
    gridcoordsTest$postMean = day1
    # coordinates(gridcoordsTest) = ~LON + LAT
    # 
    #ggplot(gridcoordsTest, aes(x = LON, y = LAT, col = postMean)) +
    # geom_point() +
    # scale_colour_gradientn(colours = bpy.colors(15)) 
    #scale_colour_gradientn(colours = bpy.colors(15), limits=c(-.7, .8))
    #ggsave(paste("plot", i, ".png", sep = ""))
    
    
    
    #png(paste("postMean", i, ".png", sep = ""), width = 800, height = 1000, res = 200)
    # png("TEST.png")
    #print(spplot(gridcoordsTest, zcol = "postMean", key.space = "right",
    #             zlim = c(-5,5), cuts = mycuts, main = paste("Day ", i)))
    #dev.off()
    #dev.copy(png, paste("TEST", i, ".png", sep = ""))
    #dev.off()
  }
  pdat1 = data.frame(pixel= pvec1, name = "P1")
  pdat2 = data.frame(pixel= pvec2, name = "P2")
  pdat3 = data.frame(pixel= pvec3, name = "P3")
  pdat4 = data.frame(pixel= pvec4, name = "P4")
  pdat5 = data.frame(pixel= pvec5, name = "P5")
  pdat6 = data.frame(pixel= pvec6, name = "P6")
  pdat7 = data.frame(pixel= pvec7, name = "P7")
  pdat8 = data.frame(pixel= pvec8, name = "P8")
  pdat9 = data.frame(pixel= pvec9, name = "P9")
  
  
  allp = rbind(pdat1, pdat2, pdat3, pdat4,pdat5,pdat6,pdat7,pdat8,pdat9)
  allp$DOY = rep(1:365, 9)
  
  
  print(ggplot(allp, aes(x = DOY, y = pixel, col = name)) + 
          geom_line())
}

gridcoordsTest = gridcoords
gridcoordsTest$postMean = day1
gridcoordsTest = gridcoordsTest[, c("LON", "LAT", "postMean")]
coordinates(gridcoordsTest) = ~LON + LAT
spplot(gridcoordsTest, "postMean", sp.layout = list("sp.points", stations2, pch = 16, cex = .9, col = "green"))
spplot(gridcoordsTest) + layer(panel.points(x, y, col="black", pch=19), data=stations2)


library(latticeExtra)
head(stations)
stations2 = dplyr::select(stations, LON, LAT)
names(stations2) = c("x", "y")
coordinates(stations2) = ~x + y
points(stations2)
plot(stations2, axes = T)

spplot(meuse.grid, "dist", scales = list(draw = TRUE),
       sp.layout = list("sp.points", stations2, pch = 16, cex = 2, col = "black"))


head(meanvec2)

par(mfrow = c(1,1))
plot(overallMean, type = "l")
plot(pvec2, col = "purple", type ='l')
plot(pvec3, col = "blue", type ='l')
plot(pvec4, col = "red", type ='l')

pdat1 = data.frame(pixel= pvec1, name = "P1")
pdat2 = data.frame(pixel= pvec2, name = "P2")
pdat3 = data.frame(pixel= pvec3, name = "P3")
pdat4 = data.frame(pixel= pvec4, name = "P4")
allp = rbind(pdat1, pdat2, pdat3, pdat4)
allp$DOY = rep(1:365, 4)


ggplot(allp, aes(x = DOY, y = pixel, col = name)) + 
  geom_line()

dat = data.frame(DOY = 1:365, Mean = meanvec)
ggplot(dat, aes(x = DOY, y = Mean)) + 
  geom_point() +
  geom_line()

head(gridcoordsTest)
ggplot(gridcoordsTest, aes(x = LON, y = LAT, col = postMean)) +
  geom_point() +
  scale_colour_gradientn(colours = bpy.colors(15), limits=c(-5,5))

range(PostMeanTest)
mycuts = seq(-5,5, by = 1)
?spplot

?key.space
plot(gridcoordsTest)
mymin = min(day1)
mymax = max(day1)
myAT = seq(mymin, mymax, by = .2)
myAT
colorkey = list(height = 1, labels = list(at = seq(0.5, length(at) -0.5), labels = at))

