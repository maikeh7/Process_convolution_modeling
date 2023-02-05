library(waveslim)
?waveslim::dwt
testdat = wrf_all$tmax_norm_VAR1[1:8192]
testdat = wrf_all$tmax_norm_VAR1[1:16]

vals = seq( 0.7128906, 365, by= 0.7128906)
vals[1]=1
head(vals)
temp  = wrf_all$tmax_norm_VAR1[1:365]
i=1977
testing=approx(1:365, temp, xout = vals)
#interpolate 1 year to 512 days so it is more sensible in wavelet transform...but does not seem to make much difference to me
datavec = vector()
for (year in 1976:1991){
  yeardat = filter(wrf_all, YEAR == year)$tmax_norm_VAR1
  test_approx = approx(1:365, yeardat, xout = vals)
  datavec = c(datavec, test_approx$y)
  
}
test_approx$y
length(datavec)

y2=testdat
par(mfrow = c(1,1))
test=up.sample(testdat, 64)
plot(test, type = "l")
n <- length(testdat)
f=64
y=NA
as.vector(rbind(x, matrix(rep(y, (f - 1) * n), nrow = f - 
                            1)))
####################################################################
set.seed(168)
myt=1:96
myint=sample(48, 2)

freq1 = myint[1] / 96
freq2 = myint[2] / 96
A1 = rnorm(1,0,2)
B1 = rnorm(1, 0,2)
A2 = rnorm(1,0,3)
B2 = rnorm(1, 0, 3)
w = 2*pi*myt

y = .5 * cos(w*freq1)  + 1*sin(w*freq1) + 2.3*cos(w*freq2)  + 1.5*sin(w*freq2) + rnorm(96, 0,.5)
plot(y, type = "l")

y2=y[1:64]
y2=y2+6
plot(y2, type='l')
y2 = y2[1:64]

wrf_all = readRDS("/data/mholthui/climate_analysis/Chapter2/TMAX/wrf_all.Rds")
testdat = wrf_all$tmax_norm_VAR1[1:8192]
#testdat = wrf_all$tmax_norm_VAR1[1:5840]
y2=testdat
y2=datavec

mydwt = dwt(y2, "la8", n.levels = 13)
mydwt = dwt(y2, "la8", n.levels = 4)

nlevels= length(mydwt)
length(mydwt)

i=1
sdvec=vector()
alist = list()
reconlist = list()
for (i in 1:nlevels){
  temp = mydwt
  for (j in 1:nlevels){
    if (j==i){
      next 
    }
    len = length(temp[[j]])
    temp[[j]] = rep(0, len)
  }
  reconlist[[i]] = idwt(temp)
  alist[[i]] = temp
  sdvec[i] = sd(temp[[i]])
}
plot(y2, type="l", xlim=c(0,4000))

days = c(2,4,8,16,32,64)
par(mfrow=c(7,1), mar=c(0,0,0,0), oma = c(4,4,1,1))
plot(reconlist[[1]][1:4000] + reconlist[[14]][1:4000], type="l",xaxt="n", xlim=c(1, 4000), ylim = c(-3,3))
text(200,-2, paste("Freq: ", days[1], " days"))
for (i in 2:6){
  plot(reconlist[[i]][1:4000]  + reconlist[[14]][1:4000] , type = "l",xaxt="n",ylim=c(-3,3), xlim=c(1,4000))
  text(200,-2, paste("Freq: ", days[i], " days"))
}
plot(y2[1:4000] , type = "l", xlim = c(1,4000), ylim=c(-3,3), col = "blue")
lines(reconlist[[1]][1:4000]  + reconlist[[2]][1:4000]  + reconlist[[3]][1:4000]  + reconlist[[4]][1:4000]  +
       reconlist[[5]][1:4000]  + reconlist[[6]][1:4000]  + reconlist[[14]][1:4000] , col = "red")


round(sdvec,2)
par(mfrow=c(1,1))
png("wavelotPlot.png", width = 1500, height = 1000, res = 250)
par(mfrow=c(5,2), mar=c(0,0,0,0), oma = c(4,4,1,1))
days = c(2,4,8,16,32,64,128,256,512)
for (i in 1:9){
  if (i %%2 == 0){
    if (i == 6 || i == 8){
      plot(reconlist[[i]] + reconlist[[14]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "", xaxt="none", 
           yaxt="none")
      text(2000, -2.5, paste("Freq: ", days[i], " days *"))
    }else{
    plot(reconlist[[i]] + reconlist[[14]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "", xaxt="none", 
         yaxt="none")
    text(2000, -2.5, paste("Freq: ", days[i], " days"))
    }
  }else if (i==9){
    plot(reconlist[[i]] + reconlist[[14]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
         cex.axis=0.75)
    text(2000, -2.5, paste("Freq: ", days[i], " days *"))
  }else{
    if (i ==7){
      plot(reconlist[[i]] + reconlist[[14]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "", xaxt = "none",
           cex.axis=0.75)
      text(2000, -2.5, paste("Freq: ", days[i], " days *"))
    }else{
    plot(reconlist[[i]] + reconlist[[14]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "", xaxt = "none",
         cex.axis=0.75)
    text(2000, -2.5, paste("Freq: ", days[i], " days"))
    }
  }
}
plot(y2, type = "l", xlim = c(0,4000), ylim=c(-3,3) , yaxt="none", cex.axis=0.75)
lines(reconlist[[6]] + reconlist[[7]] + reconlist[[8]] + reconlist[[9]] + reconlist[[14]], col = "red")
dev.off()

plot(reconlist[[10]]+ reconlist[[14]], type = "l", col = "blue",ylim=c(-3,3))
plot(reconlist[[11]]+ reconlist[[14]], type = "l", col = "blue",ylim=c(-3,3))



















#5048 days
days = c(1, 2, 4, 8, 16)
for (i in 1:length(reconlist)){
  plot(reconlist[[i]], col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "", xaxt="none")
  text(2000, -2, paste("Freq: ", days[i], " days"))
}
plot(y2, type = "l", xlim = c(0,4000), ylim=c(-3,3))
lines(reconlist[[4]] + reconlist[[5]], col = "red")

length(mydwt[[5]])


npar(mfrow = c(1,1))
plot(y2, type="l")
plot(reconlist[[1]], type = "l")
lines(reconlist[[2]], col = "red")
lines(reconlist[[3]], col = 'blue')
lines(reconlist[[4]], col = "orange")
lines(reconlist[[5]], col = "red", lwd=2)


sdvec



lines(reconlist[[1]], col = "red")

lines(reconlist[[2]], col = "blue")
lines(reconlist[[3]], col = "purple")
lines(reconlist[[4]], col = "orange")
lines(reconlist[[5]], col = "green")
lines(reconlist[[6]], col = "darkgreen")
reconlist[[6]]

tmydwt$d1
mydwt$d2

mydwt$d5
mydwt[[6]]

lev32 = mydwt
lev16= mydwt
lev8= mydwt
lev4= mydwt
lev2= mydwt
lev_mean= mydwt

lev32$d1
lev32[[1]] = rep(0, 32)
lev16[[2]] = rep(0,16)
lev8[[3]] = rep(0, 8)
lev4[[4]] = rep(0,4)
lev2[[5]] = rep(0,2)
lev_mean[[6]] = 0

recon1 = idwt(lev32)
recon2 = idwt(lev16)
recon3 = idwt(lev8)
recon4 = idwt(lev4)
recon5 = idwt(lev2)
recon_mean = idwt(lev_mean)

mycols = viridis::plasma(13)
plot(y2,type = "l")
lines(recon1, col = mycols[1])
lines(s16+smean, col = mycols[2])
lines(s8+smean, col = mycols[3])
lines(s4+smean, col = mycols[4])
lines(s2+smean, col = mycols[5])
plot(recon1, type='l')
lines(y2, col = 'red')