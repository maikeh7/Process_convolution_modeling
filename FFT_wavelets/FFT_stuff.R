
myt=1:96
set.seed(168)
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
y = .5 * cos(2*pi*1:96*2/96)  + 1*sin(2*pi*1:96*2/96) + 
  2.3*cos(2*pi*1:96*10/96)  + 1.5*sin(2*pi*1:96*10/96) # same as above...

I <- abs(fft(y))^2/length(y)
P <- (4/length(y))*I[(1:floor(length(y)/2))] #scaled periodogram --> Aj^2 + Bj^2 = scaled periodogram_j

ztsq=sum((cos(2*pi*1:96*2/96))^2) #here is how you get coefficients!!!!!!!!!!!!!
top=sum(y*cos(2*pi*1:96*2/96))
top/ztsq # equals beta1 for cos w/ freq1

#ok so this basically gives you the same thing...
#I <- fft(y)^2/length(y)
#P <- (4/length(y))*I[(1:floor(length(y)/2))] 

x1= cos(2*pi*1:96*2/96)
x2=1*sin(2*pi*1:96*2/96) 
x3 =  2.3*cos(2*pi*1:96*10/96) 
x4 =  1.5*sin(2*pi*1:96*10/96)
summary(fit <- lm(y~0 + x1 + x2 + x3 + x4))


plot(y, type = "l")
library(TSA)
specto=periodogram(y)
specto$spec
f_full=fft(y) 
f_useful=f_full

dt=0:95 #this needs to be from 0, not 1 or whatever you have in signal 
myFreqs= specto$freq 


basis1 = exp(2*pi*complex(imaginary = 1)* myFreqs[1]*dt)
basis2 = exp(2*pi*complex(imaginary = 1)* myFreqs[2]*dt)
basis3 = exp(2*pi*complex(imaginary = 1)* myFreqs[3]*dt)
basis4 = exp(2*pi*complex(imaginary = 1)* myFreqs[4]*dt)
basis5 = exp(2*pi*complex(imaginary = 1)* myFreqs[5]*dt)



basis1 = exp(2*pi*complex(imaginary = 1)* myFreqs[1]*dt) # important bases
basis23 = exp(2*pi*complex(imaginary = 1)* myFreqs[23]*dt)

plot(basis23)

#reconstruct

plot(y, type="l")
plot(recon2, type = "l")
lines(y)
recon1 = Re(basis1*f_useful[1])
recon2 = Re(basis2*f_useful[2])
recon3 = Re(basis3*f_useful[3])
recon4 = Re(basis4*f_useful[4])
recon5 = Re(basis5*f_useful[5])
recon9 = Re(basis9*f_useful[9])

recon1 = Re(basis1*f_useful[2]) # need to add one to the important frequences from TSA::periodogram
recon23 = Re(basis23*f_useful[24])

plot((recon1 + recon23)/48, type = "l") #need to divide by half of length of signal

lines(y, type='l', col = 'red')

plot(recon32, type = "l", col ='red')

head(wrf_all)
obsmodel = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/LME_model_fits/OBS_1976_2005.Rds")
obs_data = obsmodel$data$Dtmax
obs_data = obs_data[1:8000]

wrf_all = readRDS("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/data/Final_datasets/Matern/wrf_all.Rds")
test_signal  = wrf_all$tmax_norm_VAR1[1:8000]

test_periodogram = periodogram(test_signal)
head(test_periodogram$spec, n=60) #large peak at 44
test_periodogram$spec[test_periodogram$spec > 30]
which(test_periodogram$spec>70)

#WRF
dt=0:7999 #this needs to be from 0, not 1 or whatever you have in signal 
#freq=8000 / index
# 43=186 days
# 61=131
# 4=2000
# 3 = 2667
# 19 = 421
# 104 = 77 days
test_periodogram = periodogram(test_signal)
freq_indices = which(test_periodogram$spec > 30)
specs = test_periodogram$spec[test_periodogram$spec > 30]
freq_indices = freq_indices[order(specs, decreasing = TRUE)]
freq_indices

myFreqs= test_periodogram$freq 
f_full=fft(test_signal) 
f_useful=f_full

basis_mat = matrix(nrow=length(freq_indices), ncol = length(test_signal))
counter=1
for (i in freq_indices){
  A_basis = exp(2*pi*complex(imaginary = 1)* myFreqs[i]*dt)
  recon = Re(A_basis*f_useful[i+1]) / (0.5*length(test_signal))
  basis_mat[counter, ] = recon
  counter = counter+1
}

# 43=186 days
# 61=131
# 4=2000
# 3 = 2667
# 19 = 421
# 104 = 77 days
days = c(186, 131, 2000, 2667, 421, 77)
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/FFTMOD.png", height = 1000, width = 1800, res=250)
par(mfrow=c(7,1), mar=c(0,0,0,0), oma = c(4,4,1,1))
plot(basis_mat[1,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[1], " days"))
plot(basis_mat[2,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[2], " days"))
plot(basis_mat[3,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[3], " days"))
plot(basis_mat[4,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[4], " days"))
plot(basis_mat[5,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[5], " days"))
plot(basis_mat[6,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[6], " days"))
#plot(basis_mat[7,][1:4000]/4000, type = "l", xaxt="n", ylim = c(-3.3, 3.3))

plot(test_signal[1:4000], type = "l", ylim = c(-3.3, 3.3), col = "blue")
lines(colSums(basis_mat[1:6, ])[1:4000], col = "red")
dev.off()
###############################################

myFreqs= test_periodogram$freq 
plot(test_periodogram)
f_full=fft(test_signal) 
f_useful=f_full
myFreqs[which(test_periodogram$spec > 40)]



basis3 = exp(2*pi*complex(imaginary = 1)* myFreqs[3]*dt)
recon3 = Re(basis3*f_useful[4])

basis4 = exp(2*pi*complex(imaginary = 1)* myFreqs[4]*dt)
recon4 = Re(basis4*f_useful[5])

basis19 = exp(2*pi*complex(imaginary = 1)* myFreqs[19]*dt)
recon19 = Re(basis19*f_useful[20])

basis44 = exp(2*pi*complex(imaginary = 1)* myFreqs[44]*dt)
recon44 = Re(basis44*f_useful[45])

basis61 = exp(2*pi*complex(imaginary = 1)* myFreqs[61]*dt)
recon61 = Re(basis61*f_useful[62])

basis66 = exp(2*pi*complex(imaginary = 1)* myFreqs[66]*dt)
recon66 = Re(basis66*f_useful[67])

basis88 = exp(2*pi*complex(imaginary = 1)* myFreqs[88]*dt)
recon88 = Re(basis88*f_useful[89])
allrecon = (recon3 + recon4 + recon19 + recon44 + recon61 + recon66 + recon88)/4000

png("FourierPlot.png", width = 1500, height = 1000, res = 250)
par(mfrow=c(4,2), mar=c(0,0,0,0), oma = c(4,4,1,1))

#days = c(90, 121, 131, 181, 421, 2000, 2667)
plot(recon88[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3),
     xaxt="none", ylab = "")
text(2000, -2.5, paste("Freq: 90 days"))

plot(recon66[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
     xaxt="none", yaxt="none")
text(2000, -2.5, paste("Freq: 121 days"))

plot(recon61[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
     xaxt="none")
text(2000, -2.5, paste("Freq: 131 days"))

plot(recon44[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
     xaxt="none", 
     yaxt="none")
text(2000, -2.5, paste("Freq: 181 days"))

plot(recon19[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
     xaxt="none")
text(2000, -2.5, paste("Freq: 421 days"))

plot(recon4[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "",
     xaxt="none", 
     yaxt="none")
text(2000, -2.5, paste("Freq: 2000 days"))

plot(recon3[1:4000]/4000, col = "blue", type= 'l', xlim = c(0,4000), ylim=c(-3,3), ylab = "")
text(2000, -2.5, paste("Freq: 2667 days"))

plot(test_signal, type = "l", xlim = c(0,4000), ylim=c(-3,3) , yaxt="none", cex.axis=0.75)
lines(allrecon, col = "red")
dev.off()



#Observed data
#just try to reconstruct the 180 ish frequency component for mod
dt=0:7999 #this needs to be from 0, not 1 or whatever you have in signal 
#69=116 days
#18=444
#41=195
#285=28
#376=21
#144=56 days
obs_periodogram = periodogram(obs_data)
freq_indices = which(obs_periodogram$spec > 30)
specs=obs_periodogram$spec[obs_periodogram$spec > 30]
freq_indices = freq_indices[order(specs, decreasing = TRUE)]


myFreqs= obs_periodogram$freq 
obs_periodogram$spec
f_full=fft(obs_data) 
f_useful=f_full

basis_mat = matrix(nrow=length(freq_indices), ncol = length(obs_data))
counter=1
for (i in freq_indices){
  A_basis = exp(2*pi*complex(imaginary = 1)* myFreqs[i]*dt)
  recon = Re(A_basis*f_useful[i+1]) / (0.5*length(obs_data))
  basis_mat[counter, ] = recon
  counter = counter+1
}

# basis_mat = matrix(nrow=2000, ncol = length(obs_data))
# counter=1
# for (i in 1:2000){
#   A_basis = exp(2*pi*complex(imaginary = 1)* myFreqs[i]*dt)
#   recon = Re(A_basis*f_useful[i+1]) / (0.5*length(obs_data))
#   basis_mat[counter, ] = recon
#   counter = counter+1
# }
# 
# recon=colSums(basis_mat)
# par(mfrow=c(1,1))
# plot(obs_data[1:4000], type = "l")
# lines(recon, col = "red")


days = c(116, 144, 195, 28, 21, 56)
png("C:/Users/Maike/Box Sync/EPSCOR/GPstuff/TMAX/manuscript_plots/FFTOBS.png", height = 1000, width = 1800, res=250)

par(mfrow=c(7,1), mar=c(0,0,0,0), oma = c(4,4,1,1))
plot(basis_mat[1,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[1], " days"))
plot(basis_mat[2,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[2], " days"))
plot(basis_mat[3,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[3], " days"))
plot(basis_mat[4,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[4], " days"))
plot(basis_mat[5,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[5], " days"))
plot(basis_mat[6,][1:4000], type = "l", xaxt="n", ylim = c(-3.3, 3.3))
text(200,-2, paste("Freq: ", days[6], " days"))
#plot(basis_mat[7,][1:4000]/4000, type = "l", xaxt="n", ylim = c(-3.3, 3.3))

plot(obs_data[1:4000], type = "l", ylim = c(-3.3, 3.3), col = "blue")
lines(colSums(basis_mat[1:6, ])[1:4000], col = "red")
dev.off()


# use f_useful[f + 1], where f is the index of frequency from FFT
basis18 = exp(2*pi*complex(imaginary = 1)* myFreqs[18]*dt)
recon18 = Re(basis18*f_useful[19])

basis41 = exp(2*pi*complex(imaginary = 1)* myFreqs[41]*dt)
recon41 = Re(basis41*f_useful[42])

basis69 = exp(2*pi*complex(imaginary = 1)* myFreqs[69]*dt)
recon69 = Re(basis69*f_useful[70])

basis285 = exp(2*pi*complex(imaginary = 1)* myFreqs[285]*dt)
recon285 = Re(basis285*f_useful[286])

par(mfrow=c(1,1))
plot(obs_data[1:1000], type = "l")
lines(recon18/4000, col = "red")
lines(recon41/4000, col = "red")
lines(recon69/4000, col = "red")
lines((recon18 + recon41 + recon69 +recon285)/4000, type = "l", col = "red") 
plot(obs_data[1:1000], type = "l")

###############################################################################################
10950 / 
maternDat = read.csv("WRF_Matern.csv")
signal = maternDat$tmax_norm_VAR1
signal = obstrain$tmax_norm_VAR1
signal = residsw
dt=0:(length(signal)-1)

specto=periodogram(signal)
spectoBand = specto$freq[(specto$freq >=.1) & (specto$freq <=.2)]
spectoBand = which((specto$freq >=.08) & (specto$freq <=.12))
spectoBand = which((specto$freq >=.1) & (specto$freq <=.2))
range(spectoBand)
specto$spec[1125:1500]

which(specto$spec >50)
f_full = fft(signal) 
f_useful = f_full

myFreqs= specto$freq 

num_bases=3000
basis_mat = matrix(nrow=num_bases, ncol = length(signal))

for (i in 1:num_bases){
  A_basis = exp(2*pi*complex(imaginary = 1)* myFreqs[i]*dt)
  recon = Re(A_basis*f_useful[i+1]) / 0.5*length(signal)
  basis_mat[i, ] = recon
}
nrow(basis_mat_sub)
basis_mat_sub = basis_mat[spectoBand, ]
myrecon = colSums(basis_mat_sub[1:nrow(basis_mat_sub), ])
plot(myrecon[1:1000], type = "l", main = "FFT reconstructions in freq = [.1, .2]")
plot(basis_mat_sub[1, ], type = "l", xlim = c(0, 1000), ylim = c(-.1, .1))
for(i in 1:20){
  lines(basis_mat_sub[i, ])
}

periodogram(rnorm(10000))






#######################
#DWT
#######################
library(wavethresh)

write.csv(obs_data, "obs_data.csv")
obs_data = obsmodel$data$Dtmax
obs_data = obs_data[1:512]
wds1 <- wd(obs_data,family="DaubExPhase", bc = "periodic")
length(wds1$C)
length(wds1$D)
t1

plot(wds1$C)
wds1$nlevels

plot(obs_data, col = "red")
par(mfrow = c(4, 2))
t1 = accessC(wds1, level=1)
t1
t1f = accessD(wds1, level =9)
plot(t1, type = "l")

plot(t1, type = "l")
plot(t1f, type = "l")
t6= accessC(wds1, level = 6)
t6f = accessD(wds1, level = 6)

t7= accessC(wds1, level = 7)
t7f = accessD(wds1, level = 7)

t5= accessC(wds1, level = 5)
t5f = accessD(wds1, level = 5)

t4= accessC(wds1, level = 4)
t4f = accessD(wds1, level = 4)

FineCoefs <- accessD(wds1, lev=8)
sigma <- mad(FineCoefs)
sigma 
utDJ <- sigma*sqrt(2*log(512))
utDJ
ywdT <- threshold(wds1, policy="manual", value=utDJ)
ywr <- wr(ywdT)

plot( ywr, type="l")
lines(obs_data, col='red') 
lines(x, v$bumps, lty=2)
 W1 <-t(GenW(n=512, filter.number=10, family="DaubExPhase"))
mytest = W1 %*% obs_data
tail(mytest)
plot(mytest, type="l")
mytest[1:450]=0 

Winv = t(W1)
recon = Winv %*% mytest 
plot(recon, type = "l")
lines(obs_data, col = "red")

W2=W1
W2[2:512, ]  =0
image(W2)
test=W2 %*% obs_data
plot(test, type = "l")
head(t1,n=50)

test[1:20]
test=wds1$C
tail(test)
plot(test, type='l')
lines(t1, col = "red")
test[588:800]
tail(test, n=50)
test2 = wds1$D
length(test2)
test2[511]
?drawwp.default
?wd
basis <- matrix(NA,length(l$level),512)
for (i in 1:length(l$level))
  basis[i,] <- drawwp.default(l$level[i],l$pkt[i],resolution=512)

#https://www.routledgehandbooks.com/doi/10.1201/9781420033397.ch3#sec3_2_4

y <- c(1,1,7,9,2,8,8,6)
plot(y, type='l')
ywd <- wd(y, filter.number=1, family="DaubExPhase")
ywd$C
ywd$D
accessC(ywd, level=0)
accessD(ywd, level=0)

accessC(ywd, level=1)
accessD(ywd, level=1)

accessC(ywd, level=2)
accessD(ywd, level=2)

plot(y, type = 'l')
W1 <- t(GenW(n=8, filter.number=1, family="DaubExPhase"))
mydwt=W1%*% y
Winv=solve(W1)
tmean = mydwt
tmean[2:8] = 0


test = Winv %*% tmean
lines(test, type ='l')
plot(y, col='red', type='l')

dets = mydwt
dets[2:5] = 0
ds = Winv %*% dets


lines(ds)

td = tt
td[1:5] = 0
test2 = Winv %*% td
lines(test2, type = "l")
plot(td + test, type  ='l', ylim = c(-5, 10))
plot(y, col='red', type="l")
lines(td)
test=Winv%*%tt
head(tt)
lines(test, type ='l')
sum(tt[2:8])
ywd$C
ywd$D

y1=y
plot(y, type = "l")
y2=y[1:64]
y2=y2+6
plot(y2, type = "l")

dw_y2 = wd(y2, filter.number = 10, family="DaubExPhase")
plot(dw_y2)

wrf_all = readRDS("wrf_all.Rds")
y2 = wrf_all$tmax_norm_VAR1
y2 = y2[1:8192]
mydwt = wd(y2)
plot(mydwt)
W1 <- t(GenW(n=64, filter.number=1, family="DaubExPhase"))

mydwt=W1%*% y2

Wprime = W1
Wprime32=Wprime[2:33, ]
Wprime16=Wprime[34:49, ]
Wprime8=Wprime[50:57, ]
Wprime4=Wprime[58:61, ]
Wprime2=Wprime[62:63, ]
WprimeMean = Wprime[1, ]
WprimeMean = matrix(WprimeMean, nrow=1, ncol=64)

newC = mydwt
newC32 = newC[2:33, ] #only finest level coeffs
newC16 = newC[34:49, ] # next finest level...etc
newC8 = newC[50:57, ]
newC4 = newC[58:61, ]
newC2 = newC[62:63, ]
newCmean = newC[1, ]

s32 = t(Wprime32) %*% newC32
s16= t(Wprime16) %*% newC16
s8= t(Wprime8) %*% newC8
s4= t(Wprime4) %*% newC4
s2= t(Wprime2) %*% newC2

smean = t(WprimeMean) %*% newCmean


mycols = viridis::plasma(7)
plot(y2,type = "l")
lines(s32, col = mycols[1])
lines(s16+smean, col = mycols[2])
lines(s8+smean, col = mycols[3])
lines(s4+smean, col = mycols[4])
lines(s2+smean, col = mycols[5])
plot(y2,type = "l")

lines(s32+s16+s8+s4 + smean, col = "red")





lines(smean, col= mycols[6])

lines(test, type = "l", col = "red")
tmean = mydwt
length(tmean)
tmean[2:64] = 0
test = Winv %*% tmean
lines(test, type ='l', col = "red")

dets = mydwt
dets[2:8] = 0
ds = Winv %*% dets
plot(y2, type = "l")
lines(ds, col = "red")
#0.000000 -1.414214 -4.242641 1.414214
library(wavelets)
# obtain the two series listed in Percival and Walden (2000), page 42
X1 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,.7,.9,0,.3)
X2 <- c(.2,-.4,-.6,-.5,-.8,-.4,-.9,0,-.2,.1,-.1,.1,-.7,.9,0,.3)
plot(X1, type='l')
# compute the LA8 wavelet filter for both DWT and MODWT
la8.dwt <- wt.filter()
la8.dwt
la8.modwt <- wt.filter(modwt=TRUE)
# compute the DWT and MODWT level one wavelet and scaling coefficients
wt.dwt <- dwt.forward(X1, la8.dwt)
wt.modwt <- modwt.forward(X2, la8.modwt, 1)
# compute the original series with the level one coefficients
newX.dwt <- dwt.backward(wt.dwt$W, wt.dwt$V, la8.dwt)
newX.modwt <- modwt.backward(wt.modwt$W, wt.modwt$V, la8.modwt, 1)
newX.modwt
X2
plot(X2)
points(newX.modwt, col="red")
