# this file contains a number of utility functions to generate
# various entities used in spatial prediction.  Conditional
# normal functions, generalized inverses, covariance formulas

# Also, examples used in the first section of the shortcourse
# are shown here...

# psuedo inverse of  a matrix.

svdinv <- function(a, tol=1e-6){
 sva <- svd(a)
 maxsv <- max(sva$d)
 dd <- ifelse(sva$d<tol*maxsv,0,1/sva$d)
 t(sva$u %*% (t(sva$v)*dd))
}

# generate multivariate Normals
rmultnorm <- function(n,mu,sigma){
   # returns n rows of multinormal mu,sigma random vectors
   # ie: returns a n x length(mu) matrix
  p <- length(mu)
  z <- matrix(rnorm(n * p),nrow=n)
  ch <- chol(sigma,pivot=T)
  piv <- attr(ch,"pivot")
  zz <- (z%*%ch)
  zzz <- 0*zz
  zzz[,piv] <- zz
  zzz + matrix(mu,nrow=n,ncol=p,byrow=T)
}

gaupost <- function(y,Cy,mx,Cx){
 # computes the posterior distribution of x given observed
 # data y; the data model is normal centered at x with var Cy;
 # the prior for x in N(mx,Cx)
 Sxi <- t(solve(Cx+Cy,Cx))
 mpost <- mx+Sxi%*%y
 Cpost <- Cx-Sxi%*%Cx
 return(list(mean=mpost,cov=Cpost))
}
?solve

# some examples of spatial surfaces for different 
# covariance/precision specifications
locs=1:10
n=100
dim(mat)
# covariance functions
covbm <- function(locs,pow=1,scale=1){
  # assume anchored at 1st component
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)

  vario <- (mat/scale)^pow
  covm <- matrix(0,n,n)
  vj <- matrix(vario[1,],n,n,byrow=T)
  vi <- matrix(vario[,1],n,n,byrow=F)
  cc <- (vi+vj-vario)/2

  return(cc)
}
covpow <- function(locs,pow=2,scale=5){
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  cc <- exp(-(mat/scale)^pow)
  return(cc)
}

covsph <- function(locs,scale=5){
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  mat <- mat/scale
  cc <- ifelse(mat>1,0,
     ifelse(mat>0,
     1-2/pi*(sqrt(1-mat^2)*mat + asin(mat)),1))
  return(cc)
}
covmat52 <- function(locs,scale=5){
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  mat <- mat/scale
  cc <- (1+sqrt(5)*mat+5/3*mat^2)*exp(-sqrt(5)*mat)
  return(cc)
}

head(cc)
  
d1 = dist(1:20)

# a random walk / brownian motion
locs=expand.grid(1:3, 1:3)
locs <- expand.grid(1:20,1:20)
locs <- cbind(locs[,1],locs[,2])
cbm <- covbm(locs)
dim(cbm)
cbm
z <- rmultnorm(1,rep(0,20*20),cbm)
z <- matrix(z,ncol=20)
persp(z)
z <- rmultnorm(1,rep(0,50),cbm)
z <- matrix(z,ncol=10)
persp(z, theta=60)
# fractal brownian motion
cbm <- covbm(locs,pow=1.8)
z <- rmultnorm(1,rep(0,20*20),cbm)
z <- matrix(z,ncol=20)
persp(z)

brownian= covbm(locs, pow=1.8)

cga2=nearPD(cga)
cga=cga2$mat
solve(cga)
# gaussian surface
cga <- covpow(locs,pow=2,scale=6)

z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)


persp(z)

# exponentail surface
cga <- covpow(locs,pow=1,scale=12)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
persp(z)

# now make some figs for the notes
# exponential
par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covpow(locs,pow=1,scale=15)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - exponential covariogram",side=3,line=0,outer=T)

# Gaussian
par(mfrow = c(1,1))
par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covpow(locs,pow=2,scale=6)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
#par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - Gaussian covariogram",side=3,line=0,outer=T)

# p = 1.5
par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covpow(locs,pow=1.5,scale=9)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - C(d) = exp(-d^1.5)",side=3,line=0,outer=T)

# BM p = 1
par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covbm(locs,pow=1.0,scale=5)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
#par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - Brownian motion p=1",side=3,line=0,outer=T)

par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covbm(locs,pow=1.3,scale=1)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - Brownian motion p=1.3",side=3,line=0,outer=T)

par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covbm(locs,pow=1.7,scale=1)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - Brownian motion p=1.7",side=3,line=0,outer=T)

par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covbm(locs,pow=1.9,scale=1)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - Brownian motion p=1.9",side=3,line=0,outer=T)

# spherical

par(mfrow=c(1,2),oma=c(18,0,2,0))
cga <- covsph(locs,scale=15)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mar=c(0,4,0,0))
persp(z)
par(mar=c(7,4,7,2))
plot(z[,1],ylab="z",xlab="x",ylim=range(z),type="l")
mtext("random realization - spherical",side=3,line=0,outer=T)



# Some interpolation stuff:
# compute the distribution of x2 given x1

#source("/moses/local/higdon/Sfiles/svdinv.s")
condgau <- function(x,mu,cov,n1){
# computes the conditional distribution of x2 = x[n1+1..n1+n2] given
# the observed values of x1 = x[1..n1].  It is assumed that the
# first n1 values of x hold the observed values.  The function
# returns a vector of length n1+n2, where the first n1 values
# are as were entered in the function.
# browser()
 N <- length(x); n2 <- N-n1
 x1 <- x[1:n1]
 mu1 <- mu[1:n1]
 mu2 <- mu[(n1+1):N]
 s11 <- cov[1:n1,1:n1]
 s22 <- cov[(n1+1):N,(n1+1):N]
 s12 <- cov[1:n1,(n1+1):N]
 s21 <- t(s12)
 s11i <- svdinv(s11)
 #svds11$values <- ifelse(svds11$values < 1e-8,0,svds11$values)
 s22o1 <- s22-s21%*%s11i%*%s12
 mu22o1 <- mu2+s21%*%s11i%*%(x1-mu1)
 return(list(mean=mu22o1,cov=(s22o1 + t(s22o1))/2))
}

st <- covbm(cbind(1:5,rep(0,5)))
mu <- rep(0,5)
x <- c(0,1,0,0,0)
n1 <- 2

# Show conditioning on an edge
# Gaussian
cga <- covpow(locs,pow=2,scale=6)
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mfrow=c(2,2),mar=c(0,4,0,0))
persp(z)
mtext("a realization",side=3,line=-2)

a <- condgau(as.vector(z),rep(0,20*20),cga,20)
cm <- c(as.vector(z)[1:20],a$mean)
cm <- matrix(cm,ncol=20)
persp(cm,zlim=range(z))
mtext("mean conditional on Y=1 points",side=3,line=-2)
cz <- rmultnorm(1,a$mean,a$cov)
cz <- c(z[,1],cz)
cz <- matrix(cz,ncol=20)
persp(cz,zlim=range(z))
mtext("realization conditional on Y=1 points",side=3,line=-2)

# Brownian motion
locs <- expand.grid(1:20,1:20)
locs <- cbind(locs[,1],locs[,2])
cga <- covbm(locs,pow=1.5,scale=3)+1
z <- rmultnorm(1,rep(0,20*20),cga)
z <- matrix(z,ncol=20)
par(mfrow=c(2,2),mar=c(0,4,0,0))
persp(z,zlim=range(z)*1.5)
mtext("a realization",side=3,line=-2)
a <- condgau(as.vector(z),rep(0,20*20),cga,20)
cm <- c(as.vector(z)[1:20],a$mean)
cm <- matrix(cm,ncol=20)
persp(cm,zlim=range(z)*1.5)
mtext("mean conditional on Y=1 points",side=3,line=-2)
cz <- rmultnorm(1,a$mean,a$cov)
cz <- c(z[,1],cz)
cz <- matrix(cz,ncol=20)
persp(cz,zlim=range(z)*1.5)
mtext("realization conditional on Y=1 points",side=3,line=-2)

# some simple 1-d examples
#browser()

zx <- c(3,6,9,12,seq(-5,20,length=101))
zy <- c(0,0,0,0,rep(0,101))
locs <- cbind(zx,zy)
z <- c(0,1,1,0,rep(0,101))
#z <- c(-1.5,0,0,1.5,rep(0,101))

par(mfcol=c(3,3),mar=c(5,5,3,1))
zl <- c(-3.0,3.0)
 # Gaussian
par(mfrow = c(1,1))
cv <- covpow(locs,pow=2,scale=1)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Gaussian C(r), scale = 2",side=3,line=1)

cv <- covpow(locs,pow=2,scale=3)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Gaussian C(r), scale = 3",side=3,line=1)

cv <- covpow(locs,pow=2,scale=5)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Gaussian C(r), scale = 5",side=3,line=1)

 # Exponential
cv <- covpow(locs,pow=1,scale=1)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Exponential C(r), scale = 1",side=3,line=1)

cv <- covpow(locs,pow=1,scale=10)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Exponential C(r), scale = 10",side=3,line=1)

cv <- covpow(locs,pow=1,scale=20)
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Exponential C(r), scale = 20",side=3,line=1)

 # Brownian motion p=1.5
cv <- covbm(locs,pow=1.5,scale=1)+1
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Brownian motion C(r), p = 1.5 scale = 1",side=3,line=1)

cv <- covbm(locs,pow=1.5,scale=3)+1
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Brownian motion C(r), p = 1.5 scale = 3",side=3,line=1)

cv <- covbm(locs,pow=1.5,scale=5)+1
a <- condgau(z,rep(0,length(z)),cv,4)
zr <- rmultnorm(4,a$mean,a$cov)
matplot(zx[-(1:4)],t(zr),type="l",xlab="x",ylab="z",ylim=zl)
points(zx[1:4],z[1:4],cex=1.6)
lines(zx[-(1:4)],a$mean,lwd=2)
mtext("Brownian motion C(r), p = 1.5 scale = 5",side=3,line=1)


 # now a plot for a ref report
f18 <- function(s) return( exp(-1.4*s)*cos(7*pi*s/2) )
spred <- seq(0,1,length=100)
idat <- c(9,21,37,49,64,78,94)
sdat <- spred[idat]
locs <- cbind(spred,0)
zpred <- f18(spred)
z <- c(zpred[idat],zpred[-idat])

cvpred <- covbm(locs,pow=1.75,scale=1)+1
cv <- rbind( cbind(cvpred[idat,idat],cvpred[idat,-idat]),
             cbind(cvpred[-idat,idat],cvpred[-idat,-idat]))
a <- condgau(z,rep(0,length(z)),cv,length(sdat))

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(5,5,3,1))
cvpred <- covbm(locs,pow=1.75,scale=1)+1
cv <- rbind( cbind(cvpred[idat,idat],cvpred[idat,-idat]),
             cbind(cvpred[-idat,idat],cvpred[-idat,-idat]))
a <- condgau(z,rep(0,length(z)),cv,length(sdat))
plot(spred,f18(spred),type='l',xlab='x',ylab='z(x)')
points(sdat,f18(sdat),pch=16,cex=1.4)
lines(spred[-idat],a$mean,col='green',lty=2)
mtext('variogram: |d|^(1.75)',side=3,line=.5)
lines(spred,spred^1.75,col='blue',lty=3)

cvpred <- covbm(locs,pow=1.75,scale=.2)+1
cv <- rbind( cbind(cvpred[idat,idat],cvpred[idat,-idat]),
             cbind(cvpred[-idat,idat],cvpred[-idat,-idat]))
a <- condgau(z,rep(0,length(z)),cv,length(sdat))
plot(spred,f18(spred),type='l',xlab='x',ylab='z(x)')
points(sdat,f18(sdat),pch=16,cex=1.4)
lines(spred[-idat],a$mean,col='green',lty=2)
mtext('variogram: |d*5|^(1.75)',side=3,line=.5)
lines(spred,(spred*5)^1.75,col='blue',lty=3)

cvpred <- covbm(locs,pow=1.75,scale=.02)+1
cv <- rbind( cbind(cvpred[idat,idat],cvpred[idat,-idat]),
             cbind(cvpred[-idat,idat],cvpred[-idat,-idat]))
a <- condgau(z,rep(0,length(z)),cv,length(sdat))
plot(spred,f18(spred),type='l',xlab='x',ylab='z(x)')
points(sdat,f18(sdat),pch=16,cex=1.4)
lines(spred[-idat],a$mean,col='green',lty=2)
mtext('variogram: |d*50|^(1.75)',side=3,line=.5)
lines(spred,(spred*50)^1.75,col='blue',lty=3)

# a plot for the gpmsa manual
 # now a plot for a ref report
f18 <- function(s) return( exp(-1.4*s)*cos(7*pi*s/2) )
spred <- seq(0,1,length=100)
idat <- c(11,28,48,67,92)
sdat <- spred[idat]
locs <- cbind(spred,0)
zpred <- f18(spred)
#cvpred <- covpow(locs,pow=2.0,scale=.3)
#zpred <- t(rmultnorm(1,0*spred,cvpred))
z <- c(zpred[idat],zpred[-idat])

par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(5,5,3,1))
cvpred <- covpow(locs,pow=2.0,scale=.3)
cv <- rbind( cbind(cvpred[idat,idat],cvpred[idat,-idat]),
             cbind(cvpred[-idat,idat],cvpred[-idat,-idat]))
a <- condgau(z,rep(0,length(z)),cv,length(sdat))
plot(spred,f18(spred),type='l',xlab='x',ylab='z(x)')
points(sdat,f18(sdat),pch=16,cex=1.4)
lines(spred[-idat],a$mean,col='gray',lty=2)

postscript('gpeg.ps',width=5,height=4,horizontal=F,paper='special')
zr <- rmultnorm(3,a$mean,a$cov)
matplot(spred[-idat],t(zr),col=c('cyan'),lty=1,type='l',xlab='x',
  ylab='z(x)')
lines(spred,f18(spred),type='l',lwd=2)
lines(spred[-idat],a$mean,col='red',lty=2)
points(sdat,f18(sdat),pch=16,cex=1.4)
title('GP emulation of f(x) = exp(-1.4x)*cos(7*pi*x/2)')
dev.off()





