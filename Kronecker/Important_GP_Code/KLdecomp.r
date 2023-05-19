# a look at the KL decomposition...
npred <- 20    # of data points
m <- 20    # num of bumps
spred <- seq(0,10,length=npred)
locs <- cbind(spred,0)
Ssim <- covpow(locs,pow=2,scale=5) #make cov matrix w/ Gaussian cov function w/ scale =5

#pdf('truecov.pdf',height=3.5,width=3.5)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Ssim,zlim=c(0,1),xlab='location s',
      ylab='location s')
#dev.off(); # system('gv truecov.ps')


par(mfrow=c(4,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
# KL/svd decomposition
Ssvd <- svd(Ssim) # do svd decomp on our covariance matrix
SqrtSsim <- Ssvd$u %*% diag(sqrt(Ssvd$d)) #get the 'sqr root' of our cov matrix

#postscript('klbasis.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,SqrtSsim,type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
#dev.off(); # system('gv klbasis.ps')
#postscript('klbasisd.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,SqrtSsim,type='p',pch=20,cex=1.5,ylim=c(-1,1),ylab='basis',xlab='s')
#dev.off(); # system('gv klbasisd.ps')

#postscript('klcov.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
MyCov = SqrtSsim[,1:5]%*%t(SqrtSsim[,1:5]) #only take the first five cols--which have largest values, approx w/ those
image(spred,spred,SqrtSsim[,1:5]%*%t(SqrtSsim[,1:5]),zlim=c(0,1),
     xlab='location s', ylab='location s')
#dev.off(); # system('gv klcov.ps')

# chol (without pivoting)
ChS0 <- chol(Ssim+diag(rep(.000001,npred)))
matplot(spred,t(ChS0),type='l',ylim=c(-1,1),ylab='cholesky (no pivoting)')
#postscript('chbasis.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,t(ChS0),type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
#dev.off(); # system('gv chbasis.ps')
Ch0 = t(ChS0)

#postscript('chcov.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Ch0[,1:5]%*%t(Ch0[,1:5]),zlim=c(0,1),
     xlab='location s', ylab='location s')
#dev.off(); # system('gv chcov.ps')

# chol with pivoting
?order
ChS <- chol(Ssim+diag(rep(.000001,npred)),pivot=T)
ipivot <- attr(ChS,'pivot')
oo <- order(ipivot)
Chpiv = t(ChS[,oo])
matplot(spred,t(ChS[,oo]),type='l',ylim=c(-1,1),ylab='cholesky w/pivoting')
#postscript('chpivbasis.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,t(ChS[,oo]),type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
#dev.off(); # system('gv chpivbasis.ps')

#postscript('chpivcov.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Chpiv[,1:5]%*%t(Chpiv[,1:5]),zlim=c(0,1),
     xlab='location s', ylab='location s')
#dev.off(); # system('gv chpivcov.ps')

# Fourier basis - here 2pi*10 is the period
phipred <- 2*pi*spred[-npred]/10
nbasis <- 49;  #(nbasis+1 cos terms; nbasis sin terms)
mbasis <- 2*nbasis+1;
Fpred <- matrix(NA,ncol=mbasis,nrow=npred-1)
dim(Fpred)
#odd cols are cos's , evens are sin's

for(k in 0:nbasis){
  ii <- 2*k+1
  #ii <- 2*k
  Fpred[,ii] <- cos(phipred*k)
}
for(k in 1:nbasis){
  ii <- 2*k
  Fpred[,ii] <- sin(phipred*k)
}


matplot(spred[-npred],Fpred[,1:10],type='l',ylim=c(-1,1),ylab='Fourier')

# now use this basis to approximate the sqrt of the cov matrix
dim(Fpredinv)
dim(Fpred)
Fpredinv <- solve(Fpred[,(1:npred-1)]) #Fpred has only n-1 rows...as required in DFT
Df <- Fpredinv%*%Ssim[-npred,-npred]%*%t(Fpredinv)  #K^TCK  = D, from my notes

dDsqrt <- diag(sqrt(diag(Df)))
Fpred2 <- Fpred[,(1:npred-1)]%*%dDsqrt

matplot(spred[-npred],Fpred2,type='l',ylim=c(-1,1),ylab='Fourier')
testcov= (Fpred2 %*% t(Fpred2))
xvals=rnorm(19)
r0=Fpred2%*%xvals
plot(as.vector(r0), type = "l")

test=complex(imaginary = 1)
n=4
Fourier_matrix = function(n){
  mymatrix = matrix(NA, nrow = n, ncol= n)
  omega_j = function(j,k){
    exp((2*pi*complex(imaginary = 1) / n)*j*k)
  }
  for (j in 1:n){
    for (k in 1:n){
      mymatrix[j, k] = omega_j((j-1), (k-1))
    }
  }
  return(mymatrix)
}


A= c(2, -1, 0, 0, 0, 0, -1, -1, 2, -1, 0, 0, 0, 0,0, -1, 2, -1, 0, 0, 0,   0, 0, -1, 2, -1, 0, 0,
    0, 0, 0, -1, 2, -1, 0,0 ,0, 0, 0, -1, 2, -1,-1, 0, 0, 0, 0, -1, 2)
A=matrix(A, ncol = 7, nrow=7)

Fmat = Fourier_matrix(4)
result=round(Fmat %*% Lambda %*% solve(Fmat),2)

Lambda_real=diag(abs(diag(Lambda)))

Asim = Fmat %*% Lambda_real

j=1:20
n=20
k=0
k=3
test=cos((2*pi/n)*j*k)
plot(test, ylim=c(-1,1))
lines(test)
round(solve(Fmat) %*% A %*% Fmat, 3) 
Lambda=(solve(Fmat) %*% A %*% Fmat) 
my_eigens=eigen(A)
my_eigens$values
plot(A)
plot(Ssim)
# some plots of the square root matricies
postscript('cholfig.ps',height=5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),oma=c(1,4,4,1),mar=c(0,0,0,0),pty='s')
plotbasismat(t(ChS0)); mtext('rows',side=2,line=3,outer=T)
mtext('columns',side=3,line=3,outer=T)
dev.off()   # system('gv cholfig.ps')

postscript('cholpivfig.ps',height=5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),oma=c(1,4,4,1),mar=c(0,0,0,0),pty='s')
plotbasismat(t(ChS[,oo])); mtext('rows',side=2,line=3,outer=T)
mtext('columns',side=3,line=3,outer=T)
dev.off()   # system('gv cholpivfig.ps')

postscript('klfig.ps',height=5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),oma=c(1,4,4,1),mar=c(0,0,0,0),pty='s')
plotbasismat(SqrtSsim); mtext('rows',side=2,line=3,outer=T)
mtext('columns',side=3,line=3,outer=T)
dev.off()   # system('gv klfig.ps')

plotbasismat(SqrtSsim)

# normal kernel basis - 20 knot locations
sdbump <- 5/2
sbump <- seq(0-2,10+2,length=m)
Xpred <- matrix(NA,ncol=m,nrow=npred)
for(i in 1:m){
  Xpred[,i] <- dnorm(spred,mean=sbump[i],sd=sdbump)
}
intK2 <- max(Xpred%*%t(Xpred))
Xpred <- Xpred/sqrt(intK2)
postscript('kernbasis20.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,Xpred,type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
dev.off(); # system('gv kernbasis.ps')
postscript('kerncov20.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Xpred%*%t(Xpred),zlim=c(0,1),xlab='location s',
      ylab='location s')
dev.off(); # system('gv kerncov.ps')

# normal kernel basis - 10 knot locations
m <- 10; sdbump <- 5/2
sbump <- seq(0-2,10+2,length=m)
Xpred <- matrix(NA,ncol=m,nrow=npred)
for(i in 1:m){
  Xpred[,i] <- dnorm(spred,mean=sbump[i],sd=sdbump)
}
intK2 <- max(Xpred%*%t(Xpred))
Xpred <- Xpred/sqrt(intK2)
postscript('kernbasis10.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,Xpred,type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
dev.off(); # system('gv kernbasis10.ps')
postscript('kerncov10.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Xpred%*%t(Xpred),zlim=c(0,1),xlab='location s',
      ylab='location s')
dev.off(); # system('gv kerncov10.ps')

# normal kernel basis - 6 knot locations
m <- 6; sdbump <- 5/2
sbump <- seq(0-2,10+2,length=m)
Xpred <- matrix(NA,ncol=m,nrow=npred)
for(i in 1:m){
  Xpred[,i] <- dnorm(spred,mean=sbump[i],sd=sdbump)
}
intK2 <- max(Xpred%*%t(Xpred))
Xpred <- Xpred/sqrt(intK2)
postscript('kernbasis6.ps',height=3.5,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(spred,Xpred,type='l',lty=1,lwd=2,ylim=c(-1,1),ylab='basis',xlab='s')
dev.off(); # system('gv kernbasis6.ps')
postscript('kerncov6.ps',height=3.5,width=3.5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='s')
image(spred,spred,Xpred%*%t(Xpred),zlim=c(0,1),xlab='location s',
      ylab='location s')
dev.off(); # system('gv kerncov6.ps')

#% a plot to show kernel locations and form
sshow = seq(-2.5,12.5,length=100)
m <- 6; sdbump <- 5/2
sbump <- seq(0-2,10+2,length=m)
Xpred <- matrix(NA,ncol=m,nrow=length(sshow))
for(i in 1:m){
  Xpred[,i] <- dnorm(sshow,mean=sbump[i],sd=sdbump)
}
intK2 <- max(Xpred%*%t(Xpred))
Xpred <- Xpred/sqrt(intK2)
postscript('kernbasis6s.ps',height=3.0,width=5,paper='special',horizontal=F)
par(mfrow=c(1,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(sshow,Xpred,type='l',lty=1,lwd=2,ylab='basis', xlab='s')
matpoints(rbind(sbump,sbump),rbind(sbump,sbump)*0,pch=20,cex=2.2)
matlines(rbind(sbump,sbump),rbind(0*sbump,sbump*0+max(Xpred)),lty=1,lwd=2)
dev.off(); # system('gv kernbasis6s.ps')

dim(Xpred)
dim(xbasis)
# now show a single realization...
par(mfrow = c(1,2))
x = rnorm(m)
xbasis = Xpred*matrix(x,ncol=m,nrow=length(Xpred[,1]),byrow=T)

#par(mfrow=c(1,2),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
matplot(sshow,xbasis,type='l',lty=1,lwd=2, ylab='basis',xlab='s')
matpoints(rbind(sbump,sbump),rbind(sbump,sbump)*0,pch=20,cex=2.2)
kjmax = apply(abs(xbasis),2,max)*sign(x);
matlines(rbind(sbump,sbump),rbind(0*sbump,sbump*0+kjmax),lty=1,lwd=2)
#x = rnorm(m)
plot(sshow,Xpred%*%x,type='l',lty=1,lwd=2, ylab='z(t)',xlab='s')
dev.off(); # system('gv kernbasis6z.ps')




npred
# moving ave spec...
spred <- seq(0,10,length=npred)
sdbump <- 5/2
sbump <- seq(0-2,10+2,length=m)
Xpred <- matrix(NA,ncol=m,nrow=npred)
for(i in 1:m){
  Xpred[,i] <- dnorm(spred,mean=sbump[i],sd=sdbump)
}
intK2 <- max(Xpred%*%t(Xpred))
Xpred <- Xpred/sqrt(intK2)
matplot(spred,Xpred,type='l',ylim=c(-1,1),ylab='moving average basis')

# now look at the number of components required to fit the data
# given in y10...
?scan
xxxdat <- matrix(scan('statmodeg.dat'),ncol=3,byrow=T)
s <- xxxdat[,1]; y0 <- xxxdat[,2]; y <- xxxdat[,3]; e1 <- y - y0;
par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,0,0),pty='m')
plot(s,y,pch=1); lines(s,y0)
# model y = Lx + e    lamx=1; lamy=1/.2^2  
# posterior for x: Px = lamyLL' + I; mux = solve(Px,lamyL'y)
idat <- seq(1,100,by=9) #?? obv will not work to subset SqrtSsim below

ly <- 1/.2^2

# svd decomp
K = SqrtSsim
K <- SqrtSsim[idat,] #sqrt of the cov matrix
Psvd <- t(K)%*%K*ly + diag(rep(1,npred))
msvd <- solve(Psvd,ly*t(K)%*%y)

lines(spred,SqrtSsim%*%msvd,col='red')

# symmetric chol
K <- t(ChS[,oo])
K <- K[idat,]
Pch <- t(K)%*%K*ly + diag(rep(1,npred))
mch <- solve(Pch,ly*t(K)%*%y)
lines(spred,t(ChS[,oo])%*%mch,col='green')

# moving ave 
K <- Xpred[idat,]
Pma <- t(K)%*%K*ly + diag(rep(1,m))
mma <- solve(Pma,ly*t(K)%*%y)
lines(spred,Xpred%*%mma,col='blue')

# Show how the overall variance approximation changes with the
# number of bisis elements used in thesvd decomp
K <- SqrtSsim[idat,]
Psvd <- t(K)%*%K*ly + diag(rep(1,npred))
msvd <- solve(Psvd,ly*t(K)%*%y)
lines(spred,SqrtSsim%*%msvd,col='red')
# lower dim representation...
par(mfrow=c(4,5),pty='s',mar=c(.2,.2,.2,.2))
image(vtrue,zlim=zlim)
for(kk in c(1:4,seq(5,by=1,length=14),100)){
bdim <- kk; Pld <- t(K[,1:bdim])%*%K[,1:bdim]*ly + diag(rep(1,bdim))
vld <- SqrtSsim[,1:bdim]%*%solve(Pld)%*%t(SqrtSsim[,1:bdim])
vtrue <- SqrtSsim%*%solve(Psvd)%*%t(SqrtSsim)
zlim <- range(cbind(vtrue,vld))
#image(vld,zlim=zlim); text(.1,.9,kk,cex=1)
image(vld-vtrue,zlim=c(-.02,.02)); text(.1,.9,kk,cex=1)
}

# a function to plot representations of chol/basis matricies...
plotbasismat <- function(X,tol=.01,...){
 dimX <- dim(X)
 n <- dimX[1]; p <- dimX[2]
 plot(c(1,p+1),c(1,n),axes=F,type='n',xlab='',ylab=''); 
 box(); axis(2,at=c(16,11,6,1),labels=c('5','10','15','20')); axis(3);
 ym <- rbind(1:n,1:n);
 for(k in 1:p){
   xm <- rbind(rep(k,n),X[,k]*.85+rep(k,n))
   matlines(xm[,n:1],ym,lty=1,col='black')
 }
}

