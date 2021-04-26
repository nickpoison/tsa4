library('astsa')
source('grid.r')


##################################################################
##################################################################
pdf(file="bone.pdf", width=7.6, height=4)
par(mfrow=c(3,1), mar = c(0, 3.5, 0, 2), oma=c(3,0,2,0), mgp=c(2,.6,0), cex.lab=1.5, las=1, tcl=-.3)
plot.ts(blood[,1], ylab='WBC', xaxt="no", type='n')
grid(lty=1); 
lines(blood[,1], type='o', pch=19)
plot.ts(blood[,2], ylab='PLT', yaxt='no', xaxt="no", type='n')
grid(lty=1); axis(4) 
lines(blood[,2], type='o', pch=19)
plot.ts(blood[,3], ylab='HCT', type='n', ylim=c(21.5,36.5))
grid(lty=1); 
lines(blood[,3], type='o', pch=19)
title(xlab="Time", outer=TRUE)
dev.off()


#################################################################
###################################################################
pdf(file="gtemp2.pdf", width=7.6, height=3)
par(mar=c(2,2.5,.5,.5)+.5, mgp=c(1.6,.6,0))
ts.plot(globtemp, globtempl, ylab="Temperature Deviations", xlab="Time", main='', type='n')
grid(lty=1, col=gray(.9))
lines(globtemp,  lwd=2, col = rgb(.9,  0, .7, .8) )
lines(globtempl, lwd=2, col = rgb( 0, .7, .9, .7) )
legend('topleft', col=c(rgb(.9, 0, .7),  rgb(0, .7, .9)), lty=1, lwd=2, 
        legend=c("Land/Ocean", "Land Only"), bg='white')  
dev.off()


###########################################################
#########################################################
pdf(file="llmall.pdf", width=7.6, height=5.25)
par(mfrow=c(3,1), mar=c(1.5,2.5,1,0)+.5, mgp=c(1.6,.6,0))
# generate data
set.seed(1);  num = 50
w = rnorm(num+1,0,1); v = rnorm(num,0,1)
mu = cumsum(w)  # state:  mu[0], mu[1],..., mu[50]
y = mu[-1] + v  #   obs:  y[1],..., y[50]
# filter and smooth
#mu0 = 0; sigma0 = 1;  phi = 1; cQ = 1; cR = 1
ks = Ksmooth0(num, y, A=1, mu0=0, Sigma0=1, Phi=1, cQ=1, cR=1)
# start figure
Time = 1:num
plot(Time, mu[-1], main="Predict", ylim=c(-5,10), ylab=expression(mu[~t]), xlab='', panel.first=grid(lty=1))
  lines(ks$xp)
  lines(ks$xp+2*sqrt(ks$Pp), lty="dashed", col="blue")
  lines(ks$xp-2*sqrt(ks$Pp), lty="dashed", col="blue")
plot(Time, mu[-1], main="Filter", ylim=c(-5,10), ylab=expression(mu[~t]),xlab='', panel.first=grid(lty=1))
  lines(ks$xf)
  lines(ks$xf+2*sqrt(ks$Pf), lty="dashed", col="blue")
  lines(ks$xf-2*sqrt(ks$Pf), lty="dashed", col="blue")
plot(Time, mu[-1],  main="Smooth", ylim=c(-5,10), ylab=expression(mu[~t]), xlab='', panel.first=grid(lty=1))
  lines(ks$xs)
  lines(ks$xs+2*sqrt(ks$Ps), lty="dashed", col="blue")
  lines(ks$xs-2*sqrt(ks$Ps), lty="dashed", col="blue")
mtext(side=1, line=1, 'Time', cex=.8)  
#mu[1]; ks$x0n; sqrt(ks$P0n)   # initial value info
dev.off()

#################################################################
##########  newton mle example (no pic but output) #############
#################################################################
# Generate Data
set.seed(999); num = 100
x = arima.sim(n=num+1, list(ar = .8), sd=1)
y = ts(x[-1] + rnorm(num,0,1))
# Initial Estimates
u = ts.intersect(y, lag(y,-1), lag(y,-2))
varu = var(u); coru = cor(u)
phi = coru[1,3]/coru[1,2]
q = (1-phi^2)*varu[1,2]/phi;  r = varu[1,1] - q/(1-phi^2)
(init.par = c(phi, sqrt(q), sqrt(r)))  # = .91, .51, 1.03
# Function to evaluate the likelihood
Linn=function(para){
  phi = para[1]; sigw = para[2]; sigv = para[3]
  Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0
  kf = Kfilter0(num, y, 1, mu0=0, Sigma0, phi, sigw, sigv)
  return(kf$like)   }
# Estimation   (partial output shown)
(est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1, REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
round(cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]),SE), 3)


################################################################
########  glob temps ########################################
y = cbind(globtemp,globtempl); num = nrow(y); input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = -.3; Sigma0 = 1;  Phi = 1
# Function to Calculate Likelihood
Linn=function(para){
 cQ = para[1]      # sigma_w
  cR1 = para[2]    # 11 element of chol(R)
  cR2 = para[3]    # 22 element of chol(R)
  cR12 = para[4]   # 12 element of chol(R)
 cR = matrix(c(cR1,0,cR12,cR2),2)  # put the matrix together
 drift = para[5]
 kf = Kfilter1(num,y,A,mu0,Sigma0,Phi=1,drift,0,cQ,cR,input)
 return(kf$like) }
# Estimation
init.par = c(.1,.1,.1,0,.05)  # initial values of parameters
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
# display estimates
u = cbind(estimate=est$par, SE)
rownames(u)=c("sigw","cR11", "cR22", "cR12", "drift"); 
round(u,3)  
# Smooth (first set parameters to their final estimates)
cQ=est$par[1]
 cR1=est$par[2]
 cR2=est$par[3]
 cR12=est$par[4]
cR = matrix(c(cR1,0,cR12,cR2), 2)
(R = t(cR)%*%cR)    #  to view the estimated R matrix
drift = est$par[5]
ks = Ksmooth1(num,y,A,mu0,Sigma0,Phi=1,drift,0,cQ,cR,input)
# Plot
xsm = ts(as.vector(ks$xs), start=1880)
	xsm=window(xsm, start=1900)
rmspe = ts(sqrt(as.vector(ks$Ps)), start=1880)
     rmspe=window(rmspe, start=1900)
 pdf(file="gtempsmoo.pdf", width=7.6, height=4)
par(mar=c(2,2,0,.5)+.5, mgp=c(1.6,.6,0))
plot(xsm,  ylim=c(-.6,1), ylab="Temperature Deviations", type='n')
grid(lty=1)
  xx=c(time(xsm), rev(time(xsm)))
  yy=c(xsm-2*rmspe, rev(xsm+2*rmspe))
  polygon(xx, yy, border=NA, col=rgb(.6,.6,.6,alpha=.3)) 
#points(gtemp, pch=15, col=4) 
#points(gtemp2, pch=19, col=3)
lines(window(globtemp, start=1900), type='o', pch=2, col=rgb(0,0,.8), lty=6)   # color helps here
lines(window(globtempl, start=1900), type='o',  pch=3, col=rgb(0,.65,0), lty=6)
lines(xsm)
 dev.off()
 
 
########################################
#########  EM no pic #########################
###############################################
library(nlme)   # loads package nlme
# Generate data (same as Example 6.6)
set.seed(999); num = 100
x = arima.sim(n=num+1, list(ar = .8), sd=1)
y = ts(x[-1] + rnorm(num,0,1))
# Initial Estimates
u = ts.intersect(y, lag(y,-1), lag(y,-2))
varu = var(u); coru = cor(u)
phi = coru[1,3]/coru[1,2]
q = (1-phi^2)*varu[1,2]/phi
r = varu[1,1] - q/(1-phi^2)
# EM procedure - output not shown
(em = EM0(num, y, A=1, mu0=0, Sigma0=2.8, Phi=phi, cQ=sqrt(q), cR=sqrt(r), max.iter=75, tol=.00001))
# Standard Errors  (this uses nlme)
phi = em$Phi; cq = sqrt(em$Q); cr = sqrt(em$R)
mu0 = em$mu0; Sigma0 = em$Sigma0
para = c(phi, cq, cr)
Linn = function(para){  # to evaluate likelihood at estimates
  kf = Kfilter0(num, y, 1, mu0, Sigma0, para[1], para[2], para[3])
  return(kf$like)  }
emhess = fdHess(para, function(para) Linn(para))
SE = sqrt(diag(solve(emhess$Hessian)))
# Display Summary of Estimation
estimate = c(para, em$mu0, em$Sigma0); SE = c(SE, NA, NA)
u = cbind(estimate, SE)
rownames(u) = c("phi","sigw","sigv","mu0","Sigma0") 
round(u,3) 


#################################
# bonesmoo #####################
################################
y = cbind(WBC, PLT, HCT)
num = nrow(y)
A = array(0, dim=c(3,3,num))
# make array of obs matrices
for(k in 1:num) { if (y[k,1] > 0) A[,,k]= diag(1,3) }
# Initial values
mu0 = matrix(0, 3, 1)
Sigma0 = diag(c(.1, .1, 1), 3)
Phi = diag(1, 3)
cQ = diag(c(.1, .1, 1), 3)
cR = diag(c(.1, .1, 1), 3)
# EM procedure - some output previously shown
(em = EM1(num, y, A, mu0, Sigma0, Phi, cQ, cR, 100, .001))
# Graph smoother
ks = Ksmooth1(num, y, A, em$mu0, em$Sigma0, em$Phi, 0, 0, chol(em$Q), chol(em$R), 0)
y1s = ks$xs[1,,]; y2s = ks$xs[2,,]; y3s = ks$xs[3,,]
p1 = 2*sqrt(ks$Ps[1,1,])
p2 = 2*sqrt(ks$Ps[2,2,])
p3 = 2*sqrt(ks$Ps[3,3,])
#
# plot  #######
pdf(file="bonesmoo.pdf", width=7.6, height=4.5)
par(mfrow=c(3,1), mar = c(0, 3.5, 0, 2), oma=c(3,0,1,0), mgp=c(2,.6,0), cex.lab=1.5,las=1, tcl=-.3)
plot.ts(blood[,1], ylab='WBC', xaxt="no", type='n', ylim=c(1.5,4.5))
grid(lty=1); 
points(blood[,1], pch=19)
lines(y1s); 
  xx=c(time(blood), rev(time(blood)))  # same for all
  yy=c(y1s-p1, rev(y1s+p1))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .2))  
#
plot.ts(blood[,2], ylab='PLT', yaxt='no', xaxt="no", type='n', ylim=c(3.7,5.7))
grid(lty=1); axis(4) 
points(blood[,2], pch=19)
lines(y2s)
  yy=c(y2s-p2, rev(y2s+p2))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .2))  
#
plot.ts(blood[,3], ylab='HCT', type='n', ylim=c(20,40))
grid(lty=1); 
points(blood[,3], pch=19)
lines(y3s)
yy=c(y3s-p3, rev(y3s+p3))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .2)) 
title(xlab="Time", outer=TRUE, cex.lab=1.1) 
miss=ifelse(y[,1]==0,20.1,0)            # try to put a tick at missing days
lines(miss, type='h', lwd=2)
dev.off()


#############################################
###{smoothspline} #############################
###########################################
 pdf(file="smoothspline.pdf", width=7.6, height=3.5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.6,.6,0))
set.seed(123)
num=50
w = rnorm(num,0,.1)
x = cumsum(cumsum(w))
y = x +rnorm(num,0,1)
plot.ts(x, ylab="", type='n',lwd=2, ylim=c(-1,8))
grid(lty=1)
lines(x, lwd=2)
lines(y, type='o', col='#808080')
#title("Smoothing Spline")
## state space ################# 
Phi = matrix(c(2,1,-1,0),2)
A = matrix(c(1,0),1)
#mu0 = matrix(c(0,0),2)
#Sigma0= diag(c(1,1))
mu0=matrix(0,2)
Sigma0=diag(1,2)
Linn=function(para){
  sigw = para[1]; sigv = para[2]  
  cQ = diag(c(sigw,0))
  kf = Kfilter0(num,y,A,mu0,Sigma0,Phi,cQ,sigv)
  return(kf$like)   
  }
#estimation  
init.par=c(.1,1)  
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1))) 
SE = sqrt(diag(solve(est$hessian))) 
# Summary of estimation  
estimate = est$par; u = cbind(estimate, SE)
rownames(u)=c("sigw","sigv"); u  
# smooth
sigw=est$par[1]
cQ = diag(c(sigw,0))
sigv=est$par[2]
ks=Ksmooth0(num,y,A,mu0,Sigma0,Phi,cQ,sigv)
xsmoo = ts(ks$xs[1,1,])
psmoo = ts(ks$Ps[1,1,])
upp=xsmoo+2*sqrt(psmoo)
low=xsmoo-2*sqrt(psmoo)
lines(xsmoo, col=4, lty=2, lwd=3)  
lines(upp, col=4,lty=2)       
lines(low, col=4, lty=2)
lines(smooth.spline(y),lty=1,col=2)
legend("topleft", c("Observations","State"), pch=c(1,-1), lty=c(1,1), lwd=c(1,2), col=c('#808080',1), bg='white')
legend("bottomright", c("Smoother", "GCV Spline"),  lty=c(2,1), lwd=c(3,1), col=c(4,2), bg='white' )
dev.off()


#############################################
#  jjtrse and jjfore ######################
##########################################
require(astsa)
num = length(jj)
A = cbind(1,1,0,0)                                
# Function to Calculate Likelihood 
Linn=function(para){
 Phi = diag(0,4); Phi[1,1] = para[1] 
 Phi[2,]=c(0,-1,-1,-1); Phi[3,]=c(0,1,0,0); Phi[4,]=c(0,0,1,0)
 cQ1 = para[2]; cQ2 = para[3]     # sqrt q11 and sqrt q22
 cQ=diag(0,4); cQ[1,1]=cQ1; cQ[2,2]=cQ2
 cR = para[4]                     # sqrt r11
 kf = Kfilter0(num,jj,A,mu0,Sigma0,Phi,cQ,cR)
 return(kf$like)  
}
# Initial Parameters 
mu0 = c(.7,0,0,0);  Sigma0 = diag(.04,4)  
init.par = c(1.03,.1,.1,.5)   # Phi[1,1], the 2 Qs and R
# Estimation
est = optim(init.par,Linn,NULL,method="BFGS",hessian=TRUE,control=list(trace=1,REPORT=1))
SE = sqrt(diag(solve(est$hessian)))
u = cbind(estimate=est$par,SE)
rownames(u)=c("Phi11","sigw1","sigw2","sigv"); u     
# Smooth
Phi = diag(0,4); Phi[1,1] = est$par[1] 
Phi[2,]=c(0,-1,-1,-1); Phi[3,]=c(0,1,0,0); Phi[4,]=c(0,0,1,0)
cQ1 = est$par[2]; cQ2 = est$par[3]      
cQ = diag(1,4); cQ[1,1]=cQ1; cQ[2,2]=cQ2   
cR = est$par[4]   
ks = Ksmooth0(num,jj,A,mu0,Sigma0,Phi,cQ,cR)   
# Plot 
Tsm = ts(ks$xs[1,,], start=1960, freq=4)
Ssm = ts(ks$xs[2,,], start=1960, freq=4)
p1 = 3*sqrt(ks$Ps[1,1,]); p2 = 3*sqrt(ks$Ps[2,2,])
##############
pdf(file="jjtrse.pdf", width=7.6, height=5.5)
par(mfrow=c(2,1), mar=c(2.5,3,1.5,1), mgp=c(1.6,.6,0), cex.main=.9, cex.lab=.9)
plot(Tsm, main="Trend Component", ylab="Trend", type='n')
grid(lty=1)
lines(Tsm)
 xx=c(time(jj), rev(time(jj)))
 yy=c(Tsm-p1, rev(Tsm+p1)) 
polygon(xx, yy, border=NA, col=gray(.5, alpha = .3))
plot(jj, main="Data & Trend+Season", ylab="J&J QE/Share", type='n', ylim=c(-.5,17)) 
grid(lty=1)
lines(jj)
 xx = c(time(jj), rev(time(jj)))
 yy=c( (Tsm+Ssm)-(p1+p2), rev((Tsm+Ssm)+(p1+p2)) )
polygon(xx, yy, border=NA, col=gray(.5, alpha = .3))
dev.off()
#########################
# Forecast 
n.ahead=12; y = ts(append(jj, rep(0,n.ahead)), start=1960, freq=4)
rmspe = rep(0,n.ahead); x00 = ks$xf[,,num]; P00 = ks$Pf[,,num]
Q=t(cQ)%*%cQ;  R=t(cR)%*%(cR)
for (m in 1:n.ahead){ 
 xp = Phi%*%x00; Pp = Phi%*%P00%*%t(Phi)+Q
 sig = A%*%Pp%*%t(A)+R; K = Pp%*%t(A)%*%(1/sig)
 x00 = xp; P00 = Pp-K%*%A%*%Pp
 y[num+m] = A%*%xp; rmspe[m] = sqrt(sig)  
}   
pdf(file="jjfore.pdf", width=7.6, height=3.5)
par(mar=c(2.5,3,1,1), mgp=c(1.6,.6,0), cex=.9)
plot(y, type="n", main="", ylab="J&J QE/Share", ylim=c(5,30), xlim=c(1975,1984))
grid(lty=1)
lines(y, type="o")
upp = ts(y[(num+1):(num+n.ahead)]+2*rmspe, start=1981, freq=4)
low = ts(y[(num+1):(num+n.ahead)]-2*rmspe, start=1981, freq=4)
#lines(upp, lty=2);  lines(low, lty=2);  
 xx = c(time(low), rev(time(upp)))
 yy = c(low, rev(upp))
polygon(xx, yy, border=8,   col=gray(.5, alpha = .3))
abline(v=1981, lty=3)  
 dev.off() 


###################################################
##  ARMAX  ##################
##################################
# intial stuff
fit1 = sarima(cmort, 2,0,0, xreg=time(cmort)) 
acf(cbind(dmort <- resid(fit1$fit), tempr, part))
lag2.plot(tempr, dmort, 8)  
lag2.plot(part, dmort, 8) 

# easy method 1: detrend cmort then do the regression
trend = time(cmort) - mean(time(cmort))   # center time
dcmort=resid(fit2<-lm(cmort~trend,na.action=NULL))
u = ts.intersect(dM=dcmort, dM1=lag(dcmort,-1), dM2=lag(dcmort,-2),  T1=lag(tempr,-1), P=part, P4=lag(part,-4))
## summary(lm(dM~., data=u, na.action=NULL)) # could use lm, but sarima also gives residual analysis
sarima(u[,1],0,0,0,xreg=u[,2:6])  
## Coefficients:
##       intercept     dM1     dM2       T1       P      P4
##          5.9884  0.3164  0.2989  -0.1826  0.1107  0.0495
## s.e.     2.6401  0.0370  0.0395   0.0309  0.0177  0.0195
## sigma^2 estimated as 25.42:  log likelihood = -1530.45,  aic = 3074.91
###
### summary(fit2)
### Coefficients:
###              Estimate Std. Error t value    
### (Intercept) 3297.6062   276.3132   11.93    
### time(cmort)   -1.6249     0.1399  -11.61

## another method uing different model where trend is exogenous
trend = time(cmort)-mean(time(cmort))
u = ts.intersect(M=cmort, M1=lag(cmort,-1), M2=lag(cmort,-2), T1=lag(tempr,-1), P=part, P4=lag(part,-4), trend)
#summary(reg<-lm(M~., data=u, na.action=NULL))  # could again use lm, but sarima also gives residual analysis
sarima(u[,1],0,0,0,xreg=u[,2:7])
## Coefficients:
##       intercept     M1      M2       T1       P      P4    trend
##         40.3838  0.315  0.2971  -0.1845  0.1113  0.0513  -0.5214
## s.e.     4.5982  0.037  0.0394   0.0309  0.0177  0.0195   0.0956
## sigma^2 estimated as 25.32:  log likelihood = -1529.55,  aic = 3075.09
##  
##-- full run using Kfilter --## 
trend = time(cmort) - mean(time(cmort))   # center time
const = time(cmort)/time(cmort)           # a ts of 1s
ded = ts.intersect(M=cmort, T1=lag(tempr,-1), P=part, P4=lag(part,-4), trend, const)
y = ded[,1]; input =ded[,2:6]
num = length(y); A = array(c(1,0), dim = c(1,2,num))
# Function to Calculate Likelihood
Linn=function(para){
  phi1=para[1]; phi2=para[2]; cR=para[3]
  b1=para[4]; b2=para[5]; b3=para[6]; b4=para[7]; alf=para[8]
  mu0 = matrix(c(0,0), 2, 1); Sigma0 = diag(100, 2)
  Phi = matrix(c(phi1, phi2, 1, 0), 2)
  Theta = matrix(c(phi1, phi2), 2)
  Ups = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
  Gam = matrix(c(0, 0, 0, b4, alf), 1, 5); cQ = cR;  S = cR^2
  kf = Kfilter2(num, y, A, mu0, Sigma0, Phi, Ups, Gam, Theta, cQ, cR, S, input)
  return(kf$like) }
# Estimation
init.par = c(phi1=.3, phi2=.3, cR=5, b1=-.2, b2=.1, b3=.05, b4=-1.6, alf=mean(cmort)) 
L = c(0,0,1,-1,0,0,-2,70); U = c(.5,.5,10,0,.5,.5,0,90) # used in optim
est = optim(init.par, Linn, NULL, method="L-BFGS-B", lower=L, upper=U, hessian=TRUE, control=list(trace=1,REPORT=1,factr=10^8))
SE = sqrt(diag(solve(est$hessian)))
# Results
u = cbind(estimate=est$par, SE)
rownames(u)=c("phi1","phi2","sigv","TL1","P","PL4","trend", 'constant')
round(u,3)
#          estimate    SE
# phi1        0.315 0.037
# phi2        0.318 0.041
# sigv        5.061 0.161
# TL1        -0.119 0.031
# P           0.119 0.018
# PL4         0.067 0.019
# trend      -1.340 0.220
# constant   88.752 7.015
#
##  residual analysis
phi1=est$par[1]; phi2=est$par[2]; cR=est$par[3]
b1=est$par[4]; b2=est$par[5]; b3=est$par[6]; b4=est$par[7]; alf=est$par[8]
mu0 = matrix(c(0,0), 2, 1); Sigma0 = diag(100, 2)
Phi = matrix(c(phi1, phi2, 1, 0), 2)
Theta = matrix(c(phi1, phi2), 2)
Ups = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
Gam = matrix(c(0, 0, 0, b4, alf), 1, 5); cQ = cR;  S = cR^2
kf = Kfilter2(num, y, A, mu0, Sigma0, Phi, Ups, Gam, Theta, cQ, cR, S, input)
res = ts(as.vector(kf$innov), start=start(cmort), freq=frequency(cmort))
sarima(res,0,0,0,no.constant=TRUE) # gives a full resid analysis


####################################################################
#######################################################################
######################## bootstrapping section #################
####################################################################
pdf(file="newbold.pdf", width=7.6, height=3.5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(qinfl, type='n', ylab='Rate (%)')
grid(lty=1)
lines(qinfl)
lines(qintr, lty=5, col=4, lwd=2)
legend("topleft", c("Inflation","Interest"),  lty=c(1,5), lwd=c(1,2), col=c(1,4),  bg='white')
dev.off()
  
#  bootstrap
tol = sqrt(.Machine$double.eps)  # determines convergence of optimizer     
nboot = 500                      # number of bootstrap replicates     
library(plyr)
 
y = window(qinfl, c(1953,1), c(1965,2))  # inflation   
z = window(qintr, c(1953,1), c(1965,2))  # interest   
num = length(y) 
A = array(z, dim=c(1,1,num))
input = matrix(1,num,1) 
    
# Function to Calculate Likelihood   
Linn=function(para){
  phi=para[1]; alpha=para[2]; b=para[3]; Ups=(1-phi)*b
  cQ=para[4]; cR=para[5]  
  kf=Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
  return(kf$like)    
  }

# Parameter Estimation   
mu0=1; Sigma0=.01; phi=.84; alpha=-.77; b=.85; cQ=.12; cR=1.1
init.par = c(phi,alpha,b,cQ,cR)    # initial parameters   
est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1,reltol=tol))  
SE = sqrt(diag(solve(est$hessian)))                     
phi = est$par[1]; alpha=est$par[2]; b=est$par[3]; Ups=(1-phi)*b         
cQ=est$par[4]; cR=est$par[5] 
cbind(estimate=est$par, SE)  

# BEGIN BOOTSTRAP   
# Likelihood for the bootstrapped data   
Linn2=function(para){
  phi=para[1]; alpha=para[2]; b=para[3]; Ups=(1-phi)*b
  cQ=para[4]; cR=para[5]  
  kf=Kfilter2(num,y.star,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
  return(kf$like)  }

# Run the filter at the estimates   
kf = Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)          

# Pull out necessary values from the filter and initialize  
xp=kf$xp; innov=kf$innov; sig=kf$sig; K=kf$K; e=innov/sqrt(sig)
e.star=e; y.star=y; xp.star=xp; k=4:50 
para.star = matrix(0, nboot, 5)  # to store estimates
init.par=c(.84,-.77,.85,.12,1.1)    

pr <- progress_text()            # displays progress
pr$init(nboot)

for (i in 1:nboot){
 # cat("iteration:", i, "\n")  # old counter
 pr$step()                     # new counter
 e.star[k] = sample(e[k], replace=TRUE)   
 for (j in k){
  xp.star[j] = phi*xp.star[j-1]+Ups+K[j]*sqrt(sig[j])*e.star[j] 
  }   
 y.star[k] = z[k]*xp.star[k]+alpha+sqrt(sig[k])*e.star[k]  
 est.star = optim(init.par, Linn2, NULL, method="BFGS", control=list(reltol=tol))     
 para.star[i,] = cbind(est.star$par[1], est.star$par[2], est.star$par[3], abs(est.star$par[4]), 
                      abs(est.star$par[5]))
} 

write(para.star, file = "bootstrap_results.txt")
                      
# Some summary statistics  
rmse = rep(NA,5)  # SEs from the bootstrap
for(i in 1:5){rmse[i]=sqrt(sum((para.star[,i]-est$par[i])^2)/nboot)
              cat(i, rmse[i],"\n") 
}
              
# phi and sigw  
phi = para.star[,1]; sigw = abs(para.star[,4]) 
phi = ifelse(phi<0, NA, phi)  # any phi< 0 not plotted

panel.hist <- function(x, ...){    # scatterplot with histograms
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

 pdf(file="bootexample1.pdf", width=7, height=5)
 u = cbind(phi, sigw)
 colnames(u) = c("f","s")
 pairs(u, cex=1.5, pch=1, diag.panel=panel.hist, cex.labels=1.5, font.labels=5, oma=c(0,0,0,0)+1.5)
 dev.off()


#######  scatterplot with histograms ########
scatterhist <- function(x, y, xlab='' , ylab='' , plottitle="", 
                        xsize=1, cleanup=TRUE, hist.col=8,...){
  # save the old graphics settings-- they may be needed
  def.par <- par(no.readonly = TRUE)
  
  zones <- matrix(c(1,1,1, 0,5,0, 2,6,4, 0,3,0), ncol = 3, byrow = TRUE)
  layout(zones, widths=c(0.3,4,1), heights = c(1,3,10,.75))
  
  # tuning to plot histograms nicely
  xhist <- hist(x, plot = FALSE)
  yhist <- hist(y, plot = FALSE)
  top <- max(c(xhist$counts, yhist$counts))

  # for all three titles: 
  #   drop the axis titles and omit boxes, set up margins
  par(xaxt="n", yaxt="n", bty="n",  mar = c(.3,2,.3,0) +.05)
  # fig 1 from the layout
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,paste(plottitle), cex=2)
  # fig 2
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,ylab, cex=1.5, srt=90)
  # fig 3
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,xlab, cex=1.5)
  # fig 4, the first histogram, needs different margins
  # no margin on the left
  par(mar = c(2,0,1,1))
  barplot(yhist$counts, axes = FALSE, xlim = c(0, top),
          space = 0, horiz = TRUE, col=hist.col)
  # fig 5, other histogram needs no margin on the bottom
  par(mar = c(0,2,1,1))
  barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0, col=hist.col)
  # fig 6, finally, the scatterplot-- needs regular axes, different margins
  par(mar = c(2,2,.5,.5), xaxt="s", yaxt="s", bty="n")
  # this color allows traparency & overplotting-- useful if a lot of points
  plot(x, y , pch=19, col=gray(.1,alpha=.25), cex=xsize, panel.first=grid(lty=1), ...)
  
  
  # reset the graphics, if desired 
  if(cleanup) {par(def.par)}
}
pdf(file="bootexample.pdf", width=7, height=6)
scatterhist(sigw, phi, ylab=expression(phi), xlab=expression(sigma[~w]), xsize=1.5, hist.col='lightblue')
dev.off()



#############################################################################
##############################  HMM Section ####################################
#
pdf(file="EQcount.pdf", width=7, height=4)
y = EQcount;
mean(y); var(y)
U = 2/sqrt(length(y))
lmax = 20
LAG = 1:lmax
ACF = acf(y, plot=FALSE)$acf[-1]
PACF = pacf(y, plot=FALSE)$acf
par(mar=c(2,2,0,0)+.5, mgp=c(1.4,.5,0))
layout(matrix(c(1,2, 1,3), nc=2))
plot(EQcount, type='n')
grid(lty=1); lines(EQcount, type='h')
plot(LAG, ACF, type='h', ylim=c(-.2,.6), panel.first=grid())
  abline(h=c(0,-U,U), lty=c(1,2,2), col=c(1,4,4))
plot(LAG, PACF, type='h', ylim=c(-.2,.6), panel.first=grid())
  abline(h=c(0,-U,U), lty=c(1,2,2), col=c(1,4,4))  
dev.off()
#####################
#### estimation #####
#####################
require(depmixS4)
model <- depmix(EQcount ~1, nstates=2, data=data.frame(EQcount), family=poisson())
set.seed(90210)
(fm <- fit(model))   # fm contains results from estimation 
summary(fm)   #not shown
#-- get parameters --#
# note- can't control the state names in depmix, so here we
#       maintain which is state 1 [min lam] and which is state 2 [max lam]
u = as.vector(getpars(fm)) 
 if (u[7] <= u[8]) { para.mle = c(u[3:6], exp(u[7]), exp(u[8])) 
  } else {  para.mle = c(u[6:3], exp(u[8]), exp(u[7])) 
 }
mtrans = matrix(para.mle[1:4], byrow=TRUE, nrow=2)
lams = para.mle[5:6]   
pi1 = mtrans[2,1]/(2 - mtrans[1,1] - mtrans[2,2])
pi2 = 1 - pi1

## graphics ##
pdf(file="EQcount_est.pdf", width=7.6, height=4.5)
layout(matrix(c(1,2,1,3), 2))
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.6,.6,0))
 # data and states
 plot(EQcount, type='n', ylim=c(4,42))
 grid(lty=1); lines(EQcount, type='h', col=gray(.7))
 #plot(EQcount, main="", ylab='EQcount', type='h', col='#808080')
 text(EQcount, col=6*posterior(fm)[,1]-2, labels=posterior(fm)[,1], cex=.9)
 # prob of state 2
 plot(ts(posterior(fm)[,3], start=1900), type='n', ylab=expression(hat(pi)[~2]*' (t | n)'))
 grid(lty=1); abline(h=.5, lty=2)
 lines(ts(posterior(fm)[,3], start=1900))
 # histogram
 hist(EQcount, breaks=30, prob=TRUE, main="")
 xvals = seq(1,45)
 u1 = pi1*dpois(xvals, lams[1])
 u2 = pi2*dpois(xvals, lams[2])
 lines(xvals, u1, col=4)
 lines(xvals, u2, col=10)
dev.off()
## end graphics

###  Bootstap  ####
# function to generate data
pois.HMM.generate_sample = function(n,m,lambda ,Mtrans,StatDist=NULL){
  # n = data length, m = # of states, Mtrans = transition matrix, StatDist = stationary distn 
  if(is.null(StatDist))StatDist =solve(t(diag(m)-Mtrans +1),rep(1,m))
  mvect = 1:m
  state = numeric(n)
  state [1] = sample(mvect ,1,prob=StatDist)
  for (i in 2:n)
       state[i]=sample(mvect ,1,prob=Mtrans[state[i-1] ,])
  y = rpois(n,lambda=lambda[state ])
  list(y= y, state= state)
}
# start it up
set.seed(10101101)
nboot = 100
nobs = length(EQcount)
para.star = matrix(NA, nrow=nboot, ncol = 6)
for (j in 1:nboot){
 x.star = pois.HMM.generate_sample(n=nobs, m=2, lambda=lams, Mtrans=mtrans)$y
 model <- depmix(x.star ~1, nstates=2, data=data.frame(x.star), family=poisson())
 u = as.vector(getpars(fit(model, verbose=0)))
 # make sure state 1 is the one with the smaller intensity parameter 
 if (u[7] <= u[8]) { para.star[j,] = c(u[3:6], exp(u[7]), exp(u[8])) }
     else  { para.star[j,] = c(u[6:3], exp(u[8]), exp(u[7])) }
}
# bootstrapped stnd errors 
SE = sqrt(apply(para.star,2,var)+ (apply(para.star,2,mean)-para.mle)^2)[c(1,4:6)]
names(SE)=c('seM11/M12', 'seM21/M22', 'seLam1', 'seLam2')
SE
##########  end poisson  #############

#######################################
#  start s&p500 normal #############
#####################################
libary(depmixS4)
library(astsa)   # sp500w is here
y = ts(sp500w, start=2003, freq=52)  # makes data useable for depmix
mod3 <- depmix(y~1, nstates=3, data=data.frame(y))
set.seed(2)
summary(fm3 <- fit(mod3))
# graphics  
pdf(file="sp500_est.pdf", width=7.6, height=5)
#layout(matrix(c(1,1,1,2,2,  1,1,1,3,3), 5,2))
layout(matrix(c(1,2, 1,3), 2), heights=c(1,.75))
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.6,.6,0))
plot(y, main="", ylab='S&P500 Weekly Returns', type='n', ylim=c(-.11,.11))
grid(lty=1); lines(y,  col=gray(.7))
 # change green to blue
 culer=4-posterior(fm3)[,1]; 
 culer[culer==3]=4
 text(y, col=culer, labels=4-posterior(fm3)[,1], cex=1.1)
# text(y, col=4-posterior(fm3)[,1], labels=4-posterior(fm3)[,1], cex=.9)  
# NOTE:  For display purposes, we reversed the names of states 1 and 3, 
#           the graphic looks better this way. 

# the MLEs: 
para.mle = as.vector(getpars(fm3)[-(1:3)])

# for display (states 1 and 3 names switched)
permu = matrix(c(0,0,1,0,1,0,1,0,0), 3,3) 
(mtrans.mle = permu%*%round(t(matrix(para.mle[1:9],3,3)),3)%*%permu)
(norms.mle =  round(matrix(para.mle[10:15],2,3),3)%*%permu)

acf(y^2, xlab='LAG', xlim=c(.02,.5), ylim=c(-.09,.5), panel.first=grid(lty=1) ) 


#  plot a histogram  
hist(y, 25, prob=TRUE, main="", xlab='     S&P500 Weekly Returns',   col=rgb(.68,.85,.90,alpha=.3),
 border='lightblue')
u = getpars(fm3)[13:18]  # 1 and 3 need to be reversed
up =  colSums(posterior(fm3)[-1,2:4])/length(y)
#### i'm going to increase up for 1 and 2 (really 3 and 2) for the histogram plot
up[1:2]=up[1:2]+.02
##
j=0;  culer=c(4,2,1)
 for (i in seq(1,5,by=2)) {mu = u[i] ; sig = u[i+1]
  x = seq(-.15,.12, by=.001); j=j+1; if(j==2){x=seq(-.1,.01, by=.001)}
  lines(x, up[j]*dnorm(x, mean=mu, sd=sig), col=culer[j])
  }
dev.off()

# bootstrap
set.seed(666)
n.obs = length(y)
n.boot = 100       # lower this if it takes too long
para.star = matrix(NA, nrow=n.boot, ncol = 15)
respst <- para.mle[10:15]
trst <- para.mle[1:9]
for ( nb in 1:n.boot){
  mod <- simulate(mod3)
  y.star = as.vector(mod@response[[1]][[1]]@y)
  dfy = data.frame(y.star)
  mod.star <- depmix(y.star~1, data=dfy, respst=respst, trst=trst, nst=3)
  fm.star = fit(mod.star, emcontrol=em.control(maxit = 1000, tol = 1e-5), verbose=FALSE)
  para.star[nb,] = as.vector(getpars(fm.star)[-(1:3)])
}
# bootstrap stnd errors 
SE = sqrt(apply(para.star,2,var)+ (apply(para.star,2,mean)-para.mle)^2)
# for display
(SE.mtrans.mle = permu%*%round(t(matrix(SE[1:9],3,3)),3)%*%permu)
(SE.norms.mle = round(matrix(SE[10:15], 2,3),3)%*%permu)




#####################################################
## switching AR -- flu data #############################
####################################################
# flu_MSwM
require(astsa)
require(MSwM)

set.seed(90210)  #  controls labels
# fit model on diffed flu   
dflu =  diff(flu)
model = lm(dflu~1)
mod = msmFit(model, k=2, p=2, sw=rep(TRUE,4))  # 2 regimes, AR(2)s 
# summary(mod)
plotProb(mod, which=3)  
##################  this is enough here

# smoother/filter probs for state 2
sp = mod@Fit@smoProb[,2]  
fp = mod@Fit@filtProb[,2]

# graphic
pdf(file="flu_MSwM.pdf", width=7.6, height=4)
x = window(diff(flu), start= c(1968,3)) 
regime = rep(1, length(x))

regime[sp <.5] = 2        # for graphic, get regime numbers 1 or 2    
layout(matrix(1:2,2,2), heights=c(2.2,1))
par(mar=c(2, 3,.2,.5), mgp=c(1.6,.6,0), cex=.9)
plot(x, type='n',  ylab= expression(nabla~flu[~t]), xlab='')
 abline(v=1968:1978,lty=1, col=gray(.9)); abline(h=seq(-.4,.4,.2), lty=1, col=gray(.9))
 lines(x, type='c', col='#808080')
 text(x, col=regime, labels=regime, cex=.9)  
 sp = ts(1-sp, start=tsp(x)[1], freq=12)
 fp = ts(1-fp, start=tsp(x)[1]+1/12, freq=12)
plot(sp, col=1, ylab= expression(hat(pi)[~2]*' (t | n)'), axes=FALSE, xlab='')
 axis(1); axis(2, c(0,.5,1)); box()
 abline(v=1968:1978,lty=1, col=gray(.9)); abline(h=c(.25,.5,.75), lty=1, col=gray(.9))
 lines(fp, type='h', col=4)
dev.off()

####################################






######################################################################
##########    switching section ##################################
#################################################################
# flusw ##########
y = as.matrix(flu)
num = length(y)
nstate = 4;

M1 = as.matrix(cbind(1,0,0,1))  # obs matrix normal
M2 = as.matrix(cbind(1,0,1,1))  # obs matrix flu epi

prob = matrix(0,num,1); yp = y  # to store pi2(t|t-1) & y(t|t-1)
xfilter = array(0, dim=c(nstate,1,num)) # to store x(t|t)

# Function to Calculate Likelihood 
Linn = function(para){
  alpha1=para[1]; alpha2=para[2]; beta0=para[3]      
  sQ1=para[4];  sQ2=para[5];  like=0
  xf=matrix(0, nstate, 1)  # x filter
  xp=matrix(0, nstate, 1)  # x predict
  Pf=diag(.1, nstate)      # filter covar
  Pp=diag(.1, nstate)      # predict covar
  pi11 <- .75 -> pi22;  pi12 <- .25 -> pi21; pif1 <- .5 -> pif2            
  phi=matrix(0,nstate,nstate)
  phi[1,1]=alpha1; phi[1,2]=alpha2; phi[2,1]=1; phi[4,4]=1 
  Ups = as.matrix(rbind(0,0,beta0,0))
  Q = matrix(0,nstate,nstate)
  Q[1,1]=sQ1^2; Q[3,3]=sQ2^2; R=0  # R=0 in final model
  # begin filtering 
    for(i in 1:num){
    xp = phi%*%xf + Ups; Pp = phi%*%Pf%*%t(phi) + Q
    sig1 = as.numeric(M1%*%Pp%*%t(M1) + R)
    sig2 = as.numeric(M2%*%Pp%*%t(M2) + R)
    k1 = Pp%*%t(M1)/sig1; k2 = Pp%*%t(M2)/sig2 
    e1 = y[i]-M1%*%xp; e2 = y[i]-M2%*%xp
    pip1 = pif1*pi11 + pif2*pi21; pip2 = pif1*pi12 + pif2*pi22;
    den1 = (1/sqrt(sig1))*exp(-.5*e1^2/sig1); 
    den2 = (1/sqrt(sig2))*exp(-.5*e2^2/sig2);
    denom = pip1*den1 + pip2*den2;
    pif1 = pip1*den1/denom; pif2 = pip2*den2/denom;
    pif1=as.numeric(pif1); pif2=as.numeric(pif2)
    e1=as.numeric(e1); e2=as.numeric(e2)
    xf = xp + pif1*k1*e1 + pif2*k2*e2
    eye = diag(1, nstate)
    Pf = pif1*(eye-k1%*%M1)%*%Pp + pif2*(eye-k2%*%M2)%*%Pp 
    like = like - log(pip1*den1 + pip2*den2)
    prob[i]<<-pip2; xfilter[,,i]<<-xf; innov.sig<<-c(sig1,sig2)
    yp[i]<<-ifelse(pip1 > pip2, M1%*%xp, M2%*%xp)  
    }    
 return(like)   
 }
 
# Estimation
alpha1=1.4; alpha2=-.5; beta0=.3; sQ1=.1; sQ2=.1
init.par = c(alpha1, alpha2, beta0, sQ1, sQ2)

(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))    
u = cbind(estimate=est$par, SE)
rownames(u)=c("alpha1","alpha2","beta0","sQ1","sQ2"); u

# Graphics 
predepi = ifelse(prob<.5,0,1); k = 6:length(y)       
Time = time(flu)[k]
regime = predepi[k]+1

pdf(file="flusw.pdf", width=7.6, height=6)
par(mfrow=c(3,1), mar=c(2,3,1,1)+.1)     
plot(Time, y[k], type="n", ylab="")
 grid(lty=1); lines(Time, y[k],  col=gray(.7))
  text(Time, y[k], col=regime, labels=regime, cex=1.1)  
 text(1979,.95,"(a)") 
plot(Time, xfilter[1,,k], type="n", ylim=c(-.1,.4), ylab="")
  grid(lty=1); lines(Time, xfilter[1,,k]) 
 lines(Time, xfilter[3,,k]); lines(Time, xfilter[4,,k])
 text(1979,.35,"(b)")
plot(Time, y[k], type="n",   ylim=c(.1,.9),ylab="")
 grid(lty=1); # points(Time, y[k], pch=19)
 prde1 = 2*sqrt(innov.sig[1]); prde2 = 2*sqrt(innov.sig[2])
 prde = ifelse(predepi[k]<.5, prde1,prde2)
   xx = c(Time, rev(Time))
   yy = c(yp[k]-prde, rev(yp[k]+prde))
 #lines(Time, yp[k]+prde, lty=2, lwd=1.5)
 #lines(Time, yp[k]-prde, lty=2, lwd=1.5)
 polygon(xx, yy, border=8, col=gray(.6, alpha=.3)) 
  points(Time, y[k], pch=19)
 text(1979,.85,"(c)")
dev.off()
###################################################



####################################################
##### stochastic volatility #######################
###################################################
# nyse_pred
require(astsa)
y = log(nyse^2)
num = length(y)
# Initial Parameters
phi0=0; phi1=.95; sQ=.2; alpha=mean(y); sR0=1; mu1=-3; sR1=2
init.par = c(phi0, phi1, sQ, alpha, sR0, mu1, sR1)
# Innovations Likelihood
Linn = function(para){
  phi0=para[1]; phi1=para[2]; sQ=para[3]; alpha=para[4]
  sR0=para[5]; mu1=para[6]; sR1=para[7]
  sv = SVfilter(num,y,phi0,phi1,sQ,alpha,sR0,mu1,sR1)
  return(sv$like)    }
# Estimation
(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
u = cbind(estimates=est$par, SE)
rownames(u)=c('phi0','phi1','sQ','alpha','sigv0','mu1','sigv1'); u
# Graphics   (need filters at the estimated parameters)
phi0=est$par[1]; phi1=est$par[2]; sQ=est$par[3]; alpha=est$par[4]
sR0=est$par[5]; mu1=est$par[6]; sR1=est$par[7]
sv = SVfilter(num,y,phi0,phi1,sQ,alpha,sR0,mu1,sR1)
# densities plot (f is chi-sq, fm is fitted mixture)
x = seq(-15,6,by=.01)
f = exp(-.5*(exp(x)-x))/(sqrt(2*pi))
f0 = exp(-.5*(x^2)/sR0^2)/(sR0*sqrt(2*pi))
f1 = exp(-.5*(x-mu1)^2/sR1^2)/(sR1*sqrt(2*pi))
fm = (f0+f1)/2
####### start plot ##
pdf(file="nyse_density.pdf", width=7.6, height=3)
par(mar=c(2.5,2.5,1,.5), mgp=c(1.6,.6,0))
plot(x, f, type='l', ylab='density', panel.first=grid(lty=1)); 
lines(x, fm, lty=2,lwd=2)
dev.off()


pdf(file="nyse_pred.pdf", width=7.6, height=3.25)
par(mar=c(2,2,1,.5)); 
strt = c(1984,22)
freak = 252 # average number of trading days per year
tsnyse = ts(nyse, start=strt, frequency=freak)
ts1 = window(tsnyse, start=1987, end=1988.5)
tsxp = ts(sv$xp, start=strt, frequency=freak )
ts2 = window(tsxp, start=1987, end=1988.5)
tsPp = ts(sv$Pp, start=strt, frequency=freak )
ts3 = window(tsPp, start=1987, end=1988.5)
plot(ts2,  main='', ylim=c(-.18,.12), ylab='', xlab='', type='n')
grid(lty=1)
lines(ts2/10, lwd=2, col=6)
lines(ts1, lwd=2, col=4)
 dev.off()
 
 ######################################################
 ## bootstrap gnp resids
 n.boot = 500   # number of bootstrap replicates
tol = sqrt(.Machine$double.eps)  # convergence tolerance
gnpgr = diff(log(gnp))
fit = arima(gnpgr, order=c(1,0,0))
y = as.matrix(log(resid(fit)^2))
num = length(y)
plot.ts(y, ylab='')
# Initial Parameters
phi1 = .9; sQ = .5; alpha = mean(y); sR0 = 1; mu1 = -3; sR1 = 2.5
init.par = c(phi1, sQ, alpha, sR0, mu1, sR1)
# Innovations Likelihood
Linn=function(para){
  phi1 = para[1]; sQ = para[2]; alpha = para[3]
  sR0 = para[4]; mu1 = para[5]; sR1 = para[6]
  sv = SVfilter(num, y, 0, phi1, sQ, alpha, sR0, mu1, sR1)
  return(sv$like)    }
# Estimation
(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
uu = rbind(estimates=est$par, SE)
colnames(uu)=c('phi1','sQ','alpha','sig0','mu1','sig1'); uu 
# Bootstrap
para.star = matrix(0, n.boot, 6)  # to store parameter estimates
Linn2 = function(para){   # calculate likelihood
  phi1 = para[1]; sQ = para[2]; alpha = para[3]
  sR0 = para[4]; mu1 = para[5]; sR1 = para[6]
  sv = SVfilter(num, y.star, 0, phi1, sQ, alpha, sR0, mu1, sR1)
  return(sv$like)    }
for (jb in 1:n.boot){
 cat('iteration:', jb, '\n')
 phi1 = est$par[1]; sQ = est$par[2]; alpha = est$par[3]
 sR0 = est$par[4]; mu1 = est$par[5]; sR1 = est$par[6]
 Q = sQ^2; R0 = sR0^2; R1 = sR1^2
 sv = SVfilter(num, y, 0, phi1, sQ, alpha, sR0, mu1, sR1)
 sig0 = sv$Pp+R0; sig1 = sv$Pp+R1;
 K0 = sv$Pp/sig0; K1 = sv$Pp/sig1
 inn0 = y-sv$xp-alpha; inn1 = y-sv$xp-mu1-alpha
 den1 = (1/sqrt(sig1))*exp(-.5*inn1^2/sig1)
 den0 = (1/sqrt(sig0))*exp(-.5*inn0^2/sig0)
 fpi1 = den1/(den0+den1)
 # start resampling at t=4
 e0 = inn0/sqrt(sig0); e1 = inn1/sqrt(sig1)
 indx = sample(4:num, replace=TRUE)
 sinn = cbind(c(e0[1:3], e0[indx]), c(e1[1:3], e1[indx]))
 eF = matrix(c(phi1, 1, 0, 0), 2, 2)
 xi = cbind(sv$xp,y) # initialize
   for (i in 4:num){    # generate boot sample
   G = matrix(c(0, alpha+fpi1[i]*mu1), 2, 1)
   h21 = (1-fpi1[i])*sqrt(sig0[i]); h11 = h21*K0[i]
   h22 = fpi1[i]*sqrt(sig1[i]); h12 = h22*K1[i]
   H = matrix(c(h11,h21,h12,h22),2,2)
   xi[i,] = t(eF%*%as.matrix(xi[i-1,],2) + G + H%*%as.matrix(sinn[i,],2))}
# Estimates from boot data
y.star = xi[,2]
phi1=.9; sQ=.5; alpha=mean(y.star); sR0=1; mu1=-3; sR1=2.5
init.par = c(phi1, sQ, alpha, sR0, mu1, sR1)   # same as for data
est.star = optim(init.par, Linn2, NULL, method='BFGS', control=list(reltol=tol))
para.star[jb,] = cbind(est.star$par[1], abs(est.star$par[2]), est.star$par[3], abs(est.star$par[4]), est.star$par[5], abs(est.star$par[6])) }
# Some summary statistics and graphics
# rmse = rep(NA,6)  # SEs from the bootstrap
# for(i in 1:6){
#   rmse[i] = sqrt(sum((para.star[,i]-est$par[i])^2)/n.boot)
#   cat(i, rmse[i],'\n') }
# dev.new(); phi = para.star[,1]
# hist(phi, 15, prob=TRUE, main='', xlim=c(.4,1.2), xlab='')
# u = seq(.4, 1.2, by=.01)
# lines(u,dnorm(u, mean=.8790267, sd=.1061884), lty='dashed', lwd=2)

pdf(file="gnp_sv_boot.pdf", width=7.7, height=3.25)  
layout(matrix(1:2,1,2), widths=c(2.2,1))
par(mar=c(2.5, 2,.2,.5), mgp=c(1.6,.6,0), cex=.9)
y = ts(y, start=start(gnp), freq=4)
plot.ts(y, ylab='', type='n')
grid(lty=1); lines(y)
phi = para.star[,1]
hist(phi, 15, prob=TRUE, main='', xlim=c(.4,1.2), ylab='', xlab='', border=gray(.5), col ='lightblue')
xx = seq(.4, 1.2, by=.01)
lines(xx, dnorm(xx, mean=uu[1,1], sd=uu[2,1]), lty='dashed', lwd=2)
dev.off()



t(round(rbind(u,rmse),3))
#      estimates    SE  rmse
#phi1      0.884 0.109 0.057
#sQ        0.381 0.221 0.324
#alpha    -9.654 0.343 1.529
#sig0      0.835 0.204 0.527
#mu1      -2.350 0.495 0.410
#sig1      2.453 0.293 0.375

 