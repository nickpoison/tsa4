## Chapter 6 - Old Kalman Filter and Smoother Details and Code



> __Note__ The old Kalman filter and smoother scripts and the EM scripts now are marked with an `x` ... e.g., `Kfilter1` is now `xKfilter1` etc. etc. etc. BUT, older scripts can be changed to the newer ones with only slight changes.  The gain in speed is worth the effort!!

> __Warning__ Eventually, the old scripts with an `x` prefix will be removed.
<br/><br/>

---

### &#128685; this is the OLD stuff and will NOT work in version 2.0 or later

The  three levels of code  `Kfilter0/Ksmooth0`,  `Kfilter1/Ksmooth1`, and `Kfilter2/Ksmooth2`
have been superseded by the newer `Kfilter` and `Ksmooth` scripts. The new scripts are faster, easier to work with, and they remove the need for 3 different scripts.  

This page contains a description of the old code and lists the older Chapter 6 code. What you see here does NOT apply anymore. 

### StaRt Old Stuff:


+ For various models, each script provides the Kalman filter/smoother, the innovations
 and the corresponding variance-covariance matrices, and the value of the innovations likelihood at the location of the parameter values passed to the script. MLE is then accomplished by calling  the script that runs the filter. _The model is specified by passing the model parameters._


+ Level 0 is for the case of a fixed measurement matrix and no inputs; i.e., if  A<sub>t</sub> = A for all t, and there are no inputs,   then use the code at level 0.  
    
+  If the measurement matrices are time varying or there are inputs, use the  code at a higher level (1 or 2).   Many of the examples in the text can be done at level 0.
    
    
+ Level 1 allows for time varying measurement matrices and inputs, and level 2 adds
    the possibility of correlated noise processes. 
    
The models for each case are (x is state, y is observation, and t = 1, &hellip;, n):
   
&diams; **Level 0:** &nbsp; &nbsp; x<sub>t</sub> = &Phi; x<sub>t-1</sub> + w<sub>t</sub>, &nbsp; &nbsp; y<sub>t</sub> = A x<sub>t</sub> + v<sub>t</sub>, &nbsp; &nbsp; w<sub>t</sub> ~ iid N<sub>p</sub>(0, Q) &perp;   v<sub>t</sub> ~ iid N<sub>q</sub>(0, R) &perp; x<sub>0</sub> ~ N<sub>p</sub>(&mu;<sub>0</sub>, &Sigma;<sub>0</sub>)

&diams; **Level 1:**  &nbsp; &nbsp;  x<sub>t</sub> = &Phi; x<sub>t-1</sub> +  &Upsilon; u<sub>t</sub> + w<sub>t</sub>,  &nbsp; &nbsp;  y<sub>t</sub> = A<sub>t</sub> x<sub>t</sub> +  &Gamma; u<sub>t</sub> + v<sub>t</sub>,  &nbsp;   u<sub>t</sub> are r-dimensional inputs, etc.<br>
     
&diams; **Level 2:** &nbsp; &nbsp; x<sub>t</sub> = &Phi; x<sub>t-1</sub> +  &Upsilon; u<sub>t</sub> + &Theta; w<sub>t-1</sub>, &nbsp;   &nbsp;   y<sub>t</sub> = A<sub>t</sub> x<sub>t</sub> +  &Gamma; u<sub>t</sub> + v<sub>t</sub>, &nbsp;  cov(w<sub>s</sub>, v<sub>t</sub>) = S &delta;<sub>s</sub><span style="position:relative; left: -.9ex; bottom: 2pt"><sup>t</sup></span>, &nbsp; &Theta; is p &times; m, and w<sub>t</sub> is m-dimensional, etc.

---  
    
### Level 0 - Fixed Measurement Matrices and No Inputs  
   
  The call to Kfilter0 is,   `Kfilter0(n,y,A,mu0,Sigma0,Phi,cQ,cR)` in fairly obvious notation   except that `cQ` and   `cR` are the Cholesky-type decompositions of 
  `Q` and `R`.  In particular `Q = t(cQ)%*%cQ` and `R = t(cR)%*%cR` is all that is required provided `Q`and `R` are valid covariance matrices (Q can be singular and there is an example in the text).   The call to `Ksmooth0` is similar. 

 __In all three cases, the smoother also returns the filter and the likelihood.__
  
---
 	
### Level 1 - Varying Measurement Matrices and Inputs
       
 The call to the filter is `Kfilter1(n,y,A,mu0,Sigma0,Phi,Ups,Gam,cQ,cR,input)` 
   where `A` is an array with `dim=c(q,p,n)`,  `Ups` is &Upsilon; [p &times; r],  `Gam`is &Gamma; [q &times; r],    and `input` is the matrix of inputs
     that has the same row dimension as y (which  is n &times; q), `input`is n &times; r; the state dimension is p).  The call to `Ksmooth1`is similar.  Set  `Ups`, `Gam`, or `input` to 0 (zero) if you  don't  use them. 
 
---  
  		
### Level 2 - Varying Measurement Matrices,  Inputs and Correlated Noise</h3>  
   
The call to the filter is `Kfilter2(n,y,A,mu0,Sigma0,Phi,Ups,Gam,Theta,cQ,cR,S,input)`, which is similar to `Kfilter1` but that `S` must be included.  `Kfilter2` runs the filter given in Property 6.5.  The call to `Ksmooth2` is similar. Set `Ups` or `Gam` or `input` to 0 (zero) if you don't  use them.   

---
---

## OLD Chapter 6 Code

Example 6.1
```r
tsplot(blood, type='o', col=c(6,4,2), lwd=2, pch=19, cex=1) 
```

Example 6.2
```r
tsplot(cbind(globtemp, globtempl), spag=TRUE, lwd=2, col=astsa.col(c(6,4),.5), ylab="Temperature Deviations")

# or the updated version (one is land only and the other ocean only)
tsplot(cbind(gtemp_land, gtemp_ocean), spaghetti=TRUE, lwd=2, pch=20, type="o", 
        col=astsa.col(c(4,2),.5), ylab="Temperature Deviations", main="Global Warming")
legend("topleft", legend=c("Land Surface", "Sea Surface"), lty=1, pch=20, col=c(4,2), bg="white")
```

Example 6.5
```r
# generate data 
set.seed(1)  
num = 50
w = rnorm(num+1,0,1)
v = rnorm(num,0,1)
                               
mu = cumsum(w)  # states:  mu[0], mu[1], . . ., mu[50] 
y = mu[-1] + v  # obs:  y[1], . . ., y[50]

# filter and smooth (Ksmooth0 does both)
mu0 = 0; sigma0 = 1;  phi = 1; cQ = 1; cR = 1   
ks = Ksmooth0(num, y, 1, mu0, sigma0, phi, cQ, cR)   

# pictures 
par(mfrow=c(3,1))
Time = 1:num

tsplot(Time, mu[-1], type='p', main="Prediction", ylim=c(-5,10))      
  lines(ks$xp)
  lines(ks$xp+2*sqrt(ks$Pp), lty="dashed", col="blue")
  lines(ks$xp-2*sqrt(ks$Pp), lty="dashed", col="blue")

tsplot(Time, mu[-1], type='p', main="Filter", ylim=c(-5,10))
  lines(ks$xf)
  lines(ks$xf+2*sqrt(ks$Pf), lty="dashed", col="blue")
  lines(ks$xf-2*sqrt(ks$Pf), lty="dashed", col="blue")

tsplot(Time, mu[-1], type='p',  main="Smoother", ylim=c(-5,10))
  lines(ks$xs)
  lines(ks$xs+2*sqrt(ks$Ps), lty="dashed", col="blue")
  lines(ks$xs-2*sqrt(ks$Ps), lty="dashed", col="blue") 

mu[1]; ks$x0n; sqrt(ks$P0n)   # initial value info

# In case you can't see the differences in the figures...
# ... either get new glasses or ... 
# ... plot them on the same graph (not shown in text)
dev.new()
tsplot(Time, mu[-1], type='o', pch=19, cex=1)
lines(ks$xp, col=4, lwd=3)
lines(ks$xf, col=3, lwd=3) 
lines(ks$xs, col=2, lwd=3)
names = c("predictor","filter","smoother")
legend("bottomright", names, col=4:2, lwd=3, lty=1, bg="white")
```


Example 6.6
```r
# Generate Data
set.seed(999)
num = 100
N = num+1
x = sarima.sim(n=N, ar=.8)
# below used in text
# x = arima.sim(n=N, list(ar = .8))
y = ts(x[-1] + rnorm(num,0,1))     

# Initial Estimates 
u = ts.intersect(y, lag(y,-1), lag(y,-2)) 
varu = var(u)
coru = cor(u) 
phi = coru[1,3]/coru[1,2]             
q = (1-phi^2)*varu[1,2]/phi   
r = varu[1,1] - q/(1-phi^2) 
(init.par = c(phi, sqrt(q), sqrt(r))) 

# Function to evaluate the likelihood 
Linn=function(para){
  phi = para[1]; sigw = para[2]; sigv = para[3]   
  Sigma0 = (sigw^2)/(1-phi^2); Sigma0[Sigma0<0]=0   
  kf = Kfilter0(num,y,1,mu0=0,Sigma0,phi,sigw,sigv)
  return(kf$like)   
  }

# Estimation  
(est = optim(init.par, Linn, gr=NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))      
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate=c(phi=est$par[1],sigw=est$par[2],sigv=est$par[3]), SE)
```



Example 6.7
```r
##- slight change from text, the data scaled -##
##- you can remove the scales to get to the original analysis -##
# Setup 
y = cbind(globtemp/sd(globtemp), globtempl/sd(globtempl))
num = nrow(y)
input = rep(1,num)
A = array(rep(1,2), dim=c(2,1,num))
mu0 = -.35; Sigma0 = 1;  Phi = 1

# Function to Calculate Likelihood 
Linn=function(para){
 cQ = para[1]      # sigma_w
  cR1 = para[2]    # 11 element of chol(R)
  cR2 = para[3]    # 22 element of chol(R)
  cR12 = para[4]   # 12 element of chol(R)
 cR = matrix(c(cR1,0,cR12,cR2),2)  # put the matrix together
 drift = para[5]
 kf = Kfilter1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)
 return(kf$like) 
 }

# Estimation
init.par = c(.1,.1,.1,0,.05)  # initial values of parameters
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1))) 
SE = sqrt(diag(solve(est$hessian))) 

# Summary of estimation  
estimate = est$par; u = cbind(estimate, SE)
rownames(u)=c("sigw","cR11", "cR22", "cR12", "drift"); u  

# Smooth (first set parameters to their final estimates)
cQ    = est$par[1]  
 cR1  = est$par[2]   
 cR2  = est$par[3]   
 cR12 = est$par[4]  
cR    = matrix(c(cR1,0,cR12,cR2), 2)
(R    = t(cR)%*%cR)    #  to view the estimated R matrix
drift = est$par[5]  
ks    = Ksmooth1(num,y,A,mu0,Sigma0,Phi,drift,0,cQ,cR,input)  

# Plot 
tsplot(y, spag=TRUE, margins=.5, type='o', pch=2:3, col=4:3, lty=6, ylab='Temperature Deviations')
xsm  = ts(as.vector(ks$xs), start=1880)
rmse = ts(sqrt(as.vector(ks$Ps)), start=1880)
lines(xsm, lwd=2)
  xx = c(time(xsm), rev(time(xsm)))
  yy = c(xsm-2*rmse, rev(xsm+2*rmse))
polygon(xx, yy, border=NA, col=gray(.6, alpha=.25))
```


Example 6.8
```r
library(nlme)   # loads package nlme (comes with R)

# Generate data (same as Example 6.6)
set.seed(999); num = 100; N = num+1
x = sarima.sim(ar=.8, n=N)
y = ts(x[-1] + rnorm(num,0,1))     

# Initial Estimates 
u = ts.intersect(y,lag(y,-1),lag(y,-2)) 
varu = var(u); coru = cor(u) 
phi = coru[1,3]/coru[1,2]             
q = (1-phi^2)*varu[1,2]/phi   
r = varu[1,1] - q/(1-phi^2) 
cr = sqrt(r); cq = sqrt(q); mu0 = 0; Sigma0 = 2.8
(em = EM0(num, y, 1, mu0, Sigma0, phi, cq, cr, 75, .00001))   

# Standard Errors  (this uses nlme)
phi = em$Phi; cq = chol(em$Q); cr = chol(em$R)
mu0 = em$mu0; Sigma0 = em$Sigma0
para = c(phi, cq, cr)

# Evaluate likelihood at estimates 
Linn=function(para){
  kf = Kfilter0(num, y, 1, mu0, Sigma0, para[1], para[2], para[3])
  return(kf$like) 
  }
emhess = fdHess(para, function(para) Linn(para))
SE = sqrt(diag(solve(emhess$Hessian)))  

# Display summary of estimation 
estimate = c(para, em$mu0, em$Sigma0); SE = c(SE,NA,NA)
u = cbind(estimate, SE)
rownames(u) = c("phi","sigw","sigv","mu0","Sigma0")
u 
```


Example 6.9
```r
y    = cbind(WBC, PLT, HCT)
num  = nrow(y)       
A    = array(0, dim=c(3,3,num))  # creates num 3x3 zero matrices
for(k in 1:num) if (y[k,1] > 0) A[,,k]= diag(1,3) 

# Initial values 
mu0    = matrix(0,3,1) 
Sigma0 = diag(c(.1,.1,1) ,3)
Phi    = diag(1,3)
cQ     = diag(c(.1,.1,1), 3)
cR     = diag(c(.1,.1,1), 3)  
(em = EM1(num, y, A, mu0, Sigma0, Phi, cQ, cR, 100, .001))    

# Graph smoother
ks  = Ksmooth1(num, y, A, em$mu0, em$Sigma0, em$Phi, 0, 0, chol(em$Q), chol(em$R), 0)
y1s = ks$xs[1,,] 
y2s = ks$xs[2,,] 
y3s = ks$xs[3,,]
p1  = 2*sqrt(ks$Ps[1,1,]) 
p2  = 2*sqrt(ks$Ps[2,2,]) 
p3  = 2*sqrt(ks$Ps[3,3,])

par(mfrow=c(3,1))
tsplot(WBC, type='p', pch=19, ylim=c(1,5), col=6, lwd=2, cex=1)
lines(y1s) 
  xx = c(time(WBC), rev(time(WBC)))  # same for all
  yy = c(y1s-p1, rev(y1s+p1))
polygon(xx, yy, border=8, col=astsa.col(8, alpha = .1))  

tsplot(PLT, type='p', ylim=c(3,6), pch=19, col=4, lwd=2, cex=1)
lines(y2s)
  yy = c(y2s-p2, rev(y2s+p2))
polygon(xx, yy, border=8, col=astsa.col(8, alpha = .1))  

tsplot(HCT, type='p', pch=19, ylim=c(20,40), col=2, lwd=2, cex=1)
lines(y3s)
  yy = c(y3s-p3, rev(y3s+p3))
polygon(xx, yy, border=8, col=astsa.col(8, alpha = .1))  
```



Example 6.10
```r
num = length(jj)
A = cbind(1,1,0,0) 

# Function to Calculate Likelihood 
Linn=function(para){
 Phi = diag(0,4) 
 Phi[1,1] = para[1] 
 Phi[2,]=c(0,-1,-1,-1); Phi[3,]=c(0,1,0,0); Phi[4,]=c(0,0,1,0)
 cQ1 = para[2]; cQ2 = para[3]     # sqrt q11 and sqrt q22
 cQ=diag(0,4); cQ[1,1]=cQ1; cQ[2,2]=cQ2
 cR = para[4]                     # sqrt r11
 kf = Kfilter0(num,jj,A,mu0,Sigma0,Phi,cQ,cR)
 return(kf$like)  
 }

# Initial Parameters 
mu0      = c(.7,0,0,0) 
Sigma0   = diag(.04,4)  
init.par = c(1.03, .1, .1, .5)  # Phi[1,1], the 2 Qs and R

# Estimation
est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1))
SE  = sqrt(diag(solve(est$hessian)))
u   = cbind(estimate=est$par,SE)
rownames(u)=c("Phi11","sigw1","sigw2","sigv"); u     

# Smooth
Phi      = diag(0,4) 
Phi[1,1] = est$par[1] 
Phi[2,]  = c(0,-1,-1,-1) 
Phi[3,]  = c(0,1,0,0) 
Phi[4,]  = c(0,0,1,0)
cQ1      = est$par[2]
cQ2      = est$par[3]      
cQ       = diag(0,4)
cQ[1,1]  = cQ1 
cQ[2,2]  = cQ2   
cR       = est$par[4]   
ks       = Ksmooth0(num, jj, A, mu0, Sigma0, Phi, cQ, cR)   

# Plots
Tsm   = ts(ks$xs[1,,], start=1960, freq=4)
Ssm   = ts(ks$xs[2,,], start=1960, freq=4)
p1    = 3*sqrt(ks$Ps[1,1,]); p2 = 3*sqrt(ks$Ps[2,2,])
par(mfrow=c(2,1))
tsplot(Tsm, main='Trend Component', ylab='Trend')
  xx  = c(time(jj), rev(time(jj)))
  yy  = c(Tsm-p1, rev(Tsm+p1))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .3))
tsplot(jj, main='Data & Trend+Season', ylab='J&J QE/Share', ylim=c(-.5,17))
  xx  = c(time(jj), rev(time(jj)) )
  yy  = c((Tsm+Ssm)-(p1+p2), rev((Tsm+Ssm)+(p1+p2)) )
polygon(xx, yy, border=NA, col=gray(.5, alpha = .3))

# Forecast
dev.new()
n.ahead = 12
y       = ts(append(jj, rep(0,n.ahead)), start=1960, freq=4)
rmspe   = rep(0,n.ahead) 
x00     = ks$xf[,,num]
P00     = ks$Pf[,,num]
Q       = t(cQ)%*%cQ 
R       = t(cR)%*%(cR)
for (m in 1:n.ahead){
       xp = Phi%*%x00
       Pp = Phi%*%P00%*%t(Phi)+Q
      sig = A%*%Pp%*%t(A)+R
        K = Pp%*%t(A)%*%(1/sig)
      x00 = xp 
      P00 = Pp-K%*%A%*%Pp
 y[num+m] = A%*%xp
 rmspe[m] = sqrt(sig) 
}
tsplot(y, type='o', main='', ylab='J&J QE/Share', ylim=c(5,30), xlim = c(1975,1984))
upp  = ts(y[(num+1):(num+n.ahead)]+2*rmspe, start=1981, freq=4)
low  = ts(y[(num+1):(num+n.ahead)]-2*rmspe, start=1981, freq=4)
 xx  = c(time(low), rev(time(upp)))
 yy  = c(low, rev(upp))
polygon(xx, yy, border=8, col=gray(.5, alpha = .3))
abline(v=1981, lty=3)
```



Example 6.12
```r
# Preliminary analysis
fit1   = sarima(cmort, 2,0,0, xreg=time(cmort))
acf(cbind(dmort <- resid(fit1$fit), tempr, part))
lag2.plot(tempr, dmort, 8)
lag2.plot(part, dmort, 8)

# quick and dirty fit (detrend then fit ARMAX)
trend   = time(cmort) - mean(time(cmort))  
dcmort  = resid(fit2 <- lm(cmort~trend, na.action=NULL))  # detrended mort
u       = ts.intersect(dM=dcmort, dM1=lag(dcmort,-1), dM2=lag(dcmort,-2), T1=lag(tempr,-1), P=part, P4=lag(part,-4))
sarima(u[,1], 0,0,0, xreg=u[,2:6])  # ARMAX fit with residual analysis 

# all estimates at once
trend   = time(cmort) - mean(time(cmort)) # center time
const   = time(cmort)/time(cmort)         # appropriate time series of 1s
ded     = ts.intersect(M=cmort, T1=lag(tempr,-1), P=part, P4=lag(part,-4), trend, const)
y       = ded[,1]
input   = ded[,2:6]
num     = length(y)
A       = array(c(1,0), dim = c(1,2,num))

# Function to Calculate Likelihood
Linn=function(para){
 phi1   = para[1]; phi2 = para[2]; cR = para[3];  b1 = para[4]
 b2     = para[5];   b3 = para[6]; b4 = para[7]; alf = para[8]
 mu0    = matrix(c(0,0), 2, 1)
 Sigma0 = diag(100, 2)
 Phi    = matrix(c(phi1, phi2, 1, 0), 2)
 Theta  = matrix(c(phi1, phi2), 2)
 Ups    = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
 Gam    = matrix(c(0, 0, 0, b4, alf), 1, 5); cQ = cR; S = cR^2
 kf     = Kfilter2(num, y, A, mu0, Sigma0, Phi, Ups, Gam, Theta, cQ, cR, S, input)
return(kf$like) 
}

# Estimation - prelim analysis gives good starting values
init.par = c(phi1=.3, phi2=.3, cR=5, b1=-.2, b2=.1, b3=.05, b4=-1.6, alf=mean(cmort)) 
L = c( 0,  0,  1, -1,  0,  0, -2, 70)   # lower bound on parameters
U = c(.5, .5, 10,  0, .5, .5,  0, 90)   # upper bound - used in optim
est      = optim(init.par, Linn, NULL, method='L-BFGS-B', lower=L, upper=U, 
                 hessian=TRUE, control=list(trace=1, REPORT=1, factr=10^8))
SE       = sqrt(diag(solve(est$hessian)))
round(cbind(estimate=est$par, SE), 3) # results

# Residual Analysis (not shown)
phi1   = est$par[1]; phi2 = est$par[2]
cR     = est$par[3]; b1   = est$par[4]
b2     = est$par[5]; b3   = est$par[6]
b4     = est$par[7]; alf  = est$par[8]
mu0    = matrix(c(0,0), 2, 1)
Sigma0 = diag(100, 2)
Phi    = matrix(c(phi1, phi2, 1, 0), 2)
Theta  = matrix(c(phi1, phi2), 2)
Ups    = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
Gam    = matrix(c(0, 0, 0, b4, alf), 1, 5)
cQ     = cR
S      = cR^2
kf     = Kfilter2(num, y, A, mu0, Sigma0, Phi, Ups, Gam, Theta, cQ, cR, S, input)
res    = ts(as.vector(kf$innov), start=start(cmort), freq=frequency(cmort))
sarima(res, 0,0,0, no.constant=TRUE)  # gives a full residual analysis

# Similar fit with but with trend in the X of ARMAX
trend  = time(cmort) - mean(time(cmort))
u      = ts.intersect(M=cmort, M1=lag(cmort,-1), M2=lag(cmort,-2), T1=lag(tempr,-1), 
           P=part, P4=lag(part -4), trend)
sarima(u[,1], 0,0,0, xreg=u[,2:7])
```



Example 6.13
```r
################################## 
# NOTE:  If this takes a long time to run on your machine, try
#         tol   = .0001   and if you need more speed
#         nboot = 250         
tol = sqrt(.Machine$double.eps)  # determines convergence of optimizer     
nboot = 500                      # number of bootstrap replicates     
################################## 

pb = txtProgressBar(min = 0, max = nboot, initial = 0, style=3)  # progress bar

y     = window(qinfl, c(1953,1), c(1965,2))  # inflation   
z     = window(qintr, c(1953,1), c(1965,2))  # interest   
num   = length(y) 
A     = array(z, dim=c(1,1,num))
input = matrix(1,num,1)  

# Function to Calculate Likelihood   
Linn  = function(para, y.data){  # pass data also
   phi = para[1];  alpha = para[2]
   b   = para[3];  Ups   = (1-phi)*b
   cQ  = para[4];  cR    = para[5]  
   kf  = Kfilter2(num,y.data,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)
   return(kf$like)    
}

# Parameter Estimation   
mu0      = 1
Sigma0   = .01  
init.par = c(phi=.84, alpha=-.77, b=.85, cQ=.12, cR=1.1) # initial values   

est = optim(init.par,  Linn, NULL, y.data=y, method="BFGS", hessian=TRUE, 
             control=list(trace=1, REPORT=1, reltol=tol))  
SE  = sqrt(diag(solve(est$hessian)))   

phi   = est$par[1];  alpha = est$par[2]
b     = est$par[3];  Ups   = (1-phi)*b         
cQ    = est$par[4];  cR    = est$par[5] 
round(cbind(estimate=est$par, SE), 3)  


# BEGIN BOOTSTRAP   
# Run the filter at the estimates 
kf = Kfilter2(num,y,A,mu0,Sigma0,phi,Ups,alpha,1,cQ,cR,0,input)  

# Pull out necessary values from the filter and initialize  
xp      = kf$xp
innov   = kf$innov 
sig     = kf$sig 
e       = innov/sqrt(sig)
e.star  = e                      # initialize values
y.star  = y  
xp.star = xp  
k       = 4:50                   # hold first 3 observations fixed 
para.star = matrix(0, nboot, 5)  # to store estimates
init.par  =  c(.84, -.77, .85, .12, 1.1)    

for (i in 1:nboot){
 setTxtProgressBar(pb,i)                       
 e.star[k] = sample(e[k], replace=TRUE)   
 for (j in k){ 
    K        = (phi*Pp[j]*z[j])/sig[j]
  xp.star[j] = phi*xp.star[j-1] + Ups+K[j]*sqrt(sig[j])*e.star[j] }   
   y.star[k] = z[k]*xp.star[k] + alpha + sqrt(sig[k])*e.star[k]  
 est.star = optim(init.par, Linn, NULL, y.data=y.star, method="BFGS", control=list(reltol=tol))     
 para.star[i,] = cbind(est.star$par[1], est.star$par[2], est.star$par[3], 
                       abs(est.star$par[4]), abs(est.star$par[5]))   
}
close(pb) 

# Some summary statistics  
rmse = rep(NA,5)                 # SEs from the bootstrap
for(i in 1:5){rmse[i]=sqrt(sum((para.star[,i]-est$par[i])^2)/nboot)
              cat(i, rmse[i],"\n") 
             }              
# Plot phi and sigw  (scatter.hist in astsa v1.13)
phi  = para.star[,1] 
sigw = abs(para.star[,4]) 
phi  = ifelse(phi<0, NA, phi)    # any phi < 0 not plotted
scatter.hist(sigw, phi, ylab=expression(phi), xlab=expression(sigma[~w]), 
             hist.col=astsa.col(5,.4), pt.col=5, pt.size=1.5)
```

Example 6.14
```r
set.seed(123)
num   = 50
w     = rnorm(num,0,.1)
x     = cumsum(cumsum(w))
y     = x + rnorm(num,0,1)
## State Space ##
Phi   = matrix(c(2,1,-1,0),2)
A     = matrix(c(1,0),1)
mu0   = matrix(0,2); Sigma0 = diag(1,2)
Linn  = function(para){
  sigw = para[1] 
  sigv = para[2]
  cQ   = diag(c(sigw,0))
  kf   = Kfilter0(num, y, A, mu0, Sigma0, Phi, cQ, sigv)
return(kf$like) 
}
## Estimation ##
init.par = c(.1, 1)
(est  = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))
SE    = sqrt(diag(solve(est$hessian)))
# Summary of estimation
estimate = est$par; u = cbind(estimate, SE)
rownames(u) = c("sigw","sigv"); u
# Smooth
sigw  = est$par[1]
cQ    = diag(c(sigw,0))
sigv  = est$par[2]
ks    = Ksmooth0(num, y, A, mu0, Sigma0, Phi, cQ, sigv)
xsmoo = ts(ks$xs[1,1,]); psmoo = ts(ks$Ps[1,1,])
upp   = xsmoo+2*sqrt(psmoo)
low   = xsmoo-2*sqrt(psmoo)
#
tsplot(x, ylab="", ylim=c(-1,8), col=1)
lines(y, type='o', col=8)
lines(xsmoo, col=4, lty=2, lwd=3)
lines(upp, col=4, lty=2); lines(low, col=4, lty=2)
lines(smooth.spline(y), lty=1, col=2)
legend("topleft", c("Observations","State"), pch=c(1,-1), lty=1, lwd=c(1,2), col=c(8,1))
legend("bottomright", c("Smoother", "GCV Spline"), lty=c(2,1), lwd=c(3,1), col=c(4,2))
```

Example 6.16
```r
library(depmixS4)
model <- depmix(EQcount ~1, nstates=2, data=data.frame(EQcount), family=poisson('identity'), respstart=c(15,25))
set.seed(90210)
summary(fm <- fit(model))   # estimation results  
standardError(fm)           # with standard errors

##-- A little nicer display of the parameters --##
para.mle = as.vector(getpars(fm))[3:8]  
( mtrans = matrix(para.mle[1:4], byrow=TRUE, nrow=2) )
( lams   = para.mle[5:6] )  
( pi1    = mtrans[2,1]/(2 - mtrans[1,1] - mtrans[2,2]) )  
( pi2    = 1 - pi1 )

#-- Graphics --##
par(mfrow=c(3,1))
# data and states
tsplot(EQcount, main="", ylab='EQcount', type='h', col=gray(.7), ylim=c(0,50))
text(EQcount, col=6*posterior(fm)[,1]-2, labels=posterior(fm)[,1])
# prob of state 2
tsplot(ts(posterior(fm)[,3], start=1900), ylab = expression(hat(pi)[~2]*'(t|n)'));  abline(h=.5, lty=2)
# histogram
hist(EQcount, breaks=30, prob=TRUE, main="")
xvals = seq(1,45)
u1 = pi1*dpois(xvals, lams[1])  
u2 = pi2*dpois(xvals, lams[2])
lines(xvals, u1, col=4)   
lines(xvals, u2, col=2)
```

Example 6.17
```r
library(depmixS4)
y = ts(sp500w, start=2003, freq=52)  # make data depmix friendly
mod3 <- depmix(y~1, nstates=3, data=data.frame(y))
set.seed(2)
summary(fm3 <- fit(mod3))   # estimation results  

##-- a little nicer display --## 
 para.mle    = as.vector(getpars(fm3)[-(1:3)])
 permu       = matrix(c(0,0,1,0,1,0,1,0,0), 3,3)   # for the label switch 
 (mtrans.mle = permu%*%round(t(matrix(para.mle[1:9],3,3)),3)%*%permu)
 (norms.mle  = round(matrix(para.mle[10:15],2,3),3)%*%permu)

##-- Graphics --##  
layout(matrix(c(1,2, 1,3), 2), heights=c(1,.75))

tsplot(y, main="", ylab='S&P500 Weekly Returns', col=gray(.7), ylim=c(-.11,.11))
 culer = 4-posterior(fm3)[,1];  culer[culer==3]=4  # switch labels 1 and 3
 text(y, col=culer, labels=4-posterior(fm3)[,1])

acf1(y^2, 25) 

hist(y, 25, prob=TRUE, main='', col=astsa.col(8,.2))
pi.hat = colSums(posterior(fm3)[-1,2:4])/length(y)
culer = c(1,2,4)
for (i in 1:3) { 
 mu = norms.mle[1,i]; sig = norms.mle[2,i]
 x = seq(-.2,.15, by=.001)
lines(x, pi.hat[4-i]*dnorm(x, mean=mu, sd=sig), col=culer[i], lwd=2)  
 }
```

Example 6.18
```r
library(MSwM)
set.seed(90210)
dflu  = diff(flu)
model = lm(dflu~ 1)
mod   = msmFit(model, k=2, p=2, sw=rep(TRUE,4)) # 2 regimes, AR(2)s
summary(mod)
plotProb(mod, which=3)
```



Example 6.22
```r
y   = flu  
num = length(y)
nstate = 4                      # state dimenstion
M1 = as.matrix(cbind(1,0,0,1))  # obs matrix normal
M2 = as.matrix(cbind(1,0,1,1))  # obs matrix flu epi
prob = matrix(0,num,1); yp = y  # to store pi2(t|t-1) & y(t|t-1)
xfilter = array(0, dim=c(nstate,1,num)) # to store x(t|t)
# Function to Calculate Likelihood
Linn = function(para){
  alpha1 = para[1]; alpha2 = para[2]; beta0 = para[3]
  sQ1 = para[4];  sQ2 = para[5];  like=0
  xf  = matrix(0, nstate, 1)  # x filter
  xp  = matrix(0, nstate, 1)  # x pred
  Pf  = diag(.1, nstate)      # filter cov
  Pp  = diag(.1, nstate)      # pred cov
  pi11 <- .75 -> pi22;  pi12 <- .25 -> pi21; pif1 <- .5 -> pif2
  phi = matrix(0,nstate,nstate)
  phi[1,1] = alpha1; phi[1,2] = alpha2; phi[2,1]=1; phi[4,4]=1
  Ups = as.matrix(rbind(0,0,beta0,0))
  Q   = matrix(0,nstate,nstate)
  Q[1,1] = sQ1^2; Q[3,3] = sQ2^2; R=0  # R=0 in final model
  # begin filtering #
    for(i in 1:num){
    xp   = phi%*%xf + Ups; Pp = phi%*%Pf%*%t(phi) + Q
    sig1 = as.numeric(M1%*%Pp%*%t(M1) + R)
    sig2 = as.numeric(M2%*%Pp%*%t(M2) + R)
    k1   = Pp%*%t(M1)/sig1; k2 = Pp%*%t(M2)/sig2
    e1   = y[i]-M1%*%xp; e2 = y[i]-M2%*%xp
    pip1 = pif1*pi11 + pif2*pi21; pip2 = pif1*pi12 + pif2*pi22
    den1 = (1/sqrt(sig1))*exp(-.5*e1^2/sig1)
    den2 = (1/sqrt(sig2))*exp(-.5*e2^2/sig2)
    denm = pip1*den1 + pip2*den2
    pif1 = pip1*den1/denm; pif2 = pip2*den2/denm
    pif1 = as.numeric(pif1); pif2 = as.numeric(pif2)
    e1   = as.numeric(e1); e2=as.numeric(e2)
    xf   = xp + pif1*k1*e1 + pif2*k2*e2
    eye  = diag(1, nstate)
    Pf   = pif1*(eye-k1%*%M1)%*%Pp + pif2*(eye-k2%*%M2)%*%Pp
    like = like - log(pip1*den1 + pip2*den2)
    prob[i]<<-pip2; xfilter[,,i]<<-xf; innov.sig<<-c(sig1,sig2)
    yp[i]<<-ifelse(pip1 > pip2, M1%*%xp, M2%*%xp)  
    }
return(like)   
}
# Estimation
alpha1 = 1.4; alpha2 = -.5; beta0 = .3; sQ1 = .1; sQ2 = .1
init.par = c(alpha1, alpha2, beta0, sQ1, sQ2)
(est = optim(init.par, Linn, NULL, method='BFGS', hessian=TRUE, control=list(trace=1,REPORT=1)))
SE   = sqrt(diag(solve(est$hessian)))
u    = cbind(estimate=est$par, SE)
rownames(u)=c('alpha1','alpha2','beta0','sQ1','sQ2'); u

# Graphics
predepi = ifelse(prob<.5,0,1)
FLU = window(flu, start=1968.4)
Time = window(time(flu), start=1968.4)  
k = 6:num     
par(mfrow=c(3,1))     
tsplot(FLU, col=8, ylab='flu')
 text(FLU, col= predepi[k]+1, labels=predepi[k]+1, cex=1.1) 
 legend('topright', '(a)', bty='n')

filters = ts(t(xfilter[c(1,3,4),,]), start=tsp(flu)[1], frequency=tsp(flu)[3])
tsplot(window(filters, start=1968.4),  spag=TRUE, col=2:4, ylab='filter')
 legend('topright', '(b)', bty='n')

tsplot(FLU, type='p', pch=19, ylab='flu', cex=1.2)
 prde1 = 2*sqrt(innov.sig[1]); prde2 = 2*sqrt(innov.sig[2])
 prde = ifelse(predepi[k]<.5, prde1, prde2)
   xx = c(Time, rev(Time))
   yy = c(yp[k]-prde, rev(yp[k]+prde))
 polygon(xx, yy, border=8, col=gray(.6, alpha=.3)) 
 legend('topright', '(c)', bty='n')
```




Example 6.23
```r
y   = log(nyse^2) 
num = length(y)

# Initial Parameters
phi0=0; phi1=.95; sQ=.2; alpha=mean(y); sR0=1; mu1=-3; sR1=2
init.par = c(phi0,phi1,sQ,alpha,sR0,mu1,sR1)

# Innovations Likelihood 
Linn = function(para){
  phi0=para[1]; phi1=para[2]; sQ=para[3]; alpha=para[4]
  sR0=para[5]; mu1=para[6]; sR1=para[7]
  sv = SVfilter(num,y,phi0,phi1,sQ,alpha,sR0,mu1,sR1)
  return(sv$like)    
}

# Estimation  
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
u = cbind(estimates=est$par, SE)  
rownames(u)=c("phi0","phi1","sQ","alpha","sigv0","mu1","sigv1"); u

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
tsplot(x, f, xlab='x')
lines(x, fm, lty=2, lwd=2)
legend('topleft', legend=c('log chi-square', 'normal mixture'), lty=1:2)

dev.new()
Time = 701:1100
tsplot(Time, nyse[Time], type='l', col=4, lwd=2, ylab='', xlab='', ylim=c(-.18,.12))
lines(Time, sv$xp[Time]/10, lwd=2, col=6)
```



Example 6.24
```r
n.boot = 500        # number of bootstrap replicates
tol = sqrt(.Machine$double.eps)  # convergence limit 

gnpgr = diff(log(gnp))
fit = arima(gnpgr, order=c(1,0,0))
y = as.matrix(log(resid(fit)^2))
num = length(y) 
tsplot(y, ylab="")

# Initial Parameters
phi1 = .9; sQ = .5; alpha = mean(y); sR0 = 1; mu1 = -3; sR1 = 2.5
init.par = c(phi1, sQ, alpha, sR0, mu1, sR1)

# Innovations Likelihood 
Linn=function(para){
  phi1 = para[1]; sQ = para[2]; alpha = para[3]
  sR0 = para[4]; mu1 = para[5]; sR1 = para[6]
  sv = SVfilter(num, y, 0, phi1, sQ, alpha, sR0, mu1, sR1)
  return(sv$like)    
  }

# Estimation  
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1)))
SE = sqrt(diag(solve(est$hessian)))
u = cbind(estimates=est$par, SE)  
rownames(u)=c("phi1","sQ","alpha","sig0","mu1","sig1"); u

# Bootstrap 
para.star = matrix(0, n.boot, 6)  # to store parameter estimates
Linn2 = function(para){   
  phi1 = para[1]; sQ = para[2]; alpha = para[3]
  sR0 = para[4]; mu1 = para[5]; sR1 = para[6]
  sv = SVfilter(num, y.star, 0, phi1, sQ, alpha, sR0, mu1, sR1)
  return(sv$like)    
  }

for (jb in 1:n.boot){  cat("iteration:", jb, "\n")
 phi1 = est$par[1]; sQ = est$par[2]; alpha = est$par[3]    
 sR0 = est$par[4]; mu1 = est$par[5]; sR1 = est$par[6] 
 Q = sQ^2; R0 = sR0^2; R1 = sR1^2
 sv = SVfilter(num, y, 0, phi1, sQ, alpha, sR0, mu1, sR1)
 
 sig0 = sv$Pp+R0 
 sig1 = sv$Pp+R1 
 K0   = sv$Pp/sig0 
 K1   = sv$Pp/sig1
 inn0 = y-sv$xp-alpha; inn1 = y-sv$xp-mu1-alpha
 den1 = (1/sqrt(sig1))*exp(-.5*inn1^2/sig1)
 den0 = (1/sqrt(sig0))*exp(-.5*inn0^2/sig0)
 fpi1 = den1/(den0+den1)
 
 # (start resampling at t=4)
 e0   = inn0/sqrt(sig0)
 e1   = inn1/sqrt(sig1)  
 indx = sample(4:num, replace=TRUE)
 sinn = cbind(c(e0[1:3], e0[indx]), c(e1[1:3], e1[indx]))
 eF   = matrix(c(phi1, 1, 0, 0), 2, 2)
 xi   = cbind(sv$xp,y)    # initialize
   
   for (i in 4:num){    # generate boot sample
   G   = matrix(c(0, alpha+fpi1[i]*mu1), 2, 1)
   h21 = (1-fpi1[i])*sqrt(sig0[i]); h11 = h21*K0[i]
   h22 = fpi1[i]*sqrt(sig1[i]); h12 = h22*K1[i]
   H   = matrix(c(h11,h21,h12,h22),2,2)
   xi[i,] = t(eF%*%as.matrix(xi[i-1,],2) + G + H%*%as.matrix(sinn[i,],2))
   }
   
 # Estimates from boot data 
 y.star = xi[,2]
 phi1 = .9; sQ = .5; alpha = mean(y.star); sR0 = 1; mu1 = -3; sR1 = 2.5
 init.par = c(phi1, sQ, alpha, sR0, mu1, sR1)   #  same as for data
 est.star = optim(init.par, Linn2, NULL, method="BFGS", control=list(reltol=tol)) 
 para.star[jb,] = cbind(est.star$par[1], abs(est.star$par[2]), est.star$par[3], abs(est.star$par[4]), 
                        est.star$par[5], abs(est.star$par[6])) 
}

# Some summary statistics and graphics  
rmse = rep(NA, 6)  # SEs from the bootstrap
for(i in 1:6){rmse[i] = sqrt(sum((para.star[,i]-est$par[i])^2)/n.boot)
                cat(i, rmse[i],"\n") 
             }  
dev.new()
phi = para.star[,1]
hist(phi, 15, prob=TRUE, main="", xlim=c(0,2), xlab="", col=astsa.col(4,.3))
abline(v=mean(phi), col=4)
curve(dnorm(x, mean=u[1,1], sd=u[2,1]), 0, 2, add=TRUE)
abline(v=u[1,1])
```


Example 6.26 

&emsp; Adapted from code by: [Hedibert Freitas Lopes](http://hedibert.org/)

```r
##-- Notation --##
#           y(t) = x(t) + v(t);    v(t) ~ iid N(0,V)                     
#           x(t) = x(t-1) + w(t);  w(t) ~ iid N(0,W)                        
#  priors:  x(0) ~ N(m0,C0);  V ~ IG(a,b);  W ~ IG(c,d)
#    FFBS:  x(t|t) ~ N(m,C);  x(t|n) ~ N(mm,CC);  x(t|t+1) ~ N(a,R)  
##-- 
ffbs = function(y,V,W,m0,C0){
  n  = length(y);  a  = rep(0,n);  R  = rep(0,n)
  m  = rep(0,n);   C  = rep(0,n);  B  = rep(0,n-1)     
  H  = rep(0,n-1); mm = rep(0,n);  CC = rep(0,n)
  x  = rep(0,n); llike = 0.0
  for (t in 1:n){
    if(t==1){a[1] = m0; R[1] = C0 + W
      }else{ a[t] = m[t-1]; R[t] = C[t-1] + W }
    f      = a[t]
    Q      = R[t] + V
    A      = R[t]/Q
    m[t]   = a[t]+A*(y[t]-f)
    C[t]   = R[t]-Q*A**2
    B[t-1] = C[t-1]/R[t]
    H[t-1] = C[t-1]-R[t]*B[t-1]**2
    llike  = llike + dnorm(y[t],f,sqrt(Q),log=TRUE) }
  mm[n] = m[n]; CC[n] = C[n]
  x[n]  = rnorm(1,m[n],sqrt(C[n]))
  for (t in (n-1):1){
    mm[t] = m[t] + C[t]/R[t+1]*(mm[t+1]-a[t+1])
    CC[t] = C[t] - (C[t]^2)/(R[t+1]^2)*(R[t+1]-CC[t+1])
    x[t]  = rnorm(1,m[t]+B[t]*(x[t+1]-a[t+1]),sqrt(H[t]))  }
return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))   }

# Simulate states and data
set.seed(1); W = 0.5; V = 1.0
n  = 100; m0 = 0.0; C0 = 10.0; x0 = 0
w  = rnorm(n,0,sqrt(W))
v  = rnorm(n,0,sqrt(V))
x  = y = rep(0,n)
x[1] = x0   + w[1]
y[1] = x[1] + v[1]
for (t in 2:n){
  x[t] = x[t-1] + w[t]
  y[t] = x[t] + v[t]   }

# actual smoother (for plotting)
ks = Ksmooth0(num=n, y, A=1, m0, C0, Phi=1, cQ=sqrt(W), cR=sqrt(V))
xsmooth = as.vector(ks$xs)

# run it
run = ffbs(y,V,W,m0,C0)
m   = run$m; C = run$C; mm = run$mm
CC  = run$CC; L1 = m-2*C; U1  = m+2*C
L2  = mm-2*CC; U2 = mm+2*CC
N   = 50
Vs  = seq(0.1,2,length=N)
Ws  = seq(0.1,2,length=N)
likes = matrix(0,N,N)
for (i in 1:N){
 for (j in 1:N){
   V   = Vs[i]
   W   = Ws[j]
   run = ffbs(y,V,W,m0,C0)    
  likes[i,j] = run$llike  }  }  
# Hyperparameters
a = 0.01; b = 0.01; c = 0.01; d = 0.01
# MCMC step
set.seed(90210)
burn  = 10;  M = 1000
niter = burn + M
V1    = V;  W1 = W
draws = NULL
all_draws = NULL
for (iter in 1:niter){
  run   = ffbs(y,V1,W1,m0,C0)
  x     = run$x
  V1    = 1/rgamma(1,a+n/2,b+sum((y-x)^2)/2)
  W1    = 1/rgamma(1,c+(n-1)/2,d+sum(diff(x)^2)/2)
  draws = rbind(draws,c(V1,W1,x))    }
all_draws = draws[,1:2]
q025  = function(x){quantile(x,0.025)}
q975  = function(x){quantile(x,0.975)}
draws = draws[(burn+1):(niter),]
xs    = draws[,3:(n+2)]
lx    = apply(xs,2,q025)
mx    = apply(xs,2,mean)
ux    = apply(xs,2,q975)

##  plot data
par(mfrow=c(2,2))
tsplot(cbind(x,y), spag=TRUE,  ylab='', col=c(1,8), lwd=2)
points(y)
legend(0, 11, legend=c("x(t)","y(t)"), lty=1, col=c(1,8), lwd=2, bty="n", pch=c(-1,1))
contour(Vs, Ws, exp(likes), xlab=expression(sigma[v]^2), ylab=expression(sigma[w]^2), 
         drawlabels=FALSE, ylim=c(0,1.2))
points(draws[,1:2], pch=16, col=rgb(.9,0,0,0.3), cex=.7)
hist(draws[,1], ylab="Density",main="", xlab=expression(sigma[v]^2))
abline(v=mean(draws[,1]), col=3, lwd=3)
hist(draws[,2],main="", ylab="Density", xlab=expression(sigma[w]^2))
abline(v=mean(draws[,2]), col=3, lwd=3)

## plot states
dev.new()
tsplot(y, ylab='', type='o', col=8)
lines(xsmooth, lwd=4, col=rgb(1,0,1,alpha=.4))
lines(mx, col= 4) 
 xx=c(1:100, 100:1)
 yy=c(lx, rev(ux))
polygon(xx, yy, border=NA, col= gray(.6,alpha=.2))
legend('topleft', c('true smoother', 'data', 'posterior mean', '95% of draws'), lty=1, 
         lwd=c(3,1,1,10), pch=c(-1,1,-1,-1), col=c(6, gray(.4), 4, gray(.6, alpha=.5)), 
         bg='white' )  
```

Example 6.27

```r
y = jj
### setup - model and initial parameters
set.seed(90210)
n = length(y)
F = c(1,1,0,0)      # this is A 
G = diag(0,4)       # G is Phi 
  G[1,1] = 1.03 
  G[2,]  = c(0,-1,-1,-1); G[3,]=c(0,1,0,0); G[4,]=c(0,0,1,0)
a1 = rbind(.7,0,0,0)  # this is mu0
R1 = diag(.04,4)      # this is Sigma0
V = .1
W11 = .1
W22 = .1

##-- FFBS --##  
ffbs = function(y,F,G,V,W11,W22,a1,R1){
  n  = length(y)
  Ws = diag(c(W11,W22,1,1))  # this is Q with 1s as a device only
  iW = diag(1/diag(Ws),4)    
  a  = matrix(0,n,4)         # this is m_t
  R  = array(0,c(n,4,4))     # this is V_t
  m  = matrix(0,n,4)
  C  = array(0,c(n,4,4))
  a[1,]  = a1[,1]
  R[1,,] = R1
  f      = t(F)%*%a[1,]
  Q      = t(F)%*%R[1,,]%*%F + V
  A      = R[1,,]%*%F/Q[1,1]
  m[1,]  = a[1,]+A%*%(y[1]-f)
  C[1,,] = R[1,,]-A%*%t(A)*Q[1,1]
  for (t in 2:n){
    a[t,]  = G%*%m[t-1,]
    R[t,,] = G%*%C[t-1,,]%*%t(G) + Ws
    f      = t(F)%*%a[t,]
    Q      = t(F)%*%R[t,,]%*%F + V
    A      = R[t,,]%*%F/Q[1,1]
    m[t,]  = a[t,] + A%*%(y[t]-f)
    C[t,,] = R[t,,] - A%*%t(A)*Q[1,1]      }
  xb       = matrix(0,n,4)
  xb[n,]  = m[n,] + t(chol(C[n,,]))%*%rnorm(4)
  for (t in (n-1):1){
    iC  = solve(C[t,,])
    CCC = solve(t(G)%*%iW%*%G + iC)
    mmm = CCC%*%(t(G)%*%iW%*%xb[t+1,] + iC%*%m[t,])
    xb[t,] = mmm + t(chol(CCC))%*%rnorm(4)  }
  return(xb)                                
}

##-- Prior hyperparameters --##
# b0 = 0     # mean for beta = phi -1
# B0 = Inf   # var for  beta  (non-informative => use OLS for sampling beta)
n0 = 10      # use same for all- the prior is 1/Gamma(n0/2, n0*s20_/2)
s20v = .001  # for V
s20w =.05    # for Ws

##-- MCMC scheme --##
set.seed(90210)
burnin  = 100 
step    = 10   
M       = 1000  
niter   = burnin+step*M
pars    = matrix(0,niter,4)
xbs     = array(0,c(niter,n,4))
pb      = txtProgressBar(min=0, max=niter, initial=0, style=3)  # progress bar
            
for (iter in 1:niter){
    setTxtProgressBar(pb,iter)  
    xb = ffbs(y,F,G,V,W11,W22,a1,R1)
     u = xb[,1] 
    yu = diff(u); xu = u[-n]    # for phihat and se(phihat)
  regu = lm(yu~0+xu)                # est of beta = phi-1
 phies = as.vector(coef(summary(regu)))[1:2] + c(1,0) # phi estimate and SE
   dft = df.residual(regu)   
G[1,1] = phies[1] + rt(1,dft)*phies[2]  # use a t
    V  = 1/rgamma(1, (n0+n)/2, (n0*s20v/2) + sum((y-xb[,1]-xb[,2])^2)/2)
   W11 = 1/rgamma(1, (n0+n-1)/2, (n0*s20w/2) + sum((xb[-1,1]-phies[1]*xb[-n,1])^2)/2)
   W22 = 1/rgamma(1, (n0+ n-3)/2, (n0*s20w/2) + sum((xb[4:n,2] + xb[3:(n-1),2] + 
                  xb[2:(n-2),2] +xb[1:(n-3),2])^2)/2)
   xbs[iter,,] = xb
   pars[iter,] = c(G[1,1], sqrt(V), sqrt(W11), sqrt(W22))           
}
close(pb) 

# Plot results
ind = seq(burnin+1, niter, by=step)
names= c(expression(phi), expression(sigma[v]), expression(sigma[w~11]), expression(sigma[w~22]))
par(mfcol=c(3,4))
for (i in 1:4){
 tsplot(pars[ind,i],xlab="iterations", ylab="trace", main="")
 mtext(names[i], side=3, line=.5, cex=1) 
 acf(pars[ind,i],main="", lag.max=25, xlim=c(1,25), ylim=c(-.4,.4))
 hist(pars[ind,i],main="",xlab="")
 abline(v=mean(pars[ind,i]), lwd=2, col=3) 
}

dev.new()
par(mfrow=c(2,1))
  mxb = cbind(apply(xbs[ind,,1],2,mean), apply(xbs[,,2],2,mean))
  lxb = cbind(apply(xbs[ind,,1],2,quantile,0.005), apply(xbs[ind,,2],2,quantile,0.005))
  uxb = cbind(apply(xbs[ind,,1],2,quantile,0.995), apply(xbs[ind,,2],2,quantile,0.995))   
  mxb = ts(cbind(mxb,rowSums(mxb)), start = tsp(jj)[1], freq=4) 
  lxb = ts(cbind(lxb,rowSums(lxb)), start = tsp(jj)[1], freq=4)
  uxb = ts(cbind(uxb,rowSums(uxb)), start = tsp(jj)[1], freq=4)
  names=c('Trend', 'Season', 'Trend + Season')
  L = min(lxb[,1])-.01; U = max(uxb[,1]) +.01
tsplot(mxb[,1],  ylab=names[1], ylim=c(L,U))
  xx=c(time(jj), rev(time(jj)))
  yy=c(lxb[,1], rev(uxb[,1]))
  polygon(xx, yy, border=NA, col=gray(.4, alpha = .2)) 
  L = min(lxb[,3])-.01; U = max(uxb[,3]) +.01
tsplot(mxb[,3],  ylab=names[3], ylim=c(L,U))
  xx=c(time(jj), rev(time(jj)))
  yy=c(lxb[,3], rev(uxb[,3]))
  polygon(xx, yy, border=NA, col=gray(.4, alpha = .2))            
```

