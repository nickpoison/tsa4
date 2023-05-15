
## R Code Used in the Examples - tsa4 

<img align="left" src="https://github.com/nickpoison/astsa/blob/master/fun_with_astsa/figs/tsa4.jpg" alt="&nbsp; tsa4 &nbsp;"  height="200" />  This is an updated version of the code in [Time Series Analysis and Its Applications, 4th Edition](https://github.com/nickpoison/tsa4) 


#### &#x2728; See the [NEWS](https://github.com/nickpoison/astsa/blob/master/NEWS.md) for further details about the state of the package and the changelog.


#### &#10024; An intro to `astsa` capabilities can be found at  [FUN WITH ASTSA](https://github.com/nickpoison/astsa/blob/master/fun_with_astsa/fun_with_astsa.md)

#### &#10024; Here is [A Road Map](https://nickpoison.github.io/) if you want a broad view of what is available.

<br/><br/>

---
---

>  __Note__ when you are in a code block below, you can copy the contents of the block by moving your mouse to the upper right corner and clicking on the copy icon ( &#128203; ).




-----
------ 

### Table of Contents
  
  * [Chapter 1 - Characteristics of Time Series](#chapter-1)
  * [Chapter 2 - Time Series Regression and Exploratory Data Analysis](#chapter-2)
  * [Chapter 3 - ARIMA Models](#chapter-3)
  * [Chapter 4 - Spectral Analysis and Filtering](#chapter-4)
  * [Chapter 5 - Additional Time Domain Topics](#chapter-5)
  * [Chapter 6 - State Space Models](#chapter-6)
  * [Chapter 7 - Statistical Methods in the Frequency Domain](#chapter-7)
 
---

## Chapter 1


Example 1.1 

```r
tsplot(jj, col=4, type="o", ylab="Quarterly Earnings per Share")
```

Example 1.2  

```r
tsplot(globtemp, col=4, type="o", ylab="Global Temperature Deviations")

# or with the updated values
tsplot(gtemp_land, col=4, type="o", ylab="Global Temperature Deviations")
``` 

Example 1.3  

```r
tsplot(speech)  
``` 

Example 1.4  

```r
library(xts)         # install it if you don't have it
djiar = diff(log(djia$Close))[-1]        
plot(djiar, col=4, main="DJIA Returns") 
```

Example 1.5  

```r
par(mfrow = c(2,1))  # set up the graphics
tsplot(soi, col=4, ylab="", main="Southern Oscillation Index")
tsplot(rec, col=4, ylab="", main="Recruitment") 
```

Example 1.6

```r
par(mfrow=c(2,1))  
tsplot(fmri1[,2:5], col=1:4, ylab="BOLD", main="Cortex", spaghetti=TRUE)
tsplot(fmri1[,6:9], col=1:4, ylab="BOLD", main="Thalamus & Cerebellum", spaghetti=TRUE)

# each separately (not in text)
tsplot(fmri1[,2:9], col=1:8, lwd=2, ncol=2, ylim=c(-.6,.6))

# and another view (not in text)
x     = ts(fmri1[,4:9], start=0, freq=32)         
names = c("Cortex","Thalamus","Cerebellum")
u     = ts(rep(c(rep(.6,16), rep(-.6,16)), 4), start=0, freq=32) # stimulus signal
par(mfrow=c(3,1))
for (i in 1:3){ 
  j = 2*i - 1
  tsplot(x[,j:(j+1)], ylab="BOLD", xlab="", main=names[i], col=5:6, ylim=c(-.6,.6), 
         lwd=2, xaxt="n", spaghetti=TRUE)
  axis(seq(0,256,64), side=1, at=0:4)
  lines(u, type="s", col=gray(.3)) 
}
mtext("seconds", side=1, line=1.75, cex=.9)
```

Example 1.7

```r
par(mfrow=2:1)
tsplot(EQ5,  col=4, main="Earthquake")
tsplot(EXP6, col=4, main="Explosion")

# or try (not in text)
tsplot(cbind(EQ5,EXP6), col=4)
```
 
Example 1.9

```r
w = rnorm(500,0,1)                  # 500 N(0,1) variates
v = filter(w, sides=2, rep(1/3,3))  # moving average
par(mfrow=c(2,1))
tsplot(w, col=4, main="white noise")
tsplot(v, col=4, ylim=c(-3,3), main="moving average")
```

Example 1.10

```r
w = rnorm(550,0,1)  # 50 extra to avoid startup problems
x = filter(w, filter=c(1,-.9), method="recursive")[-(1:50)]
tsplot(x, col=4, main="autoregression")
```

Example 1.11

```r
set.seed(154) # so you can reproduce the results
w = rnorm(200); x = cumsum(w) # two commands in one line
wd = w +.2;    xd = cumsum(wd)
tsplot(xd, ylim=c(-5,55), main="random walk", ylab='')
lines(x, col=4) 
clip(0,200,0,50)
abline(h=0, col=4, lty=2)
abline(a=0, b=.2, lty=2)
```

Example 1.12

```r
cs = 2*cos(2*pi*(1:500)/50 + .6*pi)
w = rnorm(500,0,1)
par(mfrow=c(3,1))
tsplot(cs, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)))
tsplot(cs + w, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,1)))
tsplot(cs + 5*w, ylab="", main = expression(x[t]==2*cos(2*pi*t/50+.6*pi)+N(0,25)))
```

Example 1.24

```r
set.seed(2)
x = rnorm(100)
y = lag(x, -5) + rnorm(100)
ccf2(y, x, ylab='CCovF', type='covariance')
text( 9, 1.1, 'x leads')
text(-8, 1.1, 'y leads')
```

Example 1.25

```r
(r = round( acf1(soi, 6, plot=FALSE), 2)) # sample acf values
par(mfrow=c(1,2))
tsplot(lag(soi,-1), soi, col=4, type='p', xlab='lag(soi,-1)')
 legend("topleft", legend=r[1], bg="white", adj=.45, cex = 0.85)
tsplot(lag(soi,-6), soi, col=4, type='p', xlab='lag(soi,-6)')
 legend("topleft", legend=r[6], bg="white", adj=.25, cex = 0.8)
```

Example 1.26

```r
set.seed(666)
x1 = sample(c(-1,1), 11, replace=TRUE)  # simulated sequence of coin tosses
x2 = sample(c(-1,1), 101, replace=TRUE)
y1 = 5 + filter(x1, sides=1, filter=c(1,-.7))[-1]
y2 = 5 + filter(x2, sides=1, filter=c(1,-.7))[-1]
tsplot(y1, type='s')   # plot 1st series
 points(y1, pch=19)
c(mean(y1), mean(y2))  # the sample means
acf(y1, lag.max=4, plot=FALSE) 
acf(y2, lag.max=4, plot=FALSE) 

#########################################
# here's the version from the other text -
# same idea but the y values are 2-4-6-8
# like the children's cheer

set.seed(101011)
x1 = sample(c(-2,2), 11, replace=TRUE)   # simulated coin tosses
x2 = sample(c(-2,2), 101, replace=TRUE)
y1 = 5 + filter(x1, sides=1, filter=c(1,-.5))[-1]
y2 = 5 + filter(x2, sides=1, filter=c(1,-.5))[-1]
tsplot(y1, type="s", col=4, xaxt="n", yaxt="n")  # y2 not shown
 axis(1, 1:10); axis(2, seq(2,8,2), las=1)
 points(y1, pch=21, cex=1.1, bg=6)
acf(y1, lag.max=4, plot=FALSE) 
acf(y2, lag.max=4, plot=FALSE) 
```

Example 1.27

```r
acf1(speech, 250)
```

Example 1.28

```r
par(mfrow=c(3,1))
acf1(soi, main="Southern Oscillation Index")
acf1(rec, main="Recruitment")
ccf2(soi, rec, main="SOI vs Recruitment")
```

Example 1.29

```r
set.seed(1492)
num = 120 
t   = 1:num
X   = ts(2*cos(2*pi*t/12) + rnorm(num), freq=12)
Y   = ts(2*cos(2*pi*(t+5)/12) + rnorm(num), freq=12)
Yw  = resid( lm(Y~ cos(2*pi*t/12) + sin(2*pi*t/12), na.action=NULL) )
par(mfrow=c(3,2) )
tsplot(X)
tsplot(Y)
acf1(X, 48, ylab='ACF(X)')
acf1(Y, 48, ylab='ACF(Y)')
ccf2(X, Y, 24)
ccf2(X, Yw, 24, ylim=c(-.6,.6))
################################################

#  here's another example that's simpler
#  the series are trend stationary with 
#  just a hint of trend - but same result

set.seed(90210)
num = 250  
t   = 1:num
X   = .01*t + rnorm(num,0,2)
Y   = .01*t + rnorm(num)
par(mfrow=c(3,1))
tsplot(cbind(X,Y), spag=TRUE, col=astsa.col(c(4,2),.7), lwd=2, ylab='data')  
ccf2(X, Y,  ylim=c(-.3,.3), col=4, lwd=2)
Yw = detrend(Y)  # whiten Y by removing trend
ccf2(X, Yw, ylim=c(-.3,.3), col=4, lwd=2)
```

Example 1.30

```r
par(mar=rep(1,4))
persp(1:64, 1:36, soiltemp, phi=30, theta=30, scale=FALSE, expand=4, 
       ticktype="detailed", xlab="rows", ylab="cols", zlab="temperature")
dev.new()          
tsplot(rowMeans(soiltemp), xlab="row", ylab="Average Temperature")
```

Example 1.31

```r
fs = abs(fft(soiltemp-mean(soiltemp)))^2/(64*36) # see Ch 4 for info on FFT
cs = Re(fft(fs, inverse=TRUE)/sqrt(64*36))  # ACovF
rs = cs/cs[1,1]                             # ACF

rs2 = cbind(rs[1:41,21:2], rs[1:41,1:21])   #  these lines are just to center
rs3 = rbind(rs2[41:2,], rs2)                #  the 0 lag  

par(mar = c(1,2.5,0,0)+.1)
persp(-40:40, -20:20, rs3, phi=30, theta=30, expand=30, scale="FALSE",  
       ticktype="detailed", xlab="row lags", ylab="column lags", zlab="ACF")
```

[<sub>top</sub>](#table-of-contents)

---



## Chapter 2


Example 2.1

```r
# astsa now has a trend script, so Figure 2.1 can be done in one line
trend(chicken, lwd=2)    # includes a 95% CI

# in the text
summary(fit <- lm(chicken~time(chicken))) # regress price on time
tsplot(chicken, ylab="cents per pound", col=4, lwd=2)
abline(fit)      # add the fitted regression line to the plot            
```

Example 2.2

```r
##-- separate
par(mfrow=c(3,1))
tsplot(cmort, main="Cardiovascular Mortality", col=6, type="o", pch=19, ylab="")
tsplot(tempr, main="Temperature", col=4, type="o", pch=19, ylab="")
tsplot(part, main="Particulates", col=2, type="o", pch=19, ylab="")

##-- together 
dev.new()
tsplot(cbind(cmort, tempr, part), spag=TRUE, ylab="", col=c(6,4,2))
legend("topright", legend=c("Mortality", "Temperature", "Pollution"), lty=1, lwd=2, col=c(6,4,2), bg="white")

##-- scatterplot matrix
dev.new()  
panel.cor <- function(x, y, ...){
usr <- par("usr")
par(usr = c(0, 1, 0, 1))
r <- round(cor(x, y), 2)
text(0.5, 0.5, r, cex = 1.75)
}
pairs(cbind(Mortality=cmort, Temperature=tempr, Particulates=part), col=4, lower.panel=panel.cor)

#  Regression
temp  = tempr-mean(tempr)  # center temperature    
temp2 = temp^2             # square it  
trend = time(cmort)        # time

fit = lm(cmort~ trend + temp + temp2 + part, na.action=NULL)
            
summary(fit)       # regression results
summary(aov(fit))  # ANOVA table   (compare to next line)
summary(aov(lm(cmort~cbind(trend, temp, temp2, part)))) # Table 2.1

num = length(cmort)                                     # sample size
AIC(fit)/num - log(2*pi)                                # AIC 
BIC(fit)/num - log(2*pi)                                # BIC   
(AICc = log(sum(resid(fit)^2)/num) + (num+5)/(num-5-2)) # AICc
```

Examples 2.3

```r
fish = ts.intersect(rec, soiL6=lag(soi,-6), dframe=TRUE)   
summary(fit <- lm(rec~soiL6, data=fish, na.action=NULL))
## not shown in text but resids are not white
par(mfrow=2:1)
tsplot(resid(fit))
acf1(resid(fit))
```



Examples 2.4 and 2.5

```r
# astsa now has a detrend script, so Figure 2.4 can be done as
par(mfrow=2:1)
tsplot( detrend(chicken), main="detrended" )
tsplot( diff(chicken), main="first difference" )

# and Figure 2.5 as
dev.new()
par(mfrow=c(3,1))     # plot ACFs
acf1(chicken, 48, main="chicken")
acf1(detrend(chicken), 48, main="detrended")
acf1(diff(chicken), 48, main="first difference")
``` 


Example 2.6

```r
par(mfrow=c(2,1))
tsplot(diff(globtemp), type="o")
 mean(diff(globtemp))     # drift estimate = .008
acf1(diff(gtemp), 48, main="")
```

Example 2.7

```r
layout(matrix(1:4,2), widths=c(2.5,1))
tsplot(varve, main="", ylab="", col=4)
 mtext("varve", side=3, line=.5, cex=1.2, font=2, adj=0)
tsplot(log(varve), main="", ylab="", col=4)
 mtext("log(varve)", side=3, line=.5, cex=1.2, font=2, adj=0)
qqnorm(varve, main="", col=4)
 qqline(varve, col=2, lwd=2)
qqnorm(log(varve), main="", col=4)
 qqline(log(varve), col=2, lwd=2)
```

Example 2.8  

```r
lag1.plot(soi, 12, col=astsa.col(4, .3), cex=1.5, pch=20)
dev.new()
lag2.plot(soi, rec, 8, col=astsa.col(4, .3), cex=1.5, pch=20)
```

Example 2.9

```r
dummy = ifelse(soi<0, 0, 1)
fish  = ts.intersect(rec, soiL6=lag(soi,-6), dL6=lag(dummy,-6), dframe=TRUE)
summary(fit <- lm(rec~ soiL6*dL6, data=fish, na.action=NULL))
tsplot(fish$soiL6, fish$rec, type='p', col=4, ylab='rec', xlab='soiL6')
lines(lowess(fish$soiL6, fish$rec), col=4, lwd=2)
points(fish$soiL6, fitted(fit), pch='+', col=6)

dev.new()
par(mfrow=2:1)
tsplot(resid(fit)) # not shown ...
acf1(resid(fit))   # ... but obviously not noise
```


Example 2.10

```r
set.seed(1000)  # so you can reproduce these results
x  = 2*cos(2*pi*1:500/50 + .6*pi) + rnorm(500,0,5)
z1 = cos(2*pi*1:500/50)  
z2 = sin(2*pi*1:500/50)
summary(fit <- lm(x~0+z1+z2))  # zero to exclude the intercept
par(mfrow=c(2,1))
tsplot(x, col=4)
tsplot(x, col=astsa.col(4,.7), ylab=expression(hat(x)))
lines(fitted(fit), col=2, lwd=2)
```


Example 2.11

```r
wgts = c(.5, rep(1,11), .5)/12
soif = filter(soi, sides=2, filter=wgts)
tsplot(soi, col=4)
lines(soif, lwd=2, col=6)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
nwgts = c(rep(0,20), wgts, rep(0,20))
plot(nwgts, type="l", ylim = c(-.02,.1), xaxt='n', yaxt='n', ann=FALSE)
```


Example 2.12 

```r
tsplot(soi, col=4)
lines(ksmooth(time(soi), soi, "normal", bandwidth=1), lwd=2, col=6)
par(fig = c(.75, 1, .75, 1), new = TRUE) # the insert
curve(dnorm, -3, 3,  xaxt='n', yaxt='n', ann=FALSE)
```


Example 2.13
 
```r
# Figure 2.14 using the trend script
trend(soi, lowess=TRUE)
lines(lowess(soi, f=.05), lwd=2, col=6) # El Niño cycle
```



Example 2.14

```r
tsplot(soi)
lines(smooth.spline(time(soi), soi, spar=.5), lwd=2, col=4)
lines(smooth.spline(time(soi), soi, spar= 1), lty=2, lwd=2, col=2)
```

Example 2.15

```r
tsplot(tempr, cmort, type="p", xlab="Temperature", ylab="Mortality", pch=20, col=4)
lines(lowess(tempr, cmort), col=6, lwd=2)
```

[<sub>top</sub>](#table-of-contents)

---



## Chapter 3

> The way AIC, AICc, and BIC are calculated in `sarima` changed a few versions ago. The values in the text will be different than the current values, but the results of any data analysis in the text are the same.

<br/>




Example 3.2

```r
par(mfrow=c(2,1))                         
# in the expressions below, ~ is a space and == is equal
tsplot(sarima.sim(ar= .9, n=100), col=4, ylab="", main=(expression(AR(1)~~~phi==+.9))) 
tsplot(sarima.sim(ar=-.9, n=100), col=4, ylab="", main=(expression(AR(1)~~~phi==-.9))) 
```

Example 3.5

```r
par(mfrow=c(2,1))                                   
tsplot(sarima.sim(ma= .9, n=100), col=4, ylab="", main=(expression(MA(1)~~~theta==+.9)))    
tsplot(sarima.sim(ma=-.9, n=100), col=4, ylab="", main=(expression(MA(1)~~~theta==-.9)))    
```

Example 3.7

```r
set.seed(8675309)         # Jenny, I got your number
x = rnorm(150, mean=5)    # Jenerate iid N(5,1)s
arima(x, order=c(1,0,1))  # Jenstimation
```

Example 3.8

```r
ARMAtoMA(ar = .9, ma = .5, 10)   # first 10 psi-weights
ARMAtoAR(ar = .9, ma = .5, 10)   # first 10 pi-weights
```

Example 3.9
```r
# this is how Figure 3.3 was generated
seg1   =  seq( 0, 2,  by=0.1)
seg2   =  seq(-2, 2,  by=0.1)
name1  =  expression(phi[1])
name2  =  expression(phi[2])
tsplot(seg1, (1-seg1), ylim=c(-1,1), xlim=c(-2,2), ylab=name2, xlab=name1,
        main='Causal Region of an AR(2)')
 lines(-seg1, (1-seg1), ylim=c(-1,1), xlim=c(-2,2)) 
 abline(h=0, v=0, lty=2, col=8)
 lines(seg2, -(seg2^2 /4), ylim=c(-1,1))
 lines( x=c(-2,2), y=c(-1,-1), ylim=c(-1,1))
 text(0, .35, 'real roots')
 text(0, -.5, 'complex roots')
```


Example 3.11

```r
z = c(1,-1.5,.75)    # coefficients of the polynomial
(a = polyroot(z)[1]) # = 1+0.57735i,  print one root which is 1 + i 1/sqrt(3)
arg = Arg(a)/(2*pi)  # arg in cycles/pt  
1/arg                # = 12,  the period

par(mfrow=c(3,1))
set.seed(8675309)    # Jenny, it's me again
ar2 = sarima.sim(ar=c(1.5,-.75), n=144, S=12)
tsplot(ar2, xlab="Year")

ACF = ARMAacf(ar=c(1.5,-.75), ma=0, 50)[-1]
tsplot(ACF, type="h", xlab="lag")
abline(h=0, col=8)

# alternately - not in text
ACF = ARMAacf(ar=c(1.5,-.75), ma=0, 50)
tsplot(0:50, ACF, type="h", xlab="lag")
abline(h=0, col=8)


# psi-weights - not in text
psi = ts(ARMAtoMA(ar=c(1.5,-.75), ma=0, 50), start=0, freq=12)
tsplot(psi, type='o', cex=1.1, ylab=expression(psi-weights), xaxt='n', xlab='Index')
axis(1, at=0:4, labels=c('0','12','24','36','48'))

# you can play the same game with the ACF - not in text
ACF = ts(ARMAacf(ar=c(1.5,-.75), ma=0, 50), start=0, frequency=12)
tsplot(ACF, type='h', xaxt='n', xlab='LAG')
abline(h=0, col=8)
axis(1, at=0:4, labels=c('0','12','24','36','48'))

```

Example 3.12

```r
psi = ARMAtoMA(ar=.9, ma=.5, 50)       #  for a list        
tsplot(psi, type='h', ylab=expression(psi-weights), xlab='Index')    # for a graph
```


Example 3.16

```r
ar2.acf  = ARMAacf(ar=c(1.5,-.75), ma=0, 24)[-1]
ar2.pacf = ARMAacf(ar=c(1.5,-.75), ma=0, 24, pacf=TRUE)
par(mfrow=1:2)
tsplot(ar2.acf,  type="h", xlab="lag", lwd=3, nxm=5, col=c(rep(4,11), 6))
tsplot(ar2.pacf, type="h", xlab="lag", lwd=3, nxm=5, col=4)
```

Example 3.18   

```r
acf2(rec, 48)     # will produce values and a graphic 
(regr = ar.ols(rec, order=2, demean=FALSE, intercept=TRUE))  # regression
regr$asy.se.coef  # standard errors                             
```


Example 3.25  

```r
regr = ar.ols(rec, order=2, demean=FALSE, intercept=TRUE)
fore = predict(regr, n.ahead=24)
tsplot(cbind(rec, fore$pred), spag=TRUE, col=1:2, xlim=c(1980,1990), ylab="Recruitment")
lines(fore$pred, type="p", col=2)
lines(fore$pred+fore$se, lty="dashed", col=4)
lines(fore$pred-fore$se, lty="dashed", col=4)
```

Example 3.26

```r
set.seed(666)
x     = sarima.sim(ar=.9, ma=.5, n=100)
xr    = rev(x) # xr is the reversed data
pxr   = predict(arima(xr, order=c(1,0,1)), 10) # predict the reversed data
pxrp  = rev(pxr$pred) # reorder the predictors (for plotting)
pxrse = rev(pxr$se) # reorder the SEs
nx    = ts(c(pxrp, x), start=-9) # attach the backcasts to the data
tsplot(nx, ylab=expression(X[~t]), main='Backcasting', ylim=c(-7,4))
 U  = nx[1:10] + pxrse
 L  = nx[1:10] - pxrse
 xx = c(-9:0, 0:-9) 
 yy = c(L, rev(U))
polygon(xx, yy, border = 8, col = gray(0.6, alpha = 0.2))
lines(-9:0, nx[1:10], col=2, type='o')
```


Example 3.28

```r
rec.yw = ar.yw(rec, order=2)
rec.yw$x.mean    # = 62.26278 (mean estimate)
rec.yw$ar        # = 1.3315874, -.4445447  (parameter estimates)
sqrt(diag(rec.yw$asy.var.coef))  # = .04222637, .04222637  (standard errors)
rec.yw$var.pred  # = 94.79912 (error variance estimate)

rec.pr = predict(rec.yw, n.ahead=24)
 U = rec.pr$pred + rec.pr$se
 L = rec.pr$pred - rec.pr$se
tsplot(cbind(rec, rec.pr$pred), spag=TRUE, xlim=c(1980,1990),  ylab="Recruitment") 
lines(rec.pr$pred, col=2, type="o")
lines(U, col=4, lty=2)
lines(L, col=4, lty=2)
```

Example 3.29

```r
set.seed(1)
ma1 = sarima.sim(ma=.9, n=50)
acf1(ma1, 1, plot=FALSE)  # [1] .536 (lag 1 sample ACF)
```

Example 3.31

```r
rec.mle = ar.mle(rec, order=2)
rec.mle$x.mean
rec.mle$ar
sqrt(diag(rec.mle$asy.var.coef))
rec.mle$var.pred
```



Example 3.33

```r
x = diff(log(varve))       # data
r = acf1(x, 1, plot=FALSE) # acf(1)
c(0) -> z -> Sc -> Sz -> Szw -> para # initialize .. 
c(x[1]) -> w                         # .. all variables
num = length(x)            # 633

## Gauss-Newton Estimation
para[1] = (1-sqrt(1-4*(r^2)))/(2*r)  # MME to start (not very good)
niter   = 12             
for (j in 1:niter){  
 for (t in 2:num){ w[t] = x[t]   - para[j]*w[t-1]
                   z[t] = w[t-1] - para[j]*z[t-1]
 }
 Sc[j]  = sum(w^2)
 Sz[j]  = sum(z^2)
 Szw[j] = sum(z*w)
para[j+1] = para[j] + Szw[j]/Sz[j]
}
## Results
cbind(iteration=1:niter-1, thetahat=para[1:niter], Sc, Sz)

## Plot conditional SS and results
c(0) -> cSS
th = -seq(.3, .94, .01)
for (p in 1:length(th)){   
 for (t in 2:num){ w[t] = x[t] - th[p]*w[t-1] 
 }
cSS[p] = sum(w^2)
}
tsplot(th, cSS, ylab=expression(S[c](theta)), xlab=expression(theta))
abline(v=para[1:12], lty=2, col=4) # add previous results to plot
points(para[1:12], Sc[1:12], pch=16, col=4)
```



Example 3.36

```r
# generate data
set.seed(101010)   
e   = rexp(150, rate=.5) 
u   = runif(150,-1,1) 
de  = e*sign(u)
dex = 50 + sarima.sim(n=100, ar=.95, innov=de, burnin=50) 
tsplot(dex, ylab=expression(X[~t]))

# Bootstrap
set.seed(666)                 # not that 666
fit     = ar.yw(dex, order=1) # assumes the data were retained
m       = fit$x.mean          # estimate of mean
phi     = fit$ar              # estimate of phi
nboot   = 250                 # number of bootstrap replicates
resids  = fit$resid[-1]       # the 99 innovations
x.star  = dex                 # initialize x*
phi.star.yw = c()             # initialize phi* 

for (i in 1:nboot) {
  resid.star = sample(resids, replace=TRUE)
  x.star     = sarima.sim(n=99, ar=phi, innov=resid.star, burnin=0) + m 
 phi.star.yw[i] = ar.yw(x.star, order=1)$ar
}

# small sample distn
set.seed(111)
phi.yw = rep(NA, 1000)
for (i in 1:1000){
 e = rexp(150, rate=.5); u = runif(150,-1,1); de = e*sign(u)
 x = 50 + arima.sim(n=100,list(ar=.95), innov=de, n.start=50)
 phi.yw[i] = ar.yw(x, order=1)$ar 
}

# Picture
hist(phi.star.yw, 15, main="", prob=TRUE, xlim=c(.65,1.05), ylim=c(0,14), 
      col=astsa.col(4,.3), xlab=expression(hat(phi)))
lines(density(phi.yw, bw=.02), lwd=2)  
curve(dnorm(x, mean=.96, sd=.03), .75,1.1, lty=2, lwd=2, add=TRUE)
legend(.65, 14, bty='n', lty=c(1,0,2), lwd=c(2,0,2), col=1, pch=c(NA,22,NA), 
      pt.bg=c(NA,astsa.col(4,.3),NA), pt.cex=2.5, 
      legend=c('true distribution', 'bootstrap distribution', 'normal approximation'))
```


Example 3.38

```r
set.seed(666)    
x = sarima.sim(ma = -0.8, d=1, n = 100)
(x.ima = HoltWinters(x, beta=FALSE, gamma=FALSE))  # &alpha; is 1-&lambda; here
plot(x.ima)
```

Example 3.39, 3.40, and 3.43

```r
tsplot(gnp, col=4)

dev.new()
acf2(gnp, 50)               # compare to acf2(1:250, 50)       
gnpgr = diff(log(gnp))      # growth rate

dev.new()
tsplot(gnpgr, col=4)

dev.new()
acf2(gnpgr, 24)  
sarima(gnpgr, 1, 0, 0)      # AR(1)
sarima(gnpgr, 0, 0, 2)      # MA(2) 
ARMAtoMA(ar=.35, ma=0, 10)  # prints psi-weights
```

Example 3.41

```r
sarima(log(varve), 0,1,1, no.constant=TRUE, gg=TRUE, col=4)   # ARIMA(0,1,1)

dev.new()
sarima(log(varve), 1,1,1, no.constant=TRUE, gg=TRUE, col=4)   # ARIMA(1,1,1)
```

Example 3.44

```r
trend  = time(cmort) 
temp   = tempr - mean(tempr)
temp2  = temp^2
summary(fit <- lm(cmort~trend + temp + temp2 + part, na.action=NULL))
acf2(resid(fit), 52) # implies AR2
sarima(cmort, 2,0,0, xreg=cbind(trend,temp,temp2,part) )
```

Example 3.45

```r
# Note: this could benefit from a seasonal model fit, but it hasn't
#  been talked about yet - you could come back to this after the next section
dummy = ifelse(soi<0, 0, 1)
fish = ts.intersect(rec, soiL6=lag(soi,-6), dL6=lag(dummy,-6), dframe=TRUE)
summary(fit <- lm(rec ~soiL6*dL6, data=fish, na.action=NULL))
plot(resid(fit))
acf2(resid(fit))     # indicates AR(2)
intract = fish$soiL6*fish$dL6  # interaction term
sarima(fish$rec, 2,0,0, xreg = cbind(fish$soiL6, fish$dL6, intract))
```



Example 3.46

```r
set.seed(666)
SAR = sarima.sim(sar=.9, S=12, n=37) + 50
layout(matrix(c(1,2, 1,3), nc=2), heights=c(1.5,1))
tsplot(SAR, type="c", xlab="Year")
 abline(v=1:3, col=4, lty=2)
 points(SAR, pch=Months, cex=1.35, font=4, col=1:6)

phi  = c(rep(0,11),.9)
ACF  = ARMAacf(ar=phi, ma=0, 100)[-1] # [-1] removes 0 lag
PACF = ARMAacf(ar=phi, ma=0, 100, pacf=TRUE)
 LAG = 1:100/12
tsplot(LAG, ACF, type="h", xlab="LAG", ylim=c(-.04,1))
 abline(h=0, col=8)
tsplot(LAG, PACF, type="h", xlab="LAG", ylim=c(-.04,1))
 abline(h=0, col=8)
```


Example 3.47

```r
par(mfrow=1:2)
phi  = c(rep(0,11),.8)
ACF  = ARMAacf(ar=phi, ma=-.5, 50)[-1]
PACF = ARMAacf(ar=phi, ma=-.5, 50, pacf=TRUE)
 LAG = 1:50/12
tsplot(LAG, ACF, type="h", xlab="LAG", ylim=c(-.4,.8), col=4, lwd=2)
 abline(h=0, col=8)
tsplot(LAG, PACF, type="h", xlab="LAG", ylim=c(-.4,.8), col=4, lwd=2)
 abline(h=0, col=8)
```


Example 3.49

```r
x     = AirPassengers
lx    = log(x) 
dlx   = diff(lx) 
ddlx  = diff(dlx, 12)
tsplot(cbind(x, lx, dlx, ddlx), main="")

# below of interest for showing seasonal persistence (not shown here):
par(mfrow=c(2,1))
monthplot(dlx)
monthplot(ddlx)

sarima(lx, 1,1,1, 0,1,1, 12)   # model 1
sarima(lx, 0,1,1, 0,1,1, 12)   # model 2 (the winner)
sarima(lx, 1,1,0, 0,1,1, 12)   # model 3

dev.new()
sarima.for(lx, 12, 0,1,1, 0,1,1,12)  # forecasts
``` 

[<sub>top</sub>](#table-of-contents)

---

## Chapter 4


Example 4.1
```r
x1 = 2*cos(2*pi*1:100*6/100)  + 3*sin(2*pi*1:100*6/100)
x2 = 4*cos(2*pi*1:100*10/100) + 5*sin(2*pi*1:100*10/100)
x3 = 6*cos(2*pi*1:100*40/100) + 7*sin(2*pi*1:100*40/100)
x = x1 + x2 + x3 

par(mfrow=c(2,2))
tsplot(x1, ylim=c(-10,10), main = expression(omega==6/100~~~A^2==13))
tsplot(x2, ylim=c(-10,10), main = expression(omega==10/100~~~A^2==41))
tsplot(x3, ylim=c(-10,10), main = expression(omega==40/100~~~A^2==85))
tsplot(x, ylim=c(-16,16), main="sum")
```

Example 4.2
```r
# x from Example 4.1 is used here
dev.new()
P = abs(2*fft(x)/100)^2
Fr = 0:99/100                    
tsplot(Fr, P, type="o", xlab="frequency", ylab="periodogram")
abline(v=.5, lty=2)
```

Example 4.3
```r
# modulation
t = 1:200
tsplot(x <- 2*cos(2*pi*.2*t)*cos(2*pi*.01*t))     
lines(cos(2*pi*.19*t)+cos(2*pi*.21*t), col=2)     # the same
Px = Mod(fft(x))^2 
tsplot(0:199/200, Px, type='o')                   # the periodogram

# star mag analysis
n    = length(star)
par(mfrow=c(2,1))
tsplot(star, ylab="star magnitude", xlab="day")
Per   = Mod(fft(star-mean(star)))^2/n
Freq  = (1:n -1)/n
tsplot(Freq[1:50], Per[1:50], type='h', lwd=3, ylab="Periodogram", xlab="Frequency")
text(.05,  7000, "24 day cycle") 
text(.027, 9000, "29 day cycle")
#- a list to help find the peaks
round( cbind(1/Freq[1:30], Per[1:30]), 3)
```




Examples 4.5, 4.6, 4.7
```r
# default is to not plot on log scale 
# add log='y' otherwise  
par(mfrow=c(3,1))
arma.spec(main="White Noise", col=4)
arma.spec(ma=.5, main="Moving Average", col=4)
arma.spec(ar=c(1,-.9), main="Autoregression", col=4)
```


Example 4.10
```r
x      = c(1,2,3,2,1)        
c1     = cos(2*pi*1:5*1/5)
s1     = sin(2*pi*1:5*1/5) 
c2     = cos(2*pi*1:5*2/5)
s2     = sin(2*pi*1:5*2/5)
omega1 = cbind(c1, s1)
omega2 = cbind(c2, s2)
anova(lm(x~ omega1+omega2) )  # ANOVA Table
Mod(fft(x))^2/5               # the periodogram (as a check)
```

Example 4.13
```r
par(mfrow=c(2,1))      
soi.per = mvspec(soi)             
 abline(v=1/4, lty="dotted")
rec.per = mvspec(rec) 
 abline(v=1/4, lty="dotted")

soi.per$details[1:50,] 
  #       frequency  period spectrum
  #  [9,]     0.225  4.4444   0.0309
  # [10,]     0.250  4.0000   0.0537
  # [11,]     0.275  3.6364   0.0754
  # [12,]     0.300  3.3333   0.0567
  # 
  # [39,]     0.975  1.0256   0.0167
  # [40,]     1.000  1.0000   0.9722
  # [41,]     1.025  0.9756   0.0054

# conf intervals -  returned value:
U = qchisq(.025,2)    # 0.05063  
L = qchisq(.975,2)    # 7.37775
2*soi.per$spec[10]/L  # 0.01456
2*soi.per$spec[10]/U  # 2.12220
2*soi.per$spec[40]/L  # 0.26355
2*soi.per$spec[40]/U  # 38.40108

# Repeat lines above using rec in place of soi
```

Example 4.14
```r
soi.ave = mvspec(soi, kernel('daniell',4))
abline(v = c(.25,1,2,3), lty=2)
soi.ave$bandwidth      # = 0.225
df  = soi.ave$df       # df = 16.9875  
U   = qchisq(.025, df) # U = 7.555916
L   = qchisq(.975, df) # L = 30.17425
soi.ave$spec[10]       # 0.0495202
soi.ave$spec[40]       # 0.1190800
# intervals
df*soi.ave$spec[10]/L  # 0.0278789
df*soi.ave$spec[10]/U  # 0.1113333
df*soi.ave$spec[40]/L  # 0.0670396
df*soi.ave$spec[40]/U  # 0.2677201

# Repeat above commands with soi replaced by rec, for example:
rec.ave = mvspec(rec, k)
abline(v=c(.25,1,2,3), lty=2)
# and so on.
```

Example 4.15
```r
t = seq(0, 1, by=1/200) 
amps = c(1, .5, .4, .3, .2, .1)
x = matrix(0, 201, 6)
for (j in 1:6) x[,j] = amps[j]*sin(2*pi*t*2*j)
x = ts(cbind(x, rowSums(x)), start=0, deltat=1/200)               
tsplot(x, lty=c(1:6, 1), lwd=c(rep(1,6), 2), ylab="Sinusoids", col=1:6, spaghetti=TRUE)
names = c("Fundamental","2nd Harmonic","3rd Harmonic","4th Harmonic","5th Harmonic", 
          "6th Harmonic","Formed Signal")
legend("topright", names, lty=c(1:6, 1), lwd=c(rep(1,6), 2), col=1:6)
rm(t)                    #Redemption

##########################################################################
# another view of the idea, sawtooth signal periodic but not sinusoidal  #
##########################################################################
y = ts(rev(1:100 %% 20), freq=20)         # sawtooth signal
par(mfrow=2:1)
tsplot(1:100, y, ylab="sawtooth signal", col=4)
mvspec(y, main="", ylab="periodogram", col=5, xlim=c(0,7))  
```



Example 4.16
```r
kernel("modified.daniell", c(3,3))          # for a list
plot(kernel("modified.daniell", c(3,3)))    # for a graph

k        = kernel("modified.daniell", c(3,3))
soi.smo  = mvspec(soi, kernel=k, taper=.1)
abline(v = c(.25,1), lty=2)
## Repeat above lines with rec replacing soi 
soi.smo$df           # df = 17.42618
soi.smo$bandwidth    # B  = 0.2308103

# An easier way to obtain soi.smo:
soi.smo = mvspec(soi, spans=c(7,7), taper=.1, nxm=4)

# hightlight El Nino cycle
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5,.2))
mtext("1/4", side=1, line=0, at=.25, cex=.75)     
```



Example 4.17
```r
s0 = mvspec(soi, spans=c(7,7), plot=FALSE)             # no taper
s50 = mvspec(soi, spans=c(7,7), taper=.5, plot=FALSE)  # full taper
tsplot(s50$freq, s50$spec, log="y", type="l", ylab="spectrum", xlab="frequency") 
lines(s0$freq, s0$spec, lty=2) 
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.9)
legend(5,.04, legend=c('full taper', 'no taper'), lty=1:2)

text(1.42, 0.04, 'leakage', cex=.8)
arrows(1.4, .035, .75, .009, length=0.05,angle=30)   
arrows(1.4, .035, 1.21, .0075, length=0.05,angle=30)
par(fig = c(.65, 1, .65, 1),  new = TRUE, cex=.5,  mgp=c(0,-.1,0), tcl=-.2)
taper <- function(x) { .5*(1+cos(2*pi*x))  }
 x <- seq(from = -.5, to = .5, by = 0.001)
plot(x, taper(x), type = "l",  lty = 1,  yaxt='n', ann=FALSE)
```


Example 4.18
```r
# AR spectrum - AIC picks order=15
u <- spec.ic(soi,  detrend=TRUE, col=4, lwd=2, nxm=4)  
# plot AIC and BIC
dev.new()
tsplot(0:30, u[[1]][,2:3], type='o', col=2:3, xlab='ORDER', nxm=5, lwd=2, gg=TRUE)  
``` 


Example 4.21
```r
sr = mvspec(cbind(soi,rec), kernel("daniell",9), plot.type="coh")
sr$df                     # df = 35.8625
f = qf(.999, 2, sr$df-2)  # f = 8.529792
C = f/(18+f)              # C = 0.3188779
abline(h = C)
```


Example 4.22
```r
par(mfrow=c(3,1))
tsplot(soi, col=4)                         # plot data
tsplot(diff(soi), col=4)                   # plot first difference
k = kernel("modified.daniell", 6)          # filter weights
tsplot(soif <- kernapply(soi, k), col=4)   # plot 12 month filter
dev.new()
mvspec(soif, spans=9, lwd=2, col=5, nxm=4, taper=.1) # spectral analysis (not shown)
rect(1/7, -1e5, 1/3, 1e5, density=NA, col=gray(.5,.2))
mtext("1/4", side=1, line=0, at=.25, cex=.75)
dev.new()
##-- frequency responses --##
par(mfrow=c(2,1))
w = seq(0, .5, by=.01)
FRdiff = abs(1-exp(2i*pi*w))^2
tsplot(w, FRdiff, xlab='frequency')
u = cos(2*pi*w)+cos(4*pi*w)+cos(6*pi*w)+cos(8*pi*w)+cos(10*pi*w)
FRma = ((1 + cos(12*pi*w) + 2*u)/12)^2
tsplot(w, FRma, xlab='frequency')
```


Example 4.24
```r
LagReg(soi, rec, L=15, M=32, threshold=6)
 dev.new()
LagReg(rec, soi, L=15, M=32, inverse=TRUE, threshold=.01)
 dev.new()
fish = ts.intersect(R=rec, RL1=lag(rec,-1), SL5=lag(soi,-5))
(u = lm(fish[,1]~fish[,2:3], na.action=NULL))
acf2(resid(u))       # suggests ar1
sarima(fish[,1], 1, 0, 0, xreg=fish[,2:3], details=FALSE) 
```

Example 4.25
```r
SigExtract(soi, L=9, M=64, max.freq=.05) 
```

Example 4.26
```r
per = abs(fft(soiltemp-mean(soiltemp))/sqrt(64*36))^2       
per2 = cbind(per[1:32,18:2], per[1:32,1:18])   # this and line below is just rearranging
per3 = rbind(per2[32:2,], per2)                # results to get 0 frequency in the middle

par(mar=c(1,2.5,0,0)+.1)
persp(-31:31/64, -17:17/36, per3, phi=30, theta=30, expand=.6, ticktype="detailed", xlab="cycles/row", 
       ylab="cycles/column", zlab="Periodogram Ordinate")
```

[<sub>top</sub>](#table-of-contents)

---

## Chapter 5


Example 5.1

```r
# NOTE: The example in the text uses the package 'fracdiff', 
#       which is a dinosaur and gave questionable results - 
#       this uses 'arfima' but it didn't make it into the text.
library(arfima)
summary(varve.fd <- arfima(log(varve)))  # d.hat = 0.3728, se(d,hat) = 0.0273
# residual stuff
innov = resid(varve.fd)  
par(mfrow=2:1)
tsplot(innov[[1]])  
acf1(innov[[1]])  

# plot pi wgts
dev.new()
p = rep(1,31)
for (k in 1:30){ p[k+1] = (k-coef(varve.fd)[1])*p[k]/(k+1) }
tsplot(p[-1], ylab=expression(pi[j](d)), xlab="Index (j)", type="h", lwd=4, col=2:7, nxm=5)
```


Example 5.2
```r
series = log(varve)  # specify series to be analyzed
d0 = .1              # initial value of d
n.per = nextn(length(series))
m = (n.per)/2  - 1
per = abs(fft(series-mean(series))[-1])^2  # remove 0 freq
per = per/n.per      # R doesn't scale fft by sqrt(n)
g = 4*(sin(pi*((1:m)/n.per))^2)

# Function to calculate -log.likelihood
whit.like = function(d){  
 g.d=g^d
 sig2 = (sum(g.d*per[1:m])/m)
 log.like = m*log(sig2) - d*sum(log(g)) + m
 return(log.like)   
}

# Estimation (?optim for details - output not shown)
(est = optim(d0, whit.like, gr=NULL, method="L-BFGS-B", 
     hessian=TRUE, lower=-.5, upper=.5, control=list(trace=1,REPORT=1)))

# Results  [d.hat = .380, se(dhat) = .028]  
cat("d.hat =", est$par, "se(dhat) = ",1/sqrt(est$hessian),"\n")  
g.dhat = g^est$par
sig2 = sum(g.dhat*per[1:m])/m
cat("sig2hat =",sig2,"\n")  # sig2hat = .229 

# compart AR spectrum to long memory spectrum
u = spec.ic(log(varve), log='y', lty=2, xlim=c(0,.25), ylim=c(.2,20), col=4)        
g = 4*(sin(pi*((1:500)/2000))^2)
fhat = sig2*g^{-est$par}             # long memory spectral estimate          
lines(1:500/2000, fhat, col=6)
ar.mle(log(varve))                   # to get AR(8) estimates 

# 'fracdiff' has a GPH method, but I don't trust the pacakge
# library(fracdiff)
# fdGPH(log(varve), bandw=.9)   # m = n^bandw- it's supposed to be small- this is way too big 
```


Example 5.3
```r
library(tseries)
adf.test(log(varve), k=0)  # DF test
adf.test(log(varve))       # ADF test
pp.test(log(varve))        # PP test
```

Example 5.4
```r
gnpgr = diff(log(gnp))          # get the returns
u     = sarima(gnpgr, 1, 0, 0)  # fit an AR(1)
acf2(resid(u$fit), 20)          # get (p)acf of the squared residuals
 
library(fGarch)
summary(garchFit(~arma(1,0)+garch(1,0), gnpgr))
```



Example 5.5 and 5.6 
```r
library(xts)   # needed to handle djia
djiar = diff(log(djia$Close))[-1]
acf2(djiar)    # exhibits some autocorrelation (not shown)
acf2(djiar^2)  # oozes autocorrelation (not shown)
library(fGarch)
# GARCH fit
summary(djia.g <- garchFit(~arma(1,0)+garch(1,1), data=djiar, cond.dist='std'))
plot(djia.g)    # to see all plot options
# APARCH fit
summary(djia.ap <- garchFit(~arma(1,0)+aparch(1,1), data=djiar, cond.dist='std'))
plot(djia.ap)
```

Example 5.7
```r
tsplot(flu, type="c")
Months = c("J","F","M","A","M","J","J","A","S","O","N","D")
points(flu, pch=Months, cex=.8, font=2)
# Start analysis
dflu = diff(flu)
dev.new()
lag1.plot(dflu, corr=FALSE) # scatterplot with lowess fit
thrsh = .05 # threshold
Z = ts.intersect(dflu, lag(dflu,-1), lag(dflu,-2), lag(dflu,-3),
lag(dflu,-4) )
ind1 = ifelse(Z[,2] < thrsh, 1, NA) # indicator < thrsh
ind2 = ifelse(Z[,2] < thrsh, NA, 1) # indicator >= thrsh
X1 = Z[,1]*ind1
X2 = Z[,1]*ind2
summary(fit1 <- lm(X1~ Z[,2:5]) ) # case 1
summary(fit2 <- lm(X2~ Z[,2:5]) ) # case 2
D = cbind(rep(1, nrow(Z)), Z[,2:5]) # design matrix
p1 = D %*% coef(fit1) # get predictions
p2 = D %*% coef(fit2)
prd = ifelse(Z[,2] < thrsh, p1, p2)
dev.new()
tsplot(dflu, ylim=c(-.5,.5), type='p', pch=3)
lines(prd)
prde1 = sqrt(sum(resid(fit1)^2)/df.residual(fit1) )
prde2 = sqrt(sum(resid(fit2)^2)/df.residual(fit2) )
prde = ifelse(Z[,2] < thrsh, prde1, prde2)
tx = time(dflu)[-(1:4)]
xx = c(tx, rev(tx))
yy = c(prd-2*prde, rev(prd+2*prde))
polygon(xx, yy, border=8, col=gray(.6, alpha=.25) )
abline(h=.05, col=4, lty=6)

# Using tsDyn (not in text)
library(tsDyn)        
# vignette("tsDyn")   # for package details (it's quirky, so you'll need this)
dflu = diff(flu)
(u = setar(dflu, m=4, thDelay=0))  # fit model and view results (thDelay=0 is lag 1 delay)
BIC(u); AIC(u)                     # if you want to try other models ... m=3 works well too
plot(u)                            # graphics -  ?plot.setar for information
```

Example 5.8 and 5.9 
```r
soi.d   = resid(lm(soi~time(soi), na.action=NULL)) # detrended SOI
acf2(soi.d)
fit     = arima(soi.d, order=c(1,0,0))
ar1     = as.numeric(coef(fit)[1]) # = 0.5875
soi.pw  = resid(fit)
rec.fil = filter(rec, filter=c(1, -ar1), sides=1)
dev.new()
ccf2(soi.pw, rec.fil) 

fish  = ts.intersect(rec, RL1=lag(rec,-1), SL5=lag(soi.d,-5))
(u    = lm(fish[,1]~fish[,2:3], na.action=NULL))
dev.new()
acf2(resid(u)) # suggests ar1
(arx  = sarima(fish[,1], 1, 0, 0, xreg=fish[,2:3])) # final model
pred  = rec + resid(arx$fit) # 1-step-ahead predictions
dev.new()
tsplot(pred, col=astsa.col(8,.3), lwd=7, ylab='rec & prediction')
lines(rec)
```

Example 5.10 and 5.11 
```r
library(vars)
x = cbind(cmort, tempr, part)
summary(VAR(x, p=1, type="both"))  # "both" fits constant + trend

VARselect(x, lag.max=10, type="both")
summary(fit <- VAR(x, p=2, type="both"))
acf(resid(fit), 52)
serial.test(fit, lags.pt=12, type="PT.adjusted")

(fit.pr = predict(fit, n.ahead = 24, ci = 0.95))  # 4 weeks ahead
dev.new()
fanchart(fit.pr)  # plot prediction + error
```

Example 5.12 
```r
library(marima)
model   = define.model(kvar=3, ar=c(1,2), ma=c(1))
arp     = model$ar.pattern 
map     = model$ma.pattern
cmort.d = resid(detr <- lm(cmort~ time(cmort), na.action=NULL))
xdata   = matrix(cbind(cmort.d, tempr, part), ncol=3)  # strip ts attributes
fit     = marima(xdata, ar.pattern=arp, ma.pattern=map, means=c(0,1,1), penalty=1)
# resid analysis (not displayed)
innov   = t(resid(fit))
tsplot(innov) 
acfm(innov)    # since astsa v1.13.2
# acf(innov, na.action = na.pass)  # or use this

# fitted values for cmort
pred    = ts(t(fitted(fit))[,1], start=start(cmort), freq=frequency(cmort)) +
detr$coef[1] + detr$coef[2]*time(cmort)
plot(pred, ylab="Cardiovascular Mortality", lwd=2, col=4)
points(cmort)
# print estimates and corresponding t^2-statistic
short.form(fit$ar.estimates, leading=FALSE)
short.form(fit$ar.fvalues,   leading=FALSE)
short.form(fit$ma.estimates, leading=FALSE)
short.form(fit$ma.fvalues,   leading=FALSE)
fit$resid.cov # estimate of noise cov matrix
```

[<sub>top</sub>](#table-of-contents)

---

## Chapter 6

> __Warning__ The code here uses the updated scripts in `astsa` _version 2.0._   Details of the updates are in the help files of `Kfilter`, `Ksmooth`, and `EM`.  Original code (prior to version 2.0) may be found here: [Original Chapter 6 Info and Code](https://github.com/nickpoison/tsa4/blob/master/chap6.md)




Example 6.1
```r
tsplot(blood, type='o', col=c(6,4,2), lwd=2, pch=19, cex=1) 
```

Example 6.2
```r
tsplot(cbind(gtemp_land, gtemp_ocean), spaghetti=TRUE, lwd=2, pch=20, type="o", 
        col=astsa.col(c(4,2),.5), ylab="Temperature Deviations", main="Global Warming")
legend("topleft", legend=c("Land Surface", "Sea Surface"), lty=1, pch=20, col=c(4,2), bg="white")
```

Example 6.5
```r
# generate data 
set.seed(1)  
num = 50
w   = rnorm(num+1,0,1)
v   = rnorm(num,0,1)
mu  = cumsum(w)     # states:  mu[0], mu[1], . . ., mu[50] 
y   = mu[-1] + v    # obs:  y[1], . . ., y[50]

# filter and smooth (Ksmooth does both)
mu0 = 0;  sigma0 = 1;  phi = 1;  sQ = 1;  sR = 1   
ks = Ksmooth(y, A=1, mu0, sigma0, phi, sQ, sR)   

# pictures 
par(mfrow=c(3,1))

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Prediction", ylim=c(-5,10)) 
  lines(ks$Xp, col=6)
  lines(ks$Xp+2*sqrt(ks$Pp), lty=6, col=6)
  lines(ks$Xp-2*sqrt(ks$Pp), lty=6, col=6)

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Filter", ylim=c(-5,10)) 
  lines(ks$Xf, col=6)
  lines(ks$Xf+2*sqrt(ks$Pf), lty=6, col=6)
  lines(ks$Xf-2*sqrt(ks$Pf), lty=6, col=6)

tsplot(mu[-1], type='p', col=4, pch=19, ylab=expression(mu[~t]), main="Smoother", ylim=c(-5,10)) 
  lines(ks$Xs, col=6)
  lines(ks$Xs+2*sqrt(ks$Ps), lty=6, col=6)
  lines(ks$Xs-2*sqrt(ks$Ps), lty=6, col=6)

mu[1]; ks$X0n; sqrt(ks$P0n)   # initial value info
```


Example 6.6 
```r
# Generate Data
set.seed(999)
num = 100
N = num+1
x = sarima.sim(n=N, ar=.8)
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
  kf = Kfilter(y,A=1,mu0=0,Sigma0,phi,sigw,sigv)
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

##- ALSO note that regression parameters not used are set to NULL (the default) instead of 0 now. 

# Setup 
y = cbind(globtemp/sd(globtemp), globtempl/sd(globtempl))
num = nrow(y)
input = rep(1,num)
A = matrix(c(1,1), nrow=2)
mu0 = -.35; Sigma0 = 1;  Phi = 1

# Function to Calculate Likelihood 
Linn=function(para){
 sQ = para[1]      # sigma_w
  sR1 = para[2]    # 11 element of  sR
  sR2 = para[3]    # 22 element of sR
  sR21 = para[4]   # 21 element of sR
 sR = matrix(c(sR1,sR21,0,sR2), 2)  # put the matrix together
 drift = para[5]
 kf = Kfilter(y,A,mu0,Sigma0,Phi,sQ,sR,Ups=drift,Gam=NULL,input)  # NOTE Gamma is set to NULL now (instead of 0)
 return(kf$like) 
 }

# Estimation
init.par = c(.1,.1,.1,0,.05)  # initial values of parameters
(est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1))) 
SE = sqrt(diag(solve(est$hessian))) 

# Summary of estimation  
estimate = est$par; u = cbind(estimate, SE)
rownames(u)=c("sigw","sR11", "sR22", "sR21", "drift"); u  

# Smooth (first set parameters to their final estimates)
sQ    = est$par[1]  
 sR1  = est$par[2]   
 sR2  = est$par[3]   
 sR21 = est$par[4]  
sR    = matrix(c(sR1,sR21,0,sR2), 2)
(R    = sR%*%t(sR))   #  to view the estimated R matrix
drift = est$par[5]  
ks    = Ksmooth(y,A,mu0,Sigma0,Phi,sQ,sR,Ups=drift,Gam=NULL,input)  # NOTE Gamma is set to NULL now (instead of 0)

# Plot 
tsplot(y, spag=TRUE, margins=.5, type='o', pch=2:3, col=4:3, lty=6, ylab='Temperature Deviations')
xsm  = ts(as.vector(ks$Xs), start=1880)
rmse = ts(sqrt(as.vector(ks$Ps)), start=1880)
lines(xsm, lwd=2, col=6)
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
mu0 = 0; Sigma0 = 2.8
( em = EM(y, 1, mu0, Sigma0, phi, q, r) )   

# Standard Errors  (this uses nlme)
phi = em$Phi; sq = sqrt(em$Q); sr = sqrt(em$R)
mu0 = em$mu0; Sigma0 = em$Sigma0
para = c(phi, sq, sr)
 # evaluate likelihood at estimates 
Linn=function(para){
  kf = Kfilter(y, A=1, mu0, Sigma0, para[1], para[2], para[3])
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
y    = blood
num  = nrow(y)       
A    = array(0, dim=c(3,3,num))   
for(k in 1:num) if (!is.na(y[k,1])) A[,,k]= diag(1,3) 

# Initial values 
mu0    = matrix(0,3,1) 
Sigma0 = diag(c(.1,.1,1) ,3)
Phi    = diag(1,3)
Q      = diag(c(.1,.1,1), 3)
R     = diag(c(.1,.1,1), 3)  
( em = EM(y, A, mu0, Sigma0, Phi, Q, R) ) 

# Graph smoother
sQ = em$Q %^% .5
sR = sqrt(em$R)
ks  = Ksmooth(y, A, em$mu0, em$Sigma0, em$Phi, sQ , sR)
y1s = ks$Xs[1,,] 
y2s = ks$Xs[2,,] 
y3s = ks$Xs[3,,]
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
A   = cbind(1,1,0,0) 

# Function to Calculate Likelihood 
Linn = function(para){
 Phi = diag(0,4) 
 Phi[1,1] = para[1] 
 Phi[2,]=c(0,-1,-1,-1); Phi[3,]=c(0,1,0,0); Phi[4,]=c(0,0,1,0)
 sQ1 = para[2]; sQ2 = para[3]     # sqrt q11 and sqrt q22
 sQ  = diag(0,4); sQ[1,1]=sQ1; sQ[2,2]=sQ2
 sR = para[4]                     # sqrt r11
 kf = Kfilter(jj, A, mu0, Sigma0, Phi, sQ, sR)
 return(kf$like)  
 }

# Initial Parameters 
mu0      = c(.7,0,0,0) 
Sigma0   = diag(.04, 4)  
init.par = c(1.03, .1, .1, .5)   # Phi[1,1], the 2 Qs and R

# Estimation
est = optim(init.par, Linn, NULL, method="BFGS", hessian=TRUE, control=list(trace=1,REPORT=1))
SE  = sqrt(diag(solve(est$hessian)))
u   = cbind(estimate=est$par,SE)
rownames(u)=c("Phi11","sigw1","sigw2","sigv"); u 

# Smooth
Phi      = diag(0,4) 
Phi[1,1] = est$par[1]; Phi[2,]  = c(0,-1,-1,-1) 
Phi[3,]  = c(0,1,0,0); Phi[4,]  = c(0,0,1,0)
sQ       = diag(0,4)
sQ[1,1]  = est$par[2]
sQ[2,2]  = est$par[3]   
sR       = est$par[4]   
ks       = Ksmooth(jj, A, mu0, Sigma0, Phi, sQ, sR)   

# Plots
Tsm   = ts(ks$Xs[1,,], start=1960, freq=4)
Ssm   = ts(ks$Xs[2,,], start=1960, freq=4)
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
x00     = ks$Xf[,,num]
P00     = ks$Pf[,,num]
Q       = sQ%*%t(sQ) 
R       = sR%*%t(sR)
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

#################################################
##  All estimates at once                      ##
#################################################
trend   = time(cmort) - mean(time(cmort)) # center time
const   = time(cmort)/time(cmort)         # appropriate time series of 1s
ded     = ts.intersect(M=cmort, T1=lag(tempr,-1), P=part, P4=lag(part,-4), trend, const)
y       = ded[,1]
input   = ded[,2:6]
num     = length(y)
A       = matrix(c(1,0), 1, 2)  

# Function to Calculate Likelihood
Linn=function(para){
 phi1   = para[1]; phi2 = para[2]; sR = para[3];  b1 = para[4]
 b2     = para[5];   b3 = para[6]; b4 = para[7]; alf = para[8]
 mu0    = matrix(c(0,0), 2, 1)
 Sigma0 = diag(100, 2)
 Phi    = matrix(c(phi1, phi2, 1, 0), 2)
 S      = 1
 Ups    = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
 Gam    = matrix(c(0, 0, 0, b4, alf), 1, 5) 
 sQ     = matrix(c(phi1, phi2), 2)*sR
# S      = sR^2
 kf     = Kfilter(y, A, mu0, Sigma0, Phi, sQ, sR, Ups=Ups, Gam=Gam, input=input, S=S, version=2)
return(kf$like) 
}

# Estimation - prelim analysis gives good starting values
init.par = c(phi1=.3, phi2=.3, sR=4, b1=-.1, b2=.1, b3=.04, b4=-1.3, alf=mean(cmort)) 
L = c( 0,  0,  1, -1,  0,  0, -2, 70)   # lower bound on parameters
U = c(.5, .5, 10,  0, .5, .5,  0, 90)   # upper bound - used in optim
est      = optim(init.par, Linn, NULL, method='L-BFGS-B', lower=L, upper=U, 
                 hessian=TRUE, control=list(trace=1, REPORT=1, factr=10^8))
SE       = sqrt(diag(solve(est$hessian)))
round(cbind(estimate=est$par, SE), 3) # results
#################################################

# Residual Analysis (not shown)
phi1   = est$par[1]; phi2 = est$par[2]
sR     = est$par[3]; b1   = est$par[4]
b2     = est$par[5]; b3   = est$par[6]
b4     = est$par[7]; alf  = est$par[8]
mu0    = matrix(c(0,0), 2, 1)
Sigma0 = diag(100, 2)
Phi    = matrix(c(phi1, phi2, 1, 0), 2)
S      = 1
Ups    = matrix(c(b1, 0, b2, 0, b3, 0, 0, 0, 0, 0), 2, 5)
Gam    = matrix(c(0, 0, 0, b4, alf), 1, 5) 
sQ     = matrix(c(phi1, phi2), 2)*sR
kf     = Kfilter(y, A, mu0, Sigma0, Phi, sQ, sR, Ups=Ups, Gam=Gam, input=input, S=S, version=2)
res    = ts(drop(kf$innov), start=start(cmort), freq=frequency(cmort))
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
# NOTE:  If this takes a long time on your machine, 
#        increase `tol` and/or decrease `nboot`
tol   = .0001    # determines convergence of optimizer     
nboot = 500      # number of bootstrap replicates     
################################## 

y     = window(qinfl, c(1953,1), c(1965,2))  # inflation   
z     = window(qintr, c(1953,1), c(1965,2))  # interest   
num   = length(y) 
A     = array(z, dim=c(1,1,num))
input = matrix(1,num,1)  

# Function to Calculate Likelihood   
Linn  = function(para, y.data){  # pass data also
   phi = para[1];  alpha = para[2]
   b   = para[3];  Ups   = (1-phi)*b
   sQ  = para[4];  sR    = para[5]  
   kf  = Kfilter(y.data,A,mu0,Sigma0,phi,sQ ,sR ,Ups ,Gam=alpha,input)
   return(kf$like)    
}

# Parameter Estimation   
mu0      = 1
Sigma0   = .01  
init.par = c(phi=.84, alpha=-.77, b=.85, sQ=.12, sR=1.1)  # initial values   

est = optim(init.par,  Linn, NULL, y.data=y, method="BFGS", hessian=TRUE, 
             control=list(trace=1, REPORT=1, reltol=tol))  
SE  = sqrt(diag(solve(est$hessian)))   

phi   = est$par[1];  alpha = est$par[2]
b     = est$par[3];  Ups   = (1-phi)*b         
sQ    = est$par[4];  sR    = est$par[5] 
round(cbind(estimate=est$par, SE), 3)  


# BEGIN BOOTSTRAP   
# Run the filter at the estimates 
kf  = Kfilter(y, A, mu0, Sigma0, phi, sQ, sR, Ups, Gam=alpha, input) 

# Pull out necessary values from the filter and initialize  
xp      = kf$Xp
Pp      = kf$Pp
innov   = kf$innov 
sig     = kf$sig 
e       = innov/sqrt(sig)
e.star  = e                      # initialize values
y.star  = y  
xp.star = xp  
k       = 4:50                   # hold first 3 observations fixed 
para.star = matrix(0, nboot, 5)  # to store estimates
init.par  =  c(.84, -.77, .85, .12, 1.1)    

pb = txtProgressBar(min = 0, max = nboot, initial = 0, style=3)  # progress bar

for (i in 1:nboot){
 setTxtProgressBar(pb,i)                       
 e.star[k] = sample(e[k], replace=TRUE)   
 for (j in k){ 
   K  = (phi*Pp[j]*z[j])/sig[j]  
  xp.star[j] = phi*xp.star[j-1] + Ups +   K*sqrt(sig[j])*e.star[j]
  } 
   y.star[k] = z[k]*xp.star[k] + alpha + sqrt(sig[k])*e.star[k]  
 est.star  = optim(init.par, Linn, NULL, y.data=y.star, method='BFGS', control=list(reltol=tol))     
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
  sQ   = diag(c(sigw,0))
  kf   = Kfilter(y, A, mu0, Sigma0, Phi, sQ, sigv)
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
sQ    = diag(c(sigw,0))
sigv  = est$par[2]
ks    = Ksmooth(y, A, mu0, Sigma0, Phi, sQ, sigv)
xsmoo = ts(ks$Xs[1,1,])
psmoo = ts(ks$Ps[1,1,])
upp   = xsmoo + 2*sqrt(psmoo)
low   = xsmoo - 2*sqrt(psmoo)


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



<a name="6.23"></a>
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

```r
# generate some data from the model - 2 parameters
set.seed(1)
sQ = 1; sR = 3; n  = 100
mu0 = 0; Sigma0=10; x0=rnorm(1,mu0,Sigma0)
w  = rnorm(n); v  = rnorm(n)
x = c(x0   + sQ*w[1])  # initialize states
y = c(x[1] + sR*v[1])  # initialize obs
for (t in 2:n){ 
  x[t] = x[t-1] + sQ*w[t]
  y[t] = x[t] + sR*v[t]
  }

# set up the Gibbs sampler
burn   = 50;  n.iter = 1000
niter  = burn + n.iter
draws  = c()
# priors for R (a,b) and Q (c,d) IG distributions
a = 2; b = 2; c = 2; d = 1  
# (1) initialize - sample sQ and sR  
sR = sqrt(1/rgamma(1,a,b)); sQ =  sqrt(1/rgamma(1,c,d))

# progress bar
pb = txtProgressBar(min = 0, max = niter, initial = 0, style=3)  

# run it
for (iter in 1:niter){
## (2)  sample the states  
 run   = ffbs(y,1,0,10,1,sQ,sR)  # ffbs(y,A,mu0,Sigma0,Phi,Ups,Gam,sQ,sR,input)
## (1)  sample the parameters    
  Xs   = as.matrix(run$Xs)
  R    = 1/rgamma(1,a+n/2,b+sum((y-Xs)^2)/2)
 sR    = sqrt(R)
  Q    = 1/rgamma(1,c+(n-1)/2,d+sum(diff(Xs)^2)/2)
 sQ    = sqrt(Q)
## store everything 
 draws = rbind(draws,c(sQ,sR,Xs))
 setTxtProgressBar(pb,iter)  
}
close(pb)

# pull out the results for easy plotting
draws  = draws[(burn+1):(niter),]
 q025  = function(x){quantile(x,0.025)}
 q975  = function(x){quantile(x,0.975)}
xs     = draws[,3:(n+2)]
lx     = apply(xs,2,q025)
mx     = apply(xs,2,mean)
ux     = apply(xs,2,q975)

# some graphics 
tsplot(cbind(x,y,mx), spag=TRUE,  ylab='', col=c(6,8,4), lwd=c(1,1,1.5), type='o', pch=c(NA,1,NA))
legend('topleft', legend=c("x(t)","y(t)","xs(t)"), lty=1, col=c(6,8,4), lwd=1.5, bty="n", pch=c(NA,1,NA))
points(y)
 xx=c(1:100, 100:1)
 yy=c(lx, rev(ux))
polygon(xx, yy, border=NA, col=astsa.col(4,.1))
```



Example 6.27

```r
y = jj    # the data

### setup - model and initial parameters
set.seed(90210)
n   = length(y)
A   = matrix(c(1,1,0,0), 1, 4)
Phi = diag(0,4)
  Phi[1,1] = 1.03 
  Phi[2,]  = c(0,-1,-1,-1); Phi[3,]=c(0,1,0,0); Phi[4,]=c(0,0,1,0)
mu0 = rbind(.7,0,0,0)
Sigma0 = diag(.04, 4)
sR = 1                    # observation noise standard deviation
sQ = diag(c(.1,.1,0,0))   # state noise standard deviations on the diagonal

### initializing and hyperparameters
burn   = 50
n.iter = 1000
niter  = burn + n.iter
draws  = NULL
a = 2; b = 2; c = 2; d = 1   # hypers (c and d for both Qs)

pb = txtProgressBar(min = 0, max = niter, initial = 0, style=3)  # progress bar

### start Gibbs
for (iter in 1:niter){
# draw states 
  run  = ffbs(y,A,mu0,Sigma0,Phi,sQ,sR)   # initial values are given above
  xs   = run$Xs
# obs variance
  R    = 1/rgamma(1,a+n/2,b+sum((as.vector(y)-as.vector(A%*%xs[,,]))^2))
 sR    = sqrt(R)
# beta where phi = 1+beta  
 Y     = diff(xs[1,,])
 D     = as.vector(lag(xs[1,,],-1))[-1]
 regu  = lm(Y~0+D)  # est beta = phi-1
 phies = as.vector(coef(summary(regu)))[1:2] + c(1,0) # phi estimate and SE
 dft   = df.residual(regu)
 Phi[1,1]  = phies[1] + rt(1,dft)*phies[2]  # use a t to sample phi
# state variances
  u   = xs[,,2:n] - Phi%*%xs[,,1:(n-1)]
  uu  = u%*%t(u)/(n-2)
  Q1  = 1/rgamma(1,c+(n-1)/2,d+uu[1,1]/2)
  sQ1 = sqrt(Q1)
  Q2  = 1/rgamma(1,c+(n-1)/2,d+uu[2,2]/2)
  sQ2 = sqrt(Q2) 
  sQ  = diag(c(sQ1, sQ2, 0,0))
# store results
 trend = xs[1,,]
 season= xs[2,,] 
 draws = rbind(draws,c(Phi[1,1],sQ1,sQ2,sR,trend,season))
 setTxtProgressBar(pb,iter)  
}
close(pb)

##- display results -##

# set up
u     = draws[(burn+1):(niter),]
parms = u[,1:4]
q025  = function(x){quantile(x,0.025)}
q975  = function(x){quantile(x,0.975)}

##  plot parameters (display at end)
names= c(expression(phi), expression(sigma[w1]), expression(sigma[w2]), expression(sigma[v]))
par(mfrow=c(4,1), mar=c(2,1,2,1)+1)
for (i in 1:4){
hist(parms[,i], col=astsa.col(5,.4), main=names[i], xlab='', cex.main=2)
 u1 = apply(parms,2,q025); u2 = apply(parms,2,mean); u3 = apply(parms,2,q975);
abline(v=c(u1[i],u2[i],u3[i]), lwd=2, col=c(3,6,3))
}

###  plot states  (display at end)
# trend
dev.new()
par(mfrow=2:1)
tr    = ts(u[,5:(n+4)], start=1960, frequency=4)
ltr   = ts(apply(tr,2,q025), start=1960, frequency=4)
mtr   = ts(apply(tr,2,mean), start=1960, frequency=4)
utr   = ts(apply(tr,2,q975), start=1960, frequency=4)

tsplot(mtr, ylab='', col=4, main='trend')
 xx=c(time(mtr), rev(time(mtr)))
 yy=c(ltr, rev(utr))
polygon(xx, yy, border=NA, col=astsa.col(4,.1)) 

# trend + season
sea    = ts(u[,(n+5):(2*n)], start=1960, frequency=4)
lsea   = ts(apply(sea,2,q025), start=1960, frequency=4)
msea   = ts(apply(sea,2,mean), start=1960, frequency=4)
usea   = ts(apply(sea,2,q975), start=1960, frequency=4)
tsplot(msea+mtr, ylab='', col=4, main='trend + season')
 xx=c(time(msea), rev(time(msea)))
 yy=c(lsea+ltr, rev(usea+utr))
polygon(xx, yy, border=NA, col=astsa.col(4,.1)) 
```


[<sub>top</sub>](#table-of-contents)

---

## Chapter 7

Code in Introduction

```r
x = matrix(0, 128, 6)
for (i in 1:6) x[,i] = rowMeans(fmri[[i]])
colnames(x)=c("Brush", "Heat", "Shock", "Brush", "Heat", "Shock")
tsplot(x, ncol=2, byrow=FALSE, main='')
mtext('Awake', outer=T, adj=.25)
mtext('Sedated', outer=T, adj=.75)
```

```r
attach(eqexp)  # to use the names

P = 1:1024; S = P+1024
x = cbind(EQ5[P], EQ6[P], EX5[P], EX6[P], NZ[P], EQ5[S], EQ6[S], 
          EX5[S], EX6[S], NZ[S])
x.name = c("EQ5","EQ6","EX5","EX6","NZ")
colnames(x) = c(x.name, x.name)
tsplot(x, ncol=2, byrow=FALSE, main='')
mtext('P waves', outer=T, adj=.25)
mtext('S waves', outer=T, adj=.75)

detach(eqexp)  # Redemption 
```


Example 7.1
```r
tsplot(climhyd, ncol=2, col=2:7, lwd=2)    # figure 7.3
Y = climhyd         # Y to hold transformed series
Y[,6] = log(Y[,6])  # log inflow
Y[,5] = sqrt(Y[,5]) # sqrt precipitation 
                              
L = 25              # setup 
M = 100 
alpha = .001
fdr = .001
nq = 2              # number of inputs  (Temp and Precip)

# Spectral Matrix 
Yspec = mvspec(Y, spans=L, kernel="daniell", taper=.1, plot=FALSE)
n = Yspec$n.used          # effective sample size
Fr = Yspec$freq           # fundamental freqs 
n.freq = length(Fr)       # number of frequencies
Yspec$bandwidth           # = 0.05 
 
# Coherencies (see section 4.7 also)
Fq = qf(1-alpha, 2, L-2); cn = Fq/(L-1+Fq)
plt.name = c("(a)","(b)","(c)","(d)","(e)","(f)")
dev.new()
par(mfrow=c(2,3), cex.lab=1.2) 
# The coherencies are listed as 1,2,...,15=choose(6,2) 
for (i in 11:15){
 tsplot(Fr,Yspec$coh[,i], ylab="Sq Coherence", xlab="Frequency", ylim=c(0,1),
      main=paste("Inflow with", names(climhyd[i-10]), sep=' '))
 abline(h = cn); text(.45,.98, plt.name[i-10], cex=1.2)  
} 

 # Multiple Coherency 
coh.15 = stoch.reg(Y, cols.full = c(1,5), cols.red = NULL, alpha, L, M, plot.which = "coh", 
                     main="Inflow with Temp & Precip")  
text(.45,.98, plt.name[6], cex=1.2) 

# Partial F (note F-stat is called eF in the code)
numer.df = 2*nq
denom.df = Yspec$df-2*nq

dev.new()
par(mfrow=c(3,1), mar=c(3,3,2,1)+.5, mgp = c(1.5,0.4,0), cex.lab=1.2)  
out.15 = stoch.reg(Y, cols.full = c(1,5), cols.red = 5, alpha, L, M, plot.which = "F.stat")
 eF = out.15$eF 
 pvals = pf(eF, numer.df, denom.df, lower.tail = FALSE)
 pID = FDR(pvals, fdr)
abline(h=c(eF[pID]), lty=2)
title(main = "Partial F Statistic")

# Regression Coefficients
S = seq(from = -M/2+1, to = M/2 - 1, length = M-1)

tsplot(S, coh.15$Betahat[,1], type = "h", xlab = "", ylab =names(climhyd[1]), 
        ylim = c(-.025, .055), lwd=2)
abline(h=0)
title(main = "Impulse Response Functions")

tsplot(S, coh.15$Betahat[,2], type = "h", xlab = "Index", ylab = names(climhyd[5]), 
        ylim = c(-.015, .055), lwd=2)
abline(h=0)
```


Example 7.2
```r
attach(beamd)

tau    = rep(0,3) 
     u = ccf(sensor1, sensor2, plot=FALSE)
tau[1] = u$lag[which.max(u$acf)]    #  17
     u = ccf(sensor3, sensor2, plot=FALSE)
tau[3] = u$lag[which.max(u$acf)]    # -22

Y = ts.union(sensor1=lag(sensor1,tau[1]), lag(sensor2, tau[2]), lag(sensor3, tau[3]))
Y = ts.union(Y, rowMeans(Y))
colnames(Y) = c('sensor1', 'sensor2', 'sensor3', 'beam')
tsplot(Y, main="Infrasonic Signals and Beam")


detach(beamd)
```



Example 7.4
```r
attach(beamd)

L = 9
fdr = .001 
N = 3
Y = cbind(beamd, beam=rowMeans(beamd))
n = nextn(nrow(Y))

Y.fft = mvfft(as.ts(Y))/sqrt(n)
Df = Y.fft[,1:3]   # fft of the data
Bf = Y.fft[,4]     # beam fft 

ssr = N*Re(Bf*Conj(Bf))               # raw signal spectrum
sse = Re(rowSums(Df*Conj(Df))) - ssr  # raw error spectrum

# Smooth
SSE = filter(sse, sides=2, filter=rep(1/L,L), circular=TRUE)
SSR = filter(ssr, sides=2, filter=rep(1/L,L), circular=TRUE)
SST = SSE + SSR

par(mfrow=c(2,1))
Fr  = 0:(n-1)/n
nFr = 1:200   # freqs to plot

tsplot( Fr[nFr], SST[nFr], type="l", ylab="log Power", xlab="", main="Sum of Squares",log="y")
lines(Fr[nFr], SSE[nFr], type="l", lty=2)

eF = (N-1)*SSR/SSE; df1 = 2*L;  df2 = 2*L*(N-1)
pvals = pf(eF, df1, df2, lower=FALSE)  # p values for FDR
pID = FDR(pvals, fdr); Fq = qf(1-fdr, df1, df2)  

tsplot(Fr[nFr], eF[nFr], type="l", ylab="F-statistic", xlab="Frequency",  main="F Statistic")
abline(h=c(Fq, eF[pID]), lty=1:2)


detach(beamd)
```


Example 7.5
```r
attach(beamd)

L  = 9
M  = 100 
M2 = M/2
N  = 3   
Y  = cbind(beamd, beam <- rowMeans(beamd))
n  = nextn(nrow(Y))
n.freq = n/2
 
Y[,1:3] = Y[,1:3]-Y[,4]  # center each series  

Y.fft = mvfft(as.ts(Y))/sqrt(n)
Ef    = Y.fft[,1:3]              # fft of the error
Bf    = Y.fft[,4]                # beam fft
ssr   = N*Re(Bf*Conj(Bf))        # Raw Signal Spectrum 
sse   = Re(rowSums(Ef*Conj(Ef))) # Raw Error Spectrum

# Smooth
SSE = filter(sse, sides=2, filter=rep(1/L,L), circular=TRUE)
SSR = filter(ssr, sides=2, filter=rep(1/L,L), circular=TRUE)

# Estimate Signal and Noise Spectra
fv = SSE/(L*(N-1))          # Equation (7.77)
fb = (SSR-SSE/(N-1))/(L*N)  # Equation (7.78)
fb[fb<0] = 0 

H0 = N*fb/(fv+N*fb)     
H0[ceiling(.04*n):n] = 0    # zero out H0 beyond frequency .04

# Extend components to make it a valid transform
H0 = c(H0[1:n.freq], rev(H0[2:(n.freq+1)]))  
h0 = Re(fft(H0, inverse = TRUE))            # Impulse Response 
h0 = c(rev(h0[2:(M2+1)]), h0[1:(M2+1)])     # center it
h1 = spec.taper(h0, p = .5)                 # taper it
k1 = h1/sum(h1)                             # normalize it
f.beam = filter(Y$beam, filter=k1, sides=2) # filter it

# Graphics  
nFr = 1:50      # freqs to display
Fr = (nFr-1)/n  # frequencies 

layout(matrix(c(1, 2, 4, 1, 3, 4), nc=2))
tsplot(10*Fr, fb[nFr], type="l", ylab="Power", xlab="Frequency (Hz)")
 lines(10*Fr, fv[nFr], lty=2); text(.24, 5, "(a)", cex=1.2)
tsplot(10*Fr, H0[nFr], type="l", ylab="Frequency Response", xlab="Frequency(Hz)")
 text(.23, .84, "(b)", cex=1.2)
tsplot(-M2:M2, k1, type="l", ylab="Impulse Response", xlab="Index", lwd=1.5)
 text(45, .022, "(c)", cex=1.2)
tsplot(cbind(f.beam,beam), spag=TRUE, lty=1:2, ylab="beam")
 text(2040, 2, "(d)", cex=1.2)

detach(beamd) 
```

Example 7.6
```r
n          = 128               # length of series
n.freq     = 1 + n/2           # number of frequencies
Fr         = (0:(n.freq-1))/n  # the frequencies 
N          = c(5,4,5,3,5,4)    # number of series for each cell
n.subject  = sum(N)            # number of subjects (26)
n.trt      = 6                 # number of treatments
L          = 3                 # for smoothing
num.df     = 2*L*(n.trt-1)     # dfs for F test
den.df     = 2*L*(n.subject-n.trt)


# Design Matrix (Z): 
Z1 = outer(rep(1,N[1]), c(1,1,0,0,0,0))
Z2 = outer(rep(1,N[2]), c(1,0,1,0,0,0))
Z3 = outer(rep(1,N[3]), c(1,0,0,1,0,0))
Z4 = outer(rep(1,N[4]), c(1,0,0,0,1,0)) 
Z5 = outer(rep(1,N[5]), c(1,0,0,0,0,1)) 
Z6 = outer(rep(1,N[6]), c(1,-1,-1,-1,-1,-1)) 

Z  = rbind(Z1, Z2, Z3, Z4, Z5, Z6)
ZZ = t(Z)%*%Z 

SSEF <- rep(NA, n) -> SSER   

HatF = Z%*%solve(ZZ, t(Z))
HatR = Z[,1]%*%t(Z[,1])/ZZ[1,1]

par(mfrow=c(3,3))
loc.name = c("Cortex 1","Cortex 2","Cortex 3","Cortex 4","Caudate","Thalamus 1",
              "Thalamus 2", "Cerebellum 1","Cerebellum 2")

for(Loc in 1:9) {   
 i = n.trt*(Loc-1)   
 Y = cbind(fmri[[i+1]], fmri[[i+2]], fmri[[i+3]], fmri[[i+4]], fmri[[i+5]], fmri[[i+6]])
 Y = mvfft(spec.taper(Y, p=.5))/sqrt(n)	
 Y = t(Y)      # Y is now 26 x 128 FFTs

 # Calculation of Error Spectra 
 for (k in 1:n) {   
  SSY = Re(Conj(t(Y[,k]))%*%Y[,k])
  SSReg = Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
 SSEF[k] = SSY - SSReg
  SSReg = Re(Conj(t(Y[,k]))%*%HatR%*%Y[,k])
 SSER[k] = SSY - SSReg  
 }

# Smooth 
sSSEF = filter(SSEF, rep(1/L, L), circular = TRUE)
sSSER = filter(SSER, rep(1/L, L), circular = TRUE)

eF =(den.df/num.df)*(sSSER-sSSEF)/sSSEF

tsplot(Fr, eF[1:n.freq], xlab="Frequency", ylab="F Statistic", ylim=c(0,7), main=loc.name[Loc])
abline(h=qf(.999, num.df, den.df),lty=2) 
}
```

Example 7.7
```r
n      = 128
n.freq = 1 + n/2  
Fr     = (0:(n.freq-1))/n  
nFr    = 1:(n.freq/2) 
N      = c(5,4,5,3,5,4)  
n.para = 6          # number of parameters
n.subject = sum(N)  # total number of subjects

L = 3   
df.stm = 2*L*(3-1)              # stimulus (3 levels: Brush,Heat,Shock)
df.con = 2*L*(2-1)              # conscious (2 levels: Awake,Sedated)  
df.int = 2*L*(3-1)*(2-1)        # interaction 
den.df = 2*L*(n.subject-n.para) # df for full model 

# Design Matrix:          mu  a1  a2   b  g1  g2   
 Z1 = outer(rep(1,N[1]), c(1,  1,  0,  1,  1,  0)) 
 Z2 = outer(rep(1,N[2]), c(1,  0,  1,  1,  0,  1)) 
 Z3 = outer(rep(1,N[3]), c(1, -1, -1,  1, -1, -1)) 
 Z4 = outer(rep(1,N[4]), c(1,  1,  0, -1, -1,  0)) 
 Z5 = outer(rep(1,N[5]), c(1,  0,  1, -1,  0, -1)) 
 Z6 = outer(rep(1,N[6]), c(1, -1, -1, -1,  1,  1))
 
 Z = rbind(Z1, Z2, Z3, Z4, Z5, Z6)  
 ZZ = t(Z)%*%Z

rep(NA, n)-> SSEF -> SSE.stm -> SSE.con -> SSE.int              
HatF    = Z%*%solve(ZZ,t(Z))     
Hat.stm = Z[,-(2:3)]%*%solve(ZZ[-(2:3),-(2:3)], t(Z[,-(2:3)]))  
Hat.con = Z[,-4]%*%solve(ZZ[-4,-4], t(Z[,-4]))   
Hat.int = Z[,-(5:6)]%*%solve(ZZ[-(5:6),-(5:6)], t(Z[,-(5:6)]))                                                           

par(mfrow=c(5,3), oma=c(0,2,0,0))
loc.name = c("Cortex 1","Cortex 2","Cortex 3","Cortex 4","Caudate", "Thalamus 1",
              "Thalamus 2", "Cerebellum 1","Cerebellum 2")
for(Loc in c(1:4,9)) {   # only Loc 1 to 4 and 9 used                
  i = 6*(Loc-1)                                                                 
  Y = cbind(fmri[[i+1]], fmri[[i+2]], fmri[[i+3]], fmri[[i+4]], fmri[[i+5]], fmri[[i+6]])                     
  Y = mvfft(spec.taper(Y, p=.5))/sqrt(n)  
  Y = t(Y)  
  for (k in 1:n) {    
    SSY=Re(Conj(t(Y[,k]))%*%Y[,k])
    SSReg= Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
  SSEF[k]=SSY-SSReg 
    SSReg=Re(Conj(t(Y[,k]))%*%Hat.stm%*%Y[,k])
  SSE.stm[k] = SSY-SSReg
    SSReg=Re(Conj(t(Y[,k]))%*%Hat.con%*%Y[,k])
  SSE.con[k]=SSY-SSReg  
    SSReg=Re(Conj(t(Y[,k]))%*%Hat.int%*%Y[,k])
  SSE.int[k]=SSY-SSReg   
  }
 # Smooth    
  sSSEF    = filter(SSEF, rep(1/L, L), circular = TRUE)
  sSSE.stm = filter(SSE.stm, rep(1/L, L), circular = TRUE)
  sSSE.con = filter(SSE.con, rep(1/L, L), circular = TRUE)
  sSSE.int = filter(SSE.int, rep(1/L, L), circular = TRUE)  
  eF.stm   = (den.df/df.stm)*(sSSE.stm-sSSEF)/sSSEF
  eF.con   = (den.df/df.con)*(sSSE.con-sSSEF)/sSSEF
  eF.int   = (den.df/df.int)*(sSSE.int-sSSEF)/sSSEF

tsplot(Fr[nFr],eF.stm[nFr],xlab="Frequency", ylab="F Statistic", ylim=c(0,12))
   abline(h=qf(.999, df.stm, den.df),lty=2)       
  if(Loc==1) mtext("Stimulus", side=3, line=.3, cex=.8)
  mtext(loc.name[Loc], side=2, line=3, cex=.8)
 tsplot(Fr[nFr],eF.con[nFr], xlab="Frequency", ylab="F Statistic", ylim=c(0,12))
  abline(h=qf(.999, df.con, den.df),lty=2)
  if(Loc==1)  mtext("Consciousness", side=3, line=.3, cex=.8)   
 tsplot(Fr[nFr],eF.int[nFr], xlab="Frequency",ylab="F Statistic", ylim=c(0,12))
  abline(h=qf(.999, df.int, den.df),lty=2)
  if(Loc==1) mtext("Interaction", side=3, line= .3, cex=.8)   
} 
```



Example 7.8
```r
n      = 128
n.freq = 1 + n/2  
Fr     = (0:(n.freq-1))/n
nFr    = 1:(n.freq/2) 
N      = c(5,4,5,3,5,4)
L      = 3
n.subject = sum(N)
   
# Design Matrix            
Z1 = outer(rep(1,N[1]), c(1,0,0,0,0,0)) 
Z2 = outer(rep(1,N[2]), c(0,1,0,0,0,0)) 
Z3 = outer(rep(1,N[3]), c(0,0,1,0,0,0)) 
Z4 = outer(rep(1,N[4]), c(0,0,0,1,0,0)) 
Z5 = outer(rep(1,N[5]), c(0,0,0,0,1,0)) 
Z6 = outer(rep(1,N[6]), c(0,0,0,0,0,1))
Z  = rbind(Z1, Z2, Z3, Z4, Z5, Z6)
ZZ = t(Z)%*%Z 

A      = rbind(diag(1,3), diag(1,3))   # Contrasts:  6 x 3 
nq     = nrow(A)
num.df = 2*L*nq
den.df = 2*L*(n.subject-nq)                   
HatF   = Z%*%solve(ZZ, t(Z))           # full model hat matrix   

rep(NA, n)-> SSEF -> SSER 
eF = matrix(0,n,3)

par(mfrow=c(5,3), oma=c(0,2,0,0)) 

loc.name = c("Cortex 1","Cortex 2","Cortex 3","Cortex 4","Caudate","Thalamus 1","Thalamus 2",
             "Cerebellum 1","Cerebellum 2")
cond.name = c("Brush", "Heat", "Shock")             

for(Loc in c(1:4,9)) {
 i = 6*(Loc-1)                                                                 
 Y = cbind(fmri[[i+1]],fmri[[i+2]],fmri[[i+3]],fmri[[i+4]], fmri[[i+5]],fmri[[i+6]])                       
 Y = mvfft(spec.taper(Y, p=.5))/sqrt(n); Y = t(Y)  
 for (cond in 1:3){
  Q = t(A[,cond])%*%solve(ZZ, A[,cond])
  HR = A[,cond]%*%solve(ZZ, t(Z))      
  for (k in 1:n){   
    SSY = Re(Conj(t(Y[,k]))%*%Y[,k])
    SSReg= Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
   SSEF[k]= (SSY-SSReg)*Q 
    SSReg= HR%*%Y[,k]
   SSER[k] = Re(SSReg*Conj(SSReg))  
  }     

# Smooth      
 sSSEF = filter(SSEF, rep(1/L, L), circular = TRUE)
 sSSER = filter(SSER, rep(1/L, L), circular = TRUE)       
 eF[,cond]= (den.df/num.df)*(sSSER/sSSEF)   }                                    
 tsplot(Fr[nFr], eF[nFr,1], xlab="Frequency", ylab="F Statistic", ylim=c(0,5), main='')
  abline(h=qf(.999, num.df, den.df),lty=2)       
  if(Loc==1) mtext("Brush", side=3, line=.3, cex=1)
  mtext(loc.name[Loc], side=2, line=3, cex=.9)
 tsplot(Fr[nFr], eF[nFr,2], xlab="Frequency", ylab="F Statistic", ylim=c(0,5), main='')
  abline(h=qf(.999, num.df, den.df),lty=2)
  if(Loc==1)  mtext("Heat", side=3, line=.3, cex=1)   
 tsplot(Fr[nFr], eF[nFr,3], xlab="Frequency", ylab="F Statistic", ylim=c(0,5), main='')
  abline(h = qf(.999, num.df, den.df) ,lty=2)
  if(Loc==1) mtext("Shock", side=3, line=.3, cex=1)  
}  
```


Example 7.9
```r
P = 1:1024 
S = P+1024
N = 8
n = 1024
p.dim = 2
m = 10
L = 2*m+1

eq.P   = as.ts(eqexp[P,1:8])
eq.S   = as.ts(eqexp[S,1:8]) 
eq.m   = cbind(rowMeans(eq.P), rowMeans(eq.S))
ex.P   = as.ts(eqexp[P,9:16])
ex.S   = as.ts(eqexp[S,9:16]) 
ex.m   = cbind(rowMeans(ex.P), rowMeans(ex.S)) 
m.diff = mvfft(eq.m - ex.m)/sqrt(n)

eq.Pf = mvfft(eq.P-eq.m[,1])/sqrt(n)
eq.Sf = mvfft(eq.S-eq.m[,2])/sqrt(n)
ex.Pf = mvfft(ex.P-ex.m[,1])/sqrt(n)
ex.Sf = mvfft(ex.S-ex.m[,2])/sqrt(n) 

fv11 = rowSums(eq.Pf*Conj(eq.Pf)) + rowSums(ex.Pf*Conj(ex.Pf))/(2*(N-1))
fv12 = rowSums(eq.Pf*Conj(eq.Sf)) + rowSums(ex.Pf*Conj(ex.Sf))/(2*(N-1))
fv22 = rowSums(eq.Sf*Conj(eq.Sf)) + rowSums(ex.Sf*Conj(ex.Sf))/(2*(N-1))
fv21 = Conj(fv12)

# Equal Means  
T2 = rep(NA, 512)
for (k  in 1:512){
 fvk = matrix(c(fv11[k], fv21[k], fv12[k], fv22[k]), 2, 2)
 dk = as.matrix(m.diff[k,])
 T2[k] = Re((N/2)*Conj(t(dk))%*%solve(fvk,dk))  }
eF = T2*(2*p.dim*(N-1))/(2*N-p.dim-1)

par(mfrow=c(2,2))

freq = 40*(0:511)/n  # in Hz (cycles per second)  
tsplot(freq, eF, xlab="Frequency (Hz)", ylab="F Statistic", main="Equal Means")
abline(h=qf(.999, 2*p.dim, 2*(2*N-p.dim-1)))

# Equal P   
kd    = kernel("daniell",m);  
u     = Re(rowSums(eq.Pf*Conj(eq.Pf))/(N-1))
feq.P = kernapply(u, kd, circular=TRUE)
u     = Re(rowSums(ex.Pf*Conj(ex.Pf))/(N-1))
fex.P =	kernapply(u, kd, circular=TRUE)

tsplot(freq, feq.P[1:512]/fex.P[1:512], xlab="Frequency (Hz)", ylab="F Statistic", 
      main="Equal P-Spectra")
abline(h = qf(.999, 2*L*(N-1),  2*L*(N-1))) 

# Equal S   
u     = Re(rowSums(eq.Sf*Conj(eq.Sf))/(N-1))
feq.S = kernapply(u, kd, circular=TRUE)
u     = Re(rowSums(ex.Sf*Conj(ex.Sf))/(N-1))
fex.S =	kernapply(u, kd, circular=TRUE)

tsplot(freq, feq.S[1:512]/fex.S[1:512], xlab="Frequency (Hz)", ylab="F Statistic", 
      main="Equal S-Spectra")
abline(h=qf(.999, 2*L*(N-1),  2*L*(N-1))) 

# Equal Spectra  
u      = rowSums(eq.Pf*Conj(eq.Sf))/(N-1)
feq.PS = kernapply(u, kd, circular=TRUE)
u      = rowSums(ex.Pf*Conj(ex.Sf)/(N-1))
fex.PS = kernapply(u, kd, circular=TRUE)
fv11   = kernapply(fv11, kd, circular=TRUE)
fv22   = kernapply(fv22, kd, circular=TRUE)
fv12  = kernapply(fv12, kd, circular=TRUE)

Mi = L*(N-1) 
M  = 2*Mi
TS = rep(NA,512)

for (k  in 1:512){
 det.feq.k = Re(feq.P[k]*feq.S[k] - feq.PS[k]*Conj(feq.PS[k]))
 det.fex.k = Re(fex.P[k]*fex.S[k] - fex.PS[k]*Conj(fex.PS[k]))
 det.fv.k  = Re(fv11[k]*fv22[k] - fv12[k]*Conj(fv12[k]))
 
 log.n1 = log(M)*(M*p.dim)
 log.d1 = log(Mi)*(2*Mi*p.dim)
 log.n2 = log(Mi)*2 +log(det.feq.k)*Mi + log(det.fex.k)*Mi  
 log.d2 = (log(M)+log(det.fv.k))*M
 r = 1 - ((p.dim+1)*(p.dim-1)/6*p.dim*(2-1))*(2/Mi - 1/M)
 TS[k] = -2*r*(log.n1+log.n2-log.d1-log.d2)   
}

tsplot(freq, TS, xlab="Frequency (Hz)", ylab="Chi-Sq Statistic", main="Equal Spectral Matrices")
abline(h = qchisq(.9999, p.dim^2))   # about 23.5, so not on the plot
```




Example 7.10
```r
P     = 1:1024
S     = P+1024
mag.P = log10(apply(eqexp[P,],2,max) - apply(eqexp[P,],2,min))
mag.S = log10(apply(eqexp[S,],2,max) - apply(eqexp[S,],2,min))
eq.P  = mag.P[1:8]
eq.S  = mag.S[1:8]
ex.P  = mag.P[9:16]
ex.S =  mag.S[9:16]
NZ.P  = mag.P[17]
NZ.S  = mag.S[17] 

# Compute linear discriminant function 
cov.eq     = var(cbind(eq.P, eq.S))
cov.ex     = var(cbind(ex.P, ex.S))
cov.pooled = (cov.ex + cov.eq)/2

means.eq   =  colMeans(cbind(eq.P,eq.S)); 
means.ex   =  colMeans(cbind(ex.P,ex.S))
slopes.eq  = solve(cov.pooled, means.eq)
inter.eq   = -sum(slopes.eq*means.eq)/2
slopes.ex  = solve(cov.pooled, means.ex)
inter.ex   = -sum(slopes.ex*means.ex)/2
d.slopes   = slopes.eq - slopes.ex
d.inter    = inter.eq - inter.ex

# Classify new observation 
new.data = cbind(NZ.P, NZ.S)

d = sum(d.slopes*new.data) + d.inter
post.eq = exp(d)/(1+exp(d))

# Print (disc function, posteriors) and plot results  
cat(d.slopes[1], "mag.P +" , d.slopes[2], "mag.S +" , d.inter,"\n")  
cat("P(EQ|data) =", post.eq,  "  P(EX|data) =", 1-post.eq, "\n" )    

tsplot(eq.P, eq.S, type='p', xlim=c(0,1.5), ylim=c(.75,1.25), 
        xlab="log mag(P)", ylab ="log mag(S)", pch = 8, cex=1.1, lwd=2, 
        main="Classification Based on Magnitude Features", col=4)
 points(ex.P, ex.S, pch = 6, cex=1.1, lwd=2, col=6)
 points(new.data, pch = 3, cex=1.1, lwd=2, col=3)
 abline(a = -d.inter/d.slopes[2], b = -d.slopes[1]/d.slopes[2])
 text(eq.P-.07,eq.S+.005, label=names(eqexp[1:8]), cex=.8)
 text(ex.P+.07,ex.S+.003, label=names(eqexp[9:16]), cex=.8)
 text(NZ.P+.05,NZ.S+.003, label=names(eqexp[17]), cex=.8)
 legend("topright",c("EQ","EX","NZ"),pch=c(8,6,3),pt.lwd=2,cex=1.1, col=c(4,6,3))

# Cross-validation 
all.data = rbind(cbind(eq.P,eq.S), cbind(ex.P,ex.S))
post.eq <- rep(NA, 8) -> post.ex 

for(j in 1:16) {
 if (j <= 8) {samp.eq = all.data[-c(j, 9:16),];  samp.ex = all.data[9:16,]}
 if (j > 8)  {samp.eq = all.data[1:8,];    samp.ex = all.data[-c(j, 1:8),]}
 
 df.eq      = nrow(samp.eq)-1;  df.ex = nrow(samp.ex)-1
 mean.eq    = colMeans(samp.eq);  mean.ex = colMeans(samp.ex)
 cov.eq     = var(samp.eq);  cov.ex = var(samp.ex)
 cov.pooled = (df.eq*cov.eq + df.ex*cov.ex)/(df.eq + df.ex)
 slopes.eq  = solve(cov.pooled, mean.eq)
 inter.eq   = -sum(slopes.eq*mean.eq)/2
 slopes.ex  = solve(cov.pooled, mean.ex)
 inter.ex   = -sum(slopes.ex*mean.ex)/2
 d.slopes   = slopes.eq - slopes.ex
 d.inter    = inter.eq - inter.ex
 
 d = sum(d.slopes*all.data[j,]) + d.inter
 if (j <= 8) post.eq[j] = exp(d)/(1+exp(d))
 if (j > 8) post.ex[j-8] = 1/(1+exp(d))  
}

Posterior = cbind(1:8, post.eq, 1:8, post.ex)
colnames(Posterior) = c("EQ","P(EQ|data)","EX","P(EX|data)")
# results from cross-validation 
round(Posterior, 3)   
```



Example 7.11

```r
P = 1:1024
S = P+1024
p.dim = 2
n =1024

eq = as.ts(eqexp[,1:8])
ex = as.ts(eqexp[,9:16])
nz = as.ts(eqexp[,17])
f.eq <- array(dim=c(8,2,2,512)) -> f.ex 
f.NZ = array(dim=c(2,2,512))

# determinant for 2x2 complex matrix 
det.c = function(mat){return(Re(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1]))}
L = c(15,13,5)  # for smoothing 
for (i in 1:8){     # compute spectral matrices
 f.eq[i,,,] = mvspec(cbind(eq[P,i],eq[S,i]), spans=L, taper=.5, plot=FALSE)$fxx
 f.ex[i,,,] = mvspec(cbind(ex[P,i],ex[S,i]), spans=L, taper=.5, plot=FALSE)$fxx
}
u = mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5) 
f.NZ = u$fxx	 
bndwidth = u$bandwidth*40  # about .75 Hz
fhat.eq = apply(f.eq, 2:4, mean)  # average spectra
fhat.ex = apply(f.ex, 2:4, mean)

# plot the average spectra
par(mfrow=c(2,2))
Fr = 40*(1:512)/n
tsplot(Fr, Re(fhat.eq[1,1,]), main="", xlab="Frequency (Hz)", ylab="")
tsplot(Fr, Re(fhat.eq[2,2,]), main="", xlab="Frequency (Hz)", ylab="")
tsplot(Fr, Re(fhat.ex[1,1,]), xlab="Frequency (Hz)", ylab="")
tsplot(Fr, Re(fhat.ex[2,2,]), xlab="Frequency (Hz)", ylab="")
mtext("Average P-spectra", side=3, line=-1.5, adj=.2, outer=TRUE)
mtext("Earthquakes", side=2, line=-1, adj=.8,  outer=TRUE)
mtext("Average S-spectra", side=3, line=-1.5, adj=.82, outer=TRUE)
mtext("Explosions", side=2, line=-1, adj=.2, outer=TRUE)
par(fig = c(.75, .995, .75, .98), new = TRUE)
ker = kernel("modified.daniell", L)$coef; ker = c(rev(ker),ker[-1])
plot((-33:33)/40,ker,type="l",ylab="",xlab="",cex.axis=.7,yaxp=c(0,.04,2))

# choose alpha
Balpha = rep(0,19) 
for (i in 1:19){  alf=i/20
for (k in 1:256) {  	
Balpha[i]= Balpha[i] + Re(log(det.c(alf*fhat.ex[,,k] + (1-alf)*fhat.eq[,,k])/ det.c(fhat.eq[,,k]))- 
                           alf*log(det.c(fhat.ex[,,k])/det.c(fhat.eq[,,k])))} }
alf = which.max(Balpha)/20   # = .4  

# calculate information criteria          
rep(0,17) -> KLDiff -> BDiff -> KLeq -> KLex -> Beq -> Bex
for (i in 1:17){
 if (i <= 8) f0 = f.eq[i,,,]
 if (i > 8 & i <= 16) f0 = f.ex[i-8,,,]
 if (i == 17) f0 = f.NZ
 for (k in 1:256) {    # only use freqs out to .25
  tr = Re(sum(diag(solve(fhat.eq[,,k],f0[,,k]))))
  KLeq[i] = KLeq[i] + tr + log(det.c(fhat.eq[,,k])) - log(det.c(f0[,,k]))
  Beq[i] =  Beq[i] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.eq[,,k])/det.c(fhat.eq[,,k])) - 
                         alf*log(det.c(f0[,,k])/det.c(fhat.eq[,,k])))
  tr = Re(sum(diag(solve(fhat.ex[,,k],f0[,,k]))))
  KLex[i] = KLex[i] + tr +  log(det.c(fhat.ex[,,k])) - log(det.c(f0[,,k]))
  Bex[i] = Bex[i] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.ex[,,k])/det.c(fhat.ex[,,k])) - 
                       alf*log(det.c(f0[,,k])/det.c(fhat.ex[,,k]))) 
  }
KLDiff[i] = (KLeq[i] - KLex[i])/n 
BDiff[i] =  (Beq[i] - Bex[i])/(2*n) 
}

x.b = max(KLDiff)+.1; x.a = min(KLDiff)-.1
y.b = max(BDiff)+.01; y.a = min(BDiff)-.01

dev.new()                                 
tsplot(KLDiff[9:16], BDiff[9:16], type="p", xlim=c(x.a,x.b), ylim=c(y.a,y.b), cex=1.1, lwd=2, 
      xlab="Kullback-Leibler Difference",ylab="Chernoff Difference", col=6,
      main="Classification Based on Chernoff and K-L Distances", pch=6)
points(KLDiff[1:8], BDiff[1:8], pch=8, cex=1.1, lwd=2, col=4)
points(KLDiff[17], BDiff[17],  pch=3, cex=1.1, lwd=2, col=3)
legend("topleft", legend=c("EQ", "EX", "NZ"), pch=c(8,6,3), pt.lwd=2, col=c(4,6,3))
abline(h=0, v=0, lty=2, col=8)
text(KLDiff[-c(1,2,3,7,14)]-.075, BDiff[-c(1,2,3,7,14)], label=names(eqexp[-c(1,2,3,7,14)]), cex=.7)
text(KLDiff[c(1,2,3,7,14)]+.075, BDiff[c(1,2,3,7,14)], label=names(eqexp[c(1,2,3,7,14)]), cex=.7)
```


Example 7.12

```r
library(cluster)
P = 1:1024
S = P+1024
p.dim = 2
n =1024

eq = as.ts(eqexp[,1:8])
ex = as.ts(eqexp[,9:16])
nz = as.ts(eqexp[,17])

f = array(dim=c(17,2,2,512))  
L = c(15,15)   # for smoothing 
for (i in 1:8){     # compute spectral matrices
 f[i,,,] = mvspec(cbind(eq[P,i],eq[S,i]), spans=L, taper=.5, plot=FALSE)$fxx
 f[i+8,,,] = mvspec(cbind(ex[P,i],ex[S,i]), spans=L, taper=.5, plot=FALSE)$fxx 
}
f[17,,,] = mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5, plot=FALSE)$fxx	

# calculate symmetric information criteria 
JD = matrix(0,17,17) 
for (i in 1:16){
 for (j in (i+1):17){	
  for (k in 1:256) {    # only use freqs out to .25
    tr1 = Re(sum(diag(solve(f[i,,,k],f[j,,,k]))))
    tr2 = Re(sum(diag(solve(f[j,,,k], f[i,,,k]))))
    JD[i,j] = JD[i,j] + (tr1 + tr2 - 2*p.dim)
  }
 }
}
 
JD = (JD + t(JD))/n
colnames(JD) = c(colnames(eq), colnames(ex), "NZ") 
rownames(JD) = colnames(JD)
cluster.2 = pam(JD, k = 2, diss = TRUE)  

summary(cluster.2)  # print results
par(mar=c(2,2,1,.5)+2, cex=3/4, cex.lab=4/3, cex.main=4/3)
clusplot(JD, cluster.2$cluster, col.clus=gray(.5), labels=3, lines=0, 
         col.p = c(rep(4,8), rep(6,8), 3),  
         main="Clustering Results for Explosions and Earthquakes")
text(-3.5,-1.5, "Group I", cex=1.1, font=2) 
text(1.5,5,"Group II", cex=1.1, font=2)
```


Example 7.13

```r
n = 128
Per = abs(mvfft(fmri1[,-1]))^2/n

par(mfrow=c(2,4), mar=c(3,2,2,1),  mgp = c(1.6,.6,0), oma=c(0,1,0,0))
for (i in 1:8){
 plot(0:20, Per[1:21,i], type="n", xaxt='n', ylim=c(0,8), main=colnames(fmri1)[i+1], xlab="Cycles", ylab="")
 axis(1, seq(0,20,by=4))
 Grid(nx=NA, ny=NULL, minor=FALSE)
 abline(v=seq(0,60,by=4), col='lightgray', lty=1)
 lines(0:20, Per[1:21,i]) }
mtext("Periodogram", side=2, line=-.3, outer=TRUE, adj=c(.2,.8))

fxx = mvspec(fmri1[,-1], kernel("daniell", c(1,1)), taper=.5, plot=FALSE)$fxx
l.val = rep(NA,64)
for (k in 1:64) {
 u = eigen(fxx[,,k], symmetric=TRUE, only.values=TRUE)
 l.val[k] = u$values[1]
}

dev.new()
par(mar=c(2.25,2,.5,.5)+.5,  mgp = c(1.6,.6,0))  
plot(l.val, type="n", xaxt='n',xlab="Cycles (Frequency x 128)", ylab="First Principal Component")
axis(1, seq(4,60,by=8))
Grid(nx=NA, ny=NULL, minor=FALSE)
abline(v=seq(4,60,by=8), col='lightgray', lty=1)
lines(l.val)

# at freq k=4
u = eigen(fxx[,,4], symmetric=TRUE)
lam = u$values 
evec = u$vectors
lam[1]/sum(lam) # % of variance explained
sig.e1 = matrix(0,8,8)
for (l in 2:5){  # last 3 evs are 0
 sig.e1 = sig.e1 + lam[l]*evec[,l]%*%Conj(t(evec[,l]))/(lam[1]-lam[l])^2
}
sig.e1 = Re(sig.e1)*lam[1]*sum(kernel("daniell", c(1,1))$coef^2)
p.val = round(pchisq(2*abs(evec[,1])^2/diag(sig.e1), 2, lower.tail=FALSE), 3)
cbind(colnames(fmri1)[-1], abs(evec[,1]), p.val) # print table values
```


Example 7.14
```r
bhat = sqrt(lam[1])*evec[,1]
Dhat = Re(diag(fxx[,,4] - bhat%*%Conj(t(bhat))))
res = Mod(fxx[,,4] - Dhat - bhat%*%Conj(t(bhat)))
```

Example 7.15
```r
gr = diff(log(ts(econ5, start=1948, frequency=4))) # growth rate
tsplot(100*gr, col=2:6, lwd=2, ncol=2, main="Growth Rates (%)")


# scale each series to have variance 1
gr = ts(apply(gr,2,scale), freq=4)   # scaling strips ts attributes
dev.new()
gr.spec = mvspec(gr, spans=c(7,7), detrend=FALSE, taper=.25, col=2:6, lwd=2)
legend("topright", colnames(econ5), lty=1:5, lwd=2, col=2:6)

dev.new()
plot.spec.coherency(gr.spec, ci=NA,  main="Squared Coherencies")

# PCs
n.freq = length(gr.spec$freq)
lam = matrix(0,n.freq,5)
for (k in 1:n.freq) lam[k,] = eigen(gr.spec$fxx[,,k], symmetric=TRUE, only.values=TRUE)$values 

dev.new()
par(mfrow=c(2,1))
tsplot(gr.spec$freq, lam[,1], ylab="", xlab="Frequency", main="First Eigenvalue")
 abline(v=.25, lty=2)
tsplot(gr.spec$freq, lam[,2], ylab="", xlab="Frequency", main="Second Eigenvalue")
 abline(v=.125, lty=2)

e.vec1 = eigen(gr.spec$fxx[,,10], symmetric=TRUE)$vectors[,1] 
e.vec2 =  eigen(gr.spec$fxx[,,5], symmetric=TRUE)$vectors[,2]
round(Mod(e.vec1), 2);  round(Mod(e.vec2), 3) 
```


Example 7.17  (there is now a script for the spectral envelope)

```r
xdata = dna2vector(bnrf1ebv)
u     = specenv(xdata, spans=c(7,7))  

# details near the peak (coefs are for A, C, G, and T)
round(u,4)[1330:1336,]
```

Example 7.18

```r
x = astsa::nyse    
# possible transformations include absolute value and squared value
xdata = cbind(x, abs(x), x^2)  
par(mfrow=2:1)
u = specenv(xdata, real=TRUE,  spans=c(3,3))
# peak at freq = .001 so let's
# plot the optimal transform 
beta = u[2, 3:5]  # scalings
b = beta/beta[2]  # makes abs(x) coef=1
gopt = function(x) { b[1]*x+b[2]*abs(x)+b[3]*x^2 }
 curve(gopt, -.2, .2, col=4, lwd=2, panel.first=Grid(nym=0))
gabs = function(x) { b[2]*abs(x) } # corresponding to |x|
 curve(gabs, -.2, .2, add=TRUE, col=6)
legend('bottomright', lty=1, col=c(4,6), legend=c('optimal', 'absolute value'), bg='white') 
```

[<sub>top</sub>](#table-of-contents)

---
