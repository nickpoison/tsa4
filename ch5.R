library('astsa')
source('grid.r')

##################
pdf(file="acflogvarve.pdf",width=7.25,height=3.25) 
ACF  =  acf(log(varve), 100, plot=FALSE)
par(mar=c(2.5,2.5,.5,0)+.5, mgp=c(1.6,.6,0))
plot(ACF[1:100], xlab="LAG", ylim=c(-.05,1), panel.first=grid(lty=1))
abline(h=0)
dev.off()

###############
pdf(file="piplot.pdf",width=7.25,height=2.75) 
d  = 0.3841688
p = rep(1,31)
par(mar=c(2,2.5,.5,0)+.5, mgp=c(1.6,.6,0))
for (k in 1:30){ p[k+1] = (k-d)*p[k]/(k+1) }
plot(1:30, p[-1], ylab=expression(pi(d)), lwd=2, xlab="Index", type="h", panel.first=grid(lty=1))
dev.off()


###############
library(fracdiff)
# I might suggest another package such as 'arfima' because 
# this package gives questionable results and it's not easy to
# pull out the residuals
pdf(file="acfvarve3.pdf",width=7.25,height=4) 
lvarve = log(varve)-mean(log(varve))
varve.fd = fracdiff(lvarve, nar=0, nma=0, M=30)
#varve.fd$d  # = 0.3841688
#varve.fd$stderror.dpq  # = 4.589514e-06 (questionable result!!)
res.fd = diffseries(log(varve), varve.fd$d)  # frac diff resids
res.arima = resid(arima(log(varve), order=c(1,1,1))) # arima resids
par(mfrow=c(2,1), mar=c(2,2.5,.5,0)+.5, mgp=c(1.4,.6,0))
acf(res.arima, 100, xlim=c(4,97), ylim=c(-.2,.2), main="", xlab="LAG", panel.first=grid(lty=1))
acf(res.fd, 100, xlim=c(4,97), ylim=c(-.2,.2), main="", xlab="LAG", panel.first=grid(lty=1))
dev.off()


######################
pdf(file="longmemspec.pdf",width=7.25,height=3.5) 
par(mar=c(2,2.5,.5,0)+.5, mgp=c(1.5,.6,0))
series = log(varve) - mean(log(varve))
d0 = .1
n.per = nextn(length(series))
m=(n.per)/2  - 1
per = abs(fft(series)[-1])^2  # remove 0 freq
per = per/n.per            # R doesn't scale fft by sqrt(n)
g = 4*(sin(pi*((1:m)/n.per))^2)
whit.like = function(d){
	g.d=g^d
    sig2 = (sum(g.d*per[1:m])/m)
    log.like = m*log(sig2) + d*sum(log(g)) + m
    return(log.like)
    }	
# -- Estimation --
est=optim(d0,whit.like,gr=NULL,method="L-BFGS-B",hessian=TRUE,
           lower=-.5,upper=.5)
#cat("d.hat =", est$par,"se(dhat) =", 1/sqrt(est$hessian), "\n")
 g.dhat=g^est$par;  sig2=sum(g.dhat*per[1:m])/m
#cat("sig2hat=",sig2, "\n")  
#  spectral approach
 u = spec.ar(log(varve), plot=FALSE)  #produces ar(8)
    g= 4*(sin(pi*((1:500)/2000))^2)
    fhat = sig2*g^{-est$par}
   plot(1:500/2000, log(fhat), type="l", ylab="log(spectrum)", xlab="frequency",  panel.first=grid(lty=1))
   lines(u$freq[1:250], log(u$spec[1:250]), lty="dashed")
dev.off()
   
###############
library(fracdiff)
fdGPH(log(varve), bandw = .9)  # fdGPH doesn't seem to work very well either
##########################

# {innov2acf}
x = resid(arima(diff(log(gnp)), order=c(1,0,0)))
ACF = acf(x^2, 20, plot=FALSE)$acf[-1]
PACF = pacf(x^2, 20, plot=FALSE)$acf
num = length(gnp)-1
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
LAG = 1:20/4
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.2, maxP + 0.2, 1)
pdf(file="innov2acf.pdf",width=7.5,height=3.25)
par(mfrow=c(2,1), mar=c(2,2.5,0,0)+.5, mgp=c(1.4,.6,0))
plot(LAG, ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(LAG, PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()


###############
 # library(TTR)                                          
 # djia = getYahooData("^DJI",start=20060420,end=20160420,freq="daily") 
library(xts) 
djiar = diff(log(djia$Close))[-1]   # djia is in astsa now 
 acf2(djiar)     # exhibits some acf
 acf2(djiar^2)   # exhibits a lot of autocorr
library(fGarch)
summary(djia.g <- garchFit(~arma(1,0)+garch(2,1), data=djiar, cond.dist='std'))
u = djia.g@sigma.t
pdf(file="djiapred.pdf",width=7.5,height=3.5)
par(mar=c(2.5,2,0,.5)+.5, mgp=c(1.6,.6,0), oma=rep(0,4))
plot(djiar[400:900], type='n', main='')
lines(djiar[400:900], col=gray(.6))
lines((u+djiar-djiar)[400:900], col=4)
dev.off()


################# power arch (assuming djiar still available)
summary(fit <- garchFit(~arma(1,0)+aparch(1,2), data=djiar, cond.dist='std'))
v = volatility(fit, type = "sigma")
#par(mfrow=c(2,1))
#plot(fit, which=c(1,2))    # use plot(fit) to see all options
u = fit@sigma.t
# pdf(file="djiapredaparch.pdf",width=7.5,height=3.25)  # don't want to plot this
par(mar=c(2.5,2,0,.5)+.5, mgp=c(1.6,.6,0), oma=rep(0,4))
plot(djiar[400:900], type='n', main='')
lines(djiar[400:900], col=gray(.6))
lines((u+djiar-djiar)[400:900], col=4)
# dev.off()


###################################
pdf(file="flu.pdf",width=7.7,height=3.5)  
par(mar=c(2.5,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(flu, type='n')
grid(lty=1)
lines(flu, type="c")
Months = c("J","F","M","A","M","J","J","A","S","O","N","D")
points(flu, pch=Months, cex=.8, font=2)
dev.off()


#########  setar flu ###################
dflu = diff(flu)
thrsh = .05 # threshold
Z = ts.intersect(dflu, lag(dflu,-1), lag(dflu,-2), lag(dflu,-3), lag(dflu,-4))
ind1 = ifelse(Z[,2] < thrsh, 1, NA)  # indicator < thrsh
ind2 = ifelse(Z[,2] < thrsh, NA, 1)  # indicator >= thrsh
X1 = Z[,1]*ind1
X2 = Z[,1]*ind2
summary(fit1<-lm(X1~Z[,2:5]))  # case 1
summary(fit2<-lm(X2~Z[,2:5]))  # case 2
D = cbind(rep(1, nrow(Z)), Z[,2:5])  # get predictions
b1 = fit1$coef
b2 = fit2$coef
p1 = D%*%b1
p2 = D%*%b2
prd = ifelse(Z[,2] < thrsh, p1, p2)
pdf(file="flu_setar.pdf",width=7.7,height=3.5)
par(mar=c(2.5,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(dflu, ylim=c(-.5,.5), type='n')
grid(lty=1)
points(dflu, pch=3)
lines(prd)
prde1 = sqrt(sum(resid(fit1)^2)/df.residual(fit1))
prde2 = sqrt(sum(resid(fit2)^2)/df.residual(fit2))
prde = ifelse(Z[,2] < thrsh, prde1, prde2)
    x = time(dflu)[-(1:4)]
   xx = c(x, rev(x))
   yy = c(prd - 2*prde, rev(prd + 2*prde))
polygon(xx, yy, border=8, col=rgb(.6,.6,.6,alpha=.25))
abline(h=.05, col=4, lty=6)
dev.off()

################
## using tsDyn
#require(tsDyn)              # load package - install it if you don't have it
## vignette("tsDyn")         # for package details (it's quirky, so you'll need this)
#dflu = diff(flu)
#lag1.plot(dflu, 1, corr=FALSE)     # see the nonlinearity here
#(u = setar(dflu, m=4, thDelay=0))
#plot(dflu, ylim=c(-.5,.5), type='p', pch=3)
#lines(fitted(u))
#######

######
pdf(file="dflu_scat.pdf",width=6,height=3) 
par(mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
 U= matrix(Z, ncol=5)
 plot(U[,2],U[,1], panel.first=grid(lty=1), xlab=expression(dflu[~t-1]), ylab=expression(dflu[~t]))#,pch=20, col=v)
 lines(lowess(U[,2], U[,1], f =2/3), col=2)
 abline(v=.05, lty=2, col=4)
dev.off()
 
 ################################### transfer function models ##########
# {soiacfpacf} and {transccf}
soi.d = resid(lm(soi~time(soi), na.action=NULL))
# acf2(soi.d)
fit = arima(soi.d,  order=c(1, 0, 0)) 
ar1 = as.numeric(coef(fit)[1])    # = 0.5875
soi.pw = resid(fit)
rec.fil = filter(rec, filter=c(1, -ar1), method="conv", sides=1)
#ccf(soi.pw, rec.fil, main="", ylab="CCF", na.action=na.omit)

##  plots
u = acf2(soi.d, 36); dev.off(); ACF=u[,1]; PACF=u[,2]
num = length(rec)
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.2, maxP + 0.2, 1)
pdf(file="soiacfpacf.pdf",width=7.5,height=4) 
par(mfrow=c(2,1), mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
plot(ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
dev.off()

####
pdf(file="transccf.pdf",width=7.5,height=3.25) 
par(mar=c(2.5,2.5,1,0)+.5, mgp=c(1.6,.6,0))
ccf(soi.pw, rec.fil, 24, xlab='LAG', main="", ylab="CCF", panel.first=grid(lty=1), na.action=na.omit)
mtext(side=3, 'preSOI vs filRec')
dev.off()


######## the analysis for last trans func model
#y_t= a y_{t-1}+ b x_{t-5}+u_t
soi.d = resid(lm(soi~time(soi), na.action=NULL))
fish = ts.intersect(R=rec, RL1=lag(rec,-1), SL5=lag(soi.d,-5))
(u = lm(fish[,1]~fish[,2:3], na.action=NULL))
# acf2(resid(u))  # suggests ar1
uar = sarima(fish[,1], 1, 0, 0, xreg=fish[,2:3])   # armax model

######################
#transfer plot
pdf(file="transfer.pdf",width=7.5,height=4) 
layout(matrix(c(1,1,3,3,3,  2,2,3,3,3), nc=2))
par(mar=c(2.25,2.5,.5,0)+.5, mgp=c(1.6,.5,0))
ACF = acf(resid(u),24, plot=FALSE)$acf[-1]
PACF = pacf(resid(u),24, plot=FALSE)$acf 
num = length(resid(u))
minA = min(ACF)
maxA = max(ACF)
minP = min(PACF)
maxP = max(PACF)
U = 2/sqrt(num)
L = -U
minu = min(minA, minP, L) - 0.01
maxu = min(maxA + 0.1, maxP + 0.1, 1)
plot(ACF, type="h", xlab="LAG", ylim = c(minu, maxu), panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2)  
plot(PACF, type="h", xlab="LAG",  ylim = c(minu, maxu) , panel.first=grid(lty=1)); abline(h=0)
abline(h = c(L, U), col=4, lty=2) 
pred=rec+resid(uar$fit)
ts.plot(rec, pred, ylab="Recruitment", type='n')
grid(lty=1)
lines(rec) 
lines(pred, lwd=7, col=rgb(.6,.6,.6, alpha=.25))
dev.off()

##################  multivariate #######################
library(vars)
x = cbind(cmort, tempr, part)
summary(VAR(x, p=1, type='both'))       # 'both' fits constant + trend
#
round(VARselect(x, lag.max=10, type="both")$criteria[1:3,],3)
#
summary(fit <- VAR(x, p=2, type="both"))  # partial results displayed
#
pdf(file="var2resccf.pdf")
acf(resid(fit), 52, mar=c(2,1,2,0)+.5,    mgp=c(1.5,.5,0))
dev.off()
serial.test(fit, lags.pt=12, type="PT.adjusted") 
#
(fit.pr = predict(fit, n.ahead = 24, ci = 0.95))  # 4 weeks ahead
#
pdf(file="fanchart.pdf",width=7.5,height=4) 
par(mar=c(1,1,1,0)+.5, mgp=c(1.6,.5,0), font.main=1)
fanchart(fit.pr)  # plot prediction + error
dev.off()

###########################
# spliid
require(marima)
model = define.model(kvar=3,ar=c(1,2),ma=c(1))
arp = model$ar.pattern
map = model$ma.pattern
cmort.d = resid(detr <- lm(cmort~time(cmort), na.action=NULL))
xdata = matrix(cbind(cmort.d, tempr, part), ncol=3)  # strip ts attributes
fit <- marima(xdata, ar.pattern=arp, ma.pattern=map, means=c(0,1,1), penalty=1)
# results
options(digits=3)
short.form(fit$ar.estimates, leading=FALSE) # print estimates
short.form(fit$ar.fvalues, leading=FALSE)   # print t^2-statistic
short.form(fit$ma.estimates, leading=FALSE)
short.form(fit$ma.fvalues, leading=FALSE) # print estimates
fit$resid.cov
#  resid analysis
innov = t(resid(fit))
plot.ts(innov)
acf(innov, na.action=na.pass)
#  fitted
pred = ts(t(fitted(fit))[,1], start=start(cmort), freq=frequency(cmort))+detr$coef[1]+detr$coef[2]*time(cmort)
pdf(file="spliid.pdf",width=7.5,height=3) 
par(mar=c(2,2,.5,.5)+.5, mgp=c(1.6,.5,0))
ts.plot(cmort, pred, ylab="Cardiovascular Mortality", type='n')
grid(lty=1)
points(cmort) 
lines(pred, lwd=2, col=4)#rgb(.6,.6,.6, alpha=.5))
dev.off()



