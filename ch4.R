library('astsa')
source('grid.r')

##############
pdf(file="cosines.pdf",width=7.25,height=4.25) 
x1 = 2*cos(2*pi*1:100*6/100) + 3*sin(2*pi*1:100*6/100)
x2 = 4*cos(2*pi*1:100*10/100) + 5*sin(2*pi*1:100*10/100)
x3 = 6*cos(2*pi*1:100*40/100) + 7*sin(2*pi*1:100*40/100)
x = x1 + x2 + x3
par(mfrow = c(2,2), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0))
plot.ts(x1, ylim=c(-10,10), main=expression(omega==6/100~~~A^2==13),  panel.first=grid(lty=1))
plot.ts(x2, ylim=c(-10,10), main=expression(omega==10/100~~~A^2==41),  panel.first=grid(lty=1))
plot.ts(x3, ylim=c(-10,10), main=expression(omega==40/100~~~A^2==85),  panel.first=grid(lty=1))
plot.ts(x,  ylim=c(-16,16), main="sum",  panel.first=grid(lty=1), font.main=1)
dev.off()

##############
pdf(file="period1.pdf",width=7.25,height=3.25) 
par(mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
P = Mod(2*fft(x)/100)^2;  Fr = 0:99/100
plot(Fr, P, type="o", xlab="frequency", ylab="scaled periodogram", panel.first=grid(lty=1), ylim=c(0,90) )
abline(v=.5, lty=2, col=4)
abline(v=c(.1,.3,.7,.9), lty=1, col='lightgray')
dev.off()


######################
pdf(file="star.pdf",width=7.25,height=4) 
n = length(star)
layout(matrix(c(1,2), ncol = 1), height=c(1.25,1))
par( mar=c(2,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(star, ylab="star magnitude", xlab="day", type='n' )
grid(lty=1); lines(star)
Per = Mod(fft(star-mean(star)))^2/n
Freq = (1:n -1)/n
plot(Freq[1:50], Per[1:50], type='n', ylab="Periodogram", xlab="Frequency")
grid(); lines(Freq[1:50], Per[1:50], type='h', lwd=3)
#u = which.max(Per[1:50])         # 22  freq=21/600=.035 cycles/day
#uu = which.max(Per[1:50][-u])    # 25  freq=25/600=.041 cycles/day 
#1/Freq[22]; 1/Freq[26]           # period = days/cycle
text(.05, 7000, "24 day cycle"); text(.027, 9000, "29 day cycle")
dev.off()

###########################
pdf(file="theospec.pdf",width=7.25,height=5) 
par(mfrow=c(3,1), mar=c(3,3,1.5,1), mgp=c(1.6,.6,0), cex.main=1.1)
arma.spec(log="no", main="White Noise")
arma.spec(ma=.5, log="no", main="Moving Average")
arma.spec(ar=c(1,-.9), log="no", main="Autoregression")
dev.off()


###########################
### spectral plots
pdf(file="soirecper.pdf",width=7.5, height=6) 
par(mfrow=c(2,1), mar=c(3.5,3,2.5,1),  mgp=c(1.5,.6,0), oma=rep(0,4), font.main=1)
soi.per = mvspec(soi, log="no", type='n')
grid(lty=1); par(new=TRUE)
mvspec(soi, log="no") 
abline(v=1/4, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
rec.per = mvspec(rec, log="no", type='n')
grid(lty=1); par(new=TRUE)
mvspec(rec, log="no") 
abline(v=1/4, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
dev.off()

###########################
pdf(file="allthesame.pdf",width=7, height=2.5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.5,.6,0))
u=arma.spec(ar=c(1,-.9), xlim=c(.15,.151), ylim=c(10, 100), n.freq=100000, col='white')
grid(lty=1, equilogs = FALSE )
lines(u$freq,u$spec, lwd=2)
dev.off()


###########################
pdf(file="soirecsmooth1.pdf",width=7.5, height=6) 
par(mfrow=c(2,1), mar=c(3.5,3,2.5,1),  mgp=c(1.5,.6,0), oma=rep(0,4), font.main=1)
k = kernel("daniell", 4)
soi.ave = mvspec(soi, k, log="no", type='n')
grid(lty=1, col=gray(.9)); par(new=TRUE)
mvspec(soi, k, log="no")
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
rec.ave = mvspec(rec, k, log="no", type='n')
grid(lty=1, col=gray(.9)); par(new=TRUE)
mvspec(rec, k, log="no")
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
dev.off()

##########
pdf(file="soirecsmooth1log.pdf",width=7.5, height=6) 
par(mfrow=c(2,1), mar=c(3.5,3,2.5,1),  mgp=c(1.5,.6,0), oma=rep(0,4), font.main=1)
k = kernel("daniell", 4)
soi.ave = mvspec(soi, k, type='n')
grid(lty=1); par(new=TRUE)
mvspec(soi, k)
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
rec.ave = mvspec(rec, k, type='n')
grid(lty=1); par(new=TRUE)
mvspec(rec, k)
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
dev.off()




###################
pdf(file="harmonics_col.pdf",width=7.5,height=4.5)  # works with scale=.6
par(mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
 t=seq(0,1,by=1/200)
 x1 =  ts(sin(2*pi*2*t),start=0,deltat=1/200)
 x2 = ts(.5*sin(2*pi*t*4),start=0,deltat=1/200)
 x3 = ts(.4*sin(2*pi*t*6),start=0,deltat=1/200)
 x4 = ts(.3*sin(2*pi*t*8),start=0,deltat=1/200)
 x5 = ts(.2*sin(2*pi*t*10),start=0,deltat=1/200)
 x6 = ts(.1*sin(2*pi*t*12),start=0,deltat=1/200)
 xsum = x1+x2+x3+x4+x5+x6
 dcyan = rgb(0,.6,.6)
 sgreen =  rgb(0, .7, 0)
 ts.plot(x1,x2,x3,x4,x5,x6,xsum, lty=c(1,5,2,5,2,5,1), lwd=c(rep(1,6),2), col=c(4,sgreen,2,4,dcyan,6,1), ylab="Sinusoids")
 names=c("Fundamental","2nd Harmonic","3rd Harmonic","4th Harmonic", "5th Harmonic", "6th Harmonic", "Formed Signal")
 legend("topright", names, lty=c(1,5,2,5,2,5,1),  col=c(4,sgreen,2,4,dcyan,6,1), lwd=c(rep(1,6),2), cex=.8)
dev.off()






############## note tapering is used
pdf(file="soirecsmooth2.pdf",width=7.5, height=6) 
par(mfrow=c(2,1), mar=c(3.5,3,2.5,1),  mgp=c(1.5,.6,0), oma=rep(0,4), font.main=1)
k = kernel("modified.daniell", c(3,3))
soi.ave = mvspec(soi, k, log="no", taper=.1,  type='n')
grid(lty=1); par(new=TRUE)
mvspec(soi, k, taper=.1, log="no")
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
rec.ave = mvspec(rec, k, log="no", taper=.1,  type='n')
grid(lty=1); par(new=TRUE)
mvspec(rec, k, taper=.1, log="no")
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.75)
dev.off()



###############  tapering example
w=seq(-.05,.05,.0001); n=480; u=0
for (i in -4:4){ k=i/n
 u=u+sin(n*pi*(w+k))^2/sin(pi*(w+k))^2}
fk=u/(9*480)
u=0; wp=w+1/n; wm=w-1/n
for (i in -4:4){
 k=i/n; wk=w+k; wpk=wp+k; wmk=wm+k
 z =  complex(real=0,imag=2*pi*wk)
 zp = complex(real=0,imag=2*pi*wpk)
 zm = complex(real=0,imag=2*pi*wmk)
 d =  exp(z)*(1-exp(z*n))/(1-exp(z))
 dp = exp(zp)*(1-exp(zp*n))/(1-exp(zp))
 dm = exp(zm)*(1-exp(zm*n))/(1-exp(zm))
 D = .5*d - .25*dm*exp(pi*w/n)-.25*dp*exp(-pi*w/n)
 D2=abs(D)^2
 u=u+D2 }
sfk=u/(480*9)
#
pdf(file="fejertaper.pdf",width=7.5, height=6) 
par(mfrow=c(2,2), mar=c(3,2,2,1), mgp = c(1.6,0.6,0))
plot(w,fk, type="l", ylab="",xlab="frequency", main="Smoothed Fejer")
mtext(expression({}+{}), side=1, line=-.42, at=c(-0.009375, .009375), cex=1.5, col=4)
segments(-4.5/480, -2.2, 4.5/480, -2.2 , lty=1, lwd=2, col=4)
plot(w,log(fk), type="l", ylab="",xlab="frequency",  main="Smoothed Fejer - logged")
plot(w,sfk, type="l", ylab="",xlab="frequency", main="Cosine Taper")
mtext(expression({}+{}), side=1, line=-.38, at=c(-0.009375, .009375), cex=1.5, col=4)
segments(-4.5/480, -.8, 4.5/480, -.8 , lty=1, lwd=2, col=4)
plot(w,log(sfk), type="l", ylab="",xlab="frequency", main="Cosine Taper - logged")
dev.off()


######################
pdf(file="taper2.pdf",width=7.5, height=3.5) 
par(mar=c(2.5,2.5,1,1),  mgp=c(1.5,.6,0))
s0 = mvspec(soi, spans=c(7,7), plot=FALSE)             # no taper
s50 = mvspec(soi, spans=c(7,7), taper=.5, plot=FALSE)  # full taper
plot(s50$freq , s50$spec,   log="y", type="l", lty=1, ylab="spectrum", xlab="frequency", panel.first=grid(lty=1))     # solid line
lines(s0$freq, s0$spec, lty=2)  # dashed line
text(1.42, 0.04, 'leakage', cex=.8)
arrows(1.4, .035, .75, .009, length=0.05,angle=30)   
arrows(1.4, .035, 1.21, .0075, length=0.05,angle=30)
abline(v=.25, lty=2, col=8)
mtext('1/4',side=1, line=0, at=.25, cex=.9)
par(fig = c(.65, 1, .65, 1),  new = TRUE, cex=.5,  mgp=c(0,-.1,0), tcl=-.2)
taper <- function(x) { .5*(1+cos(2*pi*x))  }
 x <- seq(from = -.5, to = .5, by = 0.001)
plot(x, taper(x), type = "l",  lty = 1,  yaxt='n', ann=FALSE)
dev.off()


######################
pdf(file="aicbic.pdf",width=7.5, height=3.5) 
par(mar=c(2.75,2.5,1,1), mgp=c(1.5,.6,0))
n = length(soi)
AIC = rep(0, 30) -> AICc -> BIC
for (k in 1:30){
fit = ar(soi, order=k, aic=FALSE)
sigma2 = fit$var.pred
BIC[k] = log(sigma2) + (k*log(n)/n)
AICc[k] = log(sigma2) + ((n+k)/(n-k-2))
AIC[k] = log(sigma2) + ((n+2*k)/n)
}
IC = cbind(AIC, BIC+1)
ts.plot(IC,  type='n',  xlab="p", ylab="AIC / BIC")
grid(lty=1); par(new=TRUE)
ts.plot(IC, type="o", xlab='', ylab='')
text(15.2, -1.48, "AIC")
text(15, -1.35, "BIC")
dev.off()



######################
pdf(file="soiarper2.pdf",width=7.5, height=3.5) 
par(mar=c(2.75,2.5,1,1), mgp=c(1.5,.6,0))
spaic = spec.ar(soi, log="no", plot=FALSE) 
plot(spaic$freq , spaic$spec, type="n", ylab="spectrum", xlab="frequency")  
grid(lty=1)
mtext('1/4',side=1, line=0, at=.25, cex=.9)
abline(v=.25, lty=2)
lines(spaic$freq , spaic$spec, type="l",  ylab="spectrum", xlab="frequency")   
dev.off()


##################################
pdf(file="soireccoh.pdf",width=7.5, height=3.5)
par(mar=c(2.75,2.5,1,1), mgp=c(1.5,.6,0), font.main=1, cex.main=1.1)
sr = mvspec(cbind(soi,rec), kernel("daniell",9), plot=FALSE) 
f = qf(.999, 2, sr$df-2)  
C = f/(18+f)  
plot(sr, plot.type = "coh", ci.lty = 2, panel.first=grid(lty=1), main='SOI & Recruitment')
abline(h = C)
dev.off()

######################
pdf(file="soifilter.pdf",width=7.5, height=5)
par(mfrow=c(3,1), mar=c(3,2,1.5,1), mgp=c(1.6,.6,0))
plot(soi, type='n', ylab='') # plot data
grid(lty=1)
lines(soi)
mtext(side=3, 'SOI')
plot(diff(soi), type='n', ylab='') # plot first difference
grid(lty=1)
lines(diff(soi))
mtext(side=3, 'First Difference')
k = kernel("modified.daniell", 6) # filter weights
plot(soif <- kernapply(soi, k), type='n', ylab='') # plot 12 month filter
grid(lty=1)
lines(soif)
mtext(side=3, 'Seasonal Moving Average')
dev.off()


######################
pdf(file="freqresp.pdf",width=7.5, height=4.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,.5,.5), mgp=c(1.4,.4,0), cex=.9)
##-- frequency responses --##
w = seq(0, .5, by=.001)
FRdiff = abs(1-exp(2i*pi*w))^2
plot(w, FRdiff, type='l', xlab='frequency', panel.first=grid(lty=1), ylab='First Difference')
u = cos(2*pi*w)+cos(4*pi*w)+cos(6*pi*w)+cos(8*pi*w)+cos(10*pi*w)
FRma = ((1 + cos(12*pi*w) + 2*u)/12)^2
plot(w, FRma, type='l', xlab='frequency', panel.first=grid(lty=1), ylab='Seasonal Moving Average')
dev.off()


###########################
# lagreg1 and lagreg2
u1 = LagReg(soi, rec, L=15, M=32, threshold=6)
dev.off()
u2 = LagReg(rec, soi, L=15, M=32, inverse=TRUE, threshold=.01)
dev.off()
pdf(file="lagreg.pdf",width=7.5, height=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,1,.5), mgp=c(1.4,.4,0), cex=.9)
plot(u1$betas, type = "h", xlab = "s", ylab = "beta(s)", panel.first=grid(lty=1))
abline(h=0)
mtext(side=3,'Input: SOI' )
plot(u2$betas, type = "h", xlab = "s", ylab = "beta(s)", panel.first=grid(lty=1))
abline(h=0)
mtext(side=3,'Input: Recruitment' )
dev.off()



############################
## SigExtract(soi, L=9, M=64, max.freq=.05)
L = 9;  M = 64; max.freq=.05
series = ts(soi, frequency = 1)  
spectra = stats::spec.pgram(series, spans=L, plot = FALSE)
A <- function(nu) {
qwe = ifelse((nu > .01 && nu < max.freq), 1, 0) 
qwe  # Sets A(nu) to be qwe
}
N = 2*length(spectra$freq)  # This will be T'.
sampled.indices = (N/M)*(1:(M/2))  # These are the indices of the frequencies we want
fr.N = spectra$freq
fr.M = fr.N[sampled.indices]  # These will be the frequencies we want
spec.N = spectra$spec
spec.M = spec.N[sampled.indices] # Power at these frequencies
A.desired = vector(length = length(fr.M))
for(k in 1:length(fr.M)) A.desired[k] = A(fr.M[k])
# Invert A.desired, by discretizing the defining integral, to get the coefficients a:
delta = 1/M
Omega = seq(from = 1/M, to = .5, length = M/2)
aa = function(s) 2*delta*sum(exp(2i*pi*Omega*s)*A.desired)
S = ((-M/2+1):(M/2-1))
a = vector(length = length(S))
for(k in 1:length(S)) a[k] = aa(S[k])
a = Re(a)  # The filter coefficients 
# Apply a cosine taper
h = .5*(1+cos(2*pi*S/length(S)))
a = a*h    # Comment out this line, to see the effect of NOT tapering
# Compute the realized frequency response function, and the filtered series:
A.M = function(nu) Re(sum(exp(-2i*pi*nu*S)*a))
A.attained = vector(length = length(fr.N))
A.theoretical = vector(length = length(fr.N))
for(k in 1:length(fr.N)) {
A.attained[k] = A.M(fr.N[k]) # The attained freq. resp.
A.theoretical[k] = A(fr.N[k])
}
series.filt = stats::filter(series, a, sides = 2) # The filtered series

###########
pdf(file="sigextract.pdf",width=7.5, height=4)
par(mfrow=c(2,1), mar=c(2.5,2.5,1,.5), mgp=c(1.25,.6,0), cex.lab=.8, cex.axis=.8, font.main=1, cex.main=.9)
plot.ts(series, type='n') 
grid(lty=1)
lines(series)
mtext(side=3, "Original series")
plot.ts(series.filt, type='n') 
grid(lty=1)
lines(series.filt)
mtext(side=3, "Filtered series")
dev.off()

###########
pdf(file="sigextract_coef.pdf",width=7.5, height=4)
par(mfrow=c(2,1), mar=c(2.5,2.5,1,.5), mgp=c(1.25,.6,0), cex.lab=.8, cex.axis=.8, font.main=1, cex.main=.9)
plot(S, a, xlab = "s", ylab = "a(s)", main = "Filter coefficients", panel.first=grid(lty=1))
plot(fr.N, A.theoretical, type = "l", lty = 6, col=4, xlab = "freq", ylab = "freq. response", 
    main = "Desired and attained frequency response functions", panel.first=grid(lty=1))
lines(fr.N, A.attained, lty = 1, col = 2)
dev.off()


###############################
 pdf(file="2dfft.pdf", height=5)
 per = abs(fft(soiltemp-mean(soiltemp))/sqrt(64*36))^2       
 per2 = cbind(per[1:32,18:2], per[1:32,1:18])   
 per3 = rbind(per2[32:2,],per2)
 par(mar=c(1,1,0,0)+.5, cex.axis=.7)
 persp(-31:31/64, -17:17/36, per3, phi=30, theta=30, expand=.6, ticktype="detailed", xlab="cycles/row",  ylab="cycles/column", zlab="Periodogram Ordinate")
 dev.off()

 
#######################
pdf(file="sunspotz.pdf",width=7.5,height=3)  # works with scale=.6
par(mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0))
plot(sunspotz, type='n')
grid(lty=1)
lines(sunspotz)
dev.off()

