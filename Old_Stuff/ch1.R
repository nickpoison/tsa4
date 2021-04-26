library('astsa')
source('grid.r')  # I changed the defaults so that
#  grid() gives grid(lty=1, col = gray(.9))


################
pdf(file="jj.pdf", width=7.5,height=3.25)   
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(jj, ylab="Quarterly Earnings per Share",   type='n')
grid(lty=1); lines(jj, type="o")
dev.off()

############ this is in the EZ version only #########
pdf(file="4bros.pdf", width=6.5, height=2.75)
par(mfrow=1:2, mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5,.6,0))
t=1:(10*2)
x = 100*1.1^t
y = 150*1.1^t
z = 200*1.1^t
w = 75*1.1^t 
culer = c(1,2,'darkgreen',4)
u=t(cbind(x,y,z,w))
x=ts(c(u), freq=4)
plot(x, xlab='quarter', ylab='value', col=gray(.6))
points(x,  pch=c('1','2','3','4'), col=culer, cex=.8)
plot(log(x), xlab='quarter', ylab='log(value)', col=gray(.6))
points(log(x),  pch=c('1','2','3','4'), col=culer, cex=.8)
dev.off()


###########
pdf(file="globtemp.pdf",  width=7.5,height=3.25)  
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(globtemp, ylab="Global Temperature Deviations",    type='n')   
grid(lty=1); lines(globtemp, type='o')
dev.off()

###########
pdf(file="speech.pdf", width=7.5,height=3.25)   
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(speech,   type='n', ylim=c(0,4100))
grid(lty=1);   lines(speech); 
dev.off()

###########
# library(TTR)                                          
# djia = getYahooData("^DJI",start=20060420,end=20160420,freq="daily") 
library(xts)
djiar = diff(log(djia$Close))[-1]         
pdf(file="djiar.pdf", width=7.5,height=3.25)   
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot(djiar, ylab="DJIA Returns", main='', type='n')
lines(djiar)
dev.off()


#########
pdf(file="fish.pdf",width=7.5,height=6)  
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(soi, ylab="", xlab="", main="Southern Oscillation Index",   type='n')
grid(lty=1); lines(soi)
plot(rec, ylab="", main="Recruitment",  type='n')
grid(lty=1); lines(rec)
dev.off()


#########
pdf(file="fmri1.pdf",width=7.5,height=6)  
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
ts.plot(fmri1[,2:5], ylab="BOLD", xlab="", main="Cortex", type='n')
grid(lty=1); par(new=TRUE)
ts.plot(fmri1[,2:5], col=1:4, ylab="BOLD", xlab="", main="Cortex")
#
ts.plot(fmri1[,6:9], ylab="BOLD", xlab="", main="Thalamus & Cerebellum", type='n')
grid(lty=1); par(new=TRUE)
ts.plot(fmri1[,6:9], col=1:4, ylab="BOLD", xlab="", main="Thalamus & Cerebellum")
mtext("Time (1 pt = 2 sec)", side=1, line=1.5)
dev.off()


#####################
pdf(file="eqexp.pdf",width=7.5,height=6)   
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
plot(EQ5, main="Earthquake", xlab="", type='n')
grid(lty=1); lines(EQ5)
plot(EXP6, main="Explosion", xlab="", type='n')  
grid(lty=1); lines(EXP6)
mtext("Time", side=1, line=1.5)
dev.off()


####################
pdf(file="wn_ma.pdf",width=7.5,height=6)   
par(mfrow = c(2,1), mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
set.seed(1)
w = rnorm(500,0,1)                    # 500 N(0,1) variates
v = filter(w, sides=2, filter=rep(1/3,3))  # moving average
plot.ts(w, main="white noise", type='n')
grid(lty=1, col=gray(.9)); lines(w)
plot.ts(v, ylim=c(-3,3), main="moving average", type='n')
grid(lty=1, col=gray(.9)); lines(v)
dev.off()


####################
pdf(file="ar2.pdf",width=7.5,height=3.25)   
par(mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0))
w = rnorm(550,0,1)  # 50 extra to avoid startup problems
x = filter(w, filter=c(1,-.9), method="recursive")[-(1:50)]
plot.ts(x, main="autoregression", type='n')
grid(lty=1, col=gray(.9)); lines(x)
dev.off()


#################
pdf(file="rw.pdf",width=7.5,height=3.25)   
par(mar=c(2,2,1,0)+.5, mgp=c(1.6,.6,0), cex.main=1.05)
set.seed(154)                # so you can reproduce the results
w = rnorm(200,0,1);  x = cumsum(w)   # two commands in one line
wd = w +.2;   xd = cumsum(wd)
plot.ts(xd, ylim=c(-5,55), main="random walk", ylab='',   type='n')
grid(lty=1); lines(xd)
lines(x, col=4); abline(h=0, col=4, lty=2)
abline(a=0, b=.2, lty=2)
dev.off()


####################
pdf(file="cos.pdf",width=7.5,height=6)   
par(mfrow = c(3,1), mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
cs = 2*cos(2*pi*1:500/50 + .6*pi);  w = rnorm(500,0,1)
par(mfrow=c(3,1), mar=c(3,2,2,1), cex.main=1.05)
plot.ts(cs, ylab='',xlab='', main=expression(2*cos(2*pi*t/50+.6*pi)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs) 
plot.ts(cs+w, ylab='',xlab='',main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,1)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs+w) 
plot.ts(cs+5*w, ylab='', main=expression(2*cos(2*pi*t/50+.6*pi) + N(0,5^2)), type='n', cex.main=1.5)
 grid(lty=1, col=gray(.9)); lines(cs+5*w) 
dev.off()


###################
pdf(file="ma3.pdf",width=6,height=2.5) 
ACF = c(0,0,0,1/3,2/3,1,2/3,1/3,0,0,0)
LAG = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)
par(mar=c(2,2,0,0)+.5, mgp=c(1.6,.6,0))
plot(LAG, ACF, type='h', lwd=3,  panel.first=grid(lty=1))
abline(h=0)
points(LAG[-(4:8)],ACF[-(4:8)], pch=20)
dev.off()



###################
pdf(file="ccflagx5.pdf",width=6,height=2.5) 
par(mar=c(2.2,2,.7,0)+.5, mgp=c(1.6,.6,0))
set.seed(2)
x = rnorm(100)
y = lag(x,-5) + rnorm(100)
ccf(y,x, ylab='CCovF', xlab="LAG", type='covariance',panel.first=grid(lty=1, col=gray(.9)), main="", lwd=2)
abline(v=0, lty=2)
text(11, .9, 'x leads')
text(-9, .9, 'y leads')
title(main="y & x", cex.main=1)
dev.off()


#######################
pdf(file="acf_scatter_soi.pdf",width=6,height=3) 
r = round(acf(soi, 6, plot=FALSE)$acf[-1], 3) # first 6 sample acf values
par(mfrow=c(1,2), mar=c(2.5,2.5,0,0)+.5, mgp=c(1.6,.6,0))
plot(lag(soi,-1), soi,panel.first=grid(lty=1)); legend('topleft', legend=r[1], bg='white', adj=.25, cex = 0.85)
plot(lag(soi,-6), soi,panel.first=grid(lty=1)); legend('topleft', legend=r[6], bg='white', adj=.25, cex = 0.8)
dev.off()


######################
pdf(file="acfspeech.pdf",width=7.5,height=3.25)   
ACF = acf(speech, 250, plot = FALSE)$acf[-1]
LAG = 1:250
minA = min(ACF)
maxA = max(ACF)
U = 2/sqrt(num)
L = -U
minu = min(minA, L) - 0.01
maxu = min(maxA + 0.1, 1)
plot(LAG, ACF, type = "n", ylim = c(minu, maxu))
grid(lty = 1, col = gray(0.9))
abline(h = c(0, L, U), lty = c(1, 2, 2), col = c(1, 4, 4))
lines(LAG, ACF, type = "h")
#acf(speech, 250, panel.first=grid(lty=1))
dev.off()


###################
pdf(file="soiltempplot.pdf",  height=4.5)          
par(mar=c(0,1,0,0)+.5, mgp=c(1.6,.6,0))
persp(1:64, 1:36, soiltemp,  phi=25, theta=25, scale=FALSE, expand=4, ticktype="detailed", xlab="rows", ylab="cols", zlab="temperature")
dev.off()

##
pdf(file="soiltempave.pdf",width=7.5,height=3)  
par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
plot.ts(rowMeans(soiltemp), xlab="row", ylab="Average Temperature" , type='n')
grid(lty=1); lines(rowMeans(soiltemp))
dev.off()


#######################
pdf(file="soiltemp2dacf.pdf",height=4.5) 
par(mar=c(1,1,0,0)+.5)
fs = Mod(fft(soiltemp-mean(soiltemp)))^2/(64*36)
cs = Re(fft(fs, inverse=TRUE)/sqrt(64*36)) # ACovF
rs = cs/cs[1,1] # ACF
rs2 = cbind(rs[1:41,21:2], rs[1:41,1:21])
rs3 = rbind(rs2[41:2,], rs2)
persp(-40:40, -20:20, rs3, phi=30, theta=30, expand=30, scale="FALSE",
ticktype="detailed", xlab="row lags", ylab="column lags", zlab="ACF")
dev.off()


########
pdf(file="ccfsoirec.pdf",width=7.5,height=6) 
par(mfrow=c(3,1), mar=c(2,2.5,1,0)+.5, mgp=c(1.6,.6,0))
acf(soi, 48, xlab=' ', main="", panel.first=grid(lty=1))
mtext(side=3, 'Southern Oscillation Index', font=2)
acf(rec, 48, xlab='',main="", panel.first=grid(lty=1))
mtext(side=3, 'Recruitment', font=2)
ccf(soi, rec, 48, xlab='LAG', main="", ylab="CCF",panel.first=grid(lty=1))
mtext(side=3, 'SOI vs Recruitment', font=2)
dev.off()


########
pdf(file="prewhiten.pdf",width=7.5,height=6)
set.seed(1492)
num=120; t=1:num
X = ts(2*cos(2*pi*t/12) + rnorm(num), freq=12)
Y = ts(2*cos(2*pi*(t+5)/12) + rnorm(num), freq=12)
Yw = resid( lm(Y~ cos(2*pi*t/12) + sin(2*pi*t/12), na.action=NULL) )
par(mfrow=c(3,2), mgp=c(1.6,.6,0), mar=c(3,3,1,1) )
plot(X, type='n'); grid(lty=1, col=gray(.9)); lines(X)
plot(Y, type='n'); grid(lty=1, col=gray(.9)); lines(Y)
acf(X,48, ylab='ACF(X)', panel.first=grid(lty=1, col=gray(.9)))
acf(Y,48, ylab='ACF(Y)', panel.first=grid(lty=1, col=gray(.9)))
ccf(X,Y,24, ylab='CCF(X,Y)', panel.first=grid(lty=1, col=gray(.9)))
ccf(X,Yw,24, ylab='CCF(X,Yw)', ylim=c(-.6,.6), panel.first=grid(lty=1, col=gray(.9)))
dev.off()
