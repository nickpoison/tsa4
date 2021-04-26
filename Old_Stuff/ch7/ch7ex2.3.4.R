require('astsa')

##################
########  plot of beam data ##################
#infrasound
attach(beamd)

pdf(file="infrasound.pdf",width=7.7,height=4)
tau = rep(0,3)
u = ccf(sensor1, sensor2, plot=FALSE)
tau[1] = u$lag[which.max(u$acf)]    #  17
u = ccf(sensor3, sensor2, plot=FALSE)
tau[3] = u$lag[which.max(u$acf)]    # -22
Y = ts.union(lag(sensor1,tau[1]), lag(sensor2, tau[2]), lag(sensor3, tau[3]))
Y = ts.union(Y, rowMeans(Y))
colnames(Y) = c('sensor1', 'sensor2', 'sensor3', 'beam')
Time = time(Y)
par(mfrow=c(4,1), mar=c(0, 3.1, 0, 1.1), oma=c(2.75, 0, 2.5, 0), mgp=c(1.6,.6,0))
plot(Time, Y[,1], ylab='sensor1', xaxt="no", type='n')
grid(lty=1); lines(Y[,1])
title(main="Infrasonic Signals and Beam", outer=TRUE)
plot(Time, Y[,2], ylab='sensor2', xaxt="no",  type='n')
grid(lty=1); lines(Y[,2]) 
plot(Time,Y[,3], ylab='sensor3', xaxt="no", type='n')
grid(lty=1); lines(Y[,3])
plot(Time, Y[,4], ylab='beam', type='n')  
grid(lty=1); lines(Y[,4]) 
title(xlab="Time", outer=TRUE)
dev.off()
#####################################
###################################
#  now start the example

# attach(beamd, warn.conflicts = FALSE)
 L=9
 fdr=.001
 N=3
 Y = cbind(beamd, beam=rowMeans(beamd))
 n=nextn(nrow(Y))
 
 
 
 Y.fft= mvfft(as.ts(Y))/sqrt(n)
 Df = Y.fft[,1:3]         # fft of the data
 Bf = Y.fft[,4]           # beam fft
 
 
 # Raw Signal Spectrum 
 ssr = N*Re(Bf*Conj(Bf))
  
 # Raw Error Spectrum
 sse = Re(rowSums(Df*Conj(Df)))  - ssr
 
 # Smooth
 SSE=filter(sse, sides=2, filter=rep(1/L,L), circular=TRUE)
 SSR=filter(ssr, sides=2, filter=rep(1/L,L), circular=TRUE)
 SST=SSE+SSR
 
      ###### figure:  infpow1 
 
  pdf(file="infpow1.pdf",width=7.7,height=4)
  par( mfrow=c(2,1), mar=c(1.9,2,.5,0)+.5, mgp=c(1.35,.5,0), cex.main=1)
 
# par(mfrow=c(2,1), mar=c(4,4,2,1)+.1)
 Fr = 1:(n-1)/n
 nFr=1:200       # number of freqs to plot
 plot(Fr[nFr], log(SST[nFr]), type="l", ylab="log Power", xlab="", main="Sum of Squares", panel.first=grid(lty=1))
 lines(Fr[nFr], log(SSE[nFr]), type="l", lty=2)
 # x=c(Fr[nFr], rev(Fr[nFr]))   # maybe add this
 # y=c(SST[nFr], rev(SSE[nFr]))
 # polygon(x,y,col="gray")
 eF=(N-1)*SSR/SSE;
 df1=2*L
 df2=2*L*(N-1)
 # Compute F-value for false discovery probability of fdr
 p = pf(eF, df1, df2, lower=FALSE)
 pID=FDR(p,fdr)
 Fq=qf(1-fdr,df1,df2)
 plot(Fr[nFr], eF[nFr], type="l", ylab="F-statistic", xlab="Frequency", main="F Statistic", panel.first=grid(lty=1))
 abline(h=c(Fq, eF[pID]), lty=1:2)
dev.off()





  