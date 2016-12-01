require(astsa)


#attach(beamd)
L = 9; M = 100; M2 = M/2; N = 3   
Y = cbind(beamd, beam <- rowMeans(beamd))
n = nextn(nrow(Y));  n.freq=n/2 
Y[,1:3] = Y[,1:3]-Y[,4]  
Y.fft= mvfft(as.ts(Y))/sqrt(n)
Ef = Y.fft[,1:3]     # fft of the error
Bf = Y.fft[,4]       # beam fft
ssr = N*Re(Bf*Conj(Bf))        # Raw Signal Spectrum 
sse = Re(rowSums(Ef*Conj(Ef))) # Raw Error Spectrum
# Smooth
SSE=filter(sse, sides=2, filter=rep(1/L,L), circular=TRUE)
SSR=filter(ssr, sides=2, filter=rep(1/L,L), circular=TRUE)
# Estimate Signal and Noise Spectra
fv = SSE/(L*(N-1))           # Equation (7.77)
fb = (SSR-SSE/(N-1))/(L*N)   # Equation (7.78)
fb[fb<0] = 0 
H0 = N*fb/(fv+N*fb)     
H0[ceiling(.04*n):n]=0       # zero out H0 above .04
# Extend components to make it a valid transform
H0 = c(H0[1:n.freq],rev(H0[2:(n.freq+1)]))  
h0 = Re(fft(H0,inverse = TRUE))       # Impulse Response 
h0 = c(rev(h0[2:(M2+1)]),h0[1:(M2+1)])      # center it
h1 = spec.taper(h0, p = .5)                  # taper it
k1 = h1/sum(h1)                          # normalize it
f.beam = filter(Y[,4], filter=k1, sides=2) # filter it
# Graphics  

pdf(file="fig_ex7_6.pdf",width=7.6,height=5)
nFr = 1:50                  # number of freqs displayed
Fr = (nFr-1)/n              # frequencies 
layout(matrix(c(1,2,4, 1,3,4), nc=2)) 
par(mar=c(2,2,0,0)+.5, mgp=c(1.45,.6,0))
plot(10*Fr, fb[nFr], type="l", ylab="Power", xlab="Frequency (Hz)", panel.first=grid(lty=1)) 
 lines(10*Fr, fv[nFr], lty=2)   
 text(.24,5, "(a)", cex=1.2)
plot(10*Fr,H0[nFr],type="l",ylab="Frequency Response",xlab="Frequency (Hz)", panel.first=grid(lty=1)) 
 text(.23,.84, "(b)", cex=1.2)
plot(-M2:M2,k1,type="l",ylab="Impulse Response",xlab="Index", lwd=1.5, panel.first=grid(lty=1))
 text(45,.022, "(c)", cex=1.2)
plot(f.beam, type='n', xlab='Time', ylab="beam",  ylim=c(-5,3))  
grid(lty=1); lines(f.beam)
lines(beam, lty=1, col=rgb(0,0,.6, alpha=.2), lwd=5)
 text(2040,2, "(d)", cex=1.2)
dev.off()
