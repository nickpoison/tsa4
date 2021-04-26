# chapter 7 example 10
require(astsa)

P=1:1024
S=P+1024
N = 8 # = N1 = N2 = number of series in each grp
n = 1024 # length of each series



eq.P = as.ts(eqexp[P,1:8]); eq.S = as.ts(eqexp[S,1:8]) #EQ
eq.m = cbind(rowMeans(eq.P), rowMeans(eq.S))
ex.P = as.ts(eqexp[P,9:16]); ex.S = as.ts(eqexp[S,9:16]) #EX
ex.m = cbind(rowMeans(ex.P), rowMeans(ex.S))


m.diff = mvfft(eq.m-ex.m)/sqrt(n)
eq.Pf = mvfft(eq.P-eq.m[,1])/sqrt(n); eq.Sf = mvfft(eq.S-eq.m[,2])/sqrt(n)
ex.Pf = mvfft(ex.P-ex.m[,1])/sqrt(n); ex.Sf = mvfft(ex.S-ex.m[,2])/sqrt(n)

fv11=rowSums(eq.Pf*Conj(eq.Pf))+rowSums(ex.Pf*Conj(ex.Pf))/(2*(N-1))
fv12=rowSums(eq.Pf*Conj(eq.Sf))+rowSums(ex.Pf*Conj(ex.Sf))/(2*(N-1))
fv22=rowSums(eq.Sf*Conj(eq.Sf))+rowSums(ex.Sf*Conj(ex.Sf))/(2*(N-1))
fv21=Conj(fv12)
T2=rep(NA,512)

for (k  in 1:512){
fvk= matrix(c(fv11[k],fv21[k],fv12[k],fv22[k]),2,2)
 dk=as.matrix(m.diff[k,])
 T2[k] = Re((N/2)*Conj(t(dk))%*%solve(fvk,dk))}

p=2
eF= T2*(2*p*(N-1))/(2*N-p-1)

pdf(file="fig_ex7_10.pdf", width=7.7, height=5) 
par(mfrow=c(2,2), mar=c(3,3,2,1),  mgp = c(1.6,.6,0), cex.main=1.1)
freq=40*(0:511)/n  # Hz

plot(freq, eF, type="l", xlab="Frequency (Hz)", ylab="F Statistic", main="Equal Means", panel.first=grid(lty=1))
abline(h=qf(.999, 2*p, 2*(2*N-p-1)))


# equal P 
m=10; L=2*m+1
kd=kernel("daniell",m)
u = Re(rowSums(eq.Pf*Conj(eq.Pf))/(N-1))
feq.P = kernapply(u, kd, circular=TRUE)
u = Re(rowSums(ex.Pf*Conj(ex.Pf))/(N-1))
fex.P =	kernapply(u, kd, circular=TRUE)
plot(freq, feq.P[1:512]/fex.P[1:512], type="l", xlab="Frequency (Hz)", ylab="F Statistic", main="Equal P-Spectra", panel.first=grid(lty=1))
abline(h=qf(.9999, 2*L*(N-1),  2*L*(N-1)))

# equal S 
u = Re(rowSums(eq.Sf*Conj(eq.Sf))/(N-1))
feq.S = kernapply(u, kd, circular=TRUE)
u = Re(rowSums(ex.Sf*Conj(ex.Sf))/(N-1))
fex.S =	kernapply(u, kd, circular=TRUE)
plot(freq, feq.S[1:512]/fex.S[1:512], type="l", xlab="Frequency (Hz)", ylab="F Statistic", main="Equal S-Spectra", panel.first=grid(lty=1))
abline(h=qf(.9999, 2*L*(N-1),  2*L*(N-1)))


# equal matrices
u = rowSums(eq.Pf*Conj(eq.Sf))/(N-1)
feq.PS = kernapply(u, kd, circular=TRUE)
u = rowSums(ex.Pf*Conj(ex.Sf)/(N-1))
fex.PS =	kernapply(u, kd, circular=TRUE)

fv11 = kernapply(fv11, kd, circular=TRUE)
fv22 = kernapply(fv22, kd, circular=TRUE)
fv12 = kernapply(fv12, kd, circular=TRUE)


Mi=L*(N-1); M=2*Mi; p=2
TS=rep(NA,512)
for (k  in 1:512){
det.feq.k= Re(feq.P[k]*feq.S[k] - feq.PS[k]*Conj(feq.PS[k]))
det.fex.k= Re(fex.P[k]*fex.S[k] - fex.PS[k]*Conj(fex.PS[k]))
det.fv.k = Re(fv11[k]*fv22[k] - fv12[k]*Conj(fv12[k]))
log.n1 = log(M)*(M*p)
log.d1 = log(Mi)*(2*Mi*p)
log.n2 = log(Mi)*2 +log(det.feq.k)*Mi + log(det.fex.k)*Mi  
log.d2 = (log(M)+log(det.fv.k))*M
r= 1 - ((p+1)*(p-1)/6*p*(2-1))*(2/Mi - 1/M)
TS[k] = -2*r*(log.n1+log.n2-log.d1-log.d2)   }


plot(freq, TS, type="l", xlab="Frequency (Hz)", ylab="Chi-Sq Statistic", main="Equal Spectral Matrices", panel.first=grid(lty=1))
abline(h=qchisq(.9999, p^2))
dev.off()


