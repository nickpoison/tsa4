# fmri1 is the data col 2-9 is data, col 1 is time

n = 128
Per = abs(mvfft(fmri1[,-1]))^2/n

pdf(file="fig_ex7_14a.pdf", width=7.5, height=4) 
par(mfrow=c(2,4), mar=c(3,2,2,1),  mgp = c(1.6,.6,0), oma=c(0,1,0,0))
for (i in 1:8){
 plot(0:20, Per[1:21,i], type="n", xaxt='n', ylim=c(0,8), main=colnames(fmri1)[i+1], xlab="Cycles", ylab="",xaxp=c(0,20,5))
 axis(1, seq(0,20,by=4))
 grid(lty=1, nx=NA, ny=NULL)
 abline(v=seq(0,60,by=4), col='lightgray', lty=1)
 lines(0:20, Per[1:21,i]) }
mtext("Periodogram", side=2, line=-.3, outer=TRUE, adj=c(.2,.8))
dev.off()


fxx = mvspec(fmri1[,-1], kernel("daniell", c(1,1)), taper=.5, plot=FALSE)$fxx
l.val = rep(NA,64)
for (k in 1:64) {
u = eigen(fxx[,,k], symmetric=TRUE, only.values = TRUE)
l.val[k] = u$values[1]}


pdf(file="fig_ex7_14b.pdf", width=7, height=3) 
par(mar=c(2.25,2,.5,.5)+.5,  mgp = c(1.6,.6,0))  
plot(l.val, type="n", xaxt='n',xlab="Cycles (Frequency x 128)", ylab="First Principal Component")
axis(1, seq(4,60,by=8))
grid(lty=1, nx=NA, ny=NULL)
abline(v=seq(4,60,by=8), col='lightgray', lty=1)
lines(l.val)
dev.off()


# at freq k=4
u = eigen(fxx[,,4], symmetric=TRUE)
lam=u$values 
evec=u$vectors
lam[1]/sum(lam) # % of variance explained
sig.e1 = matrix(0,8,8)
for (l in 2:5){  # last 3 evs are 0
 sig.e1= sig.e1 + lam[l]*evec[,l]%*%Conj(t(evec[,l]))/(lam[1]-lam[l])^2}
sig.e1 = Re(sig.e1)*lam[1]*sum(kernel("daniell", c(1,1))$coef^2)
p.val = round(pchisq(2*abs(evec[,1])^2/diag(sig.e1), 2,  lower.tail = FALSE),3)
cbind(colnames(fmri1)[-1], abs(evec[,1]), p.val) # print table values



# for 7.15
bhat = sqrt(lam[1])*evec[,1]
Dhat = Re(diag(fxx[,,4] - bhat%*%Conj(t(bhat))))
res = Mod(fxx[,,4] - Dhat - bhat%*%Conj(t(bhat)))


#  0.00 & 1.73&  3.88&  3.61&  0.88& 2.04&  1.60&  2.81
#  2.41 & 0.00&  1.17&  3.77&  1.49& 5.58&  3.68&  4.21
#  8.49 & 5.34&  0.00&  2.94&  7.58&10.91&  8.36& 10.64
# 12.65 &11.84&  6.12&  0.00& 12.56&14.64& 13.34& 16.10
#  0.32 & 0.29&  2.10&  2.01&  0.00& 1.18&  2.01&  1.18
# 10.34 &16.69& 17.09& 15.94& 13.49& 0.00&  5.78& 14.74
#  5.71 & 8.51&  8.94& 10.18&  7.56& 0.97&  0.00&  8.66
#  6.25 & 8.00& 10.31& 10.69&  5.95& 8.69&  7.64&  0.00
