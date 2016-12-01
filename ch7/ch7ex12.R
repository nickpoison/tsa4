require(astsa)
#  2 graphs here


P=1:1024
S=P+1024
p.dim = 2
n =1024
eq = as.ts(eqexp[,1:8])
ex = as.ts(eqexp[,9:16])
nz = as.ts(eqexp[,17])
f.eq <- array(dim=c(8,2,2,512)) -> f.ex 
f.NZ = array(dim=c(2,2,512))
#### for determinants of 2x2 complex matrices
det.c = function(mat){return(Re(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1]))}


L = c(15,13,5)  # for smoothing 
for (i in 1:8){     # compute spectral matrices
	f.eq[i,,,] = mvspec(cbind(eq[P,i],eq[S,i]), spans=L, taper=.5)$fxx
	f.ex[i,,,] = mvspec(cbind(ex[P,i],ex[S,i]), spans=L, taper=.5)$fxx}
 f.NZ = mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5)$fxx	

bndwidth = mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5)$bandwidth*sqrt(12)*40  # ~.75 Hz

fhat.eq = apply(f.eq, 2:4, mean)  # average spectra
fhat.ex = apply(f.ex, 2:4, mean)
#plot the spectra
#layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

pdf(file="fig_ex7_12b.pdf", width=7.7, height=4) 
par(mfrow=c(2,2), mar=c(2,2,.5,.5)+.5,  mgp = c(1.4,.6,0))
Fr = 40*(1:512)/n
plot(Fr, Re(fhat.eq[1,1,]), type="l", xlab="Frequency (Hz)", ylab="", panel.first=grid(lty=1))
plot(Fr, Re(fhat.eq[2,2,]), type="l", xlab="Frequency (Hz)", ylab="", panel.first=grid(lty=1))
plot(Fr, Re(fhat.ex[1,1,]), type="l", xlab="Frequency (Hz)", ylab="", panel.first=grid(lty=1))
plot(Fr, Re(fhat.ex[2,2,]), type="l", xlab="Frequency (Hz)", ylab="", panel.first=grid(lty=1))
mtext("Average P-spectra", side=3, line=-1, adj=.2, outer=TRUE, cex=1)
mtext("Earthquakes", side=2, line=-1, adj=.8,  outer=TRUE)
mtext("Average S-spectra", side=3, line=-1, adj=.82, outer=TRUE, cex=1)
mtext("Explosions", side=2, line=-1, adj=.2, outer=TRUE)
 par(fig = c(.725, 1, .725, 1), new = TRUE)
 ker = kernel("modified.daniell", L)$coef; ker = c(rev(ker),ker[-1])
 plot((-33:33)/40, ker, type="l", ylab="",xlab="",cex.axis=.7, yaxt='no')#c(0,.04,2) )
dev.off()


# L = c(15,13,5)
# plot(kernel("modified.daniell", L))
# mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5)$bandwidth*sqrt(12)*40 


# choose alpha
 Balpha = rep(0,19) 
 for (i in 1:19){
 	alf=i/20
 for (k in 1:256) {  	
 Balpha[i] =  Balpha[i] + Re(log(det.c(alf*fhat.ex[,,k]+(1-alf)*fhat.eq[,,k])/det.c(fhat.eq[,,k])) 
               - alf*log(det.c(fhat.ex[,,k])/det.c(fhat.eq[,,k])))} }
alf = which.max(Balpha)/20   # = .4           


KLDiff <- rep(0,17) -> BDiff
KLeq <- rep(0,17) -> KLex -> Beq -> Bex


# calculate information criteria for eqs and exs
for (i in 1:17){
	if (i <= 8) f0 = f.eq[i,,,]
	if (i > 8 & i <= 16) f0 = f.ex[i-8,,,]
	if (i == 17) f0 = f.NZ
for (k in 1:256) {    # only use freqs out to .25
    tr1 = Re(sum(diag(solve(fhat.eq[,,k],f0[,,k]))))
   KLeq[i] = KLeq[i] + tr1 + log(det.c(fhat.eq[,,k])) - log(det.c(f0[,,k]))
   Beq[i] =  Beq[i] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.eq[,,k])/det.c(fhat.eq[,,k])) - alf*log(det.c(f0[,,k])/det.c(fhat.eq[,,k])))
    tr1 = Re(sum(diag(solve(fhat.ex[,,k],f0[,,k]))))
   KLex[i] = KLex[i] + tr1 +  log(det.c(fhat.ex[,,k])) - log(det.c(f0[,,k]))
   Bex[i] = Bex[i] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.ex[,,k])/det.c(fhat.ex[,,k])) - alf*log(det.c(f0[,,k])/det.c(fhat.ex[,,k])))
  }
KLDiff[i] = (KLeq[i] - KLex[i])/n 
BDiff[i] =  (Beq[i] - Bex[i])/(2*n)}

x.max = max(KLDiff)+.1
x.min = min(KLDiff)-.1
y.max = max(BDiff)+.01
y.min = min(BDiff)-.01

pdf(file="fig_ex7_12.pdf", width=7.5, height=5) 
par(mar=c(2,2,.5,.5)+.5,  mgp = c(1.4,.6,0))                               
plot(KLDiff[9:16], BDiff[9:16],  type="n",   xlim=c(x.min,x.max), ylim=c(y.min,y.max),  
       xlab="Kullback-Leibler Difference", ylab="Chernoff Difference", main="Classification Based on Chernoff and K-L Distances")
abline(h=0, v=0, lty=6, col=8)
points(KLDiff[9:16], BDiff[9:16], pch=6, cex=1.1, lwd=2, col=6)	   
points(KLDiff[1:8], BDiff[1:8], pch=8, cex=1.1, lwd=2, col=4)
points(KLDiff[17], BDiff[17],  pch=3, cex=1.1, lwd=2, col=rgb(0,.6,.2))
legend("topleft", legend=c("EQ", "EX", "NZ"), pch=c(8,6,3), pt.lwd=2, col=c(4,6,rgb(0,.6,.2)))
text(KLDiff[-c(1,2,3,7,14)]-.075, BDiff[-c(1,2,3,7,14)], label=names(eqexp[-c(1,2,3,7,14)]), cex=.7)
text(KLDiff[c(1,2,3,7,14)]+.075, BDiff[c(1,2,3,7,14)], label=names(eqexp[c(1,2,3,7,14)]), cex=.7)
dev.off()

#  cross-classify (loo)
#  basically same result as above



#
#KLDiff <- rep(0,16)  -> BDiff
#KLeq <- rep(0,16) -> KLex  -> Beq -> Bex
#
#for (los in 1:16){
#	
#if (los <= 8) {f0=f.eq[los,,,]; apply(f.eq[-los,,,], 2:4, mean); apply(f.ex, 2:4, mean)  }
#if (los > 8)  {f0=f.ex[los-8,,,]; apply(f.ex[-(los-8),,,], 2:4, mean); apply(f.eq, 2:4, mean)}
#
#
#
#for (k in 1:256) {    # only use freqs out to .25
#  tr1 = Re(sum(diag(solve(fhat.eq[,,k],f0[,,k]))))
#   KLeq[los] = KLeq[los] + tr1 + log(det.c(fhat.eq[,,k])) - log(det.c(f0[,,k]))
#    Beq[los] =  Beq[los] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.eq[,,k])/det.c(fhat.eq[,,k])) - alf*log(det.c(f0[,,k])/det.c(fhat.eq[,,k])))
#  tr1 = Re(sum(diag(solve(fhat.ex[,,k],f0[,,k]))))
#   KLex[los] = KLex[los] + tr1 +  log(det.c(fhat.ex[,,k])) - log(det.c(f0[,,k]))
#    Bex[los] = Bex[los] + Re(log(det.c(alf*f0[,,k]+(1-alf)*fhat.ex[,,k])/det.c(fhat.ex[,,k])) - alf*log(det.c(f0[,,k])/det.c(fhat.ex[,,k])))
#  }
#KLDiff[los] = (KLeq[los] - KLex[los])/n 
#BDiff[los] =  (Beq[los] - Bex[los])/(2*n)
#}
#
#plot(KLDiff[9:16], BDiff[9:16],  type="p", pch=6, xlim=c(x.min,x.max), ylim=c(y.min,y.max), cex=1.1, lwd=2,
#       xlab="Kullback-Leibler Difference", ylab="Chernoff Difference", main="Cross-Validation")
#points(KLDiff[1:8], BDiff[1:8], pch=8, cex=1.1, lwd=2)
#legend("topleft", legend=c("EQ", "EX", "NZ"), pch=c(8,6,3), pt.lwd=2, cex=1.1)
#abline(h=0,lty=2);abline(v=0,lty=2)
#text(KLDiff[-c(1,2,3,7,14)]-.075, BDiff[-c(1,2,3,7,14)], label=names(eqexp[-c(1,2,3,7,14)]), cex=.7)
#text(KLDiff[c(1,2,3,7,14)]+.075, BDiff[c(1,2,3,7,14)], label=names(eqexp[c(1,2,3,7,14)]), cex=.7)
#
#






