#attach(climhyd, warn.conflicts = FALSE)
# plot.ts(climhyd)  # figure 7.3   see ch7_intro.r for this plot


Y = climhyd         # Y will hold transformed series
Y[,6] = log(Y[,6])  # log inflow
Y[,5] = sqrt(Y[,5])  # sqrt precip                               



L = 25
M = 100        # even
alpha = .001
fdr = alpha
nq = 2   # number of inputs  (Temp and Precip)




#-- Spectral Matrix --#
 Yspec = mvspec(Y, spans=L, kernel="daniell", detrend=TRUE, demean=FALSE, taper=.1, plot=FALSE)
 #fYY = Yspec$fxx           # the spectral matrix
 n = Yspec$n.used          # effective sample size
 Fr = Yspec$freq           # fundamental freqs 
 n.freq = length(Fr)       # number of frequencies
 Yspec$bandwidth*sqrt(12)  # = 0.050.. is the bandwidth 

#-- Coherencies --#  (see sec 4.7 also)
 Fq = qf(1-alpha, 2, L-2)
 cn = Fq/(L-1+Fq)
 plt.name=c("(a)","(b)","(c)","(d)","(e)","(f)")
 
 pdf(file="fig_ex7_1a.pdf", width=7.7, height=5)
 par(mfrow=c(2,3), cex.lab=1, mgp=c(1.6,.6,0), mar=c(2.5,2.5,2.5,0)+.5) 
 # in R the coherencies are listed as 1,2,..., 15=choose(6,2) 
 for (i in 11:15){
  plot(Fr,Yspec$coh[,i], type="l", ylab="Sq Coherence", xlab="Frequency", panel.first=grid(lty=1),
         ylim=c(0,1), main=c("Inflow with", names(climhyd[i-10])))
  abline(h = cn); text(.45,.98, plt.name[i-10], cex=1.2)  } 
 # mult coherency 
 coh.15 = stoch.reg(Y, cols.full = c(1,5), cols.red = NULL, alpha, L, M, plot.which = "NULL")  
  plot(Fr,coh.15$coh, type="l", ylab="Sq Coherence", xlab="Frequency", panel.first=grid(lty=1), ylim=c(0,1))
  abline(h = cn); text(.45,.98, plt.name[6], cex=1.2) 
 title(main = c("Inflow with", "Temp and Precip"))
dev.off() 





#######  Partial F and Regr Coefs

numer.df = 2*nq
denom.df = Yspec$df-2*nq
 out.15 = stoch.reg(Y, cols.full = c(1,5), cols.red = 5, alpha, L, M, plot.which = "F.stat")
dev.off()
 pdf(file="fig_ex7_1b.pdf", width=7.7, height=5)
 par(mfrow=c(3,1), mar=c(2.5,2.5,1,0)+.5, mgp = c(1.5,0.4,0), cex.lab=1.2, mgp=c(1.6,.5,0)) 
 plot(Fr, out.15$eF,  type="l", ylab="F", xlab="Frequency", panel.first=grid(lty=1))
eF = out.15$eF
pvals = pf(eF, numer.df, denom.df, lower.tail = FALSE)
pID = FDR(pvals, fdr);  abline(h=c(eF[pID]), lty=2)
abline(h=qf(.001, numer.df, denom.df, lower.tail = FALSE) )
title(main = "Partial F Statistic")
# Regression Coefficients
S = seq(from = -M/2+1, to = M/2 - 1, length = M-1)
plot(S, coh.15$Betahat[,1], type = "h", xlab = "", ylab = names(climhyd[1]), ylim = c(-.025, .055), panel.first=grid(lty=1), lwd=2)
abline(h=0); title(main = "Impulse Response Functions")
plot(S, coh.15$Betahat[,2], type = "h", xlab = "Index", ylab = names(climhyd[5]), ylim = c(-.015, .055), panel.first=grid(lty=1), lwd=2)
abline(h=0)
dev.off()



