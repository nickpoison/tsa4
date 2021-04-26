gr = diff(log(ts(econ5, start=1948, frequency=4)))
name=colnames(econ5)

pdf(file="fig_ex7_16a.pdf", width=7.5, height=5) 
par(mfrow=c(5,1), mar = c(0, 3, 0, 2.5), oma=c(3,0,2,0), mgp=c(1.6,.6,0), cex.lab=1.5)
plot.ts(100*gr[,1], type='n', ylab=name[1], xaxt="no")
grid(lty=1);  lines(100*gr[,1])
title(main='Growth Rates (%)', outer=TRUE)
plot.ts(100*gr[,2], type='n', ylab=name[2], xaxt="no", yaxt='no')
grid(lty=1); axis(4); lines(100*gr[,2]) 
plot.ts(100*gr[,3], type='n', ylab=name[3], xaxt="no")
grid(lty=1);  lines(100*gr[,3])
plot.ts(100*gr[,4], type='n', ylab=name[4], xaxt='no',yaxt='no')
grid(lty=1); axis(4); lines(100*gr[,4])  
plot.ts(100*gr[,5], type='n', ylab=name[5])
grid(lty=1); lines(100*gr[,5]) 
title(xlab="Time", outer=TRUE)
dev.off()

#plot(100*gr, main="Growth Rates (%)")



gr = ts(apply(gr,2,scale), freq=4)   # scaling strips ts attributes

L= c(7,7)   
gr.spec = mvspec(gr, spans=L, demean=FALSE, detrend=FALSE, taper=.25, plot=FALSE)
#plot(kernel("modified.daniell", L))  # view the kernel - not shown 

pdf(file="fig_ex7_16b.pdf", width=7.5, height=3) 
par(mar=c(2.75,3,2,1), mgp=c(1.6,.6,0))
plot(gr.spec, log="no", col=c(4:1,6), main="Individual Spectra", lty=(1:6)[-3], lwd=2, xlab="frequency")
grid()
legend("topright", colnames(econ5), lty=(1:6)[-3], col=c(4:1,6), lwd=2, bg='white')  
dev.off()





pdf(file="fig_ex7_16c.pdf") 
plot.spec.coherency(gr.spec, ci=NA,  main='')# main="Squared Coherencies")
dev.off()


# PCs
n.freq = length(gr.spec$freq)
lam = matrix(0,n.freq,5)
for (k in 1:n.freq) lam[k,] = eigen(gr.spec$fxx[,,k], symmetric=TRUE, only.values=TRUE)$values 

#dev.new()  # note peaks at .25 = 1 cycle in 4 quarters and .125 = 1 cycle 8 quarters
pdf(file="fig_ex7_16d.pdf", width=7.25, height=3.5)   
par(mfrow=c(2,1), mar=c(2,2,1,.5)+.5, mgp=c(1.4,.6,0), cex=.9, cex.main=1)  
plot(gr.spec$freq, lam[,1], type="n", ylab="", xlab="Frequency", main="First Eigenvalue")
grid(lty=1); lines(gr.spec$freq, lam[,1])
abline(v=.25, lty=6, col=4)
plot(gr.spec$freq, lam[,2], type="n", ylab="", xlab="Frequency", main="Second Eigenvalue")
grid(lty=1); lines(gr.spec$freq, lam[,2])
abline(v=.125, lty=6, col=4) 
dev.off()

e.vec1 = eigen(gr.spec$fxx[,,10], symmetric=TRUE)$vectors[,1] 
e.vec2 =  eigen(gr.spec$fxx[,,5], symmetric=TRUE)$vectors[,2]
round(Mod(e.vec1), 2)
round(Mod(e.vec2), 2) 
