# Example R code to calculate the SPECTRAL ENVELOPE for a continuous time series 
# *** Must source mvspec.R first ***           
# 
# THIS IS NOW FOR NYSE

# fig_ex7_19 and b

#gnpgr=diff(log(gnp))
#u=resid(arima(gnpgr, order=c(0,0,2)))

u     = astsa::nyse
#u=scan("gnpres.dat")              # data (use the residuals from example 3.35 on gnp)
x=cbind(u, abs(u), u^2)           # possible transformations (identity, absolute value, square)
#
Var = var(x)                     # var-cov matrix 
xspec = mvspec(x, spans=c(5,3), taper=.5, plot=FALSE)  # spectral matrices are an array called fxx  was 5,3
fxxr = Re(xspec$fxx)             # fxxr is real(fxx) 
#
#-----   compute Q = Var^-1/2   -----#
ev = eigen(Var)  
Q = ev$vectors%*%diag(1/sqrt(ev$values))%*%t(ev$vectors)             
#
#--- compute spec env and scale vectors ---#    
num = xspec$n.used               # sample size used for FFT
nfreq = length(xspec$freq)       # number of freqs used   
specenv = matrix(0,nfreq,1)      # initialize the spec envelope
beta=matrix(0,nfreq,3)           # initialize the scale vectors
 for (k in 1:nfreq){ 
  ev = eigen(2*Q%*%fxxr[,,k]%*%Q/num)  # get evalues of normalized spectral matrix at freq k/n
  specenv[k] = ev$values[1]            # spec env at freq k/n is max evalue
  b = Q%*%ev$vectors[,1]               # beta at freq k/n 
  beta[k,] = b/b[1]                    # first coef is always 1
  }
#
#--- output and graphics ---# 
pdf(file="fig_ex7_19.pdf", width=7.5, height=3) 
par(mar=c(2.5,2.75,.5,.5), mgp=c(1.5,.6,0))
frequency = xspec$freq
plot(frequency, 100*specenv, type="l", ylab="Spectral Envelope (%)", panel.first=grid(lty=1))
 ## add significance threshold to plot ##
 m=xspec$kernel$m
 etainv=sqrt(sum(xspec$kernel[-m:m]^2))
thresh=100*(2/num)*exp(qnorm(.9999)*etainv)*matrix(1,nfreq,1)
lines(frequency,thresh, lty="dashed", col="blue")
 dev.off()
#--  details  --# 
output = cbind(frequency, specenv, beta)
colnames(output)=c("freq","specenv","x", "|x|", "x^2")
round(output,4) 

b = sign(b[2])*output[2,3:5]

pdf(file="fig_ex7_19b.pdf", width=7.5, height=3) 
par(mar=c(2.5,2.5,.5,.5), mgp=c(1.5,.6,0))
# plot transform
g = function(x) {b[1]*x+b[2]*abs(x)+b[3]*x^2}
curve(g, -.2, .2, panel.first=grid(lty=1))#, ylim=c(0,.012))
g2 = function(x) {b[2]*abs(x)}
curve(g2, -.2,.2, add=TRUE, lty=6, col=4)
dev.off()

# the peak is at .001 
head(output)
#       freq     specenv x       |x|        x^2
#[1,] 0.0005 0.014500410 1 201.71523 -628.98816
#[2,] 0.0010 0.014868595 1 -60.11826  201.32262  <---
#[3,] 0.0015 0.014466993 1 -37.80450  129.39620
#[4,] 0.0020 0.012120906 1 -34.69786  111.07664
#[5,] 0.0025 0.008913294 1 -38.87482   99.70634
#[6,] 0.0030 0.007044304 1 -60.10531  106.31674












                                                              # uncomment to write to screen
# write.table(output, file="output.txt", quote=FALSE, row.names=FALSE)  # uncomment to write to file
