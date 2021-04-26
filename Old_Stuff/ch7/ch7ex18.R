# fig_ex7_18a and b

require(astsa)

u = factor(bnrf1ebv)   # first, input the data as factors and then 
x = model.matrix(~u-1)[,1:3]     # make an indicator matrix
 # x = x[1:1000,]                  # select subsequence if desired
Var = var(x)                     # var-cov matrix 
xspec = mvspec(x, spans=c(7,7),plot=FALSE)  # spectral matrices are an array called fxx
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
  beta[k,] = b/sqrt(sum(b^2))          # helps to normalize beta
  }
#
#--- output and graphics ---# 

pdf(file="fig_ex7_18a.pdf", width=7.5, height=3) 
par(mar=c(2.5,2.75,.5,.5), mgp=c(1.5,.6,0))
frequency = (0:(nfreq-1))/num
plot(frequency, 100*specenv, type="l", ylab="Spectral Envelope (%)", panel.first=grid(lty=1))
 ## add significance threshold to plot ##
 m=xspec$kernel$m
 etainv=sqrt(sum(xspec$kernel[-m:m]^2))
thresh=100*(2/num)*exp(qnorm(.9999)*etainv)
abline(h=thresh, lty=6, col="blue" )
#lines(frequency,thresh, lty="dashed", col="blue")
dev.off()


#--  details  --# 
output = cbind(frequency, specenv, beta)
colnames(output)=c("freq","specenv","A", "C", "G")
round(output,3)                                                                # uncomment to write to screen
# write.table(output, file="output.txt", quote=FALSE, row.names=FALSE)  # uncomment to write to file

####################################
#  dynamic part - don't show in text
output=vector("list",4)
colnames(output)=c("freq","specenv","A", "C", "G")


pdf(file="fig_ex7_18b.pdf", width=7.5, height=5) 
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(1.6,.6,0))
for (j in 1:4){
  if (j == 1) {ind=1:1000; dog="(a)"}
  if (j == 2) {ind=1001:2000; dog="(b)"}
  if (j == 3) {ind=2001:3000; dog="(c)"}
  if (j == 4) {ind=3001:length(bnrf1ebv); dog="(d)"} 	
xx = x[ind,]                  # select subsequence if desired
Var = var(xx)                     # var-cov matrix 
xspec = mvspec(xx, spans=c(7,7),plot=FALSE)  # spectral matrices are an array called fxx
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
  beta[k,] = b/sqrt(sum(b^2))          # helps to normalize beta
  }
#
#--- output and graphics ---# 


frequency = (0:(nfreq-1))/num
plot(frequency, 100*specenv, type="l", ylab="Spectral Envelope (%)", ylim=c(0,1.3), panel.first=grid(lty=1))
 ## add significance threshold to plot ##
 m=xspec$kernel$m
 etainv=sqrt(sum(xspec$kernel[-m:m]^2))
thresh=100*(2/num)*exp(qnorm(.9999)*etainv)
abline(h=thresh, lty=6, col=4)
text(.475,1.25,dog)

 # output[[j]] = cbind(frequency, specenv, beta)

}
dev.off()



