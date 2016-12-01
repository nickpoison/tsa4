
set.seed(90210)
u = exp(3*sin(2*pi*1:500*.1) + rnorm(500,0,4))
spec.pgram(u, spans=c(5,3), taper=.5, log="dB")
dev.new()
x = cbind(u, sqrt(u), u^(1/3))           # transformation set 
Var = var(x)                     
xspec = mvspec(x, spans=c(5,3), taper=.5) 
fxxr = Re(xspec$fxx)            
ev = eigen(Var)  
Q = ev$vectors%*%diag(1/sqrt(ev$values))%*%t(ev$vectors)               
num = xspec$n.used               
nfreq = length(xspec$freq)         
specenv = matrix(0,nfreq,1)      
beta=matrix(0,nfreq,3)          
 for (k in 1:nfreq){ 
  ev = eigen(2*Q%*%fxxr[,,k]%*%Q/num) 
  specenv[k] = ev$values[1]            
  b = Q%*%ev$vectors[,1]               
  beta[k,] = b/sign(b[1]) }
#
#--- output and graphics ---# 
frequency = xspec$freq
plot(frequency, 100*specenv, type="l", ylab="Spectral Envelope (%)")
 ## add significance threshold to plot ##
 m=xspec$kernel$m
 etainv=sqrt(sum(xspec$kernel[-m:m]^2))
thresh=100*(2/num)*exp(qnorm(.999)*etainv)*rep(1,nfreq)
lines(frequency,thresh, lty="dashed", col="blue")

 
#--  details  --# 
output = cbind(frequency, specenv, beta)
colnames(output)=c("freq","specenv","x", "sqrt(x)", "x^(1/3)")
round(output,4) 

# plot transform
dev.new()
b = output[50,3:5]
g = function(x) 4.5 + b[1]*x+b[2]*sqrt(x)+b[3]*x^(1/3)
curve(g, 1, 4000)
lines(log(1:4000), lty=2)


# 
# > b
#             x       sqrt(x)       x^(1/3) 
#  6.935931e-05 -1.186615e-01  6.887149e-01 
# 








