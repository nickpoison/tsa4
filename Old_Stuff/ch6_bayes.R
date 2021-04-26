library('astsa')

####### this is the bayes part ######
#####################################

####################################
##### locallevel ###################
###################################
##-- Notation --##
#           y(t) = x(t) + v(t);    v(t) ~ iid N(0,V)                     
#           x(t) = x(t-1) + w(t);  w(t) ~ iid N(0,W)                        
#  priors:  x(0) ~ N(m0,C0); V ~ IG(a,b); W    ~ IG(c,d)
#    FFBS:  x(t|t) ~ N(m,C);  x(t|n) ~ N(mm,CC); x(t|t+1) ~ N(a,R)  
##-- 
ffbs = function(y,V,W,m0,C0){
  n  = length(y);  a  = rep(0,n);  R  = rep(0,n)
  m  = rep(0,n);   C  = rep(0,n);  B  = rep(0,n-1)     
  H  = rep(0,n-1); mm = rep(0,n);  CC = rep(0,n)
  x  = rep(0,n); llike = 0.0
  for (t in 1:n){
    if(t==1){
      a[1] = m0; R[1] = C0 + W
    }else{
      a[t] = m[t-1]; R[t] = C[t-1] + W    }
    f      = a[t]
    Q      = R[t] + V
    A      = R[t]/Q
    m[t]   = a[t]+A*(y[t]-f)
    C[t]   = R[t]-Q*A**2
    B[t-1] = C[t-1]/R[t]
    H[t-1] = C[t-1]-R[t]*B[t-1]**2
    llike  = llike + dnorm(y[t],f,sqrt(Q),log=TRUE)
  }
  mm[n] = m[n]; CC[n] = C[n]
  x[n]  = rnorm(1,m[n],sqrt(C[n]))
  for (t in (n-1):1){
    mm[t] = m[t] + C[t]/R[t+1]*(mm[t+1]-a[t+1])
    CC[t] = C[t] - (C[t]^2)/(R[t+1]^2)*(R[t+1]-CC[t+1])
    x[t]  = rnorm(1,m[t]+B[t]*(x[t+1]-a[t+1]),sqrt(H[t]))    }
return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))
}
# Simulate states and data
# ------------------------
set.seed(1); W  =  0.5; V  =  1.0
n  = 100; m0 = 0.0; C0 = 10.0; x0 = 0
w  = rnorm(n,0,sqrt(W))
v  = rnorm(n,0,sqrt(V))
x  = y = rep(0,n)
x[1] = x0   + w[1]
y[1] = x[1] + v[1]
for (t in 2:n){
  x[t] = x[t-1] + w[t]
  y[t] = x[t] + v[t]   }
### actual smoother (for plotting)
ks = Ksmooth0(num=n, y, A=1, m0, C0, Phi=1, cQ=sqrt(W), cR=sqrt(V))
xsmooth = as.vector(ks$xs)
###
run = ffbs(y,V,W,m0,C0)
m = run$m; C = run$C; mm = run$mm
CC = run$CC; L1 = m-2*C; U1  = m+2*C
L2 = mm-2*CC; U2 = mm+2*CC
N   = 50
Vs  = seq(0.1,2,length=N)
Ws  = seq(0.1,2,length=N)
likes = matrix(0,N,N)
for (i in 1:N){
  for (j in 1:N){
    V   = Vs[i]
    W   = Ws[j]
    run = ffbs(y,V,W,m0,C0)    
    likes[i,j] = run$llike  }
}  
# Hyperparameters
a = 0.01; b = 0.01; c = 0.01; d = 0.01
# MCMC step
set.seed(90210)
burn = 10; M = 1000
niter = burn + M
V1    = V; W1    = W
draws = NULL
all_draws = NULL
for (iter in 1:niter){
  run   = ffbs(y,V1,W1,m0,C0)
  x     = run$x
  V1    = 1/rgamma(1,a+n/2,b+sum((y-x)^2)/2)
  W1    = 1/rgamma(1,c+(n-1)/2,d+sum(diff(x)^2)/2)
  draws = rbind(draws,c(V1,W1,x))
}
all_draws = draws[,1:2]
q025 = function(x){quantile(x,0.025)}
q975 = function(x){quantile(x,0.975)}
draws = draws[(burn+1):(niter),]
xs    = draws[,3:(n+2)]
lx    = apply(xs,2,q025)
mx    = apply(xs,2,mean)
ux    = apply(xs,2,q975)
##  plot of the data

pdf(file="locallevel.pdf", width=7.6, height=5)
layout(matrix(1:4,2), widths=c(2,1))
par( mgp=c(1.6,.6,0), mar=c(2,2.5,.5,0)+.5)
ts.plot(ts(x), ts(y), ylab='', type='n')
grid(lty=1); lines(ts(y), lwd=2, col=8); lines(ts(x), lwd=2)
points(y)
legend(0,11,legend=c("x(t)","y(t)"),lty=1,col=c(1,8),lwd=2, bty='n',pch=c(-1,1) )
#
contour(Vs, Ws, exp(likes), xlab=expression(sigma[v]^2), ylab=expression(sigma[w]^2), drawlabels=FALSE, ylim=c(0,1.2), col=gray(.4) )
grid(lty=1)
points(draws[,1:2], pch=20, col=rgb(.9,0,0,alpha=.25),cex=1)
#
hist(draws[,2],main="", ylab="Density", xlab=expression(sigma[w]^2), col=rgb(.68,.85,.90, alpha=.4))
abline(v=mean(draws[,2]), col=4, lwd=3)
#
hist(draws[,1], ylab="Density",main="", xlab=expression(sigma[v]^2), col=rgb(.68,.85,.90, alpha=.4))
abline(v=mean(draws[,1]), col=4, lwd=3)
dev.off()

############
## plot states
pdf(file="locallevelstates.pdf", width=7.7, height=4)
par( mgp=c(1.6,.6,0), mar=c(2,1,.5,0)+.5)
plot(ts(mx), ylab='', type='n', ylim=c(min(y),max(y)))
grid()
points(y)
#lines(ts(lx), lty=2, col= rgb(0,.6,.3,alpha=.2)) 
#lines(ts(ux), lty=2, col= rgb(0,.6,.3,alpha=.2) )
lines(xsmooth, lwd=4, lty=1, col=rgb(1,0,1,alpha=.4))
lines(mx, lwd=1, lty=1 ,col= 4)# rgb(0,.3,0,alpha=.5))
 xx=c(1:100, 100:1)
 yy=c(lx, rev(ux))
polygon(xx, yy, border=NA, col= gray(.7,alpha=.2),lty=2)
lines(y,    col=gray(.4, alpha=1) )
legend('topleft', c('true smoother', 'data', 'posterior mean', '95% of draws'), 
                    lty=c(1,1,1,1), lwd=c(3,1,1,10), pch=c(-1,1,-1,-1),
col=c(6, gray(.4) ,4, gray(.6, alpha=.5)), bg='white' )  
dev.off()


###############################################################################
#########  structural model #######################################
###############################################################

require(plyr)
y = jj


# setup - model and initial parameters
set.seed(90210)
n = length(y)
F = c(1,1,0,0)      # this is A'
G = diag(0,4)       # G is Phi 
  G[1,1] = 1.03 
  G[2,]=c(0,-1,-1,-1); G[3,]=c(0,1,0,0); G[4,]=c(0,0,1,0)
a1 = rbind(.7,0,0,0) 
R1 = diag(.04,4)
V = .1
W11 = .1
W22 = .1

# FFBS scheme: sampling (x^n|phi,V,W11,W22) and (phi,V,W11,W22|x^n)
# -----------------------------------------------------
ffbs = function(y,F,G,V,W11,W22,a1,R1){
  n  = length(y)
  Ws = diag(c(W11,W22,1,1))
  iW = diag(1/diag(Ws),4)
  a  = matrix(0,n,4)
  R  = array(0,c(n,4,4))
  m  = matrix(0,n,4)
  C  = array(0,c(n,4,4))
  a[1,]  = a1[,1]
  R[1,,] = R1
  f      = t(F)%*%a[1,]
  Q      = t(F)%*%R[1,,]%*%F+V
  A      = R[1,,]%*%F/Q[1,1]
  m[1,]  = a[1,]+A%*%(y[1]-f)
  C[1,,] = R[1,,]-A%*%t(A)*Q[1,1]
  for (t in 2:n){
    a[t,]  = G%*%m[t-1,]
    R[t,,] = G%*%C[t-1,,]%*%t(G) + Ws
    f      = t(F)%*%a[t,]
    Q      = t(F)%*%R[t,,]%*%F+V
    A      = R[t,,]%*%F/Q[1,1]
    m[t,]  = a[t,]+A%*%(y[t]-f)
    C[t,,] = R[t,,]-A%*%t(A)*Q[1,1]
  }
  xb       = matrix(0,n,4)
  xb[n,]  = m[n,] + t(chol(C[n,,]))%*%rnorm(4)
  for (t in (n-1):1){
    iC  = solve(C[t,,])
    CCC = solve(t(G)%*%iW%*%G+iC)
    mmm = CCC%*%(t(G)%*%iW%*%xb[t+1,]+iC%*%m[t,])
    xb[t,] = mmm + t(chol(CCC))%*%rnorm(4)
  }
  return(xb)
}

##################################################################
#<i> Prior hyperparameters</i>

# b0 = 0     #  mean for beta = phi -1
# B0 = Inf   #  var for  beta  (non-informative prior => use OLS for sampling beta)
n0 = 10  # use same for all- the prior is 1/Gamma(n0/2, n0*s20_/2)
s20v = .001  # for V
s20w =.05    # for Ws

##########################################################


# MCMC scheme
set.seed(90210)
burnin  = 100 
step    = 10   
M       = 1000  
niter   = burnin+step*M
pars    = matrix(0,niter,4)
xbs     = array(0,c(niter,n,4))

pr <- progress_text()           # displays progress
pr$init(niter)
for (iter in 1:niter){
  xb = ffbs(y,F,G,V,W11,W22,a1,R1)
		u = xb[,1] 
		yu = diff(u); xu = u[-n]    # for phihat and se(phihat)
		regu = lm(yu~0+xu)                # est of beta = phi-1
		phies = as.vector(coef(summary(regu)))[1:2] + c(1,0)  # phi [e]stimate and [s]e
		dft = df.residual(regu)   
    G[1,1] = phies[1] + rt(1,dft)*phies[2]  # use a t
	V  = 1/rgamma(1, (n0+n)/2, (n0*s20v/2) + sum((y-xb[,1]-xb[,2])^2)/2)
    W11 = 1/rgamma(1, (n0+n-1)/2, (n0*s20w/2) + sum((xb[-1,1]-phies[1]*xb[-n,1])^2)/2)
	W22 = 1/rgamma(1, (n0+ n-3)/2, (n0*s20w/2) + sum((xb[4:n,2] + xb[3:(n-1),2]+ xb[2:(n-2),2] +xb[1:(n-3),2])^2)/2)
  xbs[iter,,] = xb
  pars[iter,] = c(G[1,1], sqrt(V), sqrt(W11), sqrt(W22))
  pr$step()
}
# graohs ###############
pdf(file="dlm-ch12-params.pdf", width=7.7, height=5.5)
ind = seq(burnin+1,niter,by=step)
names= c(expression(phi), expression(sigma[v]), expression(sigma[w~11]), expression(sigma[w~22]))
par(mfcol=c(3,4), mar=c(2,2,.25,0)+.75, mgp=c(1.6,.6,0), oma=c(0,0,1,0))
for (i in 1:4){
  plot.ts(pars[ind,i],xlab="iterations", ylab="trace", main="", col=gray(.1, alpha=.6))
  mtext(names[i], side=3, line=.5, cex=1) 
  acf(pars[ind,i],main="", lag.max=25, xlim=c(1,25), ylim=c(-.4,.4))
  hist(pars[ind,i],main="",xlab="")
  abline(v=mean(pars[ind,i]), lwd=2, col=4)
}
dev.off()


 ############## p(x^n | y^n) ###########
 #{dlm-ch12-smoothers}
 pdf(file="dlm-ch12-smoothers.pdf", width=7.7, height=5)
  par(mfrow=c(2,1), mar=c(2,2,0,0)+.7, mgp=c(1.6,.6,0))
  mxb = cbind(apply(xbs[ind,,1],2,mean), apply(xbs[,,2],2,mean))
  lxb = cbind(apply(xbs[ind,,1],2,quantile,0.005),apply(xbs[ind,,2],2,quantile,0.005))
  uxb = cbind(apply(xbs[ind,,1],2,quantile,0.995),apply(xbs[ind,,2],2,quantile,0.995))
  
  mxb=ts(cbind(mxb,rowSums(mxb)), start = tsp(jj)[1], freq=4) 
  lxb=ts(cbind(lxb,rowSums(lxb)), start = tsp(jj)[1], freq=4)
  uxb=ts(cbind(uxb,rowSums(uxb)), start = tsp(jj)[1], freq=4)
  names=c('Trend', 'Season', 'Trend + Season')
  
  #for (i in 1:2){
  #i=1
      L = min(lxb[,1])-.01; U = max(uxb[,1]) +.01
     plot(mxb[,1],  ylab=names[1], ylim=c(L,U), type='n')
       grid(lty=1); lines(mxb[,1])	 
      xx=c(time(jj), rev(time(jj)))
      yy=c(lxb[,1], rev(uxb[,1]))
     polygon(xx, yy, border=NA, col=gray(.4, alpha = .2)) 
  #i=3  
    L = min(lxb[,3])-.01; U = max(uxb[,3]) +.01
     plot(mxb[,3],  ylab=names[3], ylim=c(L,U), type='n')
       grid(lty=1); lines(mxb[,3]) 
      xx=c(time(jj), rev(time(jj)))
      yy=c(lxb[,3], rev(uxb[,3]))
     polygon(xx, yy, border=NA, col=gray(.4, alpha = .2)) 
      
 # }
  dev.off()
######################


