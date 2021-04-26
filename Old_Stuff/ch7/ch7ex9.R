# chapter 7 example 9
require(astsa)
n.obs=128;
n = nextn(n.obs)  # which is 128
n.freq = 1 + n/2  # number of frequencies
Fr = (0:(n.freq-1))/n  # the frequencies 
nFr = 1:(n.freq/2)  # number of freqs plotted

#number of vector series for each cell
N = c(5,4,5,3,5,4)
n.subject = sum(N)  # 26 
L = 3   # amt of smoothing
        

# Design Matrix:          
 Z1 = outer(rep(1,N[1]), c(1,0,0,0,0,0)) 
 Z2 = outer(rep(1,N[2]), c(0,1,0,0,0,0)) 
 Z3 = outer(rep(1,N[3]), c(0,0,1,0,0,0)) 
 Z4 = outer(rep(1,N[4]), c(0,0,0,1,0,0)) 
 Z5 = outer(rep(1,N[5]), c(0,0,0,0,1,0)) 
 Z6 = outer(rep(1,N[6]), c(0,0,0,0,0,1))
 Z = rbind(Z1, Z2, Z3, Z4, Z5, Z6) #Design matrix (n.subject x n.trt = 26 x 6)
 ZZ = t(Z)%*%Z

# Contrasts:  6 x 3
A=rbind(diag(1,3),diag(1,3)) 
nq =  nrow(A) 
num.df=2*L*nq   
den.df= 2*L*(n.subject-nq) # df for full model 
                
HatF = Z%*%solve(ZZ)%*%t(Z) # full model 
rep(NA, n)-> SSEF->SSER
eF=matrix(0,n,3)

pdf(file="fig_ex7_9.pdf", width=7.7, height=6) 
par(mfrow=c(5,3), mar=c(3.5,4,0,0), oma=c(0,0,2,2),  mgp = c(1.6,.6,0)) 
loc.name =c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate",
             "Thalamus 1", "Thalamus 2", "Cerebellum 1", "Cerebellum 2")
cond.name=c("Brush", "Heat", "Shock")             


for(Loc in c(1:4,9)) {   # Loc is location - only 1 to 4 and 9 used                                                                                                                             
      i = 6*(Loc-1)                                                                 
      Y=cbind(fmri[[i+1]],fmri[[i+2]],fmri[[i+3]],fmri[[i+4]],fmri[[i+5]],fmri[[i+6]])
      # Y is 128 x 26 matrix of observations for each Locations                       
      Y = spec.taper(Y, p=.5)	# taper the data                                        
      Y = mvfft(Y)/sqrt(n)  # now Y is the FFT of tapered centered data   
      Y = t(Y) 
  
 for (cond in 1:3){
      Q = t(A[,cond])%*%solve(ZZ, A[,cond])
     HR = A[,cond]%*%solve(ZZ, t(Z))   
     
      for (k in 1:n){
	   # Full Model    
          SSY=Re(Conj(t(Y[,k]))%*%Y[,k])
          SSReg= Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
        SSEF[k]= (SSY-SSReg)*Q
       # Reduced Models 
          SSReg= HR%*%Y[,k]
         SSER[k] = Re(SSReg*Conj(SSReg)) 
                     }  
        
        # Smooth    
       sSSEF = filter(SSEF, rep(1/L, L), circular = TRUE)
       sSSER = filter(SSER, rep(1/L, L), circular = TRUE) 
       
       eF[,cond]= (den.df/num.df)*(sSSER/sSSEF)
        }                      
             
  plot(Fr[nFr],eF[nFr,1], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,5), panel.first=grid(lty=1,ny=NA) )
       abline(h=qf(.999, num.df, den.df),lty=2)       
       if(Loc==1) mtext("Brush", side=3, line=.3, cex=1)
       mtext(loc.name[Loc], side=2, line=3, cex=.9)
   plot(Fr[nFr],eF[nFr,2], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,5), panel.first=grid(lty=1,ny=NA) )
       abline(h=qf(.999, num.df, den.df),lty=2)
       if(Loc==1)  mtext("Heat", side=3, line=.3, cex=1)   
    plot(Fr[nFr],eF[nFr,3], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,5), panel.first=grid(lty=1,ny=NA) )
       abline(h=qf(.999, num.df, den.df),lty=2)
       if(Loc==1) mtext("Shock", side=3, line= .3, cex=1)                           
   }  
 dev.off()
  
   
              
                     
  
        
        
       