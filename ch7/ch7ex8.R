# chapter 7 example 8
require(astsa)

n=128  # length of series
n.freq = 1 + n/2  # number of frequencies
Fr = (0:(n.freq-1))/n  # the frequencies 
nFr = 1:(n.freq/2)  # number of freqs plotted

#number of vector series for each cell
N = c(5,4,5,3,5,4)
n.subject = sum(N)  # 26 
n.para = 6    # number of parameters
L = 3   # amt of smoothing
df.stm=2*L*(3-1)   # stimulus (3 levels: Brush,Heat,Shock)
df.con=2*L*(2-1)   # conscious (2 levels: Awake,Sedated)  
df.int=2*L*(3-1)*(2-1)  # interaction 
den.df= 2*L*(n.subject-n.para) # df for full model 
                        




# Design Matrix:          mu  a1  a2   b  g1  g2   
 Z1 = outer(rep(1,N[1]), c(1,  1,  0,  1,  1,  0)) 
 Z2 = outer(rep(1,N[2]), c(1,  0,  1,  1,  0,  1)) 
 Z3 = outer(rep(1,N[3]), c(1, -1, -1,  1, -1, -1)) 
 Z4 = outer(rep(1,N[4]), c(1,  1,  0, -1, -1,  0)) 
 Z5 = outer(rep(1,N[5]), c(1,  0,  1, -1,  0, -1)) 
 Z6 = outer(rep(1,N[6]), c(1, -1, -1, -1,  1,  1))
 Z = rbind(Z1, Z2, Z3, Z4, Z5, Z6) #Design matrix (n.subject x n.trt = 26 x 6)
 ZZ = t(Z)%*%Z

rep(NA, n)-> SSEF->SSE.stm->SSE.con-> SSE.int            
# HatF = Z%*%solve(ZZ)%*%t(Z)     
# Hat.stm = Z[,-(2:3)]%*%solve(ZZ[-(2:3),-(2:3)])%*%t(Z[,-(2:3)])  
# Hat.con = Z[,-4]%*%solve(ZZ[-4,-4])%*%t(Z[,-4])   
# Hat.int = Z[,-(5:6)]%*%solve(ZZ[-(5:6),-(5:6)])%*%t(Z[,-(5:6)])   
 HatF = Z%*%solve(ZZ,t(Z))     
 Hat.stm = Z[,-(2:3)]%*%solve(ZZ[-(2:3),-(2:3)], t(Z[,-(2:3)]))  
 Hat.con = Z[,-4]%*%solve(ZZ[-4,-4], t(Z[,-4]))   
 Hat.int = Z[,-(5:6)]%*%solve(ZZ[-(5:6),-(5:6)], t(Z[,-(5:6)]))                                                            
                     
pdf(file="fig_ex7_8.pdf", width=7.7, height=6)
par(mfrow=c(5,3), mar=c(3.5,4,0,0), oma=c(0,0,2,2),  mgp = c(1.6,.6,0)) 
loc.name =c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate",
             "Thalamus 1", "Thalamus 2", "Cerebellum 1", "Cerebellum 2")
for(Loc in c(1:4,9)) {   # Loc is location - only 1 to 4 and 9 used                                                                                                                             
      i = 6*(Loc-1)                                                                 
      Y=cbind(fmri[[i+1]],fmri[[i+2]],fmri[[i+3]],fmri[[i+4]],fmri[[i+5]],fmri[[i+6]])
      # Y is 128 x 26 matrix of observations for each Locations                       
      Y = spec.taper(Y, p=.5)	# taper the data                                        
      Y = mvfft(Y)/sqrt(n)  # now Y is the FFT of tapered centered data   
      Y = t(Y)  
      


      for (k in 1:n) {
	   # Full Model    
          SSY=Re(Conj(t(Y[,k]))%*%Y[,k])
          SSReg= Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
        SSEF[k]=SSY-SSReg
       # Reduced Model 
          SSReg=Re(Conj(t(Y[,k]))%*%Hat.stm%*%Y[,k])
        SSE.stm[k] = SSY-SSReg
          SSReg=Re(Conj(t(Y[,k]))%*%Hat.con%*%Y[,k])
        SSE.con[k]=SSY-SSReg  
          SSReg=Re(Conj(t(Y[,k]))%*%Hat.int%*%Y[,k])
        SSE.int[k]=SSY-SSReg   
                    }
   
   # Smooth    
   sSSEF = filter(SSEF, rep(1/L, L), circular = TRUE)
   sSSE.stm = filter(SSE.stm, rep(1/L, L), circular = TRUE)
   sSSE.con = filter(SSE.con, rep(1/L, L), circular = TRUE)
   sSSE.int = filter(SSE.int, rep(1/L, L), circular = TRUE)
   
 
   eF.stm = (den.df/df.stm)*(sSSE.stm-sSSEF)/sSSEF
   eF.con = (den.df/df.con)*(sSSE.con-sSSEF)/sSSEF
   eF.int = (den.df/df.int)*(sSSE.int-sSSEF)/sSSEF
   plot(Fr[nFr],eF.stm[nFr], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,12), panel.first=grid(lty=1,ny=NA) )
       abline(h=qf(.999, df.stm, den.df),lty=2)       
       if(Loc==1) mtext("Stimulus", side=3, line=.3, cex=1)
       mtext(loc.name[Loc], side=2, line=3, cex=.9)
  plot(Fr[nFr],eF.con[nFr], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,12), panel.first=grid(lty=1,ny=NA))
       abline(h=qf(.999, df.con, den.df),lty=2)
       if(Loc==1)  mtext("Consciousness", side=3, line=.3, cex=1)   
   plot(Fr[nFr],eF.int[nFr], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,12), panel.first=grid(lty=1,ny=NA))
       abline(h=qf(.999, df.int, den.df),lty=2)
       if(Loc==1) mtext("Interaction", side=3, line= .3, cex=1)                           
   }  
dev.off()
   
   