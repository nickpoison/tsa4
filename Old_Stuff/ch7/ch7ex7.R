# chapter 7 example 7
require(astsa)

n=128   # length of series
n.freq = 1 + n/2  # number of frequencies
Fr = (0:(n.freq-1))/n  # the frequencies 
#number of vector series for each cell
N = c(5,4,5,3,5,4)
n.subject = sum(N)  # 26 
n.trt = 6    #number treatments
L = 3   # amt of smoothing
num.df = 2*L*(n.trt-1)   # df for F test
den.df = 2*L*(n.subject-n.trt)


# Design Matrix for the 1-way ANOVA - 6 treatment combinations T1,...T6
Z1 = outer(rep(1,N[1]), c(1,1,0,0,0,0)) # Mean = mu + alpha1
Z2 = outer(rep(1,N[2]), c(1,0,1,0,0,0)) # Mean = mu + alpha2
Z3 = outer(rep(1,N[3]), c(1,0,0,1,0,0)) # Mean = mu + alpha3
Z4 = outer(rep(1,N[4]), c(1,0,0,0,1,0)) # Mean = mu + alpha4
Z5 = outer(rep(1,N[5]), c(1,0,0,0,0,1)) # Mean = mu + alpha5
Z6 = outer(rep(1,N[6]), c(1,-1,-1,-1,-1,-1)) # Mean = mu + alpha6
Z = rbind(Z1, Z2, Z3, Z4, Z5, Z6) #Design matrix (n.subject x n.trt = 26 x 6)
ZZ = t(Z)%*%Z 

   
    SSEF = rep(NA, n)   # F = full
    SSER = rep(NA, n)   # R = reduced
    # HatF = Z%*%solve(ZZ)%*%t(Z)
     HatF = Z%*%solve(ZZ,t(Z))
     HatR = Z[,1]%*%t(Z[,1])/ZZ[1,1]

     
pdf(file="fig_ex7_7.pdf", width=7.7, height=5)
par(mfrow=c(3,3), mar=c(3.5,4,0,0), oma=c(0,0,2,2),  mgp = c(1.6,.6,0))
loc.name =c("Cortex 1", "Cortex 2", "Cortex 3", "Cortex 4", "Caudate",
              "Thalamus 1", "Thalamus 2", "Cerebellum 1", "Cerebellum 2")

for(Loc in 1:9) {   # Loc is location
      i = 6*(Loc-1) 
      Y=cbind(fmri[[i+1]],fmri[[i+2]],fmri[[i+3]],fmri[[i+4]],fmri[[i+5]],fmri[[i+6]])
      # Y is 128 x 26 matrix of observations for each Locations
      Y = spec.taper(Y, p=.5)	# taper the data
      Y = mvfft(Y)/sqrt(n)  # now Y is the FFT of tapered centered data
      Y = t(Y)        # Y is now 26 x 128 FFTs


    ###Calculation of Error Spectra  ###

   

    for (k in 1:n) {
	# Full Model Error Power    
          SSY=Re(Conj(t(Y[,k]))%*%Y[,k])
          SSReg= Re(Conj(t(Y[,k]))%*%HatF%*%Y[,k])
        SSEF[k]=SSY-SSReg
    # Reduced Model Error Power
          SSReg= Re(Conj(t(Y[,k]))%*%HatR%*%Y[,k])
        SSER[k]=SSY-SSReg
       }

    # Smooth 
    sSSEF = filter(SSEF, rep(1/L, L), circular = TRUE)
	sSSER = filter(SSER, rep(1/L, L), circular = TRUE)

    eF =(den.df/num.df)*(sSSER-sSSEF)/sSSEF

#
 
    plot(Fr,eF[1:n.freq], type="l", xlab="Frequency", ylab="F Statistic", ylim=c(0,7), panel.first=grid(lty=1,ny=NA))
    abline(h=qf(.999, num.df, den.df),lty=2)
    text(.25, 6.5, loc.name[Loc], cex=1.2)
    }
dev.off()