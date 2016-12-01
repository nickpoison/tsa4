require(astsa)

u=read.table('465CAT.DAT')
z=1:6
uu =  as.matrix(u)%*%as.matrix(z)

pdf(file="ch7_sleep.pdf", width=7.25, height=4)   
par(mfrow=c(2,1), mar=c(2.2,3,.25,.25)+.5, mgp=c(1.6,.6,0), cex=.9, cex.main=1)  

# plot data
#dev.new(height=2.5)
#par(mar=c(3,3,1,1), mgp=c(1.6,.6,0), oma=c(0,1,0,0))
plot.ts(uu[1:115], type='n',   axes=FALSE, ylab="", xlab='Minute' )
box()
axis(side = 1)
grid(lty=1)
lines(uu[1:115], type='s')
states = c('NR4', 'NR3', 'NR2', 'NR1', 'REM', 'AWAKE')
axis(side = 2, 1:6, labels = states, las=1)
mtext('Sleep State', side=2, line=2.5, cex=1)


# plot periodogram 1-6
#dev.new(height=2.5)
#par(mar=c(3,3,1,1), mgp=c(1.6,.6,0))
slper=mvspec(uu[1:115], log='n', plot=FALSE)
plot(slper$freq, slper$spec, type='l', ylab="Periodogram", xlab="Frequency", panel.first=grid(lty=1))
abline(v=slper$freq[2], lty=2)
mtext('1/60', side=1, adj=.04, cex=.75)
dev.off()




####
####sleep = read.table("/mydata/ch7sleep.dat")
####x = as.ts(sleep[,2])
####n = nextn(length(x)); n2 = n/2
####Fr = 0:n2/n
####x.per = abs(fft(x-mean(x)))^2/n
####
####par(mfrow=c(2,1),  mar=c(3.5,3,1,1), mgp=c(1.6,.6,0))
####plot(x, ylab="Sleep-State", xlab="Minute of Sleep")
####plot(Fr, x.per[1:(n2+1)], type="l", ylab="Periodogram", xlab="Frequency")
