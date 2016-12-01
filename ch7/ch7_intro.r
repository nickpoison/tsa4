require('astsa')

################fmri#######################
my.panel <- function(x, y=NULL, ..., pf = parent.frame()) {
  abline(v=seq(20,120,20), h=seq(-.4,.4,.2),  col = gray(.9), lty = 1, lwd = par("lwd"))
 lines(x, y, ...) }

pdf(file="mbrain1.pdf", width=7.7, height=5)
#dev.new(height=6,width=7.7)
x = matrix(0, 128, 6)
for (i in 1:6) x[,i] = rowMeans(fmri[[i]])
colnames(x)=c("Brush", "Heat", "Shock", "Brush", "Heat", "Shock")
plot.ts(x,    main="", panel=my.panel,  oma.multi = c(4, 0, 2, 0), mar.multi = c(0, 4.1, 0, 2.1))
mtext("Awake", side=3, line=3, adj=.1, cex=1.1 )
mtext("Sedated", side=3, line=3, adj=.85, cex=1.1 )
dev.off()

########### Eq exp ############################
my.panel <- function(x, y=NULL, ..., pf = parent.frame()) {
  abline(v=seq(200,1000,200), h=seq(-7,7,2),  col = gray(.9), lty = 1, lwd = par("lwd"))
 lines(x, y=NULL, ...) }

# eqexplot
pdf(file="eqexplot.pdf", width=7.7, height=5)
#dev.new(height=6,width=7.7)
attach(eqexp)
P = 1:1024; S = P+1024
x = cbind(EQ5[P], EQ6[P], EX5[P], EX6[P], NZ[P], EQ5[S], EQ6[S], EX5[S],EX6[S], NZ[S])
x.name = c("EQ5","EQ6","EX5","EX6","NZ")
colnames(x) = c(x.name, x.name)
plot.ts(x, main="",  panel=my.panel,  oma.multi = c(4, 0, 2, 0), mar.multi = c(0, 4.2, 0, 2.1), cex.lab=.9)
mtext("P waves", side=3, line=3, adj=.1, cex=1.1)
mtext("S waves", side=3, line=3, adj=.84, cex=1.1)
dev.off() 

#######  alternate plot without dividing P and S wave  - forget it

my.panel <- function(x, y=NULL, ..., pf = parent.frame()) {
  abline(v=seq(0,2048,2^8), h=seq(-6,6,2),  col = gray(.9), lty = 1, lwd = par("lwd"))
  abline(v=1024, col=4,lty=5)
 lines(x, y=NULL, ...) }

# eqexplot
###pdf(file="eqexplot2.pdf", width=7.7, height=5)
attach(eqexp)
x = cbind(EQ5, EQ6, EX5, EX6, NZ)
colnames(x) =  c("EQ5","EQ6","EX5","EX6","NZ")

dev.new(height=5,width=7.7)
plot.ts(x, main="",  panel=my.panel, nc=1,  oma.multi = c(4, 0, 2, 0), mar.multi = c(0, 4.2, 0, 2.1), cex.lab=1, xaxt = "n")
axis(1, at=c(0,1024,2048), line=1.5, cex.axis=.8, lty=0)


#mtext("P waves", side=3, line=3, adj=.1, cex=1.1)
#mtext("S waves", side=3, line=3, adj=.84, cex=1.1)
### dev.off() 

################

 
 
##############  climhyd ##################### 
my.panel <- function(x, y=NULL, ..., pf = parent.frame()) {
  abline(v=seq(100,400,100,), col = gray(.9), lty = 1, lwd = par("lwd"))
#  abline( h=seq(0,15,5), col = "lightgray", lty = 2, lwd = par("lwd"))
#  abline( h=seq(0,1.4,.2), col = "lightgray", lty = 2, lwd = par("lwd"))
#  abline( h=seq(0,1200,200), col = "lightgray", lty = 2, lwd = par("lwd"))
 lines(x, y=NULL, ...) }

 
  pdf(file="climhyd.pdf", width=7.7, height=5)
plot.ts(climhyd, main='', panel=my.panel, oma.multi = c(4, 0, 2, 0),   mar.multi = c(0, 4.2, 0, 2.1), cex.lab=.9)
 dev.off()

