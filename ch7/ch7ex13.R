library(cluster)
#  figure is clusterPAM
P=1:1024
S=P+1024
p.dim = 2
n =1024
eq = as.ts(eqexp[,1:8])
ex = as.ts(eqexp[,9:16])
nz = as.ts(eqexp[,17])
f <- array(dim=c(17,2,2,512)) 



L = c(15,15)  # for smoothing 
for (i in 1:8){     # compute spectral matrices
	f[i,,,] = mvspec(cbind(eq[P,i],eq[S,i]), spans=L, taper=.5,plot=FALSE)$fxx
	f[i+8,,,] = mvspec(cbind(ex[P,i],ex[S,i]), spans=L, taper=.5,plot=FALSE)$fxx}
    f[17,,,] = mvspec(cbind(nz[P],nz[S]), spans=L, taper=.5,plot=FALSE)$fxx	


JD <- matrix(0,17,17) 
# calculate sym information criteria 
for (i in 1:16){
 for (j in (i+1):17){	
  for (k in 1:256) {    # only use freqs out to .25
    tr1 = Re(sum(diag(solve(f[i,,,k],f[j,,,k]))))
    tr2 = Re(sum(diag(solve(f[j,,,k], f[i,,,k]))))
    JD[i,j] = JD[i,j] + (tr1 + tr2 - 2*p.dim)}
}}
 JD = (JD + t(JD))/n

colnames(JD) = c(colnames(eq), colnames(ex), "NZ") 
rownames(JD) = colnames(JD)

cluster.2 = pam(JD, k = 2, diss = TRUE)  
summary(cluster.2)  # not shown

pdf(file="clusterPAM.pdf", width=7.5, height=5) 
par(mar=c(2,2,1,.5)+.5,  mgp = c(1.4,.6,0), cex=3/4, cex.lab=4/3, cex.main=4/3)
clusplot(JD, cluster.2$cluster, col.clus=gray(.5), labels=3, lines=0, col.p=c(rep(4,8),rep(rgb(.9,0,.9),8), rgb(0,.6,0)),   
main="Clustering Results for Explosions and Earthquakes")
text(-7,-.4, "Group I", cex=1.1, font=2); text(1,5,"Group II", cex=1.1, font=2)
dev.off()





# dev.new()
# cluster.3 = pam(JD, k = 3, diss = TRUE)  
# summary(cluster.3)  # not shown
# par(mgp = c(1.6,.6,0), cex=3/4, cex.lab=4/3, cex.main=4/3)
# clusplot(JD, cluster.3$cluster, col.clus=1, labels=3, lines=0, col.p=1,
# main="Clustering Results for Explosions and Earthquakes - 3 Groups")

#pch=c(rep(8,8),rep(6,8),3),

   