require(astsa)


# set up data 
P=1:1024
S=P+1024
mag.P = log10(apply(eqexp[P,],2,max) - apply(eqexp[P,],2,min))
mag.S = log10(apply(eqexp[S,],2,max) - apply(eqexp[S,],2,min))
eq.P = mag.P[1:8];  eq.S = mag.S[1:8]
ex.P = mag.P[9:16]; ex.S = mag.S[9:16]
NZ.P = mag.P[17];   NZ.S = mag.S[17]


## Compute the linear discriminant function:
cov.eq = var(cbind(eq.P, eq.S))
cov.ex = var(cbind(ex.P, ex.S))
cov.pooled = (cov.ex + cov.eq)/2
means.eq =  colMeans(cbind(eq.P,eq.S))
means.ex =  colMeans(cbind(ex.P,ex.S))
slopes.eq = solve(cov.pooled, means.eq)
inter.eq = -sum(slopes.eq*means.eq)/2
slopes.ex = solve(cov.pooled, means.ex)
inter.ex = -sum(slopes.ex*means.ex)/2
d.slopes = slopes.eq - slopes.ex
d.inter = inter.eq - inter.ex
## Classify new observation:
new.data = cbind(NZ.P, NZ.S)
d = sum(d.slopes*new.data) + d.inter
post.eq = exp(d)/(1+exp(d))

## print results
cat(d.slopes[1], "mag.P +" , d.slopes[2], "mag.S +" , d.inter,"\n") # Discriminant function
cat("P(EQ|data) =", post.eq,  "  P(EX|data) =", 1-post.eq, "\n" ) # Classification of NZ event

pdf(file="fig_ex7_11.pdf", width=7.7, height=4) 
par(mar=c(3,3,2,1),  mgp = c(1.6,.6,0), cex.main=1.1)
plot(eq.P, eq.S, xlim = c(0,1.5), ylim = c(.75,1.25), panel.first=grid(lty=1),
	xlab = "log mag(P)", ylab = "log mag(S)",  pch = 8, cex=1.1, lwd=2, col=4,
	main="Classification Based on Magnitude Features")
points(ex.P, ex.S, pch = 6, cex=1.1, lwd=2, col=6)
points(new.data, pch = 3, cex=1.1, lwd=2, col=rgb(0,.6,.2))
abline(a = -d.inter/d.slopes[2], b = -d.slopes[1]/d.slopes[2])
text(eq.P-.07,eq.S+.005, label=names(eqexp[1:8]), cex=.8)
text(ex.P+.07,ex.S+.003, label=names(eqexp[9:16]), cex=.8)
text(NZ.P+.05,NZ.S+.003, label=names(eqexp[17]), cex=.8)
legend("topright", legend=c("EQ", "EX", "NZ"), pch=c(8,6,3), pt.lwd=2, cex=1.1, bg='white', col=c(4,6,rgb(0,.6,.2)))
dev.off()



################# Cross-validation (loo)  ############ (leave this out of example ... maybe mention it)#####
 
 all.data = rbind(cbind(eq.P,eq.S), cbind(ex.P,ex.S))
 colnames(all.data) = c("P", "S")
 
 post.eq <- rep(NA,8) -> post.ex 
 
 for(j in 1:16) {
 if(j<=8) {
 	sample.eq = all.data[-c(j, 9:16),]
 	sample.ex = all.data[9:16,]
 	}
 if(j>8){
 	sample.eq = all.data[1:8,]
 	sample.ex = all.data[-c(j, 1:8),]
 	}
 # Compute the discriminant function; evaluate at the hold out:
 df.eq = nrow(sample.eq)-1
 df.ex = nrow(sample.ex)-1
 mean.eq = colMeans(sample.eq)
 mean.ex = colMeans(sample.ex)
 cov.eq = var(sample.eq)
 cov.ex = var(sample.ex)
 cov.pooled = (df.eq*cov.eq + df.ex*cov.ex)/(df.eq + df.ex)
 slopes.eq = solve(cov.pooled, mean.eq)
 inter.eq = -sum(slopes.eq*mean.eq)/2
 slopes.ex = solve(cov.pooled, mean.ex)
 inter.ex = -sum(slopes.ex*mean.ex)/2
 d.slopes = slopes.eq - slopes.ex
 d.inter = inter.eq - inter.ex
 d = sum(d.slopes*all.data[j,]) + d.inter
 
 if (j<= 8) post.eq[j] = exp(d)/(1+exp(d))
 if (j > 8) post.ex[j-8] = 1/(1+exp(d))
 }
 
 #
 Posterior = cbind(1:8, post.eq, 1:8, post.ex)
 colnames(Posterior) = c("EQ", "P(EQ|data)", "EX", "P(EX|data)")
 round(Posterior,3)  # Results from Cross-Validation
  


