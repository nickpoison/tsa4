## Edition 4 Errata  

_For some reason, equations are not being displayed properly on this page ... we don't know why, but it's  not that big of a problem ... so here we go:_

<br/>

### Chapter 1 

- it's perfect

### Chapter 2 

- __Eq (2.38)-(2.39):__ To be more general, the time subscript should have been $t_i$ so the equations would read $$m_t = \sum_{i=1}^n  w_i(t) x_{t_i}  \tag{2.38}$$ where 
$$w_i(t)= K\bigl(\tfrac{t-{t_i}}{b}\bigr) \Bigm/   \sum_{j=1}^n K\bigl(\tfrac{t-{t_j}}{b}\bigr) \tag{2.39}$$
... and typically $K(z)=\exp(-z^2/2)$ is used (no need for constants because the weights are normalized).


So these two are the same:
```r
par(mfrow=c(2,1))
tsplot(soi)      # monthly data; frequency=12 and t = 1/12, 2/12, ...
lines(ksmooth(time(soi), soi, kernel="normal", bandwidth=1), lwd=2, col=4)
# and
SOI = ts(soi, frequency=1)   # change to t = 1,2,... 
tsplot(SOI)   
lines(ksmooth(time(SOI), SOI, kernel="normal", bandwidth=12), lwd=2, col=4)
```

 

###  Chapter 3 

- __Eq (3.10):__ The sum should be to $k$ (and not $k-1$):
$$x_t = \phi^{-k} x_{t+k} - \sum_{j=1}^{k} \phi^{-j} w_{t+j}\,. $$


### Chapter 4 

- as if


### Chapter 5 

 - __Example 5.1:__ I put this note on the R code page, but I thought I'd repeat it here. 
In Example 5.1, we used <kbd>fracdiff</kbd>, but it's not a very good package. 
We should have used another package such as <kbd>arfima</kbd>, but unfortunately it didn't make it into the revision. This is changed in Edition 5: 

```r
library(arfima)
summary(varve.fd <- arfima(log(varve)))  # d.hat = 0.3728, se(d,hat) = 0.0273
# residual stuff
innov = resid(varve.fd)  
plot.ts(innov[[1]])  
acf(innov[[1]])
```

### Chapter 6 

 - __Property 6.7, equation (6.137):__ Left off the conditioning arguments ... the $\pi_j(t)$  in the numerator and in the denominator  should be $\pi_j(t \mid t-1)$ . 


### Chapter 7 

- not enough people read this chapter to find the bloopers ... but we're fairly certain there are a few


###  Elsewhere

 - FYI: In Edition 5, Appendix R has been removed and put online  here: [dsstoffer.github.io/Rtoot](https://dsstoffer.github.io/Rtoot)


<br/><br/><br/>
