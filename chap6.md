## Chapter 6 - Kalman Filter and Smoother Details

	

+ There are  three levels of code depending on the complexity of the model.
They  are  `Kfilter0/Ksmooth0`,  `Kfilter1/Ksmooth1`, and `Kfilter2/Ksmooth2`.
For various models, each script provides the Kalman filter/smoother, the innovations
 and the corresponding variance-covariance matrices, and the value of the innovations likelihood at the location of the parameter values passed to the script. MLE is then accomplished by calling  the script that runs the filter. _The model is specified by passing the model parameters._
	  
+ Level 0 is for the case of a fixed measurement matrix and no inputs; i.e., if  A<sub>t</sub> = A for all t, and there are no inputs,   then use the code at level 0.  
    
+  If the measurement matrices are time varying or there are inputs, use the  code at a higher level (1 or 2).   Many of the examples in the text can be done at level 0.
    
    
+ Level 1 allows for time varying measurement matrices and inputs, and level 2 adds
    the possibility of correlated noise processes. 
    
The models for each case are (x is state, y is observation, and t = 1, &hellip;, n):
   
&diams; **Level 0:** &nbsp; &nbsp; x<sub>t</sub> = &Phi; x<sub>t-1</sub> + w<sub>t</sub>, &nbsp; &nbsp; y<sub>t</sub> = A x<sub>t</sub> + v<sub>t</sub>, &nbsp; &nbsp; w<sub>t</sub> ~ iid N<sub>p</sub>(0, Q) &perp;   v<sub>t</sub> ~ iid N<sub>q</sub>(0, R) &perp; x<sub>0</sub> ~ N<sub>p</sub>(&mu;<sub>0</sub>, &Sigma;<sub>0</sub>)

&diams; **Level 1:**  &nbsp; &nbsp;  x<sub>t</sub> = &Phi; x<sub>t-1</sub> +  &Upsilon; u<sub>t</sub> + w<sub>t</sub>,  &nbsp; &nbsp;  y<sub>t</sub> = A<sub>t</sub> x<sub>t</sub> +  &Gamma; u<sub>t</sub> + v<sub>t</sub>,  &nbsp;   u<sub>t</sub> are r-dimensional inputs, etc.<br>
     
&diams; **Level 2:** &nbsp; &nbsp; x<sub>t+1</sub> = &Phi; x<sub>t</sub> +  &Upsilon; u<sub>t+1</sub> + &Theta; w<sub>t</sub>, &nbsp;   &nbsp;   y<sub>t</sub> = A<sub>t</sub> x<sub>t</sub> +  &Gamma; u<sub>t</sub> + v<sub>t</sub>, &nbsp;  cov(w<sub>s</sub>, v<sub>t</sub>) = S &delta;<sub>s</sub><span style="position:relative; left: -.9ex; bottom: 2pt"><sup>t</sup></span>, &nbsp; &Theta; is p &times; m, and w<sub>t</sub> is m-dimensional, etc.

---  
    
### Level 0 - Fixed Measurement Matrices and No Inputs  
   
  The call to Kfilter0 is,   `Kfilter0(n,y,A,mu0,Sigma0,Phi,cQ,cR)` in fairly obvious notation   except that `cQ` and   `cR` are the Cholesky-type decompositions of 
  `Q` and `R`.  In particular `Q = t(cQ)%*%cQ` and `R = t(cR)%*%cR` is all that is required provided `Q`and `R` are valid covariance matrices (Q can be singular and there is an example in the text).   The call to `Ksmooth0` is similar. 

 __In all three cases, the smoother also returns the filter and the likelihood.__
  
---
 	
### Level 1 - Varying Measurement Matrices and Inputs
       
 The call to the filter is `Kfilter1(n,y,A,mu0,Sigma0,Phi,Ups,Gam,cQ,cR,input)` 
   where `A` is an array with `dim=c(q,p,n)`,  `Ups` is &Upsilon; [p &times; r],  `Gam`is &Gamma; [q &times; r],    and `input` is the matrix of inputs
     that has the same row dimension as y (which  is n &times; q), `input`is n &times; r; the state dimension is p).  The call to `Ksmooth1`is similar.  Set  `Ups`, `Gam`, or `input` to 0 (zero) if you  don't  use them. 
 
---  
  		
### Level 2 - Varying Measurement Matrices,  Inputs and Correlated Noise</h3>  
   
The call to the filter is `Kfilter2(n,y,A,mu0,Sigma0,Phi,Ups,Gam,Theta,cQ,cR,S,input)`, which is similar to `Kfilter1` but that `S` must be included.  `Kfilter2` runs the filter given in Property 6.5.  The call to `Ksmooth2` is similar. Set `Ups` or `Gam` or `input` to 0 (zero) if you don't  use them.   

---