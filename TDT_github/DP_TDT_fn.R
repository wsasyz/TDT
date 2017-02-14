
############ to get test statistics
stats = function(D) {
	    b = D[1]
	    c = D[2]
	    if( b == 0 & c == 0){
					T_D = 0
				} else {
					T_D = (b - c)^2/(b+c) # the statistics at D
				}
}


########### to get shortest Hamming distance scores from approximation algorithm
dev = function(D) {	  
	  b = D[1]
	  c = D[2]
	  if( b == 0 & c == 0){
				dev_D = c(0, 0)
			} else {
				dev_D = c(2*(b-c)/(b+c) - (b-c)^2/(b+c)^2, -2*(b-c)/(b+c) - (b-c)^2/(b+c)^2) # the gradient at D
        }

}



direc_upp = rbind(c(1, 1), c(-1, -1), c(0,1), c(-1, 0), c(0, 2), c(-2, 0), c(-1, 2), c(-2, 1), c(-1, 1), c(-2, 2)) # the directions for c > b
direc_low = rbind(c(1, 1), c(-1, -1), c(1,0), c(0,-1), c(2,0), c(0,-2), c(2,-1), c(1,-2), c(1,-1), c(2, -2))


HM.score.grad = function(b_0, c_0, n, t, direc_upp, direc_low){
	
		if (b_0 > c_0){ # since the searching region is symmetric about b=c, we only consider the upper region
			dum = b_0
			b_0 = c_0
			c_0 = dum
		} 
		D0 = c(b_0, c_0)
		T_D0 = stats(D0)
		path_stats = T_D0
		path_direc = rep(NA, 2)
		
		D = D0
		T_D = T_D0
		num_step = 0
		
		if (T_D0 < t){
			 direc = direc_upp # the directions from insignificant to significant
			 
			 while (T_D < t) {
			    dev_D = dev(D)
			    D_cand = D + t(direc) # D_1 candidates along all directions
		        direc_sel = direc[(colSums(D_cand) <= rep(n, ncol(D_cand))) & (D_cand[1,] >= rep(0, ncol(D_cand))) & (D_cand[2,] >= rep(0, ncol(D_cand))),] # select the valid moving directions such that in D_cand, b + c <= n, b >= 0, c >= 0
		        
		        if (dev_D[1] == 0 & dev_D[2] == 0){
			        D_cand = D + t(direc_sel)
			        grad_direct = numeric(ncol(D_cand))
			        for(j in 1:ncol(D_cand)){
			        	   dev_D = dev(D_cand[,j])
			        	   grad_direct[j] = crossprod(dev_D, -D_cand[,j])
			        }
			        direc_D1 = direc_sel[which.min(grad_direct), ]
		        } else {
		        	 direc_D1 = direc_sel[which.max(crossprod(dev_D, t(direc_sel))),] # the moving direction for the next step is the one maximizing inner product of the gradient and valid moving directions
		        }
		        		        
		        D = D + direc_D1
		        T_D = stats(D)
		        num_step = num_step + 1
		        path_direc = rbind(path_direc, direc_D1)
		        path_stats = c(path_stats, T_D)
		        #print(direc_D1)
		        #print(T_D)
		    }
			 
		} else {
			 direc = direc_low # the directions from significant to insignificant
		
			 while (T_D >= t) {
			 	   dev_D = dev(D)
			 	   D_cand = D + t(direc) # D_1 candidates along all directions
		        direc_sel = direc[(colSums(D_cand) <= rep(n, ncol(D_cand))) & (D_cand[1,] >= rep(0, ncol(D_cand))) & (D_cand[2,] >= rep(0, ncol(D_cand))),] # select the valid moving directions such that in D_cand, b + c <= n, b >= 0, c >= 0
		        
		        direc_D1 = direc_sel[which.min(crossprod(dev_D, t(direc_sel))),] # the moving direction for the next step is the one maximizing inner product of the gradient and valid moving directions
		        #print(direc_D1)
		        D = D + direc_D1
		        T_D = stats(D)
		        num_step = num_step + 1
		        path_direc = rbind(path_direc, direc_D1)
		        path_stats = c(path_stats, T_D)
		        #print(T_D)
		    }
		}
		
		if (T_D0 < t) {
			HM.score = -num_step 
		} else{
			 HM.score = num_step- 1
	          }
		return(HM.score)
}

#### to get the SHD approximate scores
SHD_fn = function (dat, N, thres) {
	  n = 2*N
	  A = dat[unique(rownames(dat)), ]
	  score.grad_A = numeric(nrow(A))
	   for(i in 1:nrow(A)){
	   	    #print(i)
			score.grad_A[i] = HM.score.grad(A[i,1], A[i,2], n, thres, direc_upp, direc_low)
	   }
	   score.grad_A = data.frame(score.grad_A)
	   rownames(score.grad_A) = rownames(A)
       return(score.grad_A)
}



########### to get significant index set from exponential mechanism
expon.set.sig = function(K, eps, score, sen.score){
	M = length(score)
	set.sig = rep(NA, K)
    w.score.0 = exp(eps*score/(2*K*sen.score))
    ind.inf = (1:M)[(w.score.0 == "Inf")]
    
    if (length(ind.inf) > K) {
    	return("set a smaller epsilon")
    } else {
    	w.score = w.score.0
	   	for (i in 1:K){
	   		if ( length(ind.inf) > 0) {
		    	w.score = rep(0, M)
		    	w.score[ind.inf] = 1
		     } 
			 set.sig[i] = sample(1:M, size=1, prob=w.score)
			 w.score.0[set.sig[i]] = 0
			 ind.inf = (1:M)[(w.score.0 == "Inf")]
			 w.score = w.score.0
		}
		return(set.sig) # return significant index of score
    }
 }


# to generate Laplace noise
rLap = function(n, a, b){
	    # n is sample size
	    # a is location parameter
	    # b is scale parameter
	    
	    X = sample(c(-1, 1), n, replace = TRUE) * rexp(n) # double exponential r.v. with mean 0, scale =1 and its density is f(x) = 1/2 exp(-|x|)
        X = (X + a)/b
	    return(X) # f_{(X+a)/b}(x) = b/2 exp(-|bx-a|), mean a/b and scale = 1/b
}
      
      
