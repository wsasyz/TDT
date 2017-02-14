rm(list=ls())

source("DP_TDT_fn.R")

###### to read TDT data
# dat = read.table("DP_TDT_data.txt")
# colnames(dat) = c("b", "c", "(2,0)", "(0,2)", "(1,1)", "(1,0)", "(0,1)", "(0,0)")
# head(dat)
# N = 25 # family number
# M = dim(dat) # SNP number

# ##### to generate a simulated data
# N = 150
# M = 5000
# n = 2*N
# M.1 = 10 
# bc.sum = sample(0:n, M, replace=TRUE)
# dat.b = rbinom(M, bc.sum, 0.5)
# dat.c = bc.sum - dat.b
# ind = order(bc.sum, decreasing=TRUE)[1:M.1]
# dat.b[ind] = rbinom(M.1, bc.sum[ind], 0.65)
# dat.c[ind] = bc.sum[ind] - dat.b[ind]
# dat.stat = (dat.b - dat.c)^2/(dat.b + dat.c)
# dat.stat[is.na(dat.stat)] = 0
# dat = cbind(dat.b, dat.c)
# dat[dat[,1] < dat[,2],] = dat[dat[,1] < dat[,2], c(2,1)]
# colnames(dat) = c("b", "c")
# rownames(dat) =  as.character(dat[,1]*n + dat[,2])


########## one example to apply DP TDT to the simulated data
load("sim_data_small_cohort.RData")

K = 1 # number of significant SNPs to release
eps = 3 # privacy budget
dat.stat = (dat[,1] - dat[,2])^2/(dat[,1] + dat[,2])
dat.stat[is.na(dat.stat)] = 0 # test statistics
ind.sig.dat = order(dat.stat, decreasing=TRUE)[1:K]

### Laplace mechanism based on the test statistic
lap.scale = (2*K*8*(N-1)/N)/eps
stat.lap = dat.stat + rLap(length(dat.stat), 0, 1/lap.scale) # test statistic + laplace noise
ind.sig.lap = order(stat.lap, decreasing=TRUE)[1:K]
utility.lap.stat = length(intersect(ind.sig.dat, ind.sig.lap))/length(ind.sig.dat)
utility.lap.stat

### exponential mechanism based on test statistic
score = dat.stat
sen.score = 8*(N-1)/N
ind.sig.expon = expon.set.sig(K, eps, score, sen.score)
utility.exp.stat = length(intersect(ind.sig.dat, ind.sig.expon))/length(ind.sig.dat)
utility.exp.stat
		
		
### exponential mechanism based on the shorest hamming distance score from approximation algorithm
thres = qchisq(1 - 0.05/M, df=1)
dat.uq.score = SHD_fn(dat, N, thres) 
dat.shd = matrix(0, nrow=nrow(dat), ncol=3)
dat.shd[,1:2] = dat
rownames(dat.shd) = rownames(dat)
score.nm = rownames(dat.uq.score)[match(rownames(dat), rownames(dat.uq.score))]
dat.shd[,3] = dat.uq.score[score.nm,1]
score = dat.shd[,3]
sen.score = 1
ind.sig.expon.ham.grad = expon.set.sig(K, eps, score, sen.score)
utility.shd.apprx = length(intersect(ind.sig.dat, ind.sig.expon.ham.grad))/length(ind.sig.dat)
utility.shd.apprx

### to run the code "SHD_score_exact_alg.java" to get the SHD scores from the exact algoirithm
# in the simulated data, since we only summarize the data to (b,c) pairs, to run the exact algorithm, we randomly generate one combination of n's satisfying the equations below table 1 in the paper.
dat.exact = read.table("sim_SHD_score_exact.txt", sep=",")
score.exact = dat.exact[,ncol(dat.exact)] # the SHD scores from exact algorithm
sen.score = 1
ind.sig.expon.ham.exact = expon.set.sig(K, eps, score, sen.score)
utility.shd.exact = length(intersect(ind.sig.dat, ind.sig.expon.ham.exact))/length(ind.sig.dat)
utility.shd.exact 	
	  	


### Laplace mechanisms based on p-value and truncated p-values
dat.pval = 1 - pchisq(dat.stat, df=1) # p-values
thres.pval = 0.05/M
dat.pval.trunc = pmin(dat.pval, thres.pval) # truncated p-values	

sen = pchisq(4, df=1)
lap.scale = (2*K*sen)/eps
pval.lap = dat.pval + rLap(length(dat.pval), 0, 1/lap.scale) # test statistic + laplace noise
ind.sig.lap = order(pval.lap, decreasing=FALSE)[1:K]	
utility.lap.pval = length(intersect(ind.sig.dat, ind.sig.lap))/length(ind.sig.dat)
utility.lap.pval

thres.t = qchisq(1-thres.pval, df=1)
sen = abs(1 - pchisq((thres.t-4)^2/thres.t, df=1) -thres.pval)
lap.scale = (2*K*sen)/eps
pval.lap = dat.pval.trunc + rLap(length(dat.pval.trunc), 0, 1/lap.scale) # test statistic + laplace noise
ind.sig.lap = order(pval.lap, decreasing=FALSE)[1:K]	
utility.lap.trunc.pval = length(intersect(ind.sig.dat, ind.sig.lap))/length(ind.sig.dat)
utility.lap.trunc.pval

