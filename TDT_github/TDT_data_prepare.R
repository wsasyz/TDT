no_col <- max(count.fields("tdt_count_output.txt", sep = ""))
dat = read.delim("tdt_count_output.txt", sep="", col.names=1:no_col, fill=TRUE, header=FALSE)
dat = dat[,-(1:4)]
names.dat = dat[,seq(3,dim(dat)[2], by=2)]
count.dat = dat[,seq(4,dim(dat)[2], by=2)]
sum(dim(names.dat)[2] == dim(count.dat)[2])


names.vec = c("M_m+M_m=M_M", "M_m+M_m=m_m", "M_m+M_m=M_m", "M_M+M_m=M_M","M_m+M_M=M_M", "M_m+m_m=M_m", "m_m+M_m=M_m", "M_M+M_m=M_m","M_m+M_M=M_m", "M_m+m_m=m_m", "m_m+M_m=m_m")
dat.mx = matrix(0, nrow=dim(dat)[1], ncol=length(names.vec))
colnames(dat.mx) = names.vec

for (i in 1:dim(dat)[1]) {
	dum = rep(0, dim(names.dat)[2])
	for (k in 1:dim(names.dat)[2]){
		 dum[k] = as.vector(count.dat[i,k])        
	     names(dum)[k] = as.vector(names.dat[i,k])
	}
	dat.mx[i,] = dum[names.vec]    
}
dat.mx[is.na(dat.mx)] = 0

N = 25 # N families in total
dat1 = matrix(0, nrow=dim(dat)[1], ncol=8)
colnames(dat1) = c("b", "c", "(2,0)", "(0,2)", "(1,1)", "(1,0)", "(0,1)", "(0,0)")
#dat1[,2] = dat[,1]
#dat1[,1] = dat[,2]
dat1[,3:5] = dat.mx[,1:3]
dat1[,6] = rowSums(dat.mx[,4:7]) #rowSums(dat.mx[,4:7], na.rm=TRUE)
dat1[,7] = rowSums(dat.mx[,8:11]) #rowSums(dat.mx[,8:11], na.rm=TRUE)
#dat1[is.na(dat1)] = 0

dat1[,1] = 2*dat1[,3]+dat1[,5]+dat1[,6]
dat1[,2] = 2*dat1[,4]+dat1[,5]+dat1[,7]
dat1[,8] = N - rowSums(dat1[,-c(1,2)])

dat1

write.table(dat1, file="real_data_eg_TDT.txt")

save(dat1, file="real_data_eg_TDT.RData")
