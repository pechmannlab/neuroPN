# BiRewire

library(BiRewire)


data1 <- as.matrix(read.table("tmp.inp1"))
data2 <- as.matrix(read.table("tmp.inp2"))
data3 <- as.matrix(read.table("tmp.inp3"))

diag(data1) <- 0
diag(data2) <- 0
diag(data3) <- 0


cons3 <- (data1 + data2 + data3) == 3
result.rand3 <- matrix(0, ncol=ncol(data1), nrow=nrow(data1))


idx <- data.frame(read.table("tmp.idx", header=T))




#adjust python indexing
idx_chap <- idx[idx$cat=="chap",2] + 1
idx_ribo <- idx[idx$cat=="ribo",2] + 1
idx_ub <- idx[idx$cat=="ub",2] + 1
idx_syn <- idx[idx$cat=="syn",2] + 1

idx_rando <- matrix(0, nrow=100, ncol=1000)
for (i in 1:100){
  idx_rando[i,] <- sample(c(1:ncol(data1)), 1000)
}  


true_chap <- sum(cons3[idx_chap,])
true_ub <- sum(cons3[idx_ub,])
true_chapS <- sum(cons3[idx_chap,][,idx_syn])
true_chapchap <- sum(cons3[idx_chap,][,idx_chap])
true_ubS <- sum(cons3[idx_ub,][,idx_syn])
true_ubub <- sum(cons3[idx_ub,][,idx_ub])
true_ribo <- sum(cons3[idx_ribo,][,idx_ribo])

true_density_chap <- sum(cons3[idx_chap,][,idx_chap])/prod(dim(cons3[idx_chap,][,idx_chap]))
true_density_ub <- sum(cons3[idx_ub,][,idx_ub])/prod(dim(cons3[idx_ub,][,idx_ub]))
true_density_syn <- sum(cons3[idx_syn,][,idx_syn])/prod(dim(cons3[idx_syn,][,idx_syn]))





true_rando <- rep(NA, 100)
for (i in 1:100){
  true_rando[i] <- sum(cons3[idx_rando[i,],][,idx_rando[i,]])
}
true_rando <- median(true_rando)



systems_chap <- c("CCT/TRiC_system", "HSP70_system", "HSP90_system", "sHSP_system")   
systems_ub  <- c("HECT", "RBR", "RING_Complex", "RING_Single")           

get_system_true <- function(listIN, matIN){
  res <- matrix(0, nrow=2, ncol=length(listIN))
  for (i in 1:length(listIN)){
    current_sys <- listIN[i]
    current_idx <- idx[idx$system==listIN[i],2] + 1
    res[1,i] <- sum( matIN[current_idx,] )
    res[2,i] <- sum( matIN[current_idx,][,idx_syn] )
    }
  res
}

true_system_chap <- get_system_true(systems_chap, cons3)
true_system_ub <- get_system_true(systems_ub, cons3)



N <- 1000

result.ribo <- rep(0, N)
result.riboR <- rep(0,N)
result.chap <- rep(0, N)
result.ub <- rep(0,N)
result.chapS <- rep(0,N)
result.chapchap <- rep(0,N)
result.ubS <- rep(0,N)
result.ubub <- rep(0, N)
result.rando <- matrix(0,N, 100)
result.sys.chap <- matrix(0, N, 4)
result.sys.chapS <- matrix(0, N, 4)
result.sys.ub <- matrix(0, N, 4)
result.sys.ubS <- matrix(0, N, 4)


result.degree <- matrix(0, nrow(data1), N)
result.density.chap <- rep(0, N)
result.density.ub <- rep(0, N)
result.density.syn <- rep(0,N)
result.density.rando <- matrix(0,N, 100)

for (i in 1:N){
	data1.rewired <- birewire.rewire.undirected(data1, max.iter="n",accuracy=0.00005, verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
	data2.rewired <- birewire.rewire.undirected(data2, max.iter="n",accuracy=0.00005, verbose=TRUE,MAXITER_MUL=10,exact=FALSE)
	data3.rewired <- birewire.rewire.undirected(data3, max.iter="n",accuracy=0.00005, verbose=TRUE,MAXITER_MUL=10,exact=FALSE)

	data.rewired <- data1.rewired + data2.rewired + data3.rewired
	data.cons <- data.rewired == 3
  result.rand3 <- result.rand3 + data.cons
  result.degree[,i] <- rowSums(data.cons)

  result.density.chap[i] <- sum(data.cons[idx_chap,][,idx_chap]) / prod(dim(data.cons[idx_chap,][,idx_chap]))
  result.density.ub[i] <- sum(data.cons[idx_ub,][,idx_ub]) / prod(dim(data.cons[idx_ub,][,idx_ub]))
  result.density.syn[i] <- sum(data.cons[idx_syn,][,idx_syn]) / prod(dim(data.cons[idx_syn,][,idx_syn]))
  
  
	#ribo
	cons.ribo <- data.cons[idx_ribo,][,idx_ribo]
 	result.ribo[i] <- sum(cons.ribo)

	idx_rand = sample(ncol(data.cons), length(idx_ribo), replace=F)
	cons.riboR <- data.cons[idx_rand,][,idx_rand]
	result.riboR[i] <- sum(cons.riboR)


	#chap
	cons.chap <- data.cons[idx_chap,]
	result.chap[i] <- sum(cons.chap)
	result.chapS[i] <- sum(cons.chap[,idx_syn])
	result.chapchap[i] <- sum(cons.chap[,idx_chap])

	#ub
	cons.ub <- data.cons[idx_ub,]
	result.ub[i] <- sum(cons.ub)
	result.ubS[i] <- sum(cons.ub[,idx_syn])
	result.ubub[i] <- sum(cons.ub[,idx_ub])

   	#rando
	for (j in 1:100){
	    cons.rando <- data.cons[idx_rando[j,],][,idx_rando[j,]]
	    result.rando[i,j] <- sum(cons.rando)
	    
	    result.density.rando[i,j] <- sum(data.cons[idx_rando[j,],][,idx_rando[j,]]) / prod(dim(data.cons[idx_rando[j,],][,idx_rando[j,]]))
	    
	}
	res.sys.chap <- get_system_true(systems_chap, data.cons)
	result.sys.chap[i,] <- res.sys.chap[1,]
	result.sys.chapS[i,] <- res.sys.chap[2,]
	
	res.sys.ub <- get_system_true(systems_ub, data.cons)
	result.sys.ub[i,] <- res.sys.ub[1,]
	result.sys.ubS[i,] <- res.sys.ub[2,]

	print(i)
}

result.rand3 <- result.rand3 / N
write.table(result.rand3, "../data/processed/coevo/coevo_rand_3_20.txt", col.names=F, row.names=F)

colnames(true_system_chap) <- systems_chap
write.table(true_system_chap, "../data/processed/coevo/true_system_chap.txt", col.names=T, row.names=F)
colnames(result.sys.chap) <- systems_chap
write.table(result.sys.chap, "../data/processed/coevo/rand_system_chap.txt", col.names=T, row.names=F)
colnames(result.sys.chapS) <- systems_chap
write.table(result.sys.chapS, "../data/processed/coevo/rand_system_chapS.txt", col.names=T, row.names=F)


colnames(true_system_ub) <- systems_ub
write.table(true_system_ub, "../data/processed/coevo/true_system_ub.txt", col.names=T, row.names=F)
colnames(result.sys.ub) <- systems_ub
write.table(result.sys.ub, "../data/processed/coevo/rand_system_ub.txt", col.names=T, row.names=F)
colnames(result.sys.ubS) <- systems_ub
write.table(result.sys.ubS, "../data/processed/coevo/rand_system_ubS.txt", col.names=T, row.names=F)


coevoDF <- data.frame(cat=c('chap', 'chapS', 'chapchap', 'ub', 'ubS', 'ubub', 'ribo', 'rand'),
                      true=c(true_chap, true_chapS, true_chapchap, true_ub, true_ubS, true_ubub, true_ribo, true_rando),
                      rand_mean=c( mean(result.chap), mean(result.chapS), mean(result.chapchap), mean(result.ub), mean(result.ubS), mean(result.ubub), mean(result.ribo), mean(as.vector(result.rando)) ),
                      rand_median=c( median(result.chap), median(result.chapS), median(result.chapchap), median(result.ub), median(result.ubS), median(result.ubub), median(result.ribo), median(as.vector(result.rando)) ),
                      rand_95=c(quantile(result.chap, 0.95), quantile(result.chapS, 0.95), quantile(result.chapchap, 0.95), quantile(result.ub, 0.95), quantile(result.ubS, 0.95), quantile(result.ubub, 0.95), quantile(result.ribo, 0.95), quantile(as.vector(result.rando), 0.95))
                      )

write.table(coevoDF, "../data/processed/coevo/coevoDF.txt", col.names=T, row.names=F, quote=F)

write.table(result.degree, "../data/processed/coevo/cons_rand_degree.txt", col.names=F, row.names=F, quote=F)

write.table(result.density.chap, "../data/processed/coevo/rand_density_chap.txt", col.names=F, row.names=F, quote=F)
write.table(result.density.ub, "../data/processed/coevo/rand_density_ub.txt", col.names=F, row.names=F, quote=F)
write.table(result.density.syn, "../data/processed/coevo/rand_density_syn.txt", col.names=F, row.names=F, quote=F)
write.table(result.density.rando, "../data/processed/coevo/rand_density_rando.txt", col.names=F, row.names=F, quote=F)



q(save='no')
