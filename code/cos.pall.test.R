source("double.covar.R")
y=dat$sp.p
times <- dat$sample_times
H0 <- double.H0(y,times)
par1 <- H0[-1]

#h1<-double.est(snp=dat$s.snp[1,],ny=dat$sp.p)
library("parallel")
core.number<-detectCores()
cl<-makeCluster(getOption("cl.cores",core.number))
clusterExport(cl,c("H0","par1","dat","com.get_mu_co","COMP.f.co","double.H1",
                   "double.mlefunc","double.est","times","AR1.get_mat","AR1.get_inv_mat"))
clusterEvalQ(cl, expr = library(mvtnorm))
clusterEvalQ(cl, expr = library(deSolve))
ret<-t(parSapply(cl,501:1000,function(c) double.est(snp=dat$s.snp[c,],ny=dat$sp.p)))
stopCluster(cl)

save(ret,file = "./optim.co.s/ret5-10.RData")
