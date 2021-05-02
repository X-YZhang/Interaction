
init.est<-function(dat){
  times <- dat$sample_times
  epp<-dat$ep.p
  spp<-dat$sp.p
  ty <- cbind(epp,spp)
  cm <- as.numeric(colMeans(ty))
  sum.curve<-function(par6){
    y<-com.get_mu1(par6[c(1,2,4,5)],times,x0=cm[1]/2,y0=cm[17]/2)+
      com.get_mu2(par6[c(1,3,4,6)],times,x0=cm[1]/2,y0=cm[17]/2)
    A <- sum((cm-y)^2)
    A
  }
  init.par<-c(0.5,23,-2,0.5,20,-2)
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=1500))
  return(a$par)
}

#par0<-init.est(dat)#0.4656164 22.6689926 -1.9603047  0.4641400 20.7590245 -1.9604165

curve.mlefunc<-function( par, y, time.std,x0,y0)
{
  len.cov <- 5;
  par.covar <- par[1:len.cov];
  sig <- SAD3.get_mat(par.covar, times=1:length(time.std),traits = 2);
  
  m  <- length(y[1,]);
  n  <- length(y[,1]);
  
  par1 <- par[-c(1:len.cov)]
  mu0 <- com.get_mu1(par1[c(1,2,4,5)],times,x0,y0)+
    com.get_mu2(par1[c(1,3,4,6)],times,x0,y0)
  
  fy0<- dmvnorm(y,mu0,sig)
  A <- -sum(log(fy0));
  
  return (A);
}

H0.est<-function(dat){
  epp<-dat$ep.p
  spp<-dat$sp.p
  allpheno <- cbind(epp,spp)
  
  times<-dat$sample_times
  me <- colMeans(epp)
  ms <- colMeans(spp)
  mpheno<-as.numeric(c(me,ms))
  
  par<-c(0.4656164, 22.6689926, -1.9603047,  0.4641400, 20.7590245, -1.9604165)
  covar.par<-c(0.5,sd(as.matrix(epp)),0.5,sd(as.matrix(spp)),0.5)
  parin<-c(covar.par,par)
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  while(loop_k<max_iter && max_err>epsi){
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[6:11])
      AA <- curve.mlefunc(nnpar,y=allpheno,time.std =times,x0=mpheno[1]/2,y0=mpheno[17]/2)
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mlefunc(nnpar,y=allpheno,time.std =times,x0=mpheno[1]/2,y0=mpheno[17]/2)
      AA
    }
    r1 <- optim(c(parin[6:11]),mle.1,method = "Nelder-Mead",control=list(maxit=32000))    
    new1 <- r1$par
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    
    loop_k <- loop_k+1; 
  }
  LL <- curve.mlefunc(parin2,y=allpheno,time.std =times,x0=mpheno[1]/2,y0=mpheno[17]/2)
  return(c(parin2,LL))
}
# H0<-H0.est(dat)
# load("H0.RData")

ES.mlefunc<-function( par, y, time.std,snp.index )
{
  len.cov <- 5;
  par.covar <- par[1:len.cov];
  sig <- SAD3.get_mat(par.covar, time.std,traits=2);
  
  m  <- length(y[1,]);
  n  <- length(y[,1]);
  len <- 0
  A <- 0
  for(k in 1:length(snp.index)){
    gen_par<- par[(len.cov+1+len):(len.cov+ 6+len)];
    yy0 <- y[snp.index[[k]],]
    mu0 <- com.get_mu1(gen_par[c(1,2,4,5)],time.std,4.258597,4.258597)+
      com.get_mu2(gen_par[c(1,3,4,6)],time.std,4.258597,4.258597);
    fy0 <- dmvnorm(yy0,mu0,sig)
    A <- A -sum(log(fy0));
    len <- len + 6
  }
  return (A);
}

H1.est <- function(y11,y22,SNP1,init.par=h01,times){
  index <- table(SNP1)
  snp.type <- as.character(names(index))
  phenos<-as.matrix( cbind(y11,y22) )
  
  g.par <- c()
  snp.index <- list()
  
  for(j in 1:length(snp.type)){
    index <- which(SNP1==snp.type[j])
    yy <- as.numeric(colMeans(phenos[index,]))
    sum.curve<-function(par6){
      y<-com.get_mu1(par6[c(1,2,4,5)],times,x0=4.258597,y0=4.258597)+
        com.get_mu2(par6[c(1,3,4,6)],times,x0=4.258597,y0=4.258597)
      A <- sum((yy-y)^2)
      A
    }
    r1<-optim(H0[6:11],sum.curve,method = "Nelder-Mead",control = list(maxit=1500))
    snp.index[[j]] <- index
    g.par <- c(g.par,r1$par)
  }
  parin <- c(H0[1:5],g.par)
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-5;
  max_err <- 1;
  while(loop_k<max_iter && max_err>epsi){
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[-(1:5)])
      AA <- ES.mlefunc(nnpar,y=phenos,time.std =times,snp.index)
      AA
    }
    r1.covar <- optim(parin[1:5],mle.covar1,method = "BFGS",control=list(maxit=2000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- ES.mlefunc(nnpar,y=phenos,time.std =times,snp.index)
      AA
    }
    r1.g <- optim(c(parin[-(1:5)]),mle1.g,method = "Nelder-Mead",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  return(c(r1.g$value,newpar))
}


com.DH.est1 <- function(dat,interval=c(1,40054)){
  
  y1 <- as.matrix( dat$ep.p)
  y2 <- as.matrix( dat$sp.p)
  
  times <- dat$sample_times
  geno_table <- t(dat$fsnp)
  nm <- dim(geno_table)[1]
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nm)
    n2 <- nm
  res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)
  
  for(i in n1:n2){
    SNP <- geno_table[i,]
    # missing <- which(is.na(SNP))
    # if ( length(missing) > 0)
    # {
    #   SNP1 <- SNP[-(missing)]
    #   y11 <- y1[ -(missing), ]
    #   y22 <- y2[ -(missing), ]
    # }else{
    #   SNP1 <- SNP
    #   y11 <- y1
    #   y22 <- y2
    # }
    
    h01 <- H0#没有NA不用重复计算
    
    h02 <- try(H1.est(y1,y2,SNP,init.par=h01,times),TRUE)
    if (class(h02) == "try-error") 
      h02 <- NA
    LR <- 2*(h01[12]-h02[1])
    if(is.na(h01)||is.na(h02)){
      allpar <- c(LR,rep(NA,25))
    }else{
      allpar <- c(LR,h02)
    }
    
    cat("snp", i, "=", allpar, "\n");
    res[(i-(n1-1)),(1:length(allpar))] <- allpar
    #save(res,file = "res.RData")
  }
  return(res)
}


