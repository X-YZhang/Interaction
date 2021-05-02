

com.get_mu_mono <- function(par, times,options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2]);
  }
  
  state0 <- c(X=8.517193);
  y <- COMP.f( par0, state0, times );
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2] ) );
}

COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-X/k1) 
            
            list(c(dX))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


curve.mlefunc<-function( par,y1,time.std)
{
  len.cov <- 2
  par.covar <- par[1:len.cov]
  n  <- length(y1[,1])
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  AR1 <- AR1.get_inv_mat(par.covar,time.std)
  sigma <- solve( AR1 )
  
  curve.par <- par[(len.cov+1):(len.cov+ 2)]
  mu <- com.get_mu_mono(curve.par, time.std)
  
  yy <- y1
  for ( i in 1:dim(y1)[2] )
  {
    y1.miss <- which( is.na(y1[,i]) )
    yy[y1.miss, i]<- mu[i]
  }
  fy <- dmvnorm(yy,mu,sigma)
  #fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps
  A <- -sum(log(fy))
  #cat("LL=",A,"\n")
  return(A)
}

double.H0 <- function(ny,times){
  
  cm <- as.numeric(colMeans(ny))
  sum.curve<-function(par10){
    A <- sum((cm-com.get_mu_mono(par10,times))^2)
    A
  }
  init.par<-c(0.5,max(cm))
  a<-optim(init.par,sum.curve,method = "Nelder-Mead",control = list(maxit=1500))
  
  covar.par <- c(0.1,sd(ny[,7]))
  parin <- c(covar.par,a$par)
  time.std <- times
  
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[3:4])
      AA <- curve.mlefunc(nnpar,y1=ny,time.std)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "Nelder-Mead",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mlefunc(nnpar,y1=ny,time.std)
      AA
    }
    r1 <- optim(c(parin[3:4]),mle.1,method = "Nelder-Mead",control=list(maxit=32000))    
    new1 <- r1$par
    #cat("new1:",unlist( new1), "\n");
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mlefunc(parin,y1=ny,time.std)
  return(c(LL,parin))
}

double.est <- function(dat,interval=c(1,10)){
  
  snp <- dat$s.snp
  nn <- dim(snp)[1]
  y <- dat$sp.pi
  times <- dat$sample_times
  H0 <- double.H0(y,times)
  par1 <- H0[-1]
  
  n1 <- interval[1]
  n2 <- interval[2]
  if(n2 >nn)
    n2 <- nn
  p.res <- matrix(NA,nrow=length(c(n1:n2)),ncol=100)
  for(i in n1:n2){
    nsnp <- snp[i,]
    ny <- y
    st <- table(as.matrix(nsnp))
    H1 <- double.H1(par1,ny,nsnp,times)
    LR <- 2*(H0[1]-H1[1])
    p.tmp <- c(LR,H0[1],H1)
    cat("snp", i, "=", p.tmp, "\n");
    p.res[(i-(n1-1)),1:length(p.tmp)] <- p.tmp
  }
  return(p.res)
}


double.H1 <- function(par1,ny,nsnp,times){
  
  snp1 <- nsnp
  mean.T <- times
  par.covar <- par1[1:2]
  
  index <- table(as.matrix(snp1))
  snp.type <- as.character(names(index))
  
  g.par <- c()
  SNP.index <- list()
  for(i in 1:length(snp.type)){
    SNP.n <- which(snp1==snp.type[i])
    SNP.p <- c(colMeans(ny[SNP.n,],na.rm=T))
    sum.curve<-function(par10){
      A <- sum((SNP.p-com.get_mu_mono(par10,mean.T))^2)
      A
    }
    s.m <- optim(par1[3:4],sum.curve,method="BFGS",control=list(maxit=32000))
    g.par <- c(g.par,s.m$par)
    SNP.index[[i]] <- SNP.n
  }
  parin <- c(par.covar,g.par)
  n.par <- length(parin)
  
  loop_k <- 1;
  max_iter <- 100;
  epsi <- 10^-4;
  max_err <- 1;
  
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin[3:n.par])
      AA <- double.mlefunc(nnpar,y1=ny,time.std=mean.T,SNP.index,snp.type)
      AA
    }
    r1.covar <- optim(parin[1:2],mle.covar1,method = "Nelder-Mead",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle1.g <- function(npar){
      nnpar <- c(new.covar1,npar)
      AA <- double.mlefunc(nnpar,y1=ny,time.std=mean.T,SNP.index,snp.type)
      AA
    }
    r1.g <- optim(c(parin[3:n.par]),mle1.g,method = "BFGS",control=list(maxit=32000))
    #cat("r1.g:",unlist(r1.g$par), "\n");
    
    newpar <- c(new.covar1,r1.g$par)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin <- newpar
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  
  LL <- double.mlefunc(parin,y1=ny,time.std=mean.T,SNP.index,snp.type)
  
  return(c(LL,newpar))
}


double.mlefunc <- function(par,y1,time.std,SNP.index=SNP.index,snp.type=snp.type)
{
  n  <- length(y1[,1])
  
  len.cov <- 2
  par.covar <- par[1:len.cov]
  if(par.covar[1]>1||par.covar<0)
    return(NaN)
  ar1 <- AR1.get_inv_mat(par.covar,time.std)
  sigma <- solve( ar1)
  len.gen <- 2
  len <- 0
  A1 <- c()
  for(i in 1:length(snp.type)){
    mu.g <- par[(len.cov+len+1):(len.cov+len+len.gen)]
    yy1 <- y1[SNP.index[[i]],]
    mu <- com.get_mu_mono(mu.g, time.std)
    nyy1 <- yy1
    fy1 <- dmvnorm( nyy1, mu, sigma)
    #fy1[which(fy1<=.Machine$double.eps)] <- .Machine$double.eps
    A1 <- c(A1,-sum(log(fy1)))
    len <- len + len.gen
  }
  A <- sum(A1)
  #cat("LL=",A,"\n")
  return (A);
}



