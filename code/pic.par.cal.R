#Figure parameters calculate
mono.figure<-function(dat,pheno=dat$ep.pi){
  
  times <- c(0,dat$sample_times)
  N <- dat$sample_N
  mpheno <- c(log(5000),colMeans(pheno))
  
  init.par<-c(max(mpheno),0.5,0.5)
  a<-optim(init.par,G.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.par<-a$par
  
  init.par<-c(max(mpheno),0.5,0.5)
  b<-optim(init.par,L.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.par<-b$par
  
  init.par<-c(max(mpheno),0,0,0)
  c<-optim(init.par,R.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.par<-c$par
  
  par(mar=c(4,4,2,1))
  plot(dat$sample_times,pheno[1,],pch=1,cex=1.5,col="light green",xlab="Time(h)",
       ylab="Abundance",cex.lab=1.5,xaxt="n",yaxt="n",
       xlim=c(-3,37),ylim=c(4,33))
  axis(1,c(0,8,16,24,32),c(0,8,16,24,32),cex.axis=1.2)
  axis(2,c(5,10,15,20,25,30),c(5,10,15,20,25,30),cex.axis=1.2)
  for(i in 2:N){
    points(dat$sample_times,pheno[i,],pch=1,cex=1.5,col="light green")
  }
  points(0,log(5000),pch=1,cex=1.5,col="light green")
  lines(times,G.get_mu(g.par, times),type="l",lwd=2,col="darkgreen")
  lines(times,L.get_mu(l.par,times),type="l",lwd=2,col="red")
  lines(times,R.get_mu(r.par,times),type="l",lwd=2,col="blue")
  mtext(expression(italic("Monoculture")),side=3, line=0.4,cex=1.3,adj=0.5)
  abline(h=max(R.get_mu(r.par,times))+0.5,lty=2,col="blue",lwd=1.5)
}

mono.figure<-function(dat,pheno=dat$sp.pi){
  
  times <- c(0,dat$sample_times)
  N <- dat$sample_N
  mpheno <- c(log(5000),colMeans(pheno))
  
  init.par<-c(max(mpheno),0.5,0.5)
  a<-optim(init.par,G.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.par<-a$par
  
  init.par<-c(max(mpheno),0.5,0.5)
  b<-optim(init.par,L.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.par<-b$par
  
  init.par<-c(max(mpheno),0,0,1)
  c<-optim(init.par,R.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.par<-c$par
  
  par(mar=c(4,4,2,1))
  plot(dat$sample_times,pheno[1,],pch=1,cex=1.5,col="light green",xlab="Time(h)",
       ylab="Abundance",cex.lab=1.5,xaxt="n",yaxt="n",
       xlim=c(-3,37),ylim=c(4,33))
  axis(1,c(0,8,16,24,32),c(0,8,16,24,32),cex.axis=1.2)
  axis(2,c(5,10,15,20,25,30),c(5,10,15,20,25,30),cex.axis=1.2)
  for(i in 2:N){
    points(dat$sample_times,pheno[i,],pch=1,cex=1.5,col="light green")
  }
  points(0,log(5000),pch=1,cex=1.5,col="light green")
  lines(times,G.get_mu(g.par, times),type="l",lwd=2,col="darkgreen")
  lines(times,L.get_mu(l.par,times),type="l",lwd=2,col="red")
  lines(times,R.get_mu(r.par,times),type="l",lwd=2,col="blue")
  mtext(expression(italic("Monoculture")),side=3, line=0.4,cex=1.3,adj=0.5)
  abline(h=max(R.get_mu(r.par,times))+0.5,lty=2,col="blue",lwd=1.5)
}

co.e.figure<-function(dat,pheno=dat$ep.p){
  
  times <- c(0,dat$sample_times)
  N <- dat$sample_N
  a0<-matrix(rep(log(5000),45),ncol = 1)
  pheno1<-cbind(a0,dat$ep.p)
  pheno2<-cbind(a0,dat$sp.p)
  mpheno <- colMeans(pheno1)
  
  init.par<-c(max(mpheno),1,1)
  a<-optim(init.par,G.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.par<-a$par
  
  b<-optim(init.par,L.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.par<-b$par
  
  init.par<-c(max(mpheno),1,0,0)
  c<-optim(init.par,R.mle,y=mpheno,times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.par<-c$par
  
  ty <- cbind(pheno1,pheno2)
  cm <- as.numeric(colMeans(ty))
  init.par<-c(0.5,25,0,0.5,25,0)
  lv<-optim(init.par,s.mle,s.y=cm,s.t=times,x0=cm[1],y0=cm[17],method = "Nelder-Mead",control = list(maxit=1500))
  lv.par<-lv$par
  
  plot(dat$sample_times,pheno[1,],pch=1,cex=1.5,col="light green",xlab="Time(h)",
       ylab="Abundance",cex.lab=1.5,xaxt="n",yaxt="n",
       xlim=c(-3,37),ylim=c(4,33))
  axis(1,c(0,8,16,24,32),c(0,8,16,24,32),cex.axis=1.2)
  axis(2,c(5,10,15,20,25,30),c(5,10,15,20,25,30),cex.axis=1.2)
  abline(h=log(5000),lty=2,col="blue",lwd=1.5)
  for(i in 2:N){
    points(dat$sample_times,pheno[i,],pch=1,cex=1.5,col="light green")
  }
  points(0,log(5000),pch=1,cex=1.5,col="light green")
  lines(times,G.get_mu(g.par, times),type="l",lwd=2,col="darkgreen")
  lines(times,L.get_mu(l.par,times),type="l",lwd=2,col="red")
  lines(times,R.get_mu(r.par,times),type="l",lwd=2,col="purple")
  a1<-com.get_mu(lv.par,times,x0=cm[1],y0=cm[17])[1:length(times)]
  a2<-com.get_mu1(lv.par[c(1,2,4,5)],times,x0=cm[1],y0=cm[17])[1:length(times)]
  a3<-com.get_mu2(lv.par[c(1,3,4,6)],times,x0=cm[1],y0=cm[17])[1:length(times)]
  lines(times,a1,type="l",lwd=2,col="blue")
  lines(times,a2,type="l",lty=5,lwd=2,col="blue")
  lines(times,a3,type="l",lty=3,lwd=2,col="blue")
  
  
  mtext(expression(italic("co-culture")),side=3, line=0.4,cex=1.3,adj=0.5)
}


co.s.figure<-function(dat,pheno=dat$sp.p){
  
  times <- c(0,dat$sample_times)
  N <- dat$sample_N
  a0<-matrix(rep(log(5000),45),ncol = 1)
  pheno1<-cbind(a0,dat$ep.p)
  pheno2<-cbind(a0,dat$sp.p)
  mpheno <- colMeans(pheno2)
  
  init.par<-c(max(mpheno),1,1)
  a<-optim(init.par,G.mle,y=mpheno,times=c(0,dat$sample_times),method = "Nelder-Mead",control = list(maxit=1500))
  g.par<-a$par
  
  b<-optim(init.par,L.mle,y=mpheno,times=c(0,dat$sample_times),method = "Nelder-Mead",control = list(maxit=1500))
  l.par<-b$par
  
  init.par<-c(max(mpheno),0,0,1)
  c<-optim(init.par,R.mle,y=mpheno,times=c(0,dat$sample_times),method = "Nelder-Mead",control = list(maxit=1500))
  r.par<-c$par
  
  ty <- cbind(pheno1,pheno2)
  cm <- as.numeric(colMeans(ty))
  init.par<-c(0.5,25,0,0.5,25,0)
  lv<-optim(init.par,s.mle,s.y=cm,s.t=c(0,dat$sample_times),x0=cm[1],y0=cm[17],method = "Nelder-Mead",control = list(maxit=1500))
  lv.par<-lv$par
  
  plot(dat$sample_times,pheno[1,],pch=1,cex=1.5,col="light green",xlab="Time(h)",
       ylab="Abundance",cex.lab=1.5,xaxt="n",yaxt="n",
       xlim=c(-3,37),ylim=c(4,28))
  axis(1,c(0,8,16,24,32),c(0,8,16,24,32),cex.axis=1.2)
  axis(2,c(5,10,15,20,25,30),c(5,10,15,20,25,30),cex.axis=1.2)
  abline(h=log(5000),lty=2,col="blue",lwd=1.5)
  for(i in 2:N){
    points(dat$sample_times,pheno[i,],pch=1,cex=1.5,col="light green")
  }
  points(0,log(5000),pch=1,cex=1.5,col="light green")
  lines(times,G.get_mu(g.par, times),type="l",lwd=2,col="darkgreen")
  lines(times,L.get_mu(l.par,times),type="l",lwd=2,col="red")
  lines(times,R.get_mu(r.par,times),type="l",lwd=2,col="purple")
  a1<-com.get_mu(lv.par,times,x0=cm[1],y0=cm[17])[(length(times)+1):(2*length(times))]
  a2<-com.get_mu1(lv.par[c(1,2,4,5)],times,x0=cm[1],y0=cm[17])[(length(times)+1):(2*length(times))]
  a3<-a1-a2+cm[1]
  lines(times,a1,type="l",lwd=2,col="blue")
  lines(times,a2,type="l",lty=5,lwd=2,col="blue")
  lines(times,a3,type="l",lty=3,lwd=2,col="blue")
  
  
  mtext(expression(italic("co-culture")),side=3, line=0.4,cex=1.3,adj=0.5)
}
