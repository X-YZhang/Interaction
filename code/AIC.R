load("dat(fsnp).RData")
source("ES.ode.R")

a0<-matrix(rep(log(5000),45),ncol = 1)
pheno1<-cbind(a0,dat$ep.p)
pheno2<-cbind(a0,dat$sp.p)
times <- c(0,dat$sample_times)
mpheno1 <- colMeans(pheno1)
mpheno2 <- colMeans(pheno2)
cm<-as.numeric(c(mpheno1,mpheno2))


#Gompertz
init.par<-c(max(mpheno1),1,1)
mean.a<-optim(init.par,G.mle,y=mpheno1,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par
init.par<-c(max(mpheno2),1,1)
mean.a<-optim(init.par,G.mle,y=mpheno2,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par

g.SSR.E<-0
g.e.par<-c()
for (i in 1:45) {
  init.par<-c(24.18455757,-0.06086857,0.19754101)
  a<-optim(init.par,G.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.e.par<-rbind(g.e.par,a$par)
  g.SSR.E<-g.SSR.E+a$value
}

g.SSR.S<-0
g.s.par<-c()
for (i in 1:45) {
  init.par<-c(20.2613855,-0.3561804 ,0.2030552)
  a<-optim(init.par,G.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.s.par<-rbind(g.s.par,a$par)
  g.SSR.S<-g.SSR.S+a$value
}
g.L<-(-45/2)*(1+log(2*pi)+log(g.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(g.SSR.S/45))
g.AIC<- (12-2*g.L)/45  
g.BIC<-(-2*g.L+log(45)*6)/45
g.HQ<-(-2*g.L+log(log(45))*6)/45




#logistic
init.par<-c(max(mpheno1),1,1)
mean.a<-optim(init.par,L.mle,y=mpheno1,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par
init.par<-c(max(mpheno2),1,1)
mean.a<-optim(init.par,L.mle,y=mpheno2,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par


l.SSR.E<-0
l.e.par<-c()
for (i in 1:45) {
  init.par<-c(24.0553574,0.3695013,0.2418854)
  a<-optim(init.par,L.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.e.par<-rbind(l.e.par,a$par)
  l.SSR.E<-l.SSR.E+a$value
}

l.SSR.S<-0
l.s.par<-c()
for (i in 1:45) {
  init.par<-c(20.18992413,-0.04997506, 0.23580494)
  a<-optim(init.par,L.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.s.par<-rbind(l.s.par,a$par)
  l.SSR.S<-l.SSR.S+a$value
}
l.L<-(-45/2)*(1+log(2*pi)+log(l.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(l.SSR.S/45))
l.AIC<- (12-2*l.L)/45  
l.BIC<-(-2*l.L+log(45)*6)/45
l.HQ<-(-2*l.L+log(log(45))*6)/45



#Richards
init.par<-c(max(mpheno1),1,0,0)
mean.a<-optim(init.par,R.mle,y=mpheno1,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par
init.par<-c(max(mpheno2),0,0,1)
mean.a<-optim(init.par,R.mle,y=mpheno2,times=times,method = "Nelder-Mead",control = list(maxit=1500))
mean.a$par

r.SSR.E<-0
r.e.par<-c()
for (i in 1:45) {
  init.par<-c(24.2531956,0.3486322, 0.1788248, 0.5555942)
  a<-optim(init.par,R.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.SSR.E<-r.SSR.E+a$value
  r.e.par<-rbind(r.e.par,a$par)
}

r.SSR.S<-0
r.s.par<-c()
for (i in 1:45) {
  init.par<-c(20.21667873,0.02032237,0.20617155,0.97109937)
  a<-optim(init.par,R.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.SSR.S<-r.SSR.S+a$value
  r.s.par<-rbind(r.s.par,a$par)
}

r.L<-(-45/2)*(1+log(2*pi)+log(r.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(r.SSR.S/45))
r.AIC<- (16-2*r.L)/45  
r.BIC<-(-2*r.L+log(45)*8)/45
r.HQ<-(-2*r.L+log(log(45))*8)/45



#DDHR
init.par<-c(0.5,25,0,0.5,25,0)
lv<-optim(init.par,s.mle,s.y=cm,s.t=times,x0=cm[1],y0=cm[17],method = "Nelder-Mead",control = list(maxit=1500))
lv$par

init.par<-c(0.3334374,26.3425300,-0.1192077,0.3532136,24.9968491,-0.1770921)
allpheno<-cbind(pheno1,pheno2)
SSR<-0
lvh.par<-c()
for (i in 1:45) {
  s<-optim(init.par,s.mle,s.y=as.numeric(allpheno[i,]),s.t=times,x0=log(5000),
           y0=log(5000),method = "Nelder-Mead",control = list(maxit=1500))
  SSR<-SSR+s$value
  lvh.par<-rbind(lvh.par,s$par)
}

lvh.L<-(-45/2)*(1+log(2*pi)+log(SSR/45))
lvh.AIC<- (12-2*lvh.L)/45  
lvh.BIC<-(-2*lvh.L+log(45)*6)/45
lvh.HQ<-(-2*lvh.L+log(log(45))*6)/45

#CoCoM
source("ES.ode2.R")
init.par<-c(0.28,83.1,-2.77,0.25,62.35,-1.73)
allpheno<-cbind(pheno1,pheno2)
SSR<-0
lv.par<-c()
for (i in 1:45) {
  s<-optim(init.par,s.mle,s.y=as.numeric(allpheno[i,]),s.t=times,x0=log(5000),
           y0=log(5000),method = "Nelder-Mead",control = list(maxit=1500))
  SSR<-SSR+s$value
  lv.par<-rbind(lv.par,s$par)
}

lv.L<-(-45/2)*(1+log(2*pi)+log(SSR/45))
lv.AIC<- (12-2*lv.L)/45  
lv.BIC<-(-2*lv.L+log(45)*6)/45
lv.HQ<-(-2*lv.L+log(log(45))*6)/45

lv<-optim(init.par,s.mle,s.y=c(mpheno1,mpheno2),s.t=times,x0=log(5000),
          y0=log(5000),method = "BFGS",control = list(maxit=1500))
lv$par

##########################################################
#µ¥¶ÀÅàÑø
a0<-matrix(rep(log(5000),45),ncol = 1)
pheno1<-cbind(a0,dat$ep.pi)
pheno2<-cbind(a0,dat$sp.pi)
times <- c(0,dat$sample_times)
mpheno1 <- colMeans(pheno1)
mpheno2 <- colMeans(pheno2)


#G
g.SSR.E<-0
g.e.par<-c()
for (i in 1:45) {
  init.par<-c(24.32603628,-0.08645619,0.19270538)
  a<-optim(init.par,G.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.e.par<-rbind(g.e.par,a$par)
  g.SSR.E<-g.SSR.E+a$value
}

g.SSR.S<-0
g.s.par<-c()
for (i in 1:45) {
  init.par<-c(22.8543247,-0.2624033,0.1494925)
  a<-optim(init.par,G.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  g.s.par<-rbind(g.s.par,a$par)
  g.SSR.S<-g.SSR.S+a$value
}
g.L<-(-45/2)*(1+log(2*pi)+log(g.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(g.SSR.S/45))
g.AIC<- (12-2*g.L)/45  
g.BIC<-(-2*g.L+log(45)*6)/45
g.HQ<-(-2*g.L+log(log(45))*6)/45

#L
l.SSR.E<-0
l.e.par<-c()
for (i in 1:45) {
  init.par<-c(24.1751390,0.3280352,0.2360378)
  a<-optim(init.par,L.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.e.par<-rbind(l.e.par,a$par)
  l.SSR.E<-l.SSR.E+a$value
}

l.SSR.S<-0
l.s.par<-c()
for (i in 1:45) {
  init.par<-c(22.72253934,0.07036654,0.17656100)
  a<-optim(init.par,L.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  l.s.par<-rbind(l.s.par,a$par)
  l.SSR.S<-l.SSR.S+a$value
}
l.L<-(-45/2)*(1+log(2*pi)+log(l.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(l.SSR.S/45))
l.AIC<- (12-2*l.L)/45  
l.BIC<-(-2*l.L+log(45)*6)/45
l.HQ<-(-2*l.L+log(log(45))*6)/45

#R
r.SSR.E<-0
r.e.par<-c()
for (i in 1:45) {
  init.par<-c(25.21550701,0.93135155, 0.09263587, -1.49180681)
  a<-optim(init.par,R.mle,y=as.numeric(pheno1[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.SSR.E<-r.SSR.E+a$value
  r.e.par<-rbind(r.e.par,a$par)
}

r.SSR.S<-0
r.s.par<-c()
for (i in 1:45) {
  init.par<-c(22.6694252775,-0.0000536573, 0.1495047723, 1.0000715003)
  a<-optim(init.par,R.mle,y=as.numeric(pheno2[i,]),times=times,method = "Nelder-Mead",control = list(maxit=1500))
  r.SSR.S<-r.SSR.S+a$value
  r.s.par<-rbind(r.s.par,a$par)
}

r.L<-(-45/2)*(1+log(2*pi)+log(r.SSR.E/45))+(-45/2)*(1+log(2*pi)+log(r.SSR.S/45))
r.AIC<- (16-2*r.L)/45  
r.BIC<-(-2*r.L+log(45)*8)/45
r.HQ<-(-2*r.L+log(log(45))*8)/45


tss.E<-sum( (pheno1 - mpheno1 )^2 )
tss.S<-sum( (pheno2 - mpheno2 )^2 )
R2.E <- 1-r.SSR.E/tss.E 
R2.S <- 1-r.SSR.S/tss.S


