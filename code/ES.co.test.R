setwd("G:/Micro-code")
library(mvtnorm)
source("ES.co.load.R")
source("ES.co.util.R")
source("ES.cov.R")
source("ES.covar.R")


# dat <- ES.co.load(e.file="./data/coculture-Ecoli.csv",s.file="./data/coculture-Saureus.csv",
#                e.snp.file="./data/gen_ecoli.txt",s.snp.file="./data/gen_saureus.txt",
#                ei.file="data/monocultrue-Ecoli.csv",si.file="./data/monocultrue-Saureus.csv")
# 
# dat1 <- ES.adjust.s(dat,eefile="Structure/E.coli/E/Eput_simple.5.meanQ",ssfile="Structure/S.aureus/S/Sput_simple.8.meanQ")
# 
# E.i <- ES.cov.test(pheno=dat1$ep.pi,snp=dat1$e.snp) #检验基因对monoculture的丰度作用
# E <- ES.cov.test(pheno=dat1$ep.p,snp=dat1$e.snp)
# S.i <- ES.cov.test(pheno=dat1$sp.pi,snp=dat1$s.snp)
# S <- ES.cov.test(pheno=dat1$sp.p,snp=dat1$s.snp)
# 
# save(E.i,file="E-i.RData"); save(E,file="E.RData");
# save(S.i,file="S-i.RData"); save(S,file="S.RData");
# 
# load("E-i.RData")
# load("E.RData")
# load("S-i.RData")
# load("S.RData")
# 
# dat <- snp.select(dat,E.i,E,S.i,S) #  combining snps
# #fsnp  significant locis of E, EI,S,Si  were combined
# save(dat,file="dat(fsnp).RData")

load("dat(fsnp).RData")

source("ES.ode.R")
source("optim.R")
# H0<-H0.est(dat)
# load("H0-1.RData")
#ret<-com.DH.est1(dat,interval=c(1,75168))
load("ret.RData")#DDHR

source("ES.ode2.R")
source("optim2.R")
# H0<-H0.est(dat)
# load("H0-2.RData")
#ret<-com.DH.est1(dat,interval=c(1,75168))
load("ret2.RData")#CoCoM

