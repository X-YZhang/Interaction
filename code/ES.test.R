setwd("C:/Users/ghy/Desktop/EcoInteraction/code")

source("ES.load.R")
source("ES.ode.R")
source("ode.part.R")
source("ES.covar.R")
source("ES.util.R")
library(mvtnorm)

dat <- SE.load(e.file="./data/coculture-Ecoli.csv",s.file="./data/coculture-Saureus.csv",
                    e.snp.file="./data/gen_ecoli.txt",s.snp.file="./data/gen_saureus.txt",
                    ei.file="data/monocultrue-Ecoli.csv",si.file="./data/monocultrue-Saureus.csv")

dat1 <- ES.adjust.s(dat,eefile="Structure/E.coli/E/Eput_simple.5.meanQ",
                    ssfile="Structure/S.aureus/S/Sput_simple.8.meanQ")

E.i <- ES.cov.test(pheno=dat1$ep.pi[,-1],snp=dat1$e.snp)
E <- ES.cov.test(pheno=dat1$ep.p[,-1],snp=dat1$e.snp)
S.i <- ES.cov.test(pheno=dat1$sp.pi[,-1],snp=dat1$s.snp)
S <- ES.cov.test(pheno=dat1$sp.p[,-1],snp=dat1$s.snp)

save(E.i,file="E-i.RData"); save(E,file="E.RData");
save(S.i,file="S-i.RData"); save(S,file="S.RData")

dat <- snp.select(dat1,E.i,E,S.i,S)
save(dat,file="dat(fsnp).RData")

load("dat(fsnp).RData")

ret<-com.DH.est1(dat,interval=c(1,5))

e<-read.csv("./ͼƬ/e.csv")[,2]
s<-read.csv("./ͼƬ/s.csv")[,2]
ddhr.e<-names(table(read.csv("./ͼƬ/ret.csv")[,3]))
ddhr.s<-names(table(read.csv("./ͼƬ/ret.csv")[,4]))
intersect(e,ddhr.e)
intersect(s,ddhr.s)


#a high percentage of QTLs are detected to be different from those expressed in monoculture
100/132
(368-32)/368
