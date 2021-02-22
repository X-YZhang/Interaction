ES.adjust.s <- function(dat,eefile="Structure/E.coli/E/Eput_simple.5.meanQ",
                        ssfile="Structure/S.aureus/S/Sput_simple.8.meanQ"){
  
  times <- dat$sample_times
  nt <- length(times)
  e1 <- read.table(eefile)
  s1 <- read.table(ssfile)
  mape1 <- unlist(apply(e1, 1, function(x) which(x == max(x, 
                                                          na.rm = TRUE))[1]))
  maps1 <- unlist(apply(s1, 1, function(x) which(x == max(x, 
                                                          na.rm = TRUE))[1]))
  npe <-c()
  eresi <- c()
  estr <- c()
  for(i in 1:nt){
    data1 <- data.frame(pheno=dat$ep.p[,i],pop=mape1)
    lm.e <- lm(pheno ~ pop,data=data1)
    e.tmp <- dat$ep.p[,i]-lm.e$coefficients[2]*mape1
    resi.tmp <- residuals(lm.e)
    npe <- cbind(npe,e.tmp)
    eresi <- cbind(eresi,resi.tmp)
    estr <- cbind(estr,mape1)
  }
  
  npei <-c()
  eresii <- c()
  estri <- c()
  for(i in 1:nt){
    data1 <- data.frame(pheno=dat$ep.pi[,i],pop=mape1)
    lm.ei <- lm(pheno ~ pop,data=data1)
    ei.tmp <- dat$ep.pi[,i]-lm.ei$coefficients[2]*mape1
    resii.tmp <- residuals(lm.ei)
    npei <- cbind(npei,ei.tmp)
    eresii <- cbind(eresii,resii.tmp)
    estri <- cbind(estri,mape1)
  }
  
  nps <-c()
  sresi <- c()
  sstr <- c()
  for(i in 1:nt){
    data1 <- data.frame(pheno=dat$sp.p[,i],pop=maps1)
    lm.s <- lm(pheno ~ pop,data=data1)
    s.tmp <- dat$sp.p[,i]-lm.s$coefficients[2]*maps1
    res.tmp <- residuals(lm.s)
    nps <- cbind(nps,s.tmp)
    sresi <- cbind(sresi,res.tmp)
    sstr <- cbind(sstr,maps1)
  }
  
  npsi <-c()
  sresii <- c()
  sstri <- c()
  for(i in 1:nt){
    data1 <- data.frame(pheno=dat$sp.pi[,i],pop=maps1)
    lm.si <- lm(pheno ~ pop,data=data1)
    si.tmp <- dat$sp.pi[,i]-lm.si$coefficients[2]*maps1
    resii.tmp <- residuals(lm.si)
    npsi <- cbind(npsi,si.tmp)
    sresii <- cbind(sresii,resii.tmp)
    sstri <- cbind(sstri,maps1)
  }
  dat1 <- dat
  rownames(npe) <- rownames(dat$ep.p);rownames(npei) <- rownames(dat$ep.pi);
  rownames(nps) <- rownames(dat$sp.p);rownames(npsi) <- rownames(dat$sp.pi);
  dat1$ep.p <- npe;dat1$sp.p <- nps;dat1$ep.pi <- npei;dat1$sp.pi <- npsi;
  dat1$sstr <- sstr;dat1$estr <- estr;
  dat1$eresi <- eresi;dat1$eresii <- eresii;dat1$sresi <- sresi;dat1$sresii <- eresii;
  return(dat1)
}


snp.select <- function(dat,E.i,E,S.i,S){
  
  
  E.tmp <- c()
  S.tmp <- c()
  for(i in 1:length(dat$sample_times)){
    e1 <- as.numeric(which(E.i$p.value[i,]<0.01))
    e2 <- as.numeric(which(E$p.value[i,]<0.01))
    E.tmp <- c(E.tmp,e1,e2)
    
    s1 <- as.numeric(which(S.i$p.value[i,]<0.01))
    s2 <- as.numeric(which(S$p.value[i,]<0.01))
    S.tmp <- c(S.tmp,s1,s2)
  }
  E.s <- sort(unique(c(E.tmp)))
  S.s <- sort(unique(c(S.tmp)))
  new.e <- dat$e.snp[E.s,]
  new.s <- dat$s.snp[S.s,]
  
  new.en <- rownames(new.e)
  new.sn <- rownames(new.s)
  snpnc <- c()
  snpc <- c()
  for(j in 1:dim(new.e)[1]){
    for(k in 1:dim(new.s)[1]){
      tmp <- paste(as.character(unlist(c(new.e[j,]))),as.character(unlist(c(new.s[k,]))),sep="")
      snpc <- cbind(snpc,tmp)
      tmp1 <- c(as.numeric(new.en[j]),as.numeric(new.sn[k]))
      snpnc <- cbind(snpnc,tmp1)
      cat("e",j,",s",k,"\n")
    }
  }
  
  del.i <- c()
  for(ii in 1:dim(snpc)[2]){
    tmpm <- min(as.numeric(table(snpc[,ii])))
    if(tmpm <4){
      del.i <- c(del.i,ii)
    }
  }
  fsnp <- snpc[,-del.i]
  fsnpn <- snpnc[,-del.i]
  colnames(fsnp) <- 1:dim(fsnp)[2]
  colnames(fsnpn) <- 1:dim(fsnp)[2]
  
  dat$fsnp <- fsnp
  dat$fsnpn <- fsnpn
  return(dat)
}



