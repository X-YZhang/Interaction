ES.cov.test <- function(pheno,snp){
  y <- pheno
  snps <- snp
  nsnp <- dim(snps)[1]
  time <- dim(y)[2]
  p.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  fdr.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  for(i in 1:time){
    for(j in 1:nsnp){
      tmpsnp <- snps[j,]
      miss <- which(is.na(tmpsnp))
      if(length(miss)>0){
        newsnp <- tmpsnp[-miss]
        yy <- y[-miss,i]
      }else{
        newsnp <- tmpsnp
        yy <- y[,i]
      }
      symbol <- names(table(as.character(unlist(c(newsnp)))))
      index1 <- which(newsnp==symbol[1])
      y1 <- yy[index1]
      index0 <- which(newsnp==symbol[2])
      y0 <- yy[index0]
      #snpdata <- data.frame(X=c(y1,y0),A=factor(c(rep(1,length(index1)),
      #rep(2,length(index0)))))
      #snp.aov <- aov(X ~ A, data=snpdata)
      var.t <- var.test(y0,y1)
      if(var.t$p.value>0.05)
        var.i <- TRUE
      else
        var.i <- FALSE
      p.value <- t.test(y0,y1,var.equal = var.i)$p.value#summary(snp.aov)[[1]][[1,"Pr(>F)"]] 
      p.matrix[i,j] <- p.value
    }
    fdr.matrix[i,] <- p.adjust(p.matrix[i,],method="fdr")      
  }
  
  colnames(p.matrix)<- rownames(snps)
  colnames(fdr.matrix)<- rownames(snps)
  
  logp <- -log10(p.matrix)
  thre1 <- -log10(0.05/nsnp)
  thre2 <- -log10(0.01/nsnp)
  return(list(p.value=p.matrix,fdr=fdr.matrix,logp=logp,thre1=thre1,thre2=thre2))
}

ES.co.cov.test <- function(pheno=dat$ep.p,snp=dat$fsnp,interval=c(1,10)){
  y <- pheno
  snps <- t(snp)
  sumsnp <- dim(snps)[1]
  time <- dim(y)[2]
  n1 <- interval[1]
  n2 <- interval[2]
  nsnp <- n2-n1+1
  p.matrix <- matrix(rep(0,nsnp*time),nrow=time) #15行,列表示检验的位点的数目
  fdr.matrix <- matrix(rep(0,nsnp*time),nrow=time)
  for(i in 1:time){
    for(j in n1:n2){
      tmpsnp <- snps[j,] #位点j的45个样本的基因
      miss <- which(is.na(tmpsnp))  #有缺失的样本
      if(length(miss)>0){
        newsnp <- tmpsnp[-miss]
        yy <- y[-miss,i]
      }else{
        newsnp <- tmpsnp
        yy <- y[,i]
      }
      wi <- data.frame(
        yy,sn=as.character(unlist(c(newsnp)))
      )
      #oneway_test(yy ~ sn, data = wi, distribution = approximate(nresample = 1000))
      #nn <- table(wi$sn)
      #symbol <- names(nn)
      #var.t <- bartlett.test(yy ~ sn,data=wi)  #检验方差齐性
      y.aov <- aov(yy ~ factor(sn),data=wi)
      p.value <- summary(y.aov)[[1]][[1,"Pr(>F)"]] #Pr <- summary(y.aov)[[1]][[5]][1]
      p.matrix[i,j-n1+1] <- p.value
      cat("time",i,";snp",j,"done","\n")
    }
    fdr.matrix[i,] <- p.adjust(p.matrix[i,],method="fdr")
    #cat("time",i,"done","\n")
  }
  
  colnames(p.matrix)<- rownames(snps)[n1:n2]
  colnames(fdr.matrix)<- rownames(snps)[n1:n2]
  
  logp <- -log10(p.matrix)
  thre1 <- -log10(0.05/sumsnp)
  thre2 <- -log10(0.01/sumsnp)
  return(list(p.value=p.matrix,fdr=fdr.matrix,logp=logp,thre1=thre1,thre2=thre2))
}


