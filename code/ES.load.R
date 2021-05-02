
SE.load <- function(e.file="./data/coculture-Ecoli.csv",s.file="./data/coculture-S1.csv",
                    e.snp.file="./data/gen_e.txt",s.snp.file="./data/gen_s.txt",
                    ei.file="data/monocultrue-Ecoli.csv",si.file="./data/monocultrue-S.csv"){
  
  SP <- read.csv(s.file)
  rownames(SP) <- as.character(SP[,1])
  SP1 <- SP[,-1]
  colnames(SP1) <- NULL
  SP1L <- log(SP1)

  SPi <- read.csv(si.file)
  rownames(SPi) <- as.character(SP[,1])
  SP1i <- SPi[,-1]
  colnames(SP1i) <- NULL
  SP1Li <- log(SP1i)

  EP <- read.csv(e.file)
  rownames(EP) <- as.character(EP[,1])
  EP1 <- EP[,-1]
  colnames(EP1) <- NULL
  EP1L <- log(EP1)

  EPi <- read.csv(ei.file)
  rownames(EPi) <- as.character(EP[,1])
  EP1i <- EPi[,-1]
  colnames(EP1i) <- NULL
  EP1Li <- log(EP1i)

  
  sample_N <- dim(SP)[1]
  sample_times <- c(0,2,4,6,8,10,12,14,16,18,20,22,24,28,32,36)
  
  ssnp <- read.table(s.snp.file,header=T)
  rownames(ssnp) <- ssnp[,1]
  ssnp1 <- ssnp[,-1]
  colnames(ssnp1) <- as.character(SP[,1])
  
  esnp <- read.table(e.snp.file,header=T)
  rownames(esnp) <- esnp[,1]
  esnp1 <- esnp[,-1]
  colnames(esnp1) <- as.character(EP[,1])
  
  del.i <- c()
  for(i in 1:dim(esnp1)[1]){
    
    index <- min(table(as.character(unlist(c(esnp1[i,]))))/sample_N)
    del.i <- c(del.i,index)
  }
  
  nesnp1 <- esnp1[-which(del.i<0.1),]
  
  del.i <- c()
  for(i in 1:dim(ssnp1)[1]){
    
    index <- min(table(as.character(unlist(c(ssnp1[i,]))))/sample_N)
    del.i <- c(del.i,index)
  }
  
  nssnp1 <- ssnp1[-which(del.i<0.1),]
  
  list(sp.p=SP1L,ep.p=EP1L,sp.pi=SP1Li,ep.pi=EP1Li,
       sample_times=sample_times,sample_N=sample_N,s.snp=nssnp1,e.snp=nesnp1)
}

