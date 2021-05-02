setwd("/home/xyzhang/micro")
require("deSolve")
require(mvtnorm)
source("double.covar.R")
source("./co.e/co.e.optim.R")

load("dat.RData")

get_con_param<-function(parm.id)
{
  for (e in commandArgs()) 
  {
    ta = strsplit(e,"=", fixed=TRUE);
    if(! is.na( ta[[1]][2])) 
    {
      temp = ta[[1]][2];
      if( ta[[1]][1] == parm.id)
        return (temp );
    }
  }
  return(NA);
}

your.no<-get_con_param("your.no")
your.no<-as.numeric(your.no)
ret <- double.est(dat,interval=c(1+(your.no-1)*1000,your.no*1000))
filename<-paste("ret-co-e",your.no,".RData",sep = "")
save(ret,file = filename)