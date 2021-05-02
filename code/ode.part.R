com.get_mu_mono <- function(par, times, x0, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2]);
  }
  
  state0 <- c(X=x0);
  y <- COMP.f.mono( par0, state0, times );
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2] ) );
}

COMP.f.mono <-function( parameters, state, times ){
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


com.get_mu_co <- function(par, times,x0,y0,options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2],
      a1 = par[3],
      r2 = par[4],
      k2 = par[5],
      a2 = par[6]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f.co( par0, state0, times);
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f.co <-function( parameters, state, times){
  
  Lorenz<-function(t, state, parameters,inter) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-X/k1) + r1*X*(a1/(1+X))*Y
            dY <- r2*Y*(1-Y/k2) + r2*Y*(a2/(1+Y))*X
            list(c(dX,dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

