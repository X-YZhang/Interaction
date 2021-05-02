
require("deSolve")
G.get_mu <- function(par, times, options=list())
{
  A <- par[1]*exp(-exp(par[2]-par[3]*times))
  return (A);
}

G.mle <- function(par,y,times){
  
  sum((y-G.get_mu(par,times))^2)
}

L.get_mu <- function(par, times, options=list())
{
  A <- par[1]/(1+exp(par[2]-par[3]*times))
  return (A);
}

L.mle <- function(par,y,times){
  
  sum((y-L.get_mu(par,times))^2)
}

R.get_mu <- function(par, times, options=list())
{
  A <- par[1]*(1-par[2]*exp(-par[3]*times))^(1/(1-par[4]))
  return (A);
}

R.mle <- function(par,y,times){
  
  sum((y-R.get_mu(par,times))^2)
}


s.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu(s.par,s.t,x0,y0))^2 )
  A
}

com.get_mu <- function(par, times, x0,y0, options=list())
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
  y <- COMP.f( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-X/k1) + r1*X*(a1/(1+X))*Y
            dY <- r2*Y*(1-Y/k2) + r2*Y*(a2/(1+Y))*X
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

com.get_mu1 <- function(par, times,x0,y0, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      k1 = par[2],
      r2 = par[3],
      k2 = par[4]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-(X/k1))
            dY <- r2*Y*(1-(Y/k2))
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


com.get_mu2 <- function(par, times,x0,y0, options=list())
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      r1 = par[1],
      a1 = par[2],
      r2 = par[3],
      a2 = par[4]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f2( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}

COMP.f2 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(a1/(1+X))*Y
            dY <- r2*Y*(a2/(1+Y))*X
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}

