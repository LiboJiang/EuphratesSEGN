







smooth.optim <- function(times,para,y,nt=seq(1,6,length=150)){
  
  
  allpar <- c()
  smooth.d <- c()
  dsmooth.d <- c()
  for(i in 1:dim(y)[1]){
    L <- optim(para,smL,DS1=y[i,],times=times,method="BFGS")
    allpar <- rbind(allpar,c(L$par,L$value))
    smooth.d <- rbind(smooth.d,Legendre.model(nt,L$par))
    dsmooth.d <- rbind(dsmooth.d,dLegendre.model(nt,L$par))
  }

  list(allpar=allpar,smooth.d=smooth.d,dsmooth.d=dsmooth.d)
}



smL <- function(times,para,DS1){
  
  sum((DS1-Legendre.model(t=times,mu=para))^2)
}





Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}



dLegendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1]*0 + 1*mu[2];
  if (np.order>=2)
    L <- L + 0.5 * (6 * ti)* mu[3] ;
  if (np.order>=3)
    L <- L +0.5 * (15 * ti ^ 2 - 3)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125 * (35 * 4 * ti ^ 3 - 60 * ti)* mu[5];
  if (np.order>=5)
    L <- L + 0.125 * (63 * 5 * ti ^ 4 - 210 * ti ^ 2 + 15)*mu[6];
  if (np.order>=6)
    L <- L + (1 / 16) * (231 * 6 * ti ^ 5 - 315 * 4 * ti ^ 3 + 105 * 2 *ti)* mu[7];
  if (np.order>=7)
    L <- L + (1 / 16) * (429 * 7 * ti ^ 6 - 693 * 5 * ti ^ 4 + 315 * 3 *ti ^ 2 - 35)* mu[8];
  return(L);
}









Legendre.model <-function( t, mu, tmin=NULL, tmax=NULL )
{
  u <- -1;
  v <- 1;
  if (is.null(tmin)) tmin<-min(t);
  if (is.null(tmax)) tmax<-max(t);
  ti    <- u + ((v-u)*(t-tmin))/(tmax - tmin);
  np.order <- length(mu)-1;
  L <- mu[1] + ti*mu[2];
  if (np.order>=2)
    L <- L + 0.5*(3*ti*ti-1)* mu[3] ;
  if (np.order>=3)
    L <- L + 0.5*(5*ti^3-3*ti)*mu[4] ;
  if (np.order>=4)
    L <- L + 0.125*(35*ti^4-30*ti^2+3)* mu[5];
  if (np.order>=5)
    L <- L + 0.125*(63*ti^5-70*ti^3+15*ti)*mu[6];
  if (np.order>=6)
    L <- L + (1/16)*(231*ti^6-315*ti^4+105*ti^2-5)* mu[7];
  if (np.order>=7)
    L <- L + (1/16)*(429*ti^7-693*ti^5+315*ti^3-35*ti)* mu[8];
  if (np.order>=8)
    L <- L + (1/128)*(6435*ti^8-12012*ti^6+6930*ti^4-1260*ti^2+35)* mu[9];
  if (np.order>=9)
    L <- L + (1/128)*(12155*ti^9-25740*ti^7+18018*ti^5-4620*ti^3+315*ti)* mu[10];
  if (np.order>=10)
    L <- L + (1/256)*(46189*ti^10-109395*ti^8+90090*ti^6-30030*ti^4+3465*ti^2-63)* mu[11];
  if (np.order>=11)
  {
    for(r in 11:(np.order))
    {
      kk <- ifelse(r%%2==0, r/2, (r-1)/2);
      for (k in c(0:kk) )
      {
        L <- L + (-1)^k*factorial(2*r-2*k)/factorial(k)/factorial(r-k)/factorial(r-2*k)/(2^r)*ti^(r-2*k)*mu[r+1];
      }
    }
  }
  return(L);
}