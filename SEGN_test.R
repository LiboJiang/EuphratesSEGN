
####load data#########
load("phy_dat.RData")

####load package and function########
library(grplasso)
library(parallel)
source("SEGN_sup1.R")
source("SEGN_sup2.R")
source("SEGN_sup3.R")



dd <- apply(phy_dat$treat1[,-1]-phy_dat$treat0[,-1],1,scale)

#################################smooth################################

stage21 <- smooth.optim(times=1:5,para=rep(.1,6),y=t(dd),nt=seq(1,5,length=30))

#################################variable selection####################

stage32 <- varsel1(X=t(stage21$smooth.d),Y=t(stage21$dsmooth.d),tt=seq(1,5,length=30),order=6)

#################################estimation############################

ST32.odee <- optim.parallel(connect=stage32$connect,effect=t(stage21$smooth.d),
                            n.cores=2,proc=ode.optim,order=6,times=seq(1,5,length=30),nstep=29)

ST32_res <- interType(con=stage32$connect,alle=ST32.odee,sme=stage21$smooth.d)

