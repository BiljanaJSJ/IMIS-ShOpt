#############################################################################
#run the IMIS-ShOpt-ABC algorithm
#############################################################################

rm(list=ls())
library(optimx)
library(mvtnorm)
library(parallel)
library(numDeriv)
library(utils)
library(matrixcalc)

source('RickerABC_functions.R')
source('IMIS_ShOpt_ABC.R')


#simulate data and get the summary statistics vector for the observed data
#K=100
#logr=0.5
#logthet=1
#sigma2=0.1^2
#phi=4
#theta.true=c(logr,phi,sigma2,logthet)
#names(theta.true)=c("logr","phi","sigma2",'logthet')
#N0=N0.true=3
#timestep=50
#obs.data = theta_ricker(N0.true,theta.true,timestep)
#obs.summ.full = summaryfun(obs.data)
#plot(obs.data$y,type='l')

# # # #
#tune the tolerance vector
#M=30
#cl=makeCluster(rep('localhost', 4))
#clusterExport(cl,varlist=ls(),envir = environment())
#clusterCall(cl,function(x) {source('RickerABC_functions.R');source('IMIS_ShOpt_ABC.R');library(numDeriv)})
#theta=c(1.115932, 3.919695,   0.3,1.095073)
#clusterExport(cl,varlist=ls(),envir = environment())
#z_eps =parLapply(cl,1:1000, function(x) {sapply(1:M,function(a) {theta_ricker(N0.true,theta,timestep)$y}) })
#clusterExport(cl,varlist=ls(),envir = environment())
#FindEps1=parLapply(cl,1:length(z_eps), function(x) {do.call(rbind,lapply(1:M,function(a){abs ( applysummaries(z_eps[[x]])[,a]-obs.summ.full ) }))})
#FindEps1=do.call(rbind,FindEps1)
#quant=0.0001
#eps1=apply(FindEps1,2,quantile,probs=quant)
#eps=eps1
#save(dataEps,file='LoadDataEps_theta_RickerData_newData_poiss_neweps_0.3_trynew.RData')
#clusterExport(cl,varlist=ls(),envir = environment())


de1=get(load('LoadDataEps_theta_RickerData_newData_poiss_neweps_0.3_trynew.RData'))
obs.data=de1$obs.data
theta.true=de1$obs.data$theta
eps=de1$eps
N0.true=de1$N0.true
timestep=length(de1$obs.data$y)
obs.summ.full = summaryfun(obs.data)

optim.fun1=logposterior_synthloglike
optim.fun2=synthloglike

#run the IMIS-ShOpt-ABC
set.seed(3588648)
samples=IMIS_ShOpt_ABC(B=1000, B.re=3000, number_T=500,D=4,Q=3,obs.data=obs.data,M=50, N0.true=N0.true,timestep=timestep,dist.metric=abs,eps=eps,optim.fun1=optim.fun1, optim.fun2=optim.fun2, other=NULL)


####################################################################
#THE END!
####################################################################


