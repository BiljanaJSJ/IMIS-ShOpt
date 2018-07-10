########################################################################################
#This file runs both IMIS-ShOpt and IMIS-Opt with the same settings and write files with results
########################################################################################
rm(list = ls(all = TRUE))

source("Two-stage-FhN-just-c-with-prior.R")                         # 2-stage functions
source("IMIS.opt.colloc.proc-3optimizers.general-no-touchups.R")    # General IMIS 3 optimizers function 
source("fhn-model-set-up-x0proc-just-c.R")                          # likelihood etc...
source("FhN-model-set-up-as-ode-x0proc-thetalik-justc.R")		        # basic FhN functions
#source('AnimationRuns.R')

# include an ode solver package
library(deSolve)
library(mvtnorm)
library(coda)
library(numDeriv)
library(fda)
library(CollocInfer)
source("makeSSElik.R")
source("makeSSEprocFHN.R")


library(doParallel)

fhn<- make.FHN()

#make some fake data
times = seq(0,20,0.2) 
print("Note that the parameter labelled 'sd' is actually a variance.  will fix this eventually")
x0	= c(-1,1) 
names(x0) = c("V","R")
pars=c(3)
parnames =names(pars)=c("c")

y	= lsoda(x0,times,fhn$fn.ode,pars) 

y	= y[,2:3] 
data = y + matrix(rnorm(dim(y)[1]*2,0,sqrt(.05^2)),length(times),2)


lambda=1000

knots = knots = seq(0,20,1) 

norder= 3									# order of the b-spline basis
nbasis = length(knots) +norder -2           # number of basis functions
range = c(0,max(times))						# time span of the observations
bbasis = create.bspline.basis(range=range, nbasis=nbasis, norder=norder, breaks=knots)
fd.data = array(data,c(nrow(data),1,ncol(data)))
DEfd = Data2fd(y=fd.data,argvals=times,basisobj=bbasis,fdnames = list(NULL,NULL, colnames(data)))
if(!exists("DEfd",where=1)){# then the data2fd line failed - this can occur when there are unobserved variables
  data.temp=array(0,dim(fd.data))
  data.temp[,,apply(is.na(fd.data),3,sum)<(dim(fd.data)[1]-2)]    =fd.data[,,apply(is.na(fd.data),3,sum)<(dim(fd.data)[1]-2)]
  DEfd.temp = Data2fd(data.temp,times,bbasis,fdnames = list(NULL,NULL, colnames(data)))	
}				
coefs = DEfd$coefs
colnames(coefs) = colnames(data)

#

lik=make.SSElik.with.prior.and.sd()
lik$bvals = eval.basis(times,bbasis)
#lik$bvals[1,] = 0###########because I need to put all of the dlik/dc(1) into dlik/dp for X(0)
lik$more	= make.id() 
lik$more$weights = array(1,dim(data)) 
lik$more$bvals  = eval.basis(times,bbasis,0)     ######## needed only if x0 is a parameter
lik$more$dbvals = eval.basis(times,bbasis,1)    ######## needed only if x0 is a parameter
lik$x0index = c(5,6)  # index for the initial conditions
proc=make.SSEproc.FHN()
proc$bvals=lik$bvals
qpts =knots
qwts = rep(diff(range(times))/length(knots),length(knots))
qwtsL = qwts %*% t(lambda)
proc$more = fhn 
proc$more$weightsNoL = qwts 
proc$more$lambda = lambda 
proc$more$weights = qwtsL 

proc$more$qpts	= qpts 
proc$more$parnames	= 'c'
proc$more$names	= c("V","R")
proc$bvals = list(bvals = eval.basis(proc$more$qpts,bbasis,0), dbvals = eval.basis(proc$more$qpts,bbasis,1))
print("Note that the lambda value ignored since proc is supplied")

# end of FDA set up stuff
# smoothing based optimizer control options

control=list()
control$trace=0
control$iter.max=1000
control$eval.max=1e3
control$rel.tol=1e-8
control.in=control
control.in$rel.tol = 1e-12
control.in$iter.max = 1000
control.out=control
control.out$trace=1


colnames(coefs) = names(x0)



#spars = c(.2,.2,14,1,-.5,.5)          # Perturbed parameters  Note that the first two should stay fixed since we are dealing with a one dimensional 
#system.  The third parameter is the one to manipulate
#names(spars)=c("c")

B=1000
B.re=10000
number_k=150
D=30
parnames = names(pars)
datamatrix=data
plots=1
likelihood.w.fit=NULL
ncoefs=0


other = list()
other$coefs        = coefs    
other$lik          = lik         
other$proc         = proc       
proc$more$parnames = names(pars) 
other$control.in   = control.in  
other$control.out  = control.out 
other$ncoefs       = ncoefs      

optim.fun1 = ident
optim.fun2 = smoother.opt
optim.fun3 = twostage.optim.fhn.justc

cl <- makeCluster(4)
registerDoParallel(cl)
clusterCall(cl,function(x) {library(deSolve);library(CollocInfer);library(numDeriv);library(lokern)})
clusterExport(cl,varlist=list('IMIS.opt.colloc.3optimizers.general.no.touch.ups','d2negnormdp2',"make.fhn","%dopar%","foreach",'make.SSEproc.FHN',"neglogprior","prior","likelihood",'times',"dnegnormdp",'make.SSElik',"dneglogpriordpar","lokerns","ksLqudratic",'simex.fun.justc','neq','der.fhn.justc','jac.fhn.justc','d2neglogpriordpar2'))
clusterExport(cl,varlist=ls())

t1=proc.time()[1]
##run the IMIS-ShOpt algorithm
output<- IMIS.opt.colloc.3optimizers.general.no.touch.ups(B, B.re, number_k, D,parnames = c("c"),data=data, plots=0, ncoefs=ncoefs, optim.fun1, optim.fun2, optim.fun3,other)
##save the results
save(output,file='FhN_1Param_IMIS_Shopt_D30.RData')
(t=proc.time()[1]-t1)/60
proc.time()


#run the IMIS-Opt
t1=proc.time()[1]
source('IMIS.R')
output_IMIS_opt<- IMIS(B, B.re, number_k, D=3,logging=TRUE,data)
#save results
save(output_IMIS_opt,file='FhN_1Param_IMIS_Opt_D3.RData')
(t=proc.time()[1]-t1)/60
stopCluster(cl)

proc.time()






#exclude Two Stage optimizer

t1=proc.time()[1]
##run the IMIS-ShOpt algorithm
output<- IMIS.opt.colloc.3optimizers.general.no.touch.ups(B, B.re, number_k,D, runIndex=1:2, parnames = c("c"),data=data, plots=0, ncoefs=ncoefs, optim.fun1, optim.fun2, optim.fun3,other)
                                                        
##save the results
save(output,file='FhN_1Param_IMIS_Shopt_D30_1_2.RData')
(t=proc.time()[1]-t1)/60
proc.time()
