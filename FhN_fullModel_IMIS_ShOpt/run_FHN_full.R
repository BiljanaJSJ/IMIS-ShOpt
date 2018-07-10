

rm(list = ls(all = TRUE))


source("makeSSElik.R")
source("IMIS.R")
source("IMIS.opt.colloc.x0proc.R")
source("IMIS.opt.colloc.proc-3optimizers.general.R")
source("Two-stage-FhN.R")# also has the full parameter functions
source("fhn-model-set-up-x0proc-Normal(c)-2SDs.R")
source("FhN-model-set-up-as-ode-x0proc-thetalik-2SD.R")
source("makeSSElik.R")
source("makeSSE.R")

# include an ode solver package
library(deSolve)
library(mvtnorm)
library(coda)
library(deSolve)
library(MASS)
library(CollocInfer)



#library(snow)

#set up the FHN functions
fhn<- make.FHN()
make.fhn=make.FHN
#make some fake data
times = seq(0,20,0.2) 
pars=c(.2,.2,3,.05,.05,-1,1)
parnames =names(pars)=c("a","b","c","sdV","sdR","V0","R0")
print("Note that the parameter labelled 'sd' is actually a variance.  will fix this eventually")
x0	= c(-1,1) 
names(x0) = c("V","R")
y	= lsoda(x0,times,fhn$fn.ode,pars) 
#data = y	= y[,2:3] 
#data[,1] = y[,1] + sqrt(pars["sdV"])*matrix(rnorm(dim(y)[1]),length(times),1)
#data[,2] = y[,2] + sqrt(pars["sdR"])*matrix(rnorm(dim(y)[1]),length(times),1)

data = y[,2:3] + matrix(rnorm(dim(y)[1]*2,0,sqrt(.05^2)),length(times),2)


#plot(times,data[,1])
#lines(times,data[,2])

#set up the parameter estimation from using colloInfer
lambda=1000

#create the basis functions using B-splines
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
#obtain initial coefficients
coefs = DEfd$coefs
colnames(coefs) = colnames(data)

fhn<- make.FHN()
#set up the lik object
lik=make.SSElik.with.prior.and.sd()
lik$bvals = eval.basis(times,bbasis)
###lik$bvals[1,] = 0###########because I need to put all of the dlik/dc(1) into dlik/dp for X(0)
lik$more	= make.id() 
lik$more$weights = array(1,dim(data)) 
lik$more$bvals  = eval.basis(times,bbasis,0)     ######## needed only if x0 is a parameter
lik$more$dbvals = eval.basis(times,bbasis,1)    ######## needed only if x0 is a parameter
lik$x0index = c(6,7)  # index for the initial conditions
proc=make.SSEproc.FHN()
proc$bvals=lik$bvals
qpts =knots
qwts = rep(diff(range(times))/length(knots),length(knots))
qwts = qwts %*% t(lambda)


#set up the proc object
proc$more = fhn 
proc$more$weights = qwts 
proc$more$qpts	= qpts 
proc$more$parnames	= names(pars[1:5]) 
proc$more$names	= c("V","R")
proc$bvals = list(bvals=eval.basis(proc$more$qpts,bbasis,0) ,dbvals = eval.basis(proc$more$qpts,bbasis,1))
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

#IMIS-ShOpt input parameters
B=1000 
B.re=10000
number_k=150
D=30
parnames = names(pars)
datamatrix=data
plots=3
likelihood.w.fit=NULL
ncoefs=2
twostage=twostage.optim.fhn

#output<- IMIS.opt.colloc.3optimizers(B, B.re, number_k, D,parnames ,lik=lik,proc=proc,coefs=coefs,datamatrix=data,plots=3,ncoefs=ncoefs, twostage=twostage,  likelihood.w.fit=likelihood.w.fit)	

#set up the other input parameter
other = list()
other$coefs        = coefs    
other$lik          = lik         
other$proc         = proc       
proc$more$parnames = names(pars) 
other$control.in   = control.in  
other$control.out  = control.out 
other$ncoefs       = ncoefs      

#set up parallel computing
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)


clusterCall(cl,function(x) {library(mvtnorm)})
clusterCall(cl,function(x) {library(coda)})
clusterCall(cl,function(x) {library(MASS)})



clusterCall(cl,function(x) {library(deSolve)})
clusterCall(cl,function(x) {library(CollocInfer)})
clusterCall(cl,function(x) {library(lokern)})


clusterExport(cl,varlist=ls())

#clusterExport(cl,varlist=list("IMIS.opt.colloc.3optimizers.general","make.fhn","%dopar%","foreach",'make.SSEproc.FHN',"neglogprior","neglogpriorx0","prior","priorx0","likelihood",'times',"dnegnormdp",'make.SSElik',"dneglogpriordpar","lokerns","ksLqudratic",'simex.fun.justc','simex.fun','sample.prior','ginv','dinvgamma','dnegIgammadp',"twostage.optim.fhn",""neq""))

#run the IMIS-ShOpt for the FhN with full parameter set to be estimated
output<- IMIS.opt.colloc.3optimizers.general(B, B.re, number_k, D,parnames ,lik=lik,proc=proc,coefs=coefs,data=data,plots=0,ncoefs=ncoefs, optim.fun1=ident,optim.fun2=smoother.opt,optim.fun3=two.stager,likelihood.w.fit=likelihood.w.fit,other)	

#save the results
save(output,file='IMIS_shopt_full_fhn_D30.RData')
stopCluster(cl)

save.image("FhN-IMIS-3-optimizers-no-touch-ups-with-prior.Rdata")