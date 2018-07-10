###################################################################################
# SIR_IMIS_Opt_functions.R
# contains all the functions needed for the algorithm to run on the SIR example
###################################################################################


###################################################################################
#install packages if they are not installed
###################################################################################
checkpackages=function(package){
  if (!package %in% installed.packages())
    install.packages(package)
}

# checkpackages("gtools")
# checkpackages("MCMCpack")
# checkpackages("mvtnorm")
# checkpackages("truncnorm")
# checkpackages("optimx")
# checkpackages("coda")
# checkpackages("deSolve")
# checkpackages("CollocInfer")
# checkpackages("lokern")
# checkpackages("matrixcalc")
# This library lets us sample from a dirichlet
library(gtools)
library(MCMCpack)
library(mvtnorm)
library(truncnorm)
library(optimx)
library(coda)
library(deSolve)
# library(CollocInfer)
# library(lokern)
 library(matrixcalc)
library(numDeriv)
###################################################################################
#logprior_alphaBeta  evaluate log prior (beta,alpha \ I0)
#input:        pars      - a vector of (beta,alpha)
#              I         - numerical solution for I
#              priorpars - prior parameters for alpha, beta and I0
#output:       evaluated log prior     
###################################################################################
logprior_alphabeta=function(pars,priorpars=c(1,1,1,1,N,5/N)){
  #n      = length(I)
  term1  = dgamma(pars[1],priorpars[1],priorpars[2],log=T)
  term2  = dgamma(pars[2],priorpars[3],priorpars[4],log=T)
  #term3  = dbinom(pars[4],priorpars[5],priorpars[6],log=T)
  #out    = term1+term2+term3
  out    = term1+term2
  return(out)
}
###################################################################################
#SIRlogpost_alphaBeta    evaluate joint log posterior  P(alpha,beta/ I0  Y)
#input:        pars      - a vector of (beta,alpha)
#              y         - data
#              R         - numerical solution for R
#              I         - numerical solution for I
#              tau       - inverse temperature 
#              priorpars - prior parameters for alpha, beta and I0
#              N         - population size
#              discrete  - when not NULL used in the optimization routine
#                          to find the solution at each function evaluation
#output:       evaluated log posterior
###################################################################################
SIRlogpost_alphaBeta = function(pars,y,R=NULL,I,I0=NULL,priorpars=c(1,1,N,5/N),N,times) {
  
  if (!is.null(I0)){
    pars=c(pars,N-I0,I0)
  }
  
  
  if (is.null(R)){
    rksol=Rhat(pars,dat=y,N,times)
    R=rksol[,'R']
    I=rksol[,'I']
    
  }
  
  
  loglik    = SIRloglik(pars,y,N=N,times=times)[[1]]
  logpriors = logprior_alphabeta(pars)
  out       = loglik+logpriors
  
  return(out)
}


####################################################################################################
# A wrapper function used for optimization called from the STstep_tau
#         optimizes alpha and beta for 8 different fixed values for I0
#         at proposed and current value of tau
#         this function is called 8 times in parallel
# input:  discrete     - when not NULL used in the optimization routine
#                        to find the solution at each function evaluation of SIRlogpost 
#         theta        - a vector of current values (beta,alpha, S0,I0)
#         y            - data
#         tau_prop     - proposed value of tau
#         te           - current value of tau
#         priorpars    - prior parameters for alpha, beta and I0
#         N            - population size
# output: max_pars_p   - maximized (beta, alpha) at proposed tau
#         eval_fun_p   - maximum value of joint posterior (beta,alpha, I0) at proposed tau
#         max_pars_i   - maximized (beta, alpha) at current tau
#         eval_fun_i   - maximum value of joint posterior (beta,alpha, I0) at current tau
####################################################################################################
wrapp_fun=function(theta,y,priorpars,N,times){
  

 
    # max_pars_p_optim=optim(theta[1:2],function(x) (-1)*SIRlogpost(pars=c(x[1],x[2],theta[3],theta[4]),y=y,
    #                                                               I=NULL,priorpars=priorpars,N=N,times=times),method='Nelder-Mead',control=list(maxit=1000))
    # 
    # hess = hessian(function(x) (-1)*SIRlogpost(pars=c(x[1],x[2],theta[3],theta[4]),y=y,
    #                                   I=NULL,priorpars=priorpars,N=N,times=times),c(max_pars_p_optim$par,theta[4]))
    max_pars_p_optim=optim(theta[1:2],function(x) (-1)*SIRlogpost_alphaBeta(pars=c(x[1],x[2]),y=y,I0=theta[4],
                                                                I=NULL,priorpars=priorpars,N=N,times=times),method='Nelder-Mead',control=list(maxit=1000))
  
  
    hess = hessian(function(x) (-1)*SIRlogpost_alphaBeta(pars=c(x[1],x[2]),y=y,I0=theta[4],N=N,times=times),max_pars_p_optim$par)

    max_pars_p       = max_pars_p_optim$par
    eval_fun_p       =(-1)*max_pars_p_optim$value
  
  return(list(max_pars_p=max_pars_p,eval_fun_p=eval_fun_p,hess=hess))
}



OptimizeSIRPost=function(theta,y,priorpars,N,times,cov_prior_global){
 
   
max_pars_p =matrix(NA,10,4)
eval_fun_p = rep(NA,10)

getlist=lapply(1:10,function(x) {wrapp_fun(theta=c(theta[1:2],N-x,x),y=y,priorpars=priorpars,N=N,times=times)})

for (i in (1:10)){
  max_pars_p[i,] =c(getlist[[i]]$max_pars_p,N-i,i)
  eval_fun_p[i]  =getlist[[i]]$eval_fun_p
}

# to find the maximum value of I0, just find the I0 that corresponds to the maximum of the 
# evaluated joint posterior SIRlogpost
#Do this for proposed and for current tau separatelly
ind=which(eval_fun_p==max(eval_fun_p))
if (length(ind)>1) {ind= sample(ind,1)}
# max_pars_prop keeps maximized (beta,alpha,S0,I0) at proposed tau
max_pars_prop=max_pars_p[ind,]


hess=getlist[[ind]]$hess
print(hess)

if (!is.positive.definite(hess)){warning('Hessian matrix is not positive definite!!') ;hess=hess+diag(1/diag(cov_prior_global))}

return(list(parameters=max_pars_prop,hessian=hess))
}

###################################################################################
#SIRloglik  evaluate tempered likelihood, which is binomial
#also untempered likelihood is evaluated which can be used to estimate the 
#marginal likelihood
#input:          y            - data
#                R            - number of removed obtained from the numerical 
#                               solver of the ODE
#                I            - number of infected obtained from the numerical 
#                               solver of the ODE
#                tau          - inverse temperature
#                N            - population size           
#output:         a vector with evaluated tempered and untempered likelihood 
###################################################################################
SIRloglik = function(thetas,y,N=261,times) {
  if (is.nan(thetas[1])){
    llik=-999999
  }else{
    rksol=Rhat(thetas,y,N,times)
    R=rksol[,4]
    I=rksol[,3]
    n       = length(I)
    
    ys      = dbinom(y[,3],N,R/N,log=T)
    xsn     = dbinom(0,N,I[n]/N,log=T)
    xsn1    = dbinom(1,N,I[n-1]/N,log=T)
    llik    = sum(ys[which(!is.nan(ys))])+ifelse(is.nan(xsn),0,xsn )+ifelse(is.nan(xsn1),0,xsn1 )
  }
  return(c(llik=llik))
}
###################################################################################


###################################################################################
#SIRlogpost    evaluate joint log posterior  P(alpha,beta,I0 / Y, tau)
#input:        pars      - a vector of (beta,alpha,S0,I0)
#              y         - data
#              R         - numerical solution for R
#              I         - numerical solution for I
#              tau       - inverse temperature 
#              priorpars - prior parameters for alpha, beta and I0
#              N         - population size
#              discrete  - when not NULL used in the optimization routine
#                          to find the solution at each function evaluation
#output:       evaluated log posterior
###################################################################################
SIRlogpost = function(pars,y,R=NULL,I,I0=NULL,priorpars=c(1,1,N,5/N),N,times) {
  
  if (!is.null(I0)){
    pars=c(pars,I0,N-I0)
  }
  
  if (is.null(R)){
  rksol=Rhat(pars,dat=y,N,times)
  R=rksol[,'R']
  I=rksol[,'I']
 
  }
  loglik    = SIRloglik(pars,y=y,N=261,times)[[1]]
  logpriors = logprior(pars,priorpars=priorpars)
  out       = loglik+logpriors
  
  return(out)
}



###################################################################################
#logprior  evaluate log prior
#input:        pars      - a vector of (beta,alpha, S0,I0)
#              I         - numerical solution for I
#              priorpars - prior parameters for alpha, beta and I0
#output:       evaluated log prior     
###################################################################################
logprior=function(pars,priorpars=c(1,1,1,1,N,5/N)){
  #n      = length(I)
  term1  = dgamma(pars[1],priorpars[1],priorpars[2],log=T)
  term2  = dgamma(pars[2],priorpars[3],priorpars[4],log=T)
  term3  = dbinom(pars[4],priorpars[5],priorpars[6],log=T)
  out    = term1+term2+term3
  return(out)
}

###################################################################################
#sampleprior  - draw initial samples from the prior
#input:        N0                - size of the initial sample
#              priorpars         - parameters of the prior
#output:       samples from the prior of dimension N0x6
###################################################################################


sampleprior=function(N0,priorpars=c(1,1000,1,10,N,5/N)){
  
  term     = matrix(NA,N0,4) 
  term[,1] = rgamma(N0,priorpars[1],priorpars[2])
  term[,2] = rgamma(N0,priorpars[3],priorpars[4])
  term[,4] = rbinom(N0,priorpars[5],priorpars[6])
  term[,3] = N-term[,4]
  
  
  return(term)
}




###################################################################################
#setup the SIR-ODE model, SIR first derivative
###################################################################################


make.SIR<-function(){
  SIR.fn <- function(times, y, p, more) {
    r = y
    r[, "S"] = -p["beta"] * y[, "S"]* y[, "I"]
    r[, "I"] =  p["beta"] * y[, "S"]* y[, "I"]-p["alpha"]*y[,"I"]
    r[, "R"] =                                 p["alpha"]*y[,"I"]        
    return(r)
  }
  SIR.fn.ode <- function(times, y, p) {
    r = y
    
    dimnames(r) = dimnames(y)
    
    r[ "S"] = -p["beta"] * y[ "S"]* y[ "I"]        
    r[ "I"] =  p["beta"] * y[ "S"]* y[ "I"]-p["alpha"]*y["I"]
    r[ "R"] =                               p["alpha"]*y["I"] 
    return(list(r))
  }
  SIR.dfdx <- function(times, y, p, more) {
    #print("SIR.dfdx")
    r = array(0, c(dim(y), dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y))
    r[, "S", "S"] = -p["beta"] * y[, "I"]
    r[, "S", "I"] = -p["beta"] * y[, "S"]
    r[, "I","S"]  =  p["beta"] * y[, "I"]
    r[, "I","I"]  =  p["beta"] * y[, "S"]-p["alpha"]
    r[, "R","I"]  =                       p["alpha"]
    return(r)
  }
  SIR.dfdp <- function(times, y, p, more) {
    #print("SIR.dfdp")
    r = array(0, c(dim(y), length(p)))
    dimnames(r) = list(NULL, colnames(y), names(p))
    
    r[, "S","beta"]  = - y[, "S"]* y[, "I"]
    r[, "I","beta"]  =   y[, "S"]* y[, "I"]
    r[, "I","alpha"] = - y[,"I"]
    r[, "R","alpha"] =   y[,"I"] 
    return(r)
  }
  SIR.d2fdx2 <- function(times, y, p, more) {
    r = array(0, c(dim(y), dim(y)[2], dim(y)[2]))
    dimnames(r) = list(NULL, colnames(y), colnames(y), colnames(y))
    r[, "S", "S","I"] = -p["beta"]
    r[, "S", "I","S"] = -p["beta"]
    r[, "I","S","I"]  =  p["beta"]
    r[, "I","I","S"]  =  p["beta"]
    return(r)
  }
  SIR.d2fdxdp <- function(times, y, p, more) {
    r = array(0, c(dim(y), dim(y)[2], length(p)))
    dimnames(r) = list(NULL, colnames(y), colnames(y), names(p))
    r[, "S", "S","beta"]  = -y[, "I"]
    r[, "S", "I","beta"]  = -y[, "S"]
    r[, "I", "S","beta"]  =  y[, "I"]
    r[, "I", "I","beta"]  =  y[, "S"]
    r[, "I", "I","alpha"] =  -1 
    r[, "R", "I","alpha"] =   1
    return(r)
  }
  return(list(fn = SIR.fn, fn.ode = SIR.fn.ode, dfdx = SIR.dfdx, 
              dfdp = SIR.dfdp, d2fdx2 = SIR.d2fdx2, d2fdxdp = SIR.d2fdxdp))
}

#############################################################
#SSE    - calculate the SSE
#############################################################



SSE <- function(pars,times, dat,N=261, rez=NULL){

   
   paradd=c(pars[1:2],N-pars[3],pars[3])
   fits=Rhat(paradd,dat=dat,N,times,rez=rez)
    SSE   = sum((fits[,ncol(fits)]-dat[,3])^2 + (fits[,ncol(fits)-1][c(length(times)-1,length(times))]-dat[c(nrow(dat)-1,nrow(dat)),2])^2) 
   
  return(SSE) 
} 


################################################################################
#NLSRezStep:        solve the SIR-ODE using the NLS, by smoothing the residuals
#input:              IniPar    - initial parameters
#                    times     - times
#                    dat       - data 
#                    cov_prior - covraince of the prior
#output:     optimized parameters and Hessian
################################################################################



NLSRezStep =function(IniPar=NULL,times,dat,cov_prior=NULL){
 
  if (is.null(IniPar)){
    (IniPar= c(b=abs(rnorm(1,0,.005)),a=abs(rnorm(1,0,.05)), I0=abs(rnorm(1,5,5*(1-5/N))) ))
  }else{
    IniPar=IniPar[c(1,2,4)]
  }
  (OptimOut         = optim(par=IniPar[1:2],I0=5,fn=minimize.f,dat=dat,method="BFGS",N=N,times=times,hessian=T,rez=1))
  optpars=OptimOut$par
  
  
  calcHess=OptimOut$hessian
  if (!is.positive.definite(calcHess)){
    warning('Hessian matrix is not positive definite!')
  }
  
  
  #get the hessian of the maxpars using hessian?
  #calculate the standard errors of the parameters estimates
  SE=(1/sqrt(diag(calcHess)))/sqrt(nrow(dat))
  
  #calculate confidence intervals
  
  #calculate confidence intervals
  #assuming assymptotical normality
  ci_lower=as.numeric(OptimOut$par)-qnorm(0.975)*SE
  ci_upper=as.numeric(OptimOut$par)+qnorm(0.975)*SE
  
  (confInt=c(ci_lower,ci_upper))
  
  
  return(list(parameters=optpars, hessian=calcHess,SE=SE,confInt=confInt,dat=dat,IniPar=IniPar))
  
}


################################################################################
#NLSStep:            solve the SIR-ODE using the NLS
#input:              IniPar   - initial parameters
#                    times    - times
#                    dat      - data 
#output:     optimized parameters and Hessian
################################################################################



NLSStep =function(IniPar=NULL,times,dat){
  
  #obtain the true ODE solution using lsoda; 0.2,0.2,3,-1,1 are the parameters of
  #the true solution
  
  
  #update the initial state for the optimizer with the
  #initial values of the observations
   

  #minimize the sum of squared error
  #where error is obtained as a difference between the data and
  #fitted curve (the solution of ODE at supplied initial parameters)

#   discrete version  
   library(numDeriv)
   parameters=matrix(NA,10,4)
   calcHess=list()
    if (is.null(IniPar)){
      #(IniPar= c(b=abs(rnorm(1,0,.005)),a=abs(rnorm(1,0,.05)), I0=abs(rnorm(1,5,5*(1-5/N))) ))
      (IniPar= c(b=abs(rnorm(1,0,.005)),a=abs(rnorm(1,0,.05)), I0=abs(rnorm(1,5,5*(1-5/N))) ))
    }else{
      IniPar=IniPar[c(1,2,4)]
    }
    for (i in (1:10)){
    
    (OptimOut         = optim(par=IniPar[1:2],I0=i,fn=minimize.f,dat=dat,method="Nelder-Mead",N=N,times=times))
    (parameters[i,]   = c(OptimOut$par[1],OptimOut$par[2],i,OptimOut$value))
     calcHess[[i]]    = hessian(func=SSE,x=parameters[i,1:3],times=times,dat=dat)
     if (!is.positive.definite(calcHess[[i]])){
       warning(paste('Hessian matrix is not positive definite!',i,sep=''))
     }
    } 
 
   
   ind=which(parameters[,4]==min(parameters[,4]))
   optpars=parameters[ind,]
   hess=calcHess[[ind]]
   for (i in (1:10)){
   parameters[i,]=optpars
   calcHess[[i]] =hess
   }
  
   
  
   return(list(parameters=parameters, hessian=calcHess,IniPar=IniPar))
  
}

  
  
#########################################################################################
#symbolicDeriv.function: symbolic derivatives of a given function
#input:   vec            - a vector of parameters of the function
#         deriv.level    - level of derivative
#         a              - alpha
#         b              - beta
#         I0             - I0
#         dat            - the data
#         times          - times
#         N              - sample size
#output:  evaluated n-th derivative
#########################################################################################

symbolicDeriv.function = function(vec,deriv.level=0,a=NULL,b=NULL,I0=NULL,
                                  dat,times,N,kernel.Deriv=NULL,rez=NULL){
  
  
  #First if the two parametres are NULL that input vector is of length 2, so we assign all 2 values from
  #the input vector to local variables for the two parametres
 if (is.null(a) && is.null(b) && is.null(I0) ){
 
    
    b =  vec[1]
    a =  vec[2]
   I0 =  vec[3] 
  }else if (is.null(a) && !is.null(b) && !is.null(I0)  ){
    #The next case is used for the grids of a. So we supply function with values for tau from the grid,
    #and the input vector should contain values for the rest of the parametars (a)
    
    a = (vec[2])
    
}else if(!is.null(a) && is.null(b) && !is.null(I0)  ){
    #Similar as previous one, here the grid is for b
    
  b = (vec[1])
    
}else if(!is.null(a) && !is.null(b) && is.null(I0)  ){
    #Similar as previous one, here the grid is for I0
    
    I0 = vec[3]

 
}else if(is.null(a) && is.null(b) && !is.null(I0)  ){
  #Similar as previous one, here the grid is for I0
  
  b =  vec[1]
  a =  vec[2]
  
}
  if (a>0 && b>0 ){
    f=expression(  SIRlogpost(pars=c(b,a,I0),y=dat,R=NULL,I,priorpars=c(1,1000,1,10,N,5/N),N=261,times=times) )
  }else{f=99999999}
  #print("NLS")
# now evaluate f,df and hessian
 
  if(deriv.level==0){# just evaluate the function
    
    output = eval(f,list(b=b, a=a,I0=I0))
    
  }
  return(output)
}
###########################################################################
#wrapper for the symbolicDeriv.function
###########################################################################
minimize.f = function(vec,deriv.level=0,a=NULL,b=NULL,I0=NULL,dat,times,N,kernel.Deriv=NULL,rez=NULL){
  symbolicDeriv.function(vec=vec,deriv.level=deriv.level,a=a,b=b,I0=I0,dat=dat,times=times,N=N,kernel.Deriv=kernel.Deriv,rez=rez)
  
}




###########################################################################
#GProfilingStep          estimate the parameters of the SIR-ODE system using
#                        the Generalized Profiling (GP)
#input:                  IniPar       - initial parameters for GP
#                        dat          - data
#                        SIR.lambda   - lambda, the smoothing parameter
#                        times        - times
#                        proc         - proc 
#output:                parameter estimates from GP
###########################################################################



GProfilingStep =function(IniPar=NULL,dat,SIR.lambda = 1,
                         times=times,proc=NULL){
  
  logdata=log(dat)
  logdata[,2]=NA
  logdata[which(logdata==-Inf)]=0
  SIR = make.SIR()
  SIR.fn=SIR$fn
  SIR.range = range(times)
  
  SIR.knots = seq(SIR.range[1],SIR.range[2],length.out = 401)
  #range(FHN.knots) 
  #FHN.norder = 3
  SIR.norder =3
  #FHN.nbasis=FHN.norder+2
  SIR.nbasis = length(SIR.knots) + SIR.norder - 2
  
  #so this is a basis system Ï
  SIR.basis = create.bspline.basis(SIR.range, SIR.nbasis,
                                   SIR.norder, SIR.knots)
     
  #evaluate b splines functions
  #eval.basis.fns=bsplineS(x=times, breaks=SIR.knots, norder=SIR.norder, nderiv=0)
  #eval.basis.fns=eval.basis(times,SIR.basis)
  #SIR.names=names(dat[,3])
  SIR.names=c('S','I','R')
  #obtain the initial values of the coefficients C by smoothing the raw noisy data 
  #using the B-spline as a 
  #basis function (which is obtained from create.bspline.basis)
  SIR.fdnames = list(NULL, NULL, SIR.names)
  
  
  #  Since we need to use richer number of basis functions we can not use smooth.basis
  SIR.knots=sort(c(SIR.knots,seq(SIR.range[1],SIR.range[2],0.01)[500]))
  
  if (is.null(IniPar)){
    IniPar= c(beta=abs(rnorm(1,0,.005)),alpha=abs(rnorm(1,0,.05)))     
  }else{
    IniPar= c(beta=IniPar[1],alpha=IniPar[2]) 
  }
  #this is initial coefficients?
  
  Defd=smooth.basis(times,logdata[,3],fdPar(SIR.basis,1,2))
  
  SIR.coefs=cbind(rep(0,length(Defd$fd$coefs)),rep(0,length(Defd$fd$coefs)),Defd$fd$coefs)
  #SIR.coefs=Defd$fd$coefs
  colnames(SIR.coefs) = SIR.names
  
  plot(logdata[,3],main='Initial basis coef')
  lines(SIR.coefs[,3],col='red')
  DEfd = fd(SIR.coefs,SIR.basis)  
  
  # extract the coeficients from the fd object
  
  #having coefficients and evaluated the basis functions
  #we can calculate the initial state paramateres (V0,R0)
  
  #now use generalized profiling to obtain 
  #the estimates of the parameters
  #check the smoothing how well it approximated

  #first set up constrol.out and control.in parameters for the 
  #Profile.LS function

  control.in = list()
  control.in$rel.tol = 1e-12
  control.in$iter.max = 100
  
  control.out = list()
  control.out$trace = 2
  control.out$iter.max=100
  
  
  # run the generalized profile step
  # set up the initial parameters 
  
 
  profile.obj  =  LS.setup(pars=IniPar,coefs=SIR.coefs,fn=make.SIR(),data=dat,
                           basisvals=SIR.basis,names=SIR.names,lambda=c(1000,100,10),
                           times=times,posproc = 1,poslik=0)
  
  SIR.lik      =  profile.obj$lik
  SIR.proc     =  profile.obj$proc
  
  #make.SSElik()$fn change this function to compare data and fdevals only at observed components
  #and also the initial coeficients should be obtained as numerical solution of the ODE
  #at the provided initial parameters
  
  fres         =  FitMatchOpt(coefs=SIR.coefs,which=c(1,2),pars= IniPar,proc=SIR.proc)
  
  SIR.coefs=fres$coefs
  
  DEfd1 = fd(fres$coefs,SIR.basis)
  plot(DEfd1,main='Coef from FitMatchOpt')
  points(times,logdata[,3])
  
 
  
  innOpt       =  inneropt(data=logdata,times=times,pars=IniPar,coefs=SIR.coefs,lik=SIR.lik,proc=SIR.proc,in.meth='nlminb',control.in=control.in)
  
  DEfd2 = fd(innOpt$coefs,SIR.basis)
  plot(DEfd2,main='Coef from inneropt')
  points(times,logdata[,3])
  save(innOpt,file='InnOptRez.RData')
  
  t=proc.time()
  resultList   = outeropt(data=dat,times=times,pars=IniPar,coefs=innOpt$coefs,lik=SIR.lik,proc=SIR.proc,control.in = control.in,control.out = control.out)
  t1=proc.time() - t
  
  save(resultList,file='resultListOuterOpt10.RData')
  DEfd3 = fd(resultList$coefs,SIR.basis)
  plot(DEfd3)
  points(times,dat[,3])
  #IniPar
  plotfit.fd(logdata[,3],times,DEfd3[3],ylab="Fit to Data")
  
  
  
  #parameters estimates
  SIR.pars    = resultList$pars
  names(SIR.pars)=c('beta','alpha')
  #coefficient C estimates
  SIR.coefs   = resultList$coefs
  
  #update I0 and plot the fit to the model (plot derivatives against the model evaluated at solutions 
  #(obtained by approximations via basis functions)) so this is model fidelity
  traj=SIR.proc$bvals$bvals%*%SIR.coefs
  I0=traj[1,2]
  dtraj = SIR.proc$bvals$dbvals%*%SIR.coefs
  
  colnames(traj) = SIR.names
  ftraj = SIR.proc$more$fn(times=SIR.proc$more$qpts,y=traj,p=SIR.pars,more=SIR.proc$more$more)
  
  X11()
  par(mfrow=c(2,2),mai=c(0.3,0.6,0.1,0.1))
  for(i in 1:3){
    plot(SIR.knots,dtraj[,i],type='l',xlab='',ylab=SIR.names[i])
    lines(SIR.knots,ftraj[,i],col='red',lty=2)
    abline(h=0)
    legend('topleft',legend=c('Smooth','Model'),lty=1:2,col=1:2)
  }
  
  matplot(SIR.knots,dtraj,type='l',lty=1,
          ylab='SIR derivatives')
  matplot(SIR.knots,ftraj,type='l',lty=2,add=TRUE)
  
  
  
  CollocInferPlots(SIR.coefs,SIR.pars,SIR.lik,SIR.proc,times=times,data=logdata,
                   cols=NULL,datacols=NULL,datanames=NULL,ObsPlot=TRUE,DerivPlot=TRUE,
                   cex.axis=1.5,cex.lab=1.5,cex=1.5,lwd=2)
  
  
  
  
  
  
  var_cov = Profile.covariance(pars=SIR.pars,active=NULL,times=times,data=dat,
                               coefs=SIR.coefs,lik=SIR.lik,proc=SIR.proc,
                               control.in = control.in)
  is.positive.definite(var_cov)
  
  SE = sqrt(diag(var_cov))
  
  
  
  ci_lower=as.numeric(SIR.pars)-qnorm(0.975)*SE
  ci_upper=as.numeric(SIR.pars)+qnorm(0.975)*SE
  
  confInt=c(c(ci_lower,rep(NA,1),c(ci_upper,rep(NA,1))))
  
  
  return(list(parameters=c(SIR.pars,I0=I0),SE=c(SE,rep(NA,1)),confInt=confInt,dat=dat,IniPar=c(IniPar,I0=NA)))
  
}

#############################################################################
#Rhat  :         obtain ODE solution using the Runge Kuta method
# input:         thetas       - parameters 
#                dat          - data
#                N=25         - data points
#                times        - times
#                rez          - residuals
#output:         time series of the states of the ODE model
#############################################################################



Rhat=function(thetas,dat,N=261,times,rez=NULL){
  
  pars = thetas[1:2]
  I0=thetas[4];
  names(pars) = c("beta","alpha");
  N = 261
  x0 = c(N-I0,I0,0)
  names(x0) = varnames = c("S","I","R")
  SIR = make.SIR()
  
  
  sol  = rk4(x0,times,SIR$fn.ode,pars)
  
  # plot(sol[,2])
  # lines(sol[,3])
  # lines(sol[,4],col='green')
  # lines(dat[,3],col='red')
  
  fits=sol
  
  if (any(is.na(fits[,2])) |  any(is.na(fits[,3]))  | any(is.na(fits[,4])) ){
    fits[which(is.nan(fits[,2]) | (fits[,2]>N)),2]  = 0
    fits[length(times)-1,2]                         = 1
    fits[which(is.nan(fits[,3]) | (fits[,3]<0)),3]  = N-1
    fits[which(is.nan(fits[,4]) | (fits[,4]<0) ),4] = 1
  }
  
  
  if (!is.null(rez)){
    if (!any(is.na(fits[,4]))){
      resid<-(dat[,3]-fits[,4])
      #ind1=which(!is.na(data[,2]))
      #resid1<-(data[ind1,2]-sol[ind1,2])
      fitres  = lokerns(times,resid,inputb=T,x.out = times,bandwidth = rep(1,length(times)) )
      #fitres1 = lokerns(times[ind1],resid1,inputb=T,bandwidth = rep(19,length(times[ind1])) ,x.out = times[ind1])
      fits[,4]<-sol[,4]+fitres$est
      #if (lambda<1e-9) fits<-sol[,3]
    }  
  }  
  # plot(sol[,4])
  # points(fits[,4],col='red')
  #lines(fits[,4],col='blue')
  return(fits)
}

#############################################################################
#ObtainWeights : obtain the weights using the likelihood
# input:         thetas       - parameters 
#                dat          - data
#                N=25         - data points
#                times        - times
#output:         calculated weigths 
#############################################################################


ObtainWeights=function(thetas,dat,N,times){
 
  rksol=Rhat(thetas,dat,N,times)
  R=rksol[,4]
  I=rksol[,3]
  w= SIRloglik(y=dat[,ncol(dat)],R=R,N=N)
  #w=ifelse(is.nan(w),0,w)
  return(w)
}

#############################################################################
#runOptimizers  : run the different optimizers
#input:           IniPars     - initial parameters
#                 times       - times of the data points
#                 dat         - data
#                 runOptim    - 
#                 D           - number of optimizers
#                 cov_prior   - 
#                 N=261
#output:       obtain optimized parameters and Hessians
#############################################################################

runOptimizers=function(IniPars=NULL,times,dat,runOptim,D=3,cov_prior=NULL,N=261){
  
  parameters=list()
  IniParam=list()
  hessian=list()
  
  i=1
  if (!is.na(runOptim[1])){
    
    NLSstep_out= NLSStep(IniPar=IniPars,times=times,dat=dat)
  
    parameters[[i]]           =  cbind(NLSstep_out$parameters[,1:2],N-NLSstep_out$parameters[,3],NLSstep_out$parameters[,3])
    IniParam[[i]]             =  c(NLSstep_out$IniPar[1:2],N-NLSstep_out$IniPar[3],NLSstep_out$IniPar[3])
    hessian[[i]]              =  NLSstep_out$hessian
  }
  output=list(parameters=parameters, IniParam= IniParam,hessian=hessian)
  return(output)
}


#############################################################################################
#runIMIS_Opt:     Run the IMIS-Opt algorithm
#input:  cl        - cluster for parallel computing 
#        N0        - initial sample size
#        D         - number of initial points for the optimization
#        B         - number of incremental samples
#        J         - number of re-sampled samples
#        niter     - max number of iterations
#        priorpars - hyperparameters of the priors 
#        N         - population size  
#        dat       - data
#        times     - times
#        par_true  - true parameters
#output:  a list of 
#         theta       - re-sampled parameters
#         theta_all   - all samples before re-sampling
#         wts         - weights
#         k           - k
#         center      - centers of the multivariate normal distributions, 
#                       obtained from either the optimization or from the importance stage
#                       as maximum weight points
#         H_k         - evaluated Gaussians used to normalize the weights
#         norm_all    - normalization of the weigths
#         log_lik_all - log likelihood estimates of all points in theta_all
#         prior_all   - prior evaluated at all points in theta_all
#############################################################################################

runIMIS_Opt=function(cl,N0=3*1000,D=3,B=1000,J=3000,niter=100,priorpars=c(1,1,1,1,N,5/N),N,dat,times,par_true=c(0.0062,0.098,256,5)){

  theta_all   = theta_k=sampleprior(N0,priorpars=priorpars)
  cov_prior_global=cov(theta_all[,c(1,2)])
  stat_all = matrix(NA, 6, niter)
  
  theta_max     = matrix(NA,niter,N0)
  theta_max[1,] = rep(0,N0)
  sigma_max     = list()
  wts           = list()
  
  log_lik_all =prior_all =NULL

  #library(deSolve)
  
  for (k in (1:niter)) {
    
  #calculate weights
    log_lik_all=c(log_lik_all,parApply(cl,theta_k,1,SIRloglik,y=dat,N=N,times=times))
    prior_all= c(prior_all, apply(theta_k,1,logprior_alphabeta,priorpars=priorpars))
    prior_all[which(is.nan(prior_all))]=-999999
    
    if (k == 1) 
      norm_all = prior_all
    if (k > 1) 
      norm_all=log(apply(rbind(exp(prior_all) * N0/B, H_k), 2, sum)/(N0/B + D + (k - 2)))
      norm_all[which(is.nan(norm_all))]=-999999
      #the initial weights are calculated only from the likelihood
      wts1 = prior_all + log_lik_all - norm_all
    
    logc=max(log_lik_all[which(!log_lik_all==0)])
    wts1=wts1-logc
    wts1=exp(wts1)
    wts1[which(is.nan(wts1))]=0
    wts1[is.na(wts1)]=0
    stat_all[1,k] = log(mean(wts1))			# the raw marginal likelihood
    wts1= wts1/sum(wts1)
    stat_all[2,k] = sum(1-(1-wts1)^J)		# the expected number of unique points
    stat_all[3,k] = max(wts1)				# the maximum weight
    stat_all[4,k] = 1/sum(wts1^2)			# the effictive sample size
    stat_all[5,k] = -sum(wts1*log(wts1), na.rm = TRUE) / log(length(wts1))	# the entropy relative to uniform
    stat_all[6,k] = var(wts1/mean(wts1))	# the variance of scaled weights
    if (k==1)	print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(c(k, round(stat_all[1:4,k], 3)))
    wts[[k]]=wts1    
    
  
 
  #Stage 1
  #run the optimizers
  
  #optimRezults=runOptimizers(IniPars=NULL,times=times,dat=dat,D=D, runOptim=c(1,NULL,NULL))
  
  #Stage 2
  #for each optimized value calculate Mahalanobis distance between theta_opt in 
  #the rest of the values in theta0. Remove the (N0/number of optimized) values
  #i.e. remove the points from th modes around the optimized thetas 
  #and then repopulate with additional B multivariate Gaussian points
  if (k==1 && D>0){
   cov_prior_global = cov(theta_all[which(log_lik_all > min(log_lik_all)),c(1,2)])
   
   center=matrix(NA,D,ncol(theta_all))
   theta_k=excluded=theta_max=NULL
   d=1
 
   for (j in (1:D)){ 
     
     if (j==1){
       indmax=which(wts1==max(wts1))  
     }else{
       indmax=which(wts1==max(wts1[-excluded]))
     }
     
     if (length(indmax)>1){indmax=sample(indmax,1)}
     theta_max=theta_all[indmax,]
     
     print(theta_max)
     optimRezults=OptimizeSIRPost(theta=theta_max,y=dat,priorpars=priorpars,N=N,times=times,cov_prior_global=cov_prior_global)    
     sigma_max[[j]]=matrix(NA,ncol(center)-1,ncol(center)-1)
     #remove them from theta0 points that are close to the mode  
     #based on Mahalanobis distance
     m=optimRezults$parameters[c(1,2)]
     sigma_max[[j]] =solve(optimRezults$hessian)
     
     center[j,] =c(m[1:2],optimRezults$parameters[3],optimRezults$parameters[4])
     
  

     A=theta_all[,c(1,2)]
     Dis=t(sapply(1:nrow(A), function(x) A[x,]-m ))  
     mahalanobis= sapply(1:nrow(Dis), function(x)  (t(Dis[x,]) %*% sigma_max[[j]] %*% Dis[x,]))
     
     excluded=c(excluded,which(mahalanobis %in% sort(mahalanobis)[1:(nrow(theta_all)/D)]))
     excluded=unique(excluded)

    
    repop=rmvnorm(B,m, sigma=sigma_max[[d]])
  
    
    theta_k=rbind(theta_k,cbind(repop[,1:2],optimRezults$parameters[3],optimRezults$parameters[4]))
   }
   theta_all=rbind(theta_all,theta_k)
  }
   

  #Stage 3
  #Importance centered Sample updating

  #set up the objects prior to running the iterations

    #evaluate likelihood for the new samples (1,..,B)
    #and then update the weigths for 1,..,Nk
  if (k > 1 | D == 0) {
 
  #step1. Choose current max weight and draw B new samples around this point
  indmax=which(wts[[k]]==max(wts[[k]]))
  if (length(indmax)>1){indmax=sample(indmax,1)}
  theta_max=theta_all[indmax,]
  
  cov_current=cov(theta_all[which(log_lik_all > min(log_lik_all[-excluded])),c(1,2)])
  
  Nk=nrow(theta_all)
  

  A=theta_all[,c(1,2)]
  m=theta_all[indmax,c(1,2)]
  center=rbind(center,c(m[1:2],theta_max[3],theta_max[4]))
  Dis=as.matrix(t(sapply(1:nrow(A), function(x) A[x,]-m ))  )
  
  
  distance_all = mahalanobis(A, m, diag(diag(cov_current)))
  
  label_dist = sort(distance_all, decreasing = FALSE, index = TRUE)
  
  indB           = setdiff(label_dist$ix[1:B],which(is.nan(A[,2])))
  
  
  
  cov_weighted=cov.wt(x=A[indB,],wt=wts[[k]][indB]+(1/length(wts[[k]][-excluded])),cor=FALSE,center=m,method='unbias')
  sigma_max[[D + k - 1]]=cov_weighted$cov
  repop_k=rmvnorm(B,m,sigma=cov_weighted$cov)
                    

  
  repop_k_4=cbind(repop_k[,1:2],theta_max[3],theta_max[4])
  
  theta_k = repop_k_4
  theta_all=rbind(theta_all,repop_k_4)
  #step2.Calculate importance weights of the new samples
  
   
  }
  
  if ((k==1)){
  
  H_k=matrix(0,D,N0 + D * B)
      for (j in (1:D)){ 
        H_k[j,]=dmvnorm(theta_all[,c(1,2)],center[j,][c(1,2)], sigma=sigma_max[[j]])
        
      }
  }else{
    H_k_new=matrix(0,D+k-1,nrow(theta_all))
    H_k_new[1:(D + k - 2),(1:(dim(theta_all)[1]-B))] =H_k
    H_k_new[D + k - 1,] =dmvnorm(theta_all[,c(1,2)],center[D+k-1,c(1,2)], sigma=sigma_max[[D + k - 1]])
    
    for (j in (1:(D + k - 2))) {
      H_k_new[j, (dim(theta_all)[1] - B + 1):dim(theta_all)[1]] = dmvnorm(theta_k[,c(1,2)], center[j,c(1,2)], sigma_max[[j]])
    }
    H_k   = H_k_new
  }
  
  if (sum(1-(1-wts[[k]])^J)>(1-exp(-1))*J) break
  }
  
  #stage 4 resample J
  nonzero     = which(wts[[k]] >0)
  
  
  indResample = sample(nonzero, J, replace = TRUE, prob = wts[[k]][nonzero])
  theta       = theta_all[indResample,]
  
  out=list(stat_all=stat_all,theta=theta,theta_all=theta_all,wts=wts[[k]],k=k,center=center,sigma_max=sigma_max,H_k=H_k,norm_all=norm_all,log_lik_all=log_lik_all,prior_all=prior_all) 

  return(out)
 
}

