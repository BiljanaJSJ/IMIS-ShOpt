


##########################################################################
#likelihood   - evaluate likelihood  A Gaussian likelihood centered around the ODE solution
#input    -  pars                a vector of parameters, or a matrix where each row is another 
#                                set of parameters
#         -  t
#         -  logs                evaluate log- likelihood
#         -  parnames            par names -needed for the numerical solver of the ODE
#         -  data                the data
#output:  evaluated likelihood or log likelihood, for each provided set of parameters
##########################################################################
likelihood <- function(pars,t=times,logs=FALSE,parnames = c("a","b","c","sdV","sdR","V0","R0"),data=datamatrix){	

	if(is.matrix(pars)){
		loglike=rep(1,dim(pars)[1])
		res = foreach(index=1:(dim(pars)[1])) %dopar%{
			if(max(pars[index,3:5]<=0)){return(-Inf)}else{
						#		for(index in 1:(dim(pars)[1])){
			params=c(pars[index,])
			names(params)=parnames
			sdV=params[4]
			sdR=params[5]                
			loglikelihood.lo=-Inf
			x0 = c(params[6],params[7]);
			names(x0) = c("V","R");

			y.fit = lsoda(x0,t,make.fhn()$fn.ode,params) 
				
			loglikel = dnorm(matrix(y.fit[,2],ncol=1),matrix(data[,1],ncol=1),sqrt(sdV),log=TRUE)
			loglikel = c(loglikel,dnorm(matrix(y.fit[,3],ncol=1),matrix(data[,2],ncol=1),sqrt(sdR),log=TRUE))                
			
			if(max(t)==max(y.fit[,1])){
			return( sum(loglikel[!is.na(matrix(data,ncol=1))]))
			}else{return( -Inf)}}
		}
# print(res)
loglike = unlist(res)
	}else{
		if(max(pars[3:5]<=0)){return(-Inf)}else{
			
		index =1
		loglike=0
		params=pars
		names(params)=parnames
		sdV=params[4]
		sdR=params[5]            
		loglikelihood.lo=-Inf
		x0 = c(params[6],params[7]);
		names(x0) = c("V","R");
		
		y.fit = lsoda(x0,times,make.fhn()$fn.ode,params) 
		
		loglikel = dnorm(matrix(y.fit[,2],ncol=1),matrix(data[,1],ncol=1),sqrt(sdV),log=TRUE)
		loglikel = c(loglikel,dnorm(matrix(y.fit[,3],ncol=1),matrix(data[,2],ncol=1),sqrt(sdR),log=TRUE))
            
		loglike =  sum(loglikel[!is.na(matrix(data,ncol=1))])
		if(max(t)==max(y.fit[,1])){
			loglike =  sum(loglikel[!is.na(matrix(data,ncol=1))])
		}else{loglike =  -Inf}}
	}
	if(logs==FALSE){
	return(exp(loglike))}else{ return(loglike)}
}















##########################################################################
#prior: Set up the prior distributions and their parameters,
#       priors are a~ N(0,.4^2),  b~ N(0,.4^2),    c~N(14,2) = Gam(n/2,.5) here n=2, 
#       sigma^2 ~ IG(2,2)  , V0~N(-1,2^2), R0~  N(1,2^2)
#input: pars      - vector of parameters at which the prior has to be evaluated
#       priorpar  - setup the parameters of the prior distribution
#       logs      - TRUE or FALSE, whether to evaluate log prior
#output: a point of evaluated prior at the pars provided
##########################################################################
prior <- function(pars,priorpar=cbind(c(0,0,14,3,3),  c(.4,.4,2,3,3)),logs=FALSE){
	if(is.matrix(pars)){
		res = foreach(index=1:(dim(pars)[1])) %dopar%{
			
			if(min(pars[index,3:5])<0 ){
			   return(-Inf)
			   }else{
			return(    sum(c(dnorm(pars[index,1:3],priorpar[1:3,1],priorpar[1:3,2],log=T),
                   dinvgamma(pars[index,4:5],priorpar[4:5,1], priorpar[4:5,2],log=T))))
			   }	
		}
		priorvec = unlist(res)
	}	else   {	
		if(min(pars[3:5])<0 ){
			return(-Inf)
		}else{priorvec=    sum(dnorm(pars[1:3],priorpar[1:3,1],priorpar[1:3,2],log=T),
            dinvgamma(pars[4:5],priorpar[4:5,1], priorpar[4:5,2],log=T))}
	}

	if(logs==FALSE){
		return(exp(priorvec))}else{ return(priorvec)}
}

##################################################################################
#sample.prior   draw samples from the prior
#input:        n            - number of samples from the prior
#              priorpar     - parameters of the prior distribution
#output:       n samples from the prior
##################################################################################	
sample.prior <- function(n=1,priorpar=cbind(c(0,0,14,3,3,-1,1),  c(.4,.4,2,3,3,.5 ,.5))){
  rparmat=matrix(0,nrow=n,ncol=7)
  colnames(rparmat)=c("a","b","c","sdV","sdR","V0","R0")
  
  
  for(index in c(1,2,3,6,7)){
    rparmat[,index] = 	rnorm(n,priorpar[index,1],priorpar[index,2])
  }
  rparmat[,4] = 	1/rgamma(n,priorpar[4,1],priorpar[4,2])	
  rparmat[,5] = 	1/rgamma(n,priorpar[5,1],priorpar[5,2])	
  return(rparmat)
}




#######################################################################################
#priorx0 prior distribution on the initial states
#input: pars      - vector of parameters at which the prior has to be evaluated
#       priorpar  - setup the parameters of the prior distribution
#       logs      - TRUE or FALSE, whether to evaluate log prior
#output: a point of evaluated prior at the pars provided
#######################################################################################
priorx0 <- function(pars,priorpar=cbind(c(-1,1),  c(.5,.5)),logs=FALSE){
	if(is.matrix(pars)){
		end=dim(pars)[2]
		res = foreach(index=1:(dim(pars)[1])) %dopar%{
			
			return(    sum( dnorm(pars[index,(end-1):end],priorpar[,1],priorpar[,2],log=T)))
			
		}
		priorvec = unlist(res)
	}	else   {
		end=length(pars)
	priorvec=    sum(dnorm(pars[(end-1):end],priorpar[,1],priorpar[,2],log=T))
	}

	if(logs==FALSE){
	return(exp(priorvec))}else{ return(priorvec)}
}



##############

# density of an inverse gamma density

dinvgamma <- function(x,alpha,beta,log=TRUE){

	if(min(x)<=0){logdensity = 0}else{
	logdensity =alpha*log(beta)-log(gamma(alpha)) - (alpha+1)*log(x) - beta/x
	}
	if(log==FALSE){return(exp(logdensity))}else{return(logdensity)}
}


################################################################
#neglogprior: Set up the negative log prior 
# input:  pars
#         priorpar
#output:  evaluated negative log prior at provided pars
################################################################

neglogprior <- function(pars,priorpar=cbind(c(0,0,14,3,3),  c(.4,.4,2,3,3))){
	return(-prior(pars,priorpar,logs=T))
}

################################################################
#neglogpriorx0: Set up the negative log prior for the initial states
# input:  pars
#         priorpar
#output:  evaluated negative log prior at provided pars
################################################################


neglogpriorx0 <- function(pars,priorpar=cbind(c(-1,1),  c(.5,.5)),logs=FALSE){
	return(-priorx0(pars,priorpar,logs=T))
}



################################################################
#dneglogpriordpar - evaluate first derivative of the prior
# input  : pars
# input:  pars
#         priorpar
#output:  evaluated first derivative of the negative log prior 
#         at provided pars 
################################################################

dneglogpriordpar <- function(pars,priorpar=cbind(c(0,0,14,3,3),  c(.4,.4,2,3,3))){

	if(is.matrix(pars)){
			res = foreach(index=1:(dim(pars)[1])) %dopar%{
				return(  c(dnegnormdp(pars[index,1:3],priorpar[1:3,1],priorpar[1:3,2]),
                dnegIgammadp(pars[index,4:5],priorpar[4,1],priorpar[4,2]))
						   )
				
			}
			dneglogpriordpar = matrix(unlist(res),dim(pars)[1],byrow=T)
			return(  dneglogpriordpar)
	}	else   {
		dneglogpriordpar=   c(dnegnormdp(pars[1:3],priorpar[1:3,1],priorpar[1:3,2]),
        dnegIgammadp(pars[4:5],priorpar[4:5,1],priorpar[4:5,2])
										)
			return(  dneglogpriordpar)
	}

}


################################################################
#dneglogpriordx0 - evaluate first derivative of the prior for the
#                  initial states
# input  : pars
# input:  pars
#         priorpar
#output:  evaluated first derivative of the negative log prior 
#         at provided pars 
################################################################
dneglogpriordx0 <- function(pars,priorpar=cbind(c(-1,1),  c(.5,.5))){
	
	if(is.matrix(pars)){
		end=dim(pars)[2]
		res = foreach(index=1:(dim(pars)[1])) %dopar%{
			return(    dnegnormdp(pars[index,(end-1):end],priorpar[,1],priorpar[,2]))
		}
		dneglogpriordpar = matrix(unlist(res),dim(pars)[1],byrow=T)
		return(dneglogpriordpar)
	}	else   {
		end=length(pars)
		dneglogpriordpar=   dnegnormdp(pars[(end-1):end],priorpar[,1],priorpar[,2])
		return(   dneglogpriordpar)
	}
	
}

################################################################
#d2neglogpriordx02 - evaluate second derivative of the prior
#                    for the initial states
# input  : pars
# input:  pars
#         priorpar
#output:  evaluated first derivative of the negative log prior 
#         at provided pars 
################################################################
d2neglogpriordx02<- function(pars,priorpar=cbind(c(-1,1),  c(.5,.5))){

	if(is.matrix(pars)){
		end=dim(pars)[2]
		res = foreach(index=1:(dim(pars)[1])) %dopar%{
			return( d2negnormdp2(pars[index,(end-1):end],priorpar1=priorpar[,1],priorpar2=priorpar[,2]))
		}
		dneglogpriordpar = matrix(unlist(res),dim(pars)[1],byrow=T)
		return(   dneglogpriordpar)
	}	else   {
		end=length(pars)
		dneglogpriordpar=   d2negnormdp2(pars[(end-1):end],priorpar1=priorpar[,1],priorpar2=priorpar[,2])
		return(   dneglogpriordpar)
	}
	
}


dnegnormdp<-function(pars,priorpar1,priorpar2){
		return((pars-priorpar1)/priorpar2^2)}
	
d2negnormdp2<-function(pars,priorpar1,priorpar2){

	return(diag(1/priorpar2^2))}

										 
dneggammadp<-function(pars,priorpar1,priorpar2){
		return(-(priorpar1-1)/pars+priorpar2)}
										 
										 
dnegIgammadp<-function(pars,priorpar1,priorpar2){
		return(
			   (priorpar1+1)/pars-1/(priorpar2*pars^2))
										 }
##########################################################################################
#ident    -- NLS method, smooth the data and then estimate parameters
#input:    pars   - a vector of parameters
#          data   - the data
#          times  - a vector of times of observations
#          other  - other input
#output: optimized parameters
##########################################################################################										 
ident <-function(pars=X_imp,data=data,times=times,other=list(method="BFGS")){
  
  post = function(pars,data=data){sum(neglogprior(pars))-likelihood(pars,logs=TRUE,data=data) }
  
  output<-optim(pars, post, data=data, method=other$method, hessian=FALSE, control=list(maxit=1000))
  return(pars)
}


##########################################################################################
#two.stager   -- Two-Stage method, smooth the data and then estimate parameters
#input:    pars   - a vector of parameters
#          data   - the data
#          times  - a vector of times of observations
#          other  - other input
#output: optimized parameters
##########################################################################################
two.stager<-function(pars=X_imp,data=data,times=times,other=NULL){
  return(twostage.optim.fhn(pars=pars,data=data,times=times))}



##########################################################################################
#smoother.opt   -- GP method, smooth the data, perform inner and outer optimization
#input:    pars   - a vector of parameters
#          data   - the data
#          times  - a vector of times of observations
#          other  - other input
#output:  optimized parameters
##########################################################################################

smoother.opt <- function(pars=X_imp,data=data,times=times,other){
  coefs       = other$coefs
  ncoefs       = other$ncoefs
  lik         = other$lik
  proc        = other$proc
  names(pars) = proc$more$parnames
  control.in  = other$control.in
  control.out = other$control.out
  ncoefs      = other$ncoefs
  
  # first perform a model based smooth of the data
  res0 = inneropt(coefs, times=times, data=data, lik=lik, proc=proc, pars=pars, in.meth="nlminb", control.in=control.in)
  res1 = outeropt(data=data, times=times, pars=pars[1:(length(pars)-ncoefs)], coefs=res0$coefs, lik=lik, proc=proc, in.meth="nlminb",out.meth = "nlminb",control.in=control.in, control.out=control.out)
  if(ncoefs>0){
    output<-c(res1$pars,res1$coefs[1,1:ncoefs])
  }else{
    output<-c(res1$pars)
  }
  
  return(output)
}



