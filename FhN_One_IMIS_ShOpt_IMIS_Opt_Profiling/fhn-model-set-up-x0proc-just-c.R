

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
likelihood <- function(pars,t=times,logs=FALSE,parnames = c("c"),data){	
# note that this requires the data.frame to be attached so that V and R are included
# note that x0.conditions[,7:8] are altered to reflect the given parameter values


	if(is.vector(pars)||is.matrix(pars)){
		loglike=rep(1,length(pars))
		res = foreach(index=1:(length(pars))) %dopar%{
			params=c(.2,.2,pars[index])
			names(params)=list("a","b","c")
			sd=.05
			loglikelihood.lo =-Inf
			x0               = c(-1,1);
			names(x0)        = c("V","R");

			y.fit            = lsoda(x0,t,make.fhn()$fn.ode,params) 
				
			loglikel         = dnorm(matrix(data,ncol=1),matrix(y.fit[,-1],ncol=1),sqrt(sd),log=TRUE)
			
			if(max(t)==max(y.fit[,1])){
			return( sum(loglikel[!is.na(matrix(data,ncol=1))]))
			}else{return( -Inf)}
      }
      loglike            = unlist(res)
	}else{
      index              = 1
	  	loglike            = 0
  		params             = c(.2,.2,pars)
	  	names(params)      = list("a","b","c")
	  	sd                 = .05
	  	loglikelihood.lo   = -Inf
		  x0                 = c(-1,1);
		  names(x0)          = c("V","R");
		  y.fit              = lsoda(x0,t,make.fhn()$fn.ode,params) 
		  loglikel           = dnorm(matrix(data,ncol=1),matrix(y.fit[,-1],ncol=1),sqrt(sd),log=TRUE)
		  loglike            = sum(loglikel[!is.na(matrix(data,ncol=1))])
		if(max(t)==max(y.fit[,1])){
			loglike            = sum(loglikel[!is.na(matrix(data,ncol=1))])
		}else{loglike        = -Inf}
	}
	if(logs==FALSE){
	return(exp(loglike))}else{ return(loglike)}
}


##########################################################################
#prior: Set up the prior distributions and their parameters
#input: pars      - vector of parameters at which the prior has to be evaluated
#       priorpar  - setup the parameters of the prior distribution
#       logs      - TRUE or FALSE, whether to evaluate log prior
#output: a point of evaluated prior at the pars provided
##########################################################################

prior <- function(pars,priorpar=c(14,2),logs=FALSE){
	if(is.vector(pars)||is.matrix(pars)){
		res = foreach(index=1:(length(pars))) %dopar%{
			
			if(pars[index]<0 ){
			   return(-Inf)
			   }else{
			return(dnorm(pars[index],priorpar[1],priorpar[2],log=T))
			   }		
		}
		priorvec = unlist(res)
	}	else   {	
		if(pars<0 ){
			return(-Inf)
		}else{priorvec= dnorm(pars,priorpar[1],priorpar[2],log=T)}
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
sample.prior <- function(n=1,priorpar=c(14,2)){
  rnorm(n,priorpar[1],priorpar[2])
}

################################################################
#neglogprior: Set up the negative log prior 
# input:  pars
#         priorpar
#output:  evaluated negative log prior at provided pars
################################################################

neglogprior <- function(pars,priorpar=c(14,2)){
	return(-prior(pars,priorpar,logs=T))
}


################################################################
#dneglogpriordpar - evaluate first derivative of the prior
# input  : pars
# input:  pars
#         priorpar
#output:  evaluated first derivative of the negative log prior 
#         at provided pars 
################################################################

dneglogpriordpar <- function(pars,priorpar=c(14,2)){

	if(is.vector(pars)||is.matrix(pars)){
			res = foreach(index=1:(length(pars))) %dopar%{
				return(  dnegnormdp(pars[index],priorpar[1],priorpar[2]))
			}
			dneglogpriordpar = matrix(unlist(res),length(pars),byrow=T)
			return(  dneglogpriordpar)
	}	else   {
		dneglogpriordpar=   c(dnegnormdp(pars,priorpar[1],priorpar[2]))
			return(  dneglogpriordpar)
	}

}


################################################################
#d2neglogpriordpar2 - evaluate second derivative of the prior
# input  : pars
# input:  pars
#         priorpar
#output:  evaluated first derivative of the negative log prior 
#         at provided pars 
################################################################

d2neglogpriordpar2 <- function(pars,priorpar=c(14,2)){

	if(is.vector(pars)||is.matrix(pars)){
			res = foreach(index=1:(length(pars))) %dopar%{
				return(  d2negnormdp2(pars[index],priorpar[1],priorpar[2]))
			}
			d2neglogpriordpar2 = matrix(unlist(res),length(pars),byrow=T)
			return(  d2neglogpriordpar2)
	}	else   {
		d2neglogpriordpar2=   c(d2negnormdp2(pars,priorpar[1],priorpar[2]))
			return(  d2neglogpriordpar2)
	}

}



dnegnormdp<-function(pars,priorpar1,priorpar2){
		return((pars-priorpar1)/priorpar2^2)}

d2negnormdp2<-function(pars,priorpar1,priorpar2){
		return(1/priorpar2^2)}
	

		





