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
ident <-function(pars=X_imp,data=data,times=times,other=list(method="BFGS")){
  
  post = function(pars,data=data){sum(neglogprior(pars))-likelihood(pars,logs=TRUE,data=data) }
  
  output<-optim(pars, post, data=data, method=other$method, hessian=FALSE, control=list(maxit=1000))
  return(pars)
}