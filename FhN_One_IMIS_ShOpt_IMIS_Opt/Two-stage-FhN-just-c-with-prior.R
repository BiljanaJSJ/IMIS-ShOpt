#####################################################
#The Two-stage method
#####################################################
library(deSolve)
library(lokern)











#-------------------------------------------------
# FitzHugh-Nagumo system
#-------------------------------------------------


fhn.fun=function(t,y,parms) 
{
    u=y[1]; v=y[2];
    du= parms[3]*(u+v-u^3/3);
    dv=-(u-parms[1]+parms[2]*v)/parms[3];
    dY=c(du,dv);
    return(list(dY));
}

der.fhn<-function(b,du,dv,u,v)
{
# d/dp [Kernel Smooth Derivative - FhN(Kernel Smooth)]^2
#parameter b=c(alpha, beta, gamma)
#du: du/dt
#dv: dv/dt

a=b[1]
beta=b[2]
g=b[3]

arg1<- u-a+beta*v
arg2<- u+v-u^3/3
dpr <- dneglogpriordpar(b,priorpar=c(14,2))
d.a <- -sum(dv)/g-sum(arg1)/(g^2)+ dpr[1]        #d(obj)/d(alpha)
d.b <- sum(v*dv)/g+sum(arg1*v)/(g^2) + dpr[2]    #d(obj)/d(beta)
d.g <- -sum(arg2*du)+g*sum(arg2^2)-sum(dv*arg1)/(g^2)-sum(arg1^2)/(g^3) + dpr[3]
                    #d(obj)/d(gamma)
c(d.a,d.b,d.g)
}

jac.fhn<-function(b,du,dv,u,v)
{
# calculate Jacobian
#parameter b=c(alpha, beta, gamma)
#du: du/dt
#dv: dv/dt

a=b[1]
beta=b[2]
g=b[3]
n<-length(u)
arg1<- u-a+beta*v
arg2<- u+v-u^3/3

J = Jpr = d2neglogpriordpar2(b,priorpar=c(14,2))

J[1,1]<- n/(g^2)
J[1,2]<-J[2,1]<- -sum(v)/(g^2)
J[1,3]<-J[3,1]<- sum(dv)/(g^2)+2*sum(arg1)/(g^3)

J[2,2]<- sum(v^2)/(g^2)
J[2,3]<-J[3,2]<- -sum(v*dv)/(g^2)-2*sum(arg1*v)/(g^3)
J[3,3]<- sum(arg2^2)+2*sum(dv*arg1)/(g^3)+3*sum(arg1^2)/(g^4)
return(J)
}

linear.fhn.fit<-function(du,dv,u,v)
{
n<-length(u)
arg2<- u+v-u^3/3

J<-matrix(NA,3,3)
J[1,1]<- n
J[1,2]<-J[2,1]<- -sum(v)
J[1,3]<-J[3,1]<- -sum(dv)

J[2,2]<- sum(v^2)
J[2,3]<-J[3,2]<- sum(v*dv)
J[3,3]<- sum(arg2^2+dv^2)

d<-rep(NA,3)
d[1]<- sum(u)
d[2]<- -sum(u*v)
d[3]<- sum(du*arg2-u*dv)

return(solve(J)%*%d)
}








der.fhn.justc<-function(b,du,dv,u,v,...)
{
    #parameter b=c(alpha, beta, gamma)
    #du: du/dt
    #dv: dv/dt
    
    a=.2
    beta=.2
    g=b
    
    arg1<- u-a+beta*v
    arg2<- u+v-u^3/3
    
    #d.a<- -sum(dv)/g-sum(arg1)/(g^2)        #d(obj)/d(alpha)
    #d.b<- sum(v*dv)/g+sum(arg1*v)/(g^2)     #d(obj)/d(beta)
    return( -sum(arg2*du)+g*sum(arg2^2)-sum(dv*arg1)/(g^2)-sum(arg1^2)/(g^3)  + dneglogpriordpar(b,priorpar=c(14,2)))
    #d(obj)/d(gamma)

}

jac.fhn.justc<-function(b,du,dv,u,v,...)
{
    # calculate Jacobian
    #parameter b=c(alpha, beta, gamma)
    #du: du/dt
    #dv: dv/dt
    
    a=.2
    beta=.2
    g=b
    n<-length(u)
    arg1<- u-a+beta*v
    arg2<- u+v-u^3/3
    
    J = d2neglogpriordpar2(b,priorpar=c(14,2))[1]
    J<- sum(arg2^2)+2*sum(dv*arg1)/(g^3)+3*sum(arg1^2)/(g^4) + J
    
    return(J)
}





#-------------------------------------------------
# end for FitzHugh-Nagumo system
#-------------------------------------------------








simex.fun<-function(par,u.obs,v.obs,s1,s2,B=100,...)
{
lambda<-seq(0,2,by=0.2)
m.lam<-length(lambda)

beta.lambda<-var.lambda<-matrix(NA,m.lam,3)
for (s in 1:m.lam)
{
  beta.B<-matrix(NA,B,3)
  for (b in 1:B)
    {

    u.obs.lam<- u.obs+sqrt(lambda[s])*rnorm(length(u.obs),sd=s1)
    v.obs.lam<- v.obs+sqrt(lambda[s])*rnorm(length(v.obs),sd=s2)
    derv.u.b=lokerns(times,u.obs.lam,x.out=times,deriv=1)$est
    derv.v.b=lokerns(times,v.obs.lam,x.out=times,deriv=1)$est
    beta.B[b,]<-neq(b=par,der.fhn,jac.fhn,du=derv.u.b,dv=derv.v.b,
                    u=u.obs.lam,v=v.obs.lam,...)$b
    } # end loop b

    beta.lambda[s,]<-apply(beta.B,2,mean,na.rm=T)
    var.lambda[s,]<-apply(beta.B,2,var,na.rm=T)

} # end loop s

############################
#  quadratic extrapolation  #
#--------------------------#
beta.simex.qua<-var.simex.qua<-rep(NA,3)
lam.arg<-cbind(lambda,lambda^2)
for(qi in 1:3)
{
Gamma.qua<-lm(beta.lambda[,qi]~lam.arg)$coef
beta.simex.qua[qi]<-Gamma.qua[1]-Gamma.qua[2]+Gamma.qua[3]

#Gamma.v.qua<-lm(var.lambda[,qi]~lam.arg)$coef
#var.simex.qua[qi]<-Gamma.v.qua[1]-Gamma.v.qua[2]+Gamma.v.qua[3]
}
#--------------------------------
# end quadratic extrapolation
#################################
return(beta.simex.qua)
}
# end of the function simex.fun()





neq <- function(b,gn,jn,gtol=1e-6,iter=50,stepf=.5,steptol=1e-8,...) {
  #  b = obtain parameters b numerically, by minimizing the model residuals
  # gn = derivative of model residual {sum(FhN via Kernel Smooth States - Kernel derivatives)^2} wrt pars:   der.fhn
  # jn = Jacobian of model residual (FhN) with respect to pars x pars
  # ... extra inputs for the functions
  
  n <- length(b)
  steph <- 0
  g0 <- gn(b,priorpar=c(14,2),...)
  f0 <- sum(g0*g0)/2
  
  for (ll in 1:iter) {
    j <- jn(b,priorpar=c(14,2),...)
    sc <- -c(solve(j,g0))
    bn <- b+sc
    g1 <- gn(bn,...)
    f1 <- sum(g1*g1)/2
    i <- 0
    lam <- -2*f0
    while (is.na(f1) || f1>f0+(1e-4)*lam) {
      i <- i+1
      steph <- steph+1
      sc <- sc*stepf
      lam <- lam*stepf
      bn <- b+sc
      g1 <- gn(bn,...)
      f1 <- sum(g1*g1)/2
      if (i>20) return(list(b=b,f=f0,comp=c(iter=ll,error=1,steph=steph)))
    }
    if (max(abs(g1))<gtol) # if true, iteration converged
      return(list(b=bn,f=f1,comp=c(iter=ll,error=0,steph=steph)))
    if (max(abs(b-bn)/pmax(abs(b),1))<steptol)
      return(list(b=bn,f=f1,comp=c(iter=ll,error=3,steph=steph)))
    b <- bn
    g0 <- g1
    f0 <- f1
  }
  list(b=bn,f=f1,comp=c(iter=ll,error=2,steph=steph))
}






simex.fun.justc<-function(par,u.obs,v.obs,s1,s2,B=100,...)
{
    lambda<-seq(0,2,by=0.2)
    m.lam<-length(lambda)
    beta.lambda<-var.lambda<-matrix(NA,m.lam,3)
    for (s in 1:m.lam)
    {
        beta.B<-matrix(NA,B,3)
        for (b in 1:B)
        {
            u.obs.lam<- u.obs+sqrt(lambda[s])*rnorm(length(u.obs),sd=s1)
            v.obs.lam<- v.obs+sqrt(lambda[s])*rnorm(length(v.obs),sd=s2)
            derv.u.b=lokerns(times,u.obs.lam,x.out=times,deriv=1)$est
            derv.v.b=lokerns(times,v.obs.lam,x.out=times,deriv=1)$est
            beta.B[b,]<-neq(b=par,der.fhn.justc,jac.fhn.justc,du=derv.u.b,dv=derv.v.b,
            u=u.obs.lam,v=v.obs.lam,...)$b
        } # end loop b
        beta.lambda[s,]<-apply(beta.B,2,mean,na.rm=T)
        var.lambda[s,]<-apply(beta.B,2,var,na.rm=T)
        
    } # end loop s
    
    ############################
    #  quadratic extrapolation #
    #--------------------------#
    beta.simex.qua<-var.simex.qua<-rep(NA,length(par))
    lam.arg<-cbind(lambda,lambda^2)
    (qi = 1)
    
        Gamma.qua<-lm(beta.lambda[,qi]~lam.arg)$coef
        beta.simex.qua[qi]<-Gamma.qua[1]-Gamma.qua[2]+Gamma.qua[3]
        
        #Gamma.v.qua<-lm(var.lambda[,qi]~lam.arg)$coef
        #var.simex.qua[qi]<-Gamma.v.qua[1]-Gamma.v.qua[2]+Gamma.v.qua[3]
    
    #--------------------------------
    # end quadratic extrapolation
    #################################
    return(beta.simex.qua)
}
# end of the function simex.fun()






# local qudratic regression

ksLqudratic<- function(x, y, xnew=x, band1=0.5,band2=.25, qud=T)
{
h1 = band1
n = length(y)
if(length(h1)==1){h1=rep(h1,n)}

ones=rep(1,n)
if(qud) {hatm = matrix(NA, n, 3);} # because of qudartic 
if(!qud){hatm = matrix(NA, n, 2);} 

for (i in 1:n)
{
      if(qud) {xi = cbind(ones,(x-x[i]),(x-x[i])^2)}
      if(!qud){xi = cbind(ones,(x-x[i]))}
      k = exp( -xi[,2]^2/(2*h1[i]*h1[i]))/sqrt(2*3.14)/h1[i]
      bfW = diag(k+1.e-10)
      arg=t(xi)%*%bfW
      hatm[i,]=solve(arg%*%xi)%*%(arg%*%y)
}
hatsigma2 = mean((y - hatm[,1])^2)

# calculate the expected function m(x) and density f(x) at xnew

if(qud) {m1 = matrix(NA, length(xnew), 3)}
if(!qud) {m1 = matrix(NA, length(xnew), 2)}

f = matrix(NA, length(xnew), 1)
m1.sig<-m1

h2 = band2
if(length(h2)==1){h2=rep(h2,length(xnew))}

for (i in 1:length(xnew))
{
      if(qud){xi    = cbind(ones,(x-xnew[i]),(x-xnew[i])^2)}
      if(!qud){xi = cbind(ones,(x-xnew[i]))}
      k     = exp( -xi[,2]^2/(2*h2[i]*h2[i]))/sqrt(2*3.14)/h2[i]
      bfW   = diag(k+1.e-10)
      arg = t(xi)%*%bfW
      arg2    = solve(arg%*%xi)
      m1[i,]= arg2%*%(arg%*%y)
      f[i]    = mean(k)
      m1.sig[i,]=c(diag(arg2))
}
m=m1[,1]
L = m - 1.96*sqrt( 0.2821*hatsigma2/(n*h2*f))
U = m + 1.96*sqrt( 0.2821*hatsigma2/(n*h2*f))

L1 = m - 1.96*sqrt( hatsigma2* m1.sig[,1]/h2)
U1 = m + 1.96*sqrt( hatsigma2* m1.sig[,1]/h2)

# pointwise CI of the derivative

L.d = m1[,2] - 1.96*sqrt( hatsigma2* m1.sig[,2] /h2)
U.d = m1[,2] + 1.96*sqrt( hatsigma2* m1.sig[,2] /h2)

return(list(mx=m,CR=cbind(L,U),CR1=cbind(L1,U1), 
    density=f,derv=m1,CR.d=cbind(L.d,U.d),sig_e=hatsigma2))

}







# NOTE: I am actually already estimating the residual measurement error as another parameter.
# It would therefore be worthwhile to choose the bandwidth for which the measurement error matches the parameter OR use the asymptotic residual error 



twostage.optim.fhn.justc <- function(pars,data=datamatrix,times=timepoints,...){
    n=dim(data)[1]

    # ideal asymptotic bandwidth
    fac.arg=max(1.2,n^(2/35)*(log(n))^(-1/4))
    h.dT=lokerns(times,data[,1],x.out=times,deriv=1)$bandwidth*fac.arg
    h.T=lokerns(times,data[,1],x.out=times,deriv=0)$bandwidth*fac.arg
    h.dV=lokerns(times,data[,2],x.out=times,deriv=1)$bandwidth*fac.arg
    h.V=lokerns(times,data[,2],x.out=times,deriv=0)$bandwidth*fac.arg
    
    du.ker.arg=ksLqudratic(times,data[,1],band1=h.dT,xnew=times,band2=h.dT,qud=T)
    u.ker.arg=ksLqudratic(times,data[,1],band1=h.T,xnew=times,band2=h.T,qud=T)
    dv.ker.arg=ksLqudratic(times,data[,2],band1=h.dV,xnew=times,band2=h.dV,qud=T)
    v.ker.arg=ksLqudratic(times,data[,2],band1=h.V,xnew=times,band2=h.V,qud=T)
    
    
    du.ker=du.ker.arg$derv[,2]
    u.ker=u.ker.arg$mx
    dv.ker=dv.ker.arg$derv[,2]
    v.ker=v.ker.arg$mx
    
    #######################Estimates for other parameters
    #resids = c(u.ker-data[,1],v.ker-data[,2])
    #par[4]=var.estimate=var(resids)
    
    #par[5]=v.ker[1]
    #par[6]=u.ker[1]
    
    sigej=c(u.ker.arg$sig_e, v.ker.arg$sig_e)
    
    
    parms=c(.2,.2,pars)
    
    du.rk= parms[3]*(u.ker+v.ker-u.ker^3/3);   # based on the equation
    dv.rk=-(u.ker-parms[1]+parms[2]*v.ker)/parms[3];
    
    #ker.est<-neq(b=rep(1,3),der.fhn,jac.fhn,du=du.ker,dv=dv.ker,u=u.ker,v=v.ker,...)$b
    
    #cal.arg[j,]<-ker.est
    
    #lin.est<-linear.fhn.fit(du=du.ker,dv=dv.ker,u=out[,2],v=out[,3])
    #lin.arg[j,]<-lin.est
    
    s1=0.1; s2=0.1
    return(simex.fun.justc(par=pars,data[,1],data[,2],s1,s2,B=100))
    #    return(c(simex.fun(data[,1],data[,2],s1,s2,B=100),par[4:6],...))
    
}




twostage.optim.fhn<-function(pars=X_imp,data=datamatrix,times=timepoints,...){
n=dim(data)[1]

# ideal asymptotic bandwidth
fac.arg=max(1.2,n^(2/35)*(log(n))^(-1/4))
h.dT=lokerns(times,data[,1],x.out=times,deriv=1)$bandwidth*fac.arg
h.T=lokerns(times,data[,1],x.out=times,deriv=0)$bandwidth*fac.arg
h.dV=lokerns(times,data[,2],x.out=times,deriv=1)$bandwidth*fac.arg
h.V=lokerns(times,data[,2],x.out=times,deriv=0)$bandwidth*fac.arg

du.ker.arg=ksLqudratic(times,data[,1],band1=h.dT,xnew=times,band2=h.dT,qud=T)
u.ker.arg=ksLqudratic(times,data[,1],band1=h.T,xnew=times,band2=h.T,qud=T)
dv.ker.arg=ksLqudratic(times,data[,2],band1=h.dV,xnew=times,band2=h.dV,qud=T)
v.ker.arg=ksLqudratic(times,data[,2],band1=h.V,xnew=times,band2=h.V,qud=T)


   
du.ker=du.ker.arg$derv[,2]
u.ker=u.ker.arg$mx
dv.ker=dv.ker.arg$derv[,2]
v.ker=v.ker.arg$mx

    
    
    u.resid=u.ker-data[,1]
    v.resid=v.ker-data[,2]
    
    var1.est=t(u.resid)%*%u.resid/length(u.resid)
    var2.est=t(v.resid)%*%u.resid/length(v.resid)
    
    
    #######################Estimates for other parameters
    #resids = c(u.ker-data[,1],v.ker-data[,2])
    #par[4]=var.estimate=var(resids)

    #par[5]=v.ker[1]
    #par[6]=u.ker[1]

sigej=c(u.ker.arg$sig_e, v.ker.arg$sig_e)


parms=pars

du.rk= parms[3]*(u.ker+v.ker-u.ker^3/3);   # based on the equation
dv.rk=-(u.ker-parms[1]+parms[2]*v.ker)/parms[3];

    #ker.est<-neq(b=rep(1,3),der.fhn,jac.fhn,du=du.ker,dv=dv.ker,u=u.ker,v=v.ker)$b

    #cal.arg[j,]<-ker.est

    #lin.est<-linear.fhn.fit(du=du.ker,dv=dv.ker,u=out[,2],v=out[,3])
    #lin.arg[j,]<-lin.est

    
    s1=0.1;
    s2=0.1
    return(c(simex.fun(par=pars[1:3],u.obs=data[,1],v.obs=data[,2],s1,s2,B=100,...), var1.est, var1.est,   u.ker[1],v.ker[1]))
    #    return(c(simex.fun(data[,1],data[,2],s1,s2,B=100),par[4:6]))

}



smoother.opt<-function(pars=X_imp,data=data,times=times,other){
    coefs       = other$coefs
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


ident<-function(pars=X_imp,data=data,times=times,other=list(method="BFGS")){

    post = function(pars,data=data){sum(neglogprior(pars)) -likelihood(pars,logs=TRUE,data=data) }

    output<-optim(pars, post, data=data,method=other$method, hessian=FALSE, control=list( maxit=1000))
    return(pars)
}


