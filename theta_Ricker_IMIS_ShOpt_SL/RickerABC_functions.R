

#####################################################################################
#functions needed to run the IMIS-ShOpt-ABC
#####################################################################################

#####################################################################################
#theta_ricker    - simulates from the theta-Ricker model
#Input:          N0         the initial popluation
#                theta      the parameter vector
#                timestep   the time steps
#output: a list of 
#                y              - simulated data 
#                x              - at which time points the data were simulated
#                Nt             - the underlying stochastic process fo each data point
#                et             - process noise for each data point
#                theta          - vector of parameters
#                N0.true        - initial population
#####################################################################################
theta_ricker = function(N0,theta,timestep = 100){
  r      = exp(theta[1])
  phi    = theta[2]
  sigma2 = theta[3]
  if (is.na(theta[4])){thet=1
  }else{
    thet   = exp(theta[4])}
  K=100
  
  Nt = rep(N0,timestep)
  et = rep(0,timestep)
  Y = rep(NA,timestep)
  for(timer in 2:timestep){
    et[timer] = rnorm(n=1,0,sd = sqrt(sigma2))
    Nt[timer] = Nt[timer-1]*exp(r*(1-((Nt[timer-1])/K)^thet) +et[timer])
    Y[timer]  = rpois(n=1,phi*Nt[timer])
  }
  return(list(y=Y,x=1:timestep,Nt=Nt,et=et,theta=theta,N0.true=N0))
}
#####################################################################################



#####################################################################################
#prior     - evaluates the prior
#input:
#           theta      - parameter vector
#           meanr      - mean for the log r
#           sdr        - sd for log r 
#           dfPhi      - degrees of reedom for dfphi
#           sdPhi      - sd of the phi
#           sh         - shape of the variance of the process noise (sigma2) 
#           sc         - scale of the variance of the process noise (sigma2)
#           meanThet   - mean of log theta_tilde
#           sdThet     - sd of log theta_tilde
#           log        - boolean: whether to calculate log prior or prior, log=T or log=F
#           reduced    - the model is full or reduced, reduced=NULL corresponds to full model
#output: evaluated prior 
#####################################################################################

prior=function(theta,meanr=0.5,sdr=1,dfPhi=4,sdPhi=1,sh=2,sc=0.05,meanThet=1,sdThet=1,log=T,reduced=NULL){
  if (theta[4]==0 | is.na(theta[4])) {reduced=1}
  if (log){
    prior=dnorm(theta[1],meanr,sdr,log=T)+dchisq(theta[2], df=dfPhi,log=T)+dgamma(1/theta[3],shape= sh, scale=1/sc,log=T)+ifelse(is.null(reduced),dnorm(theta[4],meanThet,sdThet,log=T),0)
    
  }else{
    prior=dnorm(theta[1],meanr,sdr)*dchisq(theta[2], df=4)*dgamma(1/theta[3], shape=sh, scale=1/sc)*ifelse(is.null(reduced),dnorm(theta[4],meanThet,sdThet),1)
  }
  names(prior)=NULL
  return(prior)
}

#####################################################################################
#sample.prior     - draws samples from the prior
#input:
#           N          - number of samples to simulate
#           meanr      - mean for the log r
#           sdr        - sd for log r 
#           dfPhi      - degrees of reedom for dfphi
#           sdPhi      - sd of the phi
#           sh         - shape of the variance of the process noise (sigma2) 
#           sc         - scale of the variance of the process noise (sigma2)
#           meanThet   - mean of log theta_tilde
#           sdThet     - sd of log theta_tilde
#           log        - boolean: whether to calculate log prior or prior, log=T or log=F
#           reduced    - the model is full or reduced, reduced=NULL corresponds to full model
#output: evaluated prior 
#####################################################################################
sample.prior=function(N,meanr=0.5,sdr=1,dfPhi=4,sdPhi=1,sh=2,sc=0.05,meanThet=1,sdThet=1,reduced=NULL){

  r=rnorm(N,meanr,sdr)
  Phi=rchisq(N, df=dfPhi)
  Sigma2=1/rgamma(N,shape=sh,scale= 1/sc)
  if (is.null(reduced)){
    thet=rnorm(N,meanThet,sdThet)
    model=2
  }else{
    thet=0
    model=1}
  
  ret=cbind(r,Phi,Sigma2,thet)
 
  return(ret)
}

##############################################################################################################
# synthloglike  evaluates the negative synthetic log likelihood sum(indicator(matching set defined by distance measure and tolerance level))
#input:
#                theta                - vector of parameters
#                simDs                - k replicate data sets
#                k                    - number of replicate data sets to simulate if simDs is null                   
#                N0                   - initial population for the replicate data sets to simulate if simDs is null     
#                timestep             - time steps at which to simulate replicate data sets  if simDs is null     
#                obs.data             - observed data
#                indices              - indices of the summary statistics, if entire set of summary statistics used, 
#                                       then the likelihood is the target likelihood, otherwise if any subset of 
#                                       summary statistics is used, the likelihood is an approximation to the target one.
#output:         evaluated  negative synthetic log likelihood
##############################################################################################################
synthloglike = function(theta, simDs=NULL,k,N0,timestep,obs.data,indices){
  s=summaryfun(obs.data)

  if (length(theta)==3) {theta=c(theta,0);print(paste("theta=",theta[1],',',theta[2],',',theta[3],sep=''))}
  if (length(theta)==4){print(paste("theta=",theta[1],',',theta[2],',',theta[3],',',theta[4],sep=''))}
  if (theta[2]<0 | theta[3]<0) {synth_log_lik=-999999}else{
   if (is.null(simDs)){
   oout = lapply(rep(list(N0),k),theta_ricker,theta=theta,timestep = timestep)
   }else{
   oout =simDs
   }
    summ = do.call(rbind,applysummaries(oout))
    if (length(indices)==1){
      #get muhat for the summary statistics
      muhat     = mean(summ)[indices]
      #get sigmahat for the summary statistics
      Sigmahat = var(summ)[,indices]
      synth_log_lik=dnorm(s[indices],mean=muhat,sd=Sigmahat,log=T)
    }else{
      muhat = apply(summ,2,mean)[indices]
      Sigmahat = cov(summ)[indices,indices]+10^-03
      synth_log_lik=dmvnorm(s[indices],mean=muhat,sigma=Sigmahat,log=T)  
    }
 
  }
  if (is.na(synth_log_lik) | (synth_log_lik==-Inf)  | (is.nan(synth_log_lik)))  {synth_log_lik=-999999}
  return((-1)*synth_log_lik)
}

logposterior_synthloglike=function(theta,k,N0,timestep,obs.data,indices,reduced=NULL){
  logpost=synthloglike(theta=theta,k=k,N0=N0,timestep=timestep,obs.data=obs.data,indices=indices)-prior(theta=theta,reduced=reduced,log=T)
  return(logpost)
}


##############################################################################################################
#  loglik_Indicator  evaluates the log likelihood sum(indicator(matching set defined by distance measure and tolerance level))
#input:
#                z                    - M replicate data sets
#                theta                - vector of parameters
#                dist.metric          - distance metric
#                eps                  - vector of tolerance levels for each summary statistics
#                obs.data             - observed data
#                indices              - indices of the summary statistics, if entire set of summary statistics used, 
#                                       then the likelihood is the target likelihood, otherwise if any subset of 
#                                       summary statistics is used, the likelihood is an approximation to the target one.
#output:         evaluated log likelihood
##############################################################################################################
loglik_Indicator=function(z,theta,dist.metric=abs,eps,obs.data,indices){

  if ( (theta[2]<0) | (theta[3]<0) | all(is.na(z)) ) {loglik_Indicator=(-1)*Inf
  }else{
    summary_y=summaryfun(obs.data)[indices]
    if (length(indices)==1){
      summaryz=applysummaries(z)[indices,]
      Indic=sapply(1: length(summaryz), function(x) {dist.metric(summaryz[x] -summary_y)<eps[indices]} )
      loglik_Indicator=log(length(which(Indic==T)))
    }else{
      summaryz=applysummaries(z)[indices,]
      Indic_more=lapply(1: ncol(summaryz), function(x) {dist.metric(summaryz[,x]-summary_y)<eps[indices]} )
      loglik_Indicator=log(sum(unlist(lapply(1:length(Indic_more), function(x) {if (all(Indic_more[[x]])){1}else{0} }))))
    }
  }
  
  return(loglik_Indicator)
}


##############################################################################################################
# logposterior_indicator  evaluates the posterior using the log likelihood sum(indicator(matching set defined by distance measure and tolerance level))
#input:
#                z                    - M replicate data sets
#                theta                - vector of parameters
#                dist.metric          - distance metric
#                eps                  - vector of tolerance levels for each summary statistics
#                obs.data             - observed data
#                indices              - indices of the summary statistics, if entire set of summary statistics used, 
#                                       then the likelihood is the target likelihood, otherwise if any subset of 
#                                       summary statistics is used, the likelihood is an approximation to the target one.
#output:         evaluated log posterior

##############################################################################################################
#for this function I will need distance measure, tolerance and simulated data
#to find the number of simulated data points that are within the set:
# A={z: all(dist(s(z)-s(y))<eps) } that is the likelihood
logposterior_indicator=function(theta,z,dist.metric,eps,obs.data,indices){
  loglik_Indicator=loglik_Indicator(z,theta,dist.metric=dist.metric,eps,obs.data,indices)
  logpost=loglik_Indicator +prior(theta,log=T)
  return(logpost)
}

#######################################################
#applysummaries: a wrapper function to apply summaryfun
#                to a list of datasets or a matrix of datasets, where
#                each column corresponds to a different data set
#input:       output -  list from the theta-ricker function or 
#                       a list of data vectors or a matrix of 
#                       datasets, where
#                       each column corresponds to a different data set 
#output: a list of vectors of summary statistics for each data set
#######################################################


applysummaries = function(output){
  if(is.list(output)){
    out = lapply(1:length(output), function(x) {summaryfun(output[[x]])} )
  }else{
    out = apply(output,2,summaryfun)
  }
  return(out)
}

#######################################################
#summaryfun   calculate summaries
#input:       output - list from the theta-ricker function or 
#                       a data vector
#output: a vector of summary statistics
#######################################################


summaryfun = function(output){
  
  if (is.list(output)){
    data=output$y
  }else{data=output}
  
  medi         = median(data,na.rm=T)
  meany        = mean(data,na.rm=T)
  meanover1    = mean(data[data>1],na.rm=T)
  if (is.nan(meanover1)){meanover1=0}
  
  morethan10   = sum(data[which(data>10)],na.rm=T)
  #number of non-zeros
  data[1]      = 0
  Count0       = length(which(data==0))
  q3           = quantile(data,p=.75,na.rm=T)
  maxy         = max(data,na.rm=T)
  CountOver100 = (length(which(data>100  ) )) 
  CountOver300 =(length(which(data>300  ))) 
  CountOver500 = (length(which(data>500 ))) 
 
  if (length(which(data>800))==0) {
  sumOver800 = 0
  }else{
  sumOver800  = sum(data[which(data>=800)],na.rm=T)
  }

  outputsum = c(medi, meany, meanover1, morethan10, Count0,  q3,  maxy,CountOver100,CountOver300,CountOver500,sumOver800)
  names(outputsum) = c("mediy","meany","meanover1","morethan10","Count0","q3","maxy",'CountOver100','CountOver300','CountOver500','sumOver800')
  return(outputsum)
}

