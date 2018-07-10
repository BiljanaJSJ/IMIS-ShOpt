#####################################################################################################
#IMIS_ShOpt_ABC   - implemented the IMIS-ShOpt-ABC Algorithm
#input:       B              - number of incremental points to be added to the importance sampling distribution
#             B.re           - J, the number of re-sampled point
#             number_T       - number of iterations
#             D              - number of intial points for the optimizers
#             Q              - number of optimization methods
#             obs.data       - observed data
#             M=30           - Nz - number of replicate data  
#             N0.true        - initial abundance in the population for the Ricker and Theta-Ricker model
#             timestep       - number of data points
#             dist.metric    - distance metric for the ABC
#             eps            - vector of tolerance levels
#             optim.fun1     - the first optimization criteria - negative log posterior using synthetic likelihood
#             optim.fun2     - the second optimization criteria - negative log synthetic likelihood
#output:      a list of 
#             stat           - statistics 
#             X_all          - importance sampling distribution before re-sampling, comprising particles of both models M1 and M2
#             Weights        - weights for each particle in X_all
#             sigma_all      - covariance function for each optimized point in the optimization stage
#                              and for each selected highest-weight point in the importnace stage    
#             sigma_all_lst
#             resample        - re-sampled particles from both models M1 and M2
#             resample_X_all  - resampled points divided in two lists, one for each model M1 and M2 
#             like_all        - likelihood evaluation for each point in  X_all
#             prior_all       - prior evaluation for each point in  X_all
#             pModel          - a vector of posterior probabilities of each model
#             center          - center points from optimization, and the highest-weigth points from the
#                               importance stage 
#             data            - data
#             sumStatsUsed_opt - keeps subset of randomly drawn 3 summary statistics
#####################################################################################################


IMIS_ShOpt_ABC <- function(B, B.re, number_T=1000, D,Q,obs.data=obs.data,M=30, N0.true=N0.true,timestep=timestep,dist.metric=abs,eps, optim.fun1,optim.fun2, other=NULL){

  cl=makeCluster(rep('localhost', 4))
  clusterExport(cl,varlist=ls(),envir = environment())
  clusterCall(cl,function(x) {source('RickerABC_functions.R');source('IMIS_ShOpt_ABC.R');library(numDeriv)})
  
  
  B0 = B*10
  

  X_all               =  sample.prior(B0)
  X_k                 =  X_all
  

  Sig2_global            =  cov(X_all)
  stat_all               =  matrix(0, 6, number_T)				# 6 diagnostic statistics at each iteration
  

  
  allcomb=combn(11,7)
  

 
  
  sumStatsUsed_opt=list()

#res=get(load('SaveOptimResults.RData'))
#X_all=res$X_all;X_k=res$X_k;Weights=res$Weights;center_all_lst=res$center_all_lst; gaussian_all=res$gaussian_all;
#prior_all=res$prior_all;like_all=res$like_all;like_scaling=res$like_scaling;sumStatsUsed_opt=res$sumStatsUsed_opt
#like_temp=res$like_temp;prior_temp=res$prior_temp;
#Sig2_global=res$Sig2_global;prior_all_lst=res$prior_all_lst;sigma_all_lst=res$sigma_all_lst
#which_remain=c()
  prior_all=like_temp=prior_temp=sigma_all_lst=center_all_lst=NULL

  for (k in 1:number_T){
    
    
    print(paste("k is:",k,sep=''))
   
    ptm.like = proc.time()
    
    
    if (!(is.null(X_k))){
    prior_temp  = unlist(c(prior_temp, sapply(1:nrow(X_k), function(x) {prior(X_k[x,])})))
    clusterExport(cl,varlist=ls(),envir = environment())
    #simulate M replicate data sets for the indicator likelihood function with eps and distance measure
    #z1 =parLapply(cl,1:nrow(X_k), function(x) {sapply(1:M,function(a) {theta_ricker(N0.true,X_k[x,],timestep)$y}) })
    z1 =parLapply(cl,1:nrow(X_k), function(x) {lapply(1:M,function(a) {theta_ricker(N0.true,X_k[x,],timestep)}) })

    clusterExport(cl,varlist=ls(),envir = environment())
    clusterCall(cl,function(x) {library(mvtnorm)})

    #calculate likelihood using the indicator likelihood function with eps and distance measure
    #like_temp = c(like_temp,unlist(parLapply(cl,1:nrow(X_k), function(x){loglik_Indicator(z1[[x]],X_k[x,],dist.metric=dist.metric,eps=eps,obs.data=obs.data,indices=1:11)})))
    like_temp = c(like_temp,unlist(parLapply(cl,1:nrow(X_k), function(x){(-1)*synthloglike(X_k[x,],simDs=z1[[x]], k=200,N0=N0.true,timestep=timestep,obs.data=obs.data,indices=1:11)})))
    }

    print('Non-zero Weights:')
    print(length(which(like_temp!=-Inf)))
    
    clusterExport(cl,varlist=ls(),envir = environment())

    like_all               = exp(like_temp +max(like_temp[which(like_temp!=-Inf)]) )
    prior_all              = exp(prior_temp+max(prior_temp[which(prior_temp!=-Inf)]))
    
  
    
    ptm.use = (proc.time() - ptm.like)[3]
    if (k==1){
      print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
      envelop_all = prior_all			# envelop stores the sampling densities
    }else{
      Nk=dim(X_all)[1]
      envelop_all= apply( rbind(prior_all*B0/Nk, gaussian_all*B/Nk), 2, sum)
    }
    print('before weigths')
  
    Weights = prior_all*like_all/envelop_all
    Weights[which(is.na(Weights) | is.nan(Weights) | (Weights %in% c(-Inf,Inf) )) ]=0
    Weights = Weights / sum(Weights,na.rm=TRUE)	

    
    print(Weights)
    print(log(mean(Weights,na.rm=TRUE))		)
    stat_all[1,k] = log(mean(Weights,na.rm=TRUE))			# the raw marginal likelihood
    stat_all[2,k] = sum(1-(1-Weights)^B.re,na.rm=TRUE)		# the expected number of unique points
    stat_all[3,k] = max(Weights,na.rm=TRUE)				# the maximum weight
    stat_all[4,k] = 1/sum(Weights^2,na.rm=TRUE)			# the effictive sample size
    stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))	# the entropy relative to uniform
    stat_all[6,k] = var(Weights/mean(Weights,na.rm=TRUE),na.rm=TRUE)	# the variance of scaled weights
    print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
    print(unlist(c(t, round(stat_all[1:4,k], 3))))
   
    if (k==1){
      print("Here 1")
      
      which_exclude  =label_weight=which_remain=important=size_remain=list() 
     
      X_k =which_exclude=label_weight=which_remain=important=size_remain=X_imp=NULL
      
    
      label_weight = sort(Weights, decreasing = TRUE, index=TRUE)
      which_remain = unlist(which(Weights>label_weight$x[min(B0,length(label_weight$x))]) )
      size_remain = length(which_remain)
      
      ############################################ WHILE LOOP OVER THE i ###############	
      d=1;i=1

  while(length(which_remain)>1 && i<=D/3){
        print("herei")
        print(i)
        print(paste("i=",i))


        important = which_remain[which(Weights[which_remain]==max(Weights[which_remain]))]
        if (length(important)>1)	{
           important = sample(important,1)
          }
        X_imp=X_all[important,]
          

          # Remove the selected input from candidates
          which_exclude = unlist(union( which_exclude, important ))
          which_remain  = setdiff(which_remain, which_exclude)
      
          	
          print("optimization starting from a value of: ")
          print(paste(round(X_imp[1],2),',',round(X_imp[2],2),',',round(X_imp[3],2),',',round(X_imp[4],2),sep=''))
          print("here3--matrix--")
          
          #1. run the optimizer
          res=pars=hess=list()
         
         
         
          indSamples=sample(1:ncol(allcomb),Q)
          for (ind in indSamples){ 
            print(ind)
            if (size_remain>1){
            
            sumStatsUsed_opt[[d]]=allcomb[,ind]
           
            
            
            
            res[[d]]     = optimx(X_imp, optim.fun2,k=200, N0=N0.true, timestep=timestep,obs.data=obs.data,indices=sumStatsUsed_opt[[d]], method="Nelder-Mead",control=list( maxit=100))
            #res[[d]]     = optimx(X_imp[1:4], optim.fun2,k=200, N0=N0.true, timestep=timestep,obs.data=obs.data,indices=1:10, method="Nelder-Mead",control=list( maxit=100))
            pars[[d]]    = c(res[[d]][[1]],res[[d]][[2]],res[[d]][[3]],res[[d]][[4]])
            hess         = hessian(optim.fun1,pars[[d]],k=200, N0=N0.true, timestep=timestep,obs.data=obs.data,indices=1:11)
           
            if (any(is.nan(hess)))   hess[is.nan(hess)]=0

            print(length(hess))
           
            
            if (!(any(is.na(hess)) | any(is.nan(hess)))){
            if (!(is.singular.matrix(hess))){
 
             sigma_all_lst[[d]]  = tryCatch({
                 solve(hess)
               }, warning = function(w) {
                msg('Singular matrix')
              }, error = function(e) {
                 Sig2_global
               }
               )

             # sigma_all_lst[[d]] = solve(hess)
            }else{
              sigma_all_lst[[d]]=Sig2_global 
            }
            
              if (min(eigen(sigma_all_lst[[d]])$values) <=0) {
                warning('Hessian is not positive definite!')             
                eigen.values = eigen(sigma_all_lst[[d]])$values
                eigen.values[which(eigen.values < 0)] = 0
                sigma_all_lst[[d]]=sigma_all_lst[[d]]+abs( min(eigen(sigma_all_lst[[d]])$values))*diag(dim(sigma_all_lst[[d]])[1])*1.01
              }
              
            
     
              if (any(is.complex(eigen(sigma_all_lst[[d]])$values))){
                sigma_all_lst[[d]]=Sig2_global
              }
            
              # If the hessian matrix is still not positive definite, we define the covariance as following
            if (min(eigen(sigma_all_lst[[d]])$values)<=0){			
              eig = eigen(sigma_all_lst[[d]])
              eigen.values = eig$values
              eigen.values[which(eigen.values<0)] = 0
              hess_new1 = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
              if (!(is.singular.matrix(hess_new1 + diag(1/diag(Sig2_global))))){
            
                sigma_all_lst[[d]]  = tryCatch({
                  solve(hess_new1 + diag(1/diag(Sig2_global)) )
                }, warning = function(w) {
                  msg('Singular matrix')
                }, error = function(e) {
                  Sig2_global
                }
                )
           
             }
            }
           }
            
            
            # }
            # }
           
              
              X_k                 = rbind(X_k, rmvnorm(B, pars[[d]], sigma_all_lst[[d]]))
              center_all_lst      = rbind(center_all_lst, pars[[d]])	
          
            
            
              
              
              if (length(which_remain)==1) {break}  
              
              
              distance_remain = apply(abs(X_all[which_remain,]-pars[[d]]),1,sum)
              label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
             
              
              
              which_exclude = unlist(union( which_exclude, which_remain[label_dist$ix[(1:floor(size_remain/(Q*D)))]]))
              which_remain = setdiff(which_remain, which_exclude)
              size_remain=length(which_remain)
              
            }
        
            d=d+1
            if (length(which_remain)==1) {break}  
          
          }
          i=i+1
          if (length(which_remain)==1) {break}  
          }
          print("here4----")
          
        
          if  (!(is.null(X_k))) {X_all = rbind(X_all, X_k)}
         
          
        
         print("Here 3")
        

       
        

         if (!(is.null(center_all_lst))){
             gaussian_all = matrix(NA, dim(center_all_lst)[1], dim(X_all)[1])
             for (index in 1:dim(center_all_lst)[1]){
               print(index)
                     
               gaussian_all[index,] = dmvnorm(X_all[,1:ncol(center_all_lst)], center_all_lst[index,], sigma_all_lst[[index]])
             }
          
           }
         save.image("Temporary2.results.Rdata")  
         out_ls=list(X_all=X_all,X_k=X_k,k=k,Weights=Weights,center_all_lst=center_all_lst,gaussian_all=gaussian_all,prior_all=prior_all,like_all=like_all,sumStatsUsed_opt=sumStatsUsed_opt, Sig2_global=Sig2_global,prior_temp=prior_temp,sigma_all_lst=sigma_all_lst, which_remain= which_remain,envelop_all=envelop_all)
         save(out_ls,file=paste('SaveOptimResults',k,'_FixGaussian.RData',sep=''))   
            
  }
  
      
      
      ######################################################   END OF i LOOP   ############################
      
     
      
    
    ######## 
    if (k>1 ){
      print(paste("k>1, k=",k))
     
      X_k      =  NULL
      
      important= which(Weights == max(Weights,na.rm=TRUE))
      if (length(important)>1)	important = sample(important,1)
      X_imp= X_all[important,]
      
      print('Max. weight point')
      print(X_imp)

      
     
      
      center_all_lst                = rbind(center_all_lst, X_imp)
      
      distance_all                  = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )
      label_nr                      = sort(distance_all, decreasing = FALSE, index=TRUE)		
      which_var                     = label_nr$ix[1:B]								
      W1                            = Weights
      W1[is.na(W1)]                 = 0
      
      Sig2                          = cov.wt(X_all[which_var,], wt = W1[which_var]+1/length(W1), cor = FALSE, center = X_imp, method = "unbias")$cov
      
      sigma_all_lst[[length(sigma_all_lst)+1]] = Sig2
      X_k                           = rmvnorm(B, X_imp, Sig2)				# Draw new samples
      
   
      
      X_all                         = rbind(X_all, X_k )
      print("here it is")
      
      
     
      

      
      gaussian_new = matrix(0, length(sigma_all_lst), dim(X_all)[1])
      gaussian_new[1:(dim(gaussian_all)[1]), 1:(dim(gaussian_all)[2])] = gaussian_all
      
     
   

      gaussian_new[length(sigma_all_lst), ] = dmvnorm(X_all, X_imp,   sigma_all_lst[[length(sigma_all_lst)]])
      
      for (j in 1:(dim(gaussian_new)[1]-1))	gaussian_new[j, (dim(X_all)[1]-dim(X_k)[1]+1):dim(X_all)[1] ] = dmvnorm(X_k, center_all_lst[j,], sigma_all_lst[[j]])
      
      gaussian_all = gaussian_new
      

      out_ls=list(stat_all=stat_all,X_all=X_all,X_k=X_k,k=k,Weights=Weights,center_all_lst=center_all_lst,gaussian_all=gaussian_all,prior_all=prior_all,like_all=like_all,sumStatsUsed_opt=sumStatsUsed_opt, like_temp=like_temp,Sig2_global=Sig2_global,prior_temp=prior_temp,sigma_all_lst=sigma_all_lst, which_remain= which_remain,envelop_all=envelop_all)    
     
      #save(out_ls,file=paste('SaveOptimResults',k,'.RData',sep=''))
      save(out_ls,file=paste('SaveOptimResults.RData',sep=''))

     # save.image("Temporary3.results.Rdata")
      
      if (stat_all[2,k] > (1-exp(-1))*B.re)	break	
      
    }
  }
############################################## end of k loop  ############################################## 
  
  
  nonzero = which(Weights>0)
 
  which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  
  resample_X=X_all[which_X,]

 
  saveLst=list(stat=t(stat_all),X_all=X_all,Weights=Weights,resample=resample_X,center_all_lst=center_all_lst,gaussian_all=gaussian_all,prior_all=prior_all,like_all=like_all,sumStatsUsed_opt=sumStatsUsed_opt, like_temp=like_temp,Sig2_global=Sig2_global,prior_temp=prior_temp,sigma_all_lst=sigma_all_lst, which_remain= which_remain,envelop_all=envelop_all)
  
 
  save(saveLst,file='samples.RData')
  stopCluster(cl)
  return(saveLst)
  
  
} # end of IMIS-ShOpt-ABC


  
  
  
