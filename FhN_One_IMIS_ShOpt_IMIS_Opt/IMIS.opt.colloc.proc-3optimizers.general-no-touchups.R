#########################################################################################################
#IMIS.opt.colloc.3optimizers.general.no.touch.ups function to run the IMIS-ShOpt algorithm
#input:  
#            B                      - number of incremental particles
#            B.re                   - number of re-sampled particles
#            number_k               - number of iterations
#            D                      - number of initial points for the optimizaer to start
#            parnames               - a vector of parameter names
#            data                   - data
#            plots                  - diagnostic plots of likelihood, prior, and sampling Weights, 
#                                     0 if one does not want to generate diagnostic plots
#            ncoefs                 - number of coefficients for the
#            optim.fun1             - optimization criteria 1 for the Shotgun optimization
#            optim.fun2             - optimization criteria 2 for the Shotgun optimization
#            optim.fun3             - optimization criteria 3 for the Shotgun optimization
#            other                  - other parameters
#output:     a list of 
#            stat                   - statistics like number of unique particles, marginal likelihood..
#            X_all                  - the importance sampling distribution 
#            Weights                - weigths for each point in the importance sampling ditribution
#            sigma_all              - the list of matrices for the variances at the modes (optimization stage)
#                                     or highest weigths points (in the importance stage)
#            resample               - resampled points
#            center                 - modes (optimization stage) and highest weigths points (in the importance stage)
#            optimres               - optimization results from the three optimizers
#            data                   - data
#########################################################################################################
IMIS.opt.colloc.3optimizers.general.no.touch.ups <- function(B, B.re, number_k, D,parnames = c("c"),data=data, plots=3, ncoefs=0, optim.fun1, optim.fun2, optim.fun3,other=NULL,logging=TRUE){
	B0 = B*10
	#clusterExport(cl,varlist=ls())
	X_all = X_k = sample.prior(B0)				# Draw initial samples from the prior distribution
	if(is.matrix(X_all)){
     parend=dim(X_all)[2]
    }else{parend=1
    ncoefs=0}
	like_temp = c()
    like_scaling = 0
	if (is.vector(X_all))	Sig2_global = var(X_all)	# the prior covariance
	if (is.matrix(X_all))	Sig2_global = cov(X_all)	# the prior covariance
	stat_all = matrix(0, 6, number_k)				# 6 diagnostic statistics at each iteration
	center_all = prior_all = like_all  = NULL			# centers of Gaussian components, prior densities, and likelihoods
	sigma_all = list()						# covariance matrices of Gaussian components
	if (D>=1)	option.opt = 1					# use optimizer
	if (D==0){	option.opt = 0; D=1	}			# NOT use optimizer
	
    if(!exists("neglogpriorx0", mode="function")){
        neglogpriorx0<-function(...){return(0)}
       # clusterExport(cl,varlist=ls())
    }
    posterior = function(pars,data=data){sum(	neglogprior(pars[1:(length(pars)-ncoefs)])) + sum(	neglogpriorx0(pars[(length(pars)-ncoefs+1):length(pars)]))-likelihood(pars,logs=TRUE,data=data) } 
    #clusterExport(cl,varlist=ls())
    if(plots!=0){	par(mfrow=c(4,3))}
############################################## start of k loop  ############################################## 
    #k=1
	  #while(k <=number_k & (stat_all[2,max(1,k-1)] < (1-exp(-1))*B.re)){        
	  for (k in 1:number_k){
	  print("herek")
		print(k)
		ptm.like = proc.time()
		
		#prior_all = c(prior_all, prior(X_k))		# Calculate the prior densities
		
    # In pathological cases the likelihood may need rescaling by a constant. This is to deal with the case where the parameters are really bad.
    # Note that this doesnt change the relative weights which is important.
		
    #like_temp = c(like_temp+like_scaling,likelihood(X_k,logs=TRUE,data=data))  #undo the scaling part (two lines down)
    #like_scaling = max(like_temp,na.rm=TRUE)
    #like_all = exp(like_temp-like_scaling)		# Calculate the likelihoods while removing a constant scaling part.  The like_scaling numerically stabilizes the likelihood and since the impact thereof is applied to all entries, the impact is cancelled out when computing the Weights/sum(Weights) step.
	#	if(plots!=0){
	#		if(is.matrix(X_all)){plot(X_all[,plots],like_all,main="Like")
	#			plot(X_all[,plots],prior_all,main="prior")
	#		}else{
	#			plot(X_all,like_all,main="Like")
	#				plot(X_all,prior_all,main="prior")
	#			}
	#		}
	#	ptm.use = (proc.time() - ptm.like)[3]
	#	if (k==1){
	#		print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
	#		envelop_all = prior_all			# envelop stores the sampling densities
	#	}else{
     	#envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
	#	}

	#	Weights = prior_all*like_all / envelop_all	# importance weight is determined by the posterior density divided by the sampling density
	         prior_all = c(prior_all, prior(X_k,logs=logging))		# Calculate the prior densities
        like_all = c(like_all, likelihood(X_k,logs=logging,data=data))		# Calculate the likelihoods
		
      ptm.use = (proc.time() - ptm.like)[3]


    if (k==1){
      print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
      # envelop_all = prior_all			# envelop stores the sampling densities
      #}else{
      #envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
    }

                  if (logging){
			
			if (k==1)	envelop_all = prior_all			# envelop stores the sampling densities
			if (k>1)	envelop_all = log(apply( rbind(exp(prior_all)*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2)))
			Weights = exp(prior_all+like_all-envelop_all-max(like_all[which(!like_all==0)],na.rm=T))	# importance weight is determined by the posterior density divided by the sampling density
			
		}else{
			if (k==1)	envelop_all = prior_all			# envelop stores the sampling densities
			if (k>1)	envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
			Weights = prior_all*like_all / envelop_all	# importance weight is determined by the posterior density divided by the sampling density
		}

   


   	Weights[which(is.na(Weights) | is.nan(Weights) | (Weights %in% c(-Inf,Inf) )) ]=0
		Weights = Weights / sum(Weights,na.rm=TRUE)	
		stat_all[1,k] = log(mean(Weights,na.rm=TRUE))			# the raw marginal likelihood
		stat_all[2,k] = sum(1-(1-Weights)^B.re,na.rm=TRUE)		# the expected number of unique points
		stat_all[3,k] = max(Weights,na.rm=TRUE)				# the maximum weight
		stat_all[4,k] = 1/sum(Weights^2,na.rm=TRUE)			# the effictive sample size
		stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))	# the entropy relative to uniform
		stat_all[6,k] = var(Weights/mean(Weights,na.rm=TRUE),na.rm=TRUE)	# the variance of scaled weights
		print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
		print(c(k, round(stat_all[1:4,k], 3)))
		if(plots!=0){	
			if(is.matrix(X_all)){	plot(X_all[,plots],Weights,main="Weights")
		}else{	
			plot(X_all,Weights,main="Weights")}
		}
		if (k==1){
    print("Here 1")
# set up FDA smoothing stuff  
            if(is.matrix(X_all)){
			if (is.matrix(X_all[which(like_all>min(like_all)),]))	Sig2_global = cov(X_all[which(like_all>min(like_all,na.rm=TRUE)),])
            }else{Sig2_global=1}
			X_k = which_exclude = NULL					# exclude the neighborhood of the local optima 
			label_weight = sort(Weights, decreasing = TRUE, index=TRUE)
			which_remain = which(Weights>label_weight$x[min(B0,length(label_weight$x))]) 	# the candidate inputs for the starting points
			size_remain = length(which_remain)
# print("which_remain")
# print(which_remain)
			
			i=1
############################################ WHILE LOOP OVER THE i ###############	
 while(length(which_remain)>0 && i<=D){
		#	for (i in 1:D){
				print("herei")
				print(i)
				print(paste("i=",i))
				important = which_remain[which(Weights[which_remain]==max(Weights[which_remain]))]
				if (length(important)>1)	important = sample(important,1)	
				if (is.vector(X_all))	X_imp = X_all[important]
				if (is.matrix(X_all))	X_imp = X_all[important,]
# Remove the selected input from candidates
				which_exclude = union( which_exclude, important )
				which_remain = setdiff(which_remain, which_exclude)
				
				if (is.vector(X_all)){	
				 	print(paste("optimization starting from a value of ",X_imp))

                    
# run 3 optimizers,  
				
                    res = foreach(index=1:3) %dopar%{
                    

                        if(index == 1){
#1. old optimizer
                                res1<-optim.fun1(pars=X_imp,data=data,times=times,other=list(method='BFGS'))
                                r1 = optim(res1, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=0))
                                while(r1$hessian<0){
                                     r1 = optim(r1$par, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=25))
                                }
                                return(r1)
                       }else if(index==2){
                            
                                res2<-optim.fun2(pars=X_imp,data=data,times=times,other) 
                                r2 = optim(res2, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=1000))
                                while(r2$hessian<0){
                                     r2 = optim(r2$par, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=1000))
                                }
                               return(r2)

                      }else{
                               res3<-optim.fun3(pars=X_imp,data=data,times=times,other)
                               r3 = optim(res3, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=1000))
                               while(r3$hessian<0){
                                    r3 = optim(r3$par, posterior, data=data,method="BFGS", hessian=TRUE, control=list( maxit=25))
                                }
                               return(r3)
                     }
                    }
                    optimizer  = res[[1]]
                    optimizer2 = res[[2]]
                    optimizer3 = res[[3]]
                        
                            #, control=list(parscale=sqrt(Sig2_global)/10,maxit=5000))
					
					print(paste("standard optimizer " ,optimizer$par))	
#2. new optimizer				
# We have to have parameter names.  That is the way they are dealt with in CollocInfer	
		            print(paste("next step smoothing based optimizer  " ,optimizer2$par))	    
                    print(paste("2 step based optimizer  " ,optimizer3$par))	    

                        
                        
                        
                        
#   Estimated posterior covariance using the original optimizer to ensure the correct target posterior
					
# The center and variance estimates for teh Gaussian approximations based on the output from the optimizers.
					center_all = rbind(center_all, optimizer$par, optimizer2$par, optimizer3$par)
					ls1=length(sigma_all)+1
					ls2=length(sigma_all)+2
					ls3=length(sigma_all)+3                                        
                    sigma_all[[ls1]] = solve(optimizer$hessian)
					          sigma_all[[ls2]] = solve(optimizer2$hessian)
                    sigma_all[[ls3]]     = solve(optimizer3$hessian)
					
                        
					X_k = c(X_k, rnorm(B, center_all[ls1], sqrt(sigma_all[[ls1]])) )	                        # Draw new samples based on old optimizer
          X_k = c(X_k, rnorm(B, center_all[ls2], sqrt(sigma_all[[ls2]])) )	                        # Draw new samples based on smoother optimizer
					X_k = c(X_k, rnorm(B, center_all[ls3], sqrt(sigma_all[[ls3]])) )	                                                # Draw new samples based on 2-stage optimizer
                        

					
				}# end of vector operations
				if (is.matrix(X_all)){	
					print("optimization starting from a value of ")
					print(X_imp)
					print("here3--matrix--")

#1. old optimizer
                    
                        
                        # run 3 optimizers,  
                        ###########still need res 1 but it keeps failing############
                        res = foreach(index=1:3) %dopar%{
                        
                            if(index == 1){
                            
                                #1. old optimizer
                                res1<-optim.fun1(pars=X_imp,data=data,times=times,other)
                                return(res1)
                          }else if(index==2){
                                
                                res2<-optim.fun2(pars=X_imp,data=data,times=times,other)
                                return(res2)

                            }else{
                                res3<-optim.fun3(pars=X_imp,data=data,times=times,other)
                              
                                

                                return(res3)
                                
                            }
                      
                            
                        }
                    
                    
                    r1=res[[1]]
                    r2=res[[2]]
                    r3=res[[3]]
                    
                    optimizer   = list(par=r1,hessian=hessian(posterior,r1,data=data),value = posterior(r1,data=data))
                    optimizer2  = list(par=r2,hessian=hessian(posterior,r2,data=data),value = posterior(r1,data=data))
                    optimizer3  = list(par=r3,hessian=hessian(posterior,r3,data=data),value = posterior(r1,data=data))
                    
					names(optimizer$par)=names(optimizer2$par)=names(optimizer3$par)=parnames
                    
					print(paste("maximum posterior=", round(-optimizer$value,2),
                    ", likelihood=", round(likelihood(optimizer$par,log=T,data=data),2),
                    ", prior=", round(log(prior(optimizer$par[1:(length(pars)-ncoefs)])),2),
                    ", time used=", round(ptm.use/60,2),
                    "minutes, convergence=", optimizer$convergence))
                    
                    
                    
                    
                    
                    
                    
                    center_all = rbind(center_all, optimizer$par, optimizer2$par, optimizer3$par)					# the center of new samples
                    
					
                    #2. new optimizer		(the new smoothing based optimizer)
					
					
                    # first prep the parameters by naming them
					
                    
					
                    
					
                    
					print("standard optimizer " )
					print(optimizer$par)
					print(optimizer$value)
                    
					print("hybrid GPE step2 ")
					print(optimizer2$par)
					print(optimizer2$value)
					
					print("hybrid 2stage step2 ")
					print(optimizer3$par)
					print(optimizer3$value)
					ls1=length(sigma_all)+1
					ls2=length(sigma_all)+2
					ls3=length(sigma_all)+3
					
					if (min(eigen(optimizer$hessian)$values)>0)    	sigma_all[[ls1]] = solve(optimizer$hessian)						# the covariance of new samples
                    if (min(eigen(optimizer2$hessian)$values)>0)    sigma_all[[ls2]] = solve(optimizer2$hessian)												# the covariance of new samples using the smoothing based optimizer
                    if (min(eigen(optimizer3$hessian)$values)>0)    sigma_all[[ls3]] = solve(optimizer3$hessian)												# the covariance of new samples using the smoothing based optimizer
                    
					if (min(eigen(optimizer$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
						eig = eigen(optimizer$hessian)
						eigen.values = eig$values
						eigen.values[which(eigen.values<0)] = 0
						hessian = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
						sigma_all[[ls1]] = solve(hessian + diag(1/diag(Sig2_global)) )
					}
					
					
                    if (min(eigen(optimizer2$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
                        eig = eigen(optimizer2$hessian)
                        eigen.values = eig$values
                        eigen.values[which(eigen.values<0)] = 0
                        hessian = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
                        sigma_all[[ls2]] = solve(hessian + diag(1/diag(Sig2_global)) )
                    }
                    
                    
                    
                    if (min(eigen(optimizer3$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
                        eig = eigen(optimizer2$hessian)
                        eigen.values = eig$values
                        eigen.values[which(eigen.values<0)] = 0
                        hessian = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
                        sigma_all[[ls3]] = solve(hessian + diag(1/diag(Sig2_global)) )
                    }
                    
                    
					
					
                    ############# NOTE original code used								X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[1]]) )
                    ############ meaning that the variance term was fixed even though it should have been changing.  I believe that this is an error in the original code
					
					X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[ls1]]) )			# Draw new samples
                    # print("X_k from optimizer")
                    # print(X_k)
                    ############# END OF NOTE
					
					X_k = rbind(X_k, rmvnorm(B, optimizer2$par, sigma_all[[ls2]]) )			# Draw new samples using the smoothing based optimizer.
                    
                    
                    X_k = rbind(X_k, rmvnorm(B, optimizer3$par, sigma_all[[ls3]]) )			# Draw new samples using the smoothing based optimizer.
                    
                    
                    # print("X_k from smoothing based")
                    # print(X_k)
					
					
					print("here4----")
                }#end of matrix
     
     
     if(i==1){opt123=c(optimizer$par,optimizer2$par,optimizer3$par)
     }else{
         opt123=rbind(opt123,optimizer$par,optimizer2$par,optimizer3$par)
         
     }
     
     # exclude the neighborhood of the local optima
     # first do this with the original optimizer
     if (is.matrix(X_all))	distance_remain = apply(abs(X_all[which_remain,]-optimizer$par),1,sum)
     if (is.vector(X_all))	distance_remain = abs(X_all[which_remain]-optimizer$par)
     label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
     which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
     which_remain = setdiff(which_remain, which_exclude)
     
     
     if(length(which_remain)>0){

         # next do this with the smoother optimizer
         if (is.matrix(X_all))	distance_remain_alt = apply(abs(X_all[which_remain,]-optimizer2$par),1,sum)
         if (is.vector(X_all))	distance_remain_alt = abs(X_all[which_remain]-optimizer2$par)
         label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
         which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
         which_remain = setdiff(which_remain, which_exclude)
     }
     
     

     
     save.image("Temporary2.results.Rdata")
     

     i=i+1
                    

			}  # end of while "i" loop

 ######################################################   END OF i LOOP   ############################
 
 print("Here 3")
			if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
			if (is.vector(X_all))	X_all = c(X_all, X_k)		
		}
#print(k)

######## 
		if (k>1 | option.opt==0){
			print(paste("k>1, k=",k))
			important = which(Weights == max(Weights,na.rm=TRUE))
			if (length(important)>1)	important = important[1]
			if (is.matrix(X_all))	X_imp = X_all[important,]				# X_imp is the maximum weight input
			if (is.vector(X_all))	X_imp = X_all[important]

			if (is.matrix(X_all))	center_all = rbind(center_all, X_imp)
			if (is.vector(X_all))	center_all = c(center_all, X_imp)

			if (is.matrix(X_all))	distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )
			if (is.vector(X_all))	distance_all = abs(X_all-X_imp)			# Calculate the distances to X_imp

			
			label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)		# Sort the distances
			which_var = label_nr$ix[1:B]								# Pick B inputs for covariance calculation
      W1=Weights
      W1[is.na(Weights)]=0
			if (is.matrix(X_all))	Sig2 = cov.wt(X_all[which_var,], wt = W1[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
			if (is.vector(X_all)){
				Weights_var = W1[which_var]+1/length(X_all)
				Weights_var = Weights_var/sum(Weights_var,na.rm=T)
				Sig2 = (X_all[which_var]-X_imp)^2 %*% Weights_var
			}

			sigma_all[[length(sigma_all)+1]] = Sig2
			if (is.matrix(X_all))	X_k = rmvnorm(B, X_imp, Sig2)				# Draw new samples
			if (is.vector(X_all))	X_k = rnorm(B, X_imp, sqrt(Sig2))			# Draw new samples
			if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
			if (is.vector(X_all))	X_all = c(X_all, X_k)
		}
		
		
		if (k==1){
			if (is.matrix(X_all))	{
				gaussian_all = matrix(NA, dim(center_all)[1], dim(X_all)[1])
				for (index in 1:dim(center_all)[1]){
					gaussian_all[index,] = dmvnorm(X_all, center_all[index,], sigma_all[[index]])
				}
			}
			if (is.vector(X_all))	{  
				gaussian_all = matrix(NA, length(center_all), B0+length(center_all)*B)
				for (index in 1:length(center_all)){
					gaussian_all[index,] = dnorm(X_all, center_all[index], sqrt(sigma_all[[index]]))
				}
			}
		}
		
		out_ls=list(X_all=X_all,X_k=X_k,Weights=Weights,center_all=center_all,gaussian_all=gaussian_all,prior_all=prior_all,like_temp=like_temp,like_scaling=like_scaling)
		save(out_ls,file='SaveOptimResults.RData')	
# print(paste("X_all",X_all))
 print("here it is")
		
		if (k>1){
			if (is.matrix(X_all)){
				gaussian_new = matrix(0, length(sigma_all), dim(X_all)[1] )
				gaussian_new[1:(dim(gaussian_all)[1]), 1:(dim(gaussian_all)[2])] = gaussian_all
				gaussian_new[length(sigma_all), ] = dmvnorm(X_all, X_imp, sigma_all[[length(sigma_all)]])
				for (j in 1:(length(sigma_all)-1))	gaussian_new[j, (dim(X_all)[1]-dim(X_k)[1]+1):dim(X_all)[1] ] = dmvnorm(X_k, center_all[j,], sigma_all[[j]])
			}
			if (is.vector(X_all)){
				gaussian_new = matrix(0, length(sigma_all)+1, length(X_all) )
				gaussian_new[1:(dim(gaussian_all)[1]), 1:(dim(gaussian_all)[2])] = gaussian_all
				gaussian_new[length(sigma_all)+1, ] = dnorm(X_all, X_imp, sqrt(sigma_all[[length(sigma_all)]]))
				for (j in 1:length(sigma_all))	gaussian_new[j, (length(X_all)-length(X_k)+1):length(X_all) ] = dnorm(X_k, center_all[j], sqrt(sigma_all[[j]]))
			}
			gaussian_all = gaussian_new
		}
		if (stat_all[2,k] > (1-exp(-1))*B.re)	break	
		save.image("Temporary2.results.Rdata")
  #      k=k+1
	}
    # end of k
    
############################################## end of k loop  ############################################## 

	
	nonzero = which(Weights>0)
	which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
	if (is.matrix(X_all))	resample_X = X_all[which_X,]
	if (is.vector(X_all))	resample_X = X_all[which_X]
	options("warn"=0)
	return(list(stat=t(stat_all),  X_all=X_all,Weights=Weights,sigma_all=sigma_all,resample=resample_X, center=center_all, optimres=opt123,data=data))
	
	
} # end of IMIS


