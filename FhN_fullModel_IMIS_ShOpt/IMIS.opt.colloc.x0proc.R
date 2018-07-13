#library(mvtnorm)
#library(CollocInfer)
IMIS.opt.colloc <- function(B, B.re, number_k, D,parnames = c("a","b","c"),lik=NULL,proc=NULL,coefs=NULL,datamatrix=data,plots=3,ncoefs=2,likelihood.w.fit=NULL){
# NOTE that at this time the function makes plots to show the likelihood, prior, and sampling Weights... 
# Just to make sure that it is all working properly.	
# plots provides some diagnostic plots, but could hog memory and/or break the function if Infinite values appear in weights
#  you need to have defined the following:
# likelihood
# prior
# lik
# proc
# control.in
# control.out
# knots
# norder
# nbasis
# range
# bbasis 
# fd.data
# DEfd
# coefs	
	
	options("warn"=-2)
	B0 = B*10
	X_all = X_k = sample.prior(B0)				# Draw initial samples from the prior distribution
	
parend=dim(X_all)[2]
	
	if (is.vector(X_all))	Sig2_global = var(X_all)	# the prior covariance
	if (is.matrix(X_all))	Sig2_global = cov(X_all)	# the prior covariance
	stat_all = matrix(NA, 6, number_k)				# 6 diagnostic statistics at each iteration
	center_all = prior_all = like_all = like_scaling_old  = NULL			# centers of Gaussian components, prior densities, and likelihoods
	sigma_all = list()						# covariance matrices of Gaussian components
	if (D>=1)	option.opt = 1					# use optimizer
	if (D==0){	option.opt = 0; D=1	}			# NOT use optimizer
	posterior = function(theta){sum(	neglogprior(theta))+sum(neglogpriorx0(theta))-likelihood(theta,logs=TRUE) } 
	if(plots!=0){	par(mfrow=c(4,3))}

	for (k in 1:number_k ){
	print("herek")
		print(k)
		ptm.like = proc.time()
		
		prior_all = c(prior_all, prior(X_k))		# Calculate the prior densities
		
# In pathological cases the likelihood may need rescaling by a constant. This is to deal with the case where the parameters are really bad.
# Note that this doesnt change the relative weights which is important.
#TODO:########################### set up another function likelihood.w.fit that outputs not only the likelihood value but also the data fit 
		
		if(!is.null(likelihood.w.fit)){like.and.fit.all = likelihood(X_all,logs=TRUE)
			like_temp=like.and.fit$like
			fits_temp=like.and.fit$fits   #starting points for the collocation method data fit.
		}else{like_temp = likelihood(X_all,logs=TRUE)
			fits_temp=lik$more$bvals%*%coefs   #starting points for the collocation method data fit.
		}
		like_scaling = max(like_temp)
		like_all = exp(like_temp-like_scaling)		# Calculate the likelihoods
		if(plots!=0){
			if(is.matrix(X_all)){plot(X_all[,plots],like_all,main="Like")
				plot(X_all[,plots],prior_all,main="prior")
			}else{
				plot(X_all,like_all,main="Like")
					plot(X_all,prior_all,main="prior")
				}
			}
		ptm.use = (proc.time() - ptm.like)[3]
		if (k==1){
			print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
			envelop_all = prior_all			# envelop stores the sampling densities
		}else{
#	print(dim(gaussian_all))
			envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
		}
#	if(is.matrix(X_all)){		plot(X_all[,dim(X_all)[2]],envelop_all,main="envelop")
#	}else{  		plot(X_all,envelop_all,main="envelop")}
		
		Weights = prior_all*like_all / envelop_all	# importance weight is determined by the posterior density divided by the sampling density
		stat_all[1,k] = log(mean(Weights))			# the raw marginal likelihood
		Weights = Weights / sum(Weights)			
		stat_all[2,k] = sum(1-(1-Weights)^B.re)		# the expected number of unique points
		stat_all[3,k] = max(Weights)				# the maximum weight
		stat_all[4,k] = 1/sum(Weights^2)			# the effictive sample size
		stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))	# the entropy relative to uniform
		stat_all[6,k] = var(Weights/mean(Weights))	# the variance of scaled weights
		if (k==1)	print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
		print(c(k, round(stat_all[1:4,k], 3)))
		if(plots!=0){	
			if(is.matrix(X_all)){	plot(X_all[,plots],Weights,main="Weights")
		}else{	
			plot(X_all,Weights,main="Weights")}
		}
		if (k==1){
 print("Here 1")
# set up FDA smoothing stuff  
				
			if (is.matrix(X_all[which(like_all>min(like_all)),]))	Sig2_global = cov(X_all[which(like_all>min(like_all)),])
			X_k = which_exclude = NULL					# exclude the neighborhood of the local optima 
			label_weight = sort(Weights, decreasing = TRUE, index=TRUE)
			which_remain = which(Weights>label_weight$x[B0]) 	# the candidate inputs for the starting points
			size_remain = length(which_remain)
# print("which_remain")
# print(which_remain)
			
			i=1
			while(length(which_remain)>0 && i<=D){
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
# run both optimizers,  maybe turn this off later or permit a choice of optimizer
					
#1. old optimizer
					optimizer = optim(X_imp, posterior, method="BFGS", hessian=TRUE)#, control=list(parscale=sqrt(Sig2_global)/10,maxit=5000))
					
					print(paste("standard optimizer " ,optimizer$par))	
#2. new optimizer				
# We have to have parameter names.  That is the way they are dealt with in CollocInfer	
					pars = X_imp
					names(pars) = proc$more$parnames
# first perform a model based smooth of the data
					if(!is.null(fits_temp)){
						coeftemp=ginv(t(lik$more$bvals)%*%lik$more$bvals)%*%t(lik$more$bvals)%*%fits_temp
					}else{coeftemp=ginv(t(lik$more$bvals)%*%lik$more$bvals)%*%t(lik$more$bvals)%*%data}
					
					print("here3--vector--")
					res0 = inneropt(coefs, times=times, data=datamatrix, lik=lik, proc=proc, pars=pars, in.meth="nlminb", control.in=control.in)
					
					res1 = outeropt(data=datamatrix, times=times, pars=pars[1:(parend - ncoefs)], coefs=res0$coefs, lik=lik, proc=proc, in.meth="nlminb", out.meth="nlminb",
										control.in=control.in, control.out=control.out)
					
if(plots!=0){
#par(mfrow=c(2,2))
plot(times,data[,1]);
lines(times,lik$bvals%*%res1$coefs[,1],col=2);lines(times,lik$bvals%*%coefs[,1],col=3)
lines(times,lik$bvals%*%solve(t(lik$bvals)%*%lik$bvals)%*%t(lik$bvals)%*%data[,1],col=4)

plot(times,data[,2]);

lines(times,lik$bvals%*%res1$coefs[,2],col=2);lines(times,lik$bvals%*%coefs[,2],col=3)
lines(times,lik$bvals%*%solve(t(lik$bvals)%*%lik$bvals)%*%t(lik$bvals)%*%data[,2],col=4)

legend(2,1,lty=c(1,1,1),col=c(2,3,4),list("1","c","d"))

}
					
					print("here2-----")

					if(is.null(res1$par)){res1 = outeropt(data=datamatrix, times=times, pars=pars, coefs=res0$coefs, lik=lik, proc=proc, in.meth="optim", 
														  out.meth="optim",	control.in=control.in, control.out=control.out)}
					print(paste("smoothing based optimizer  " ,res1$par))	
					 optimizer2 =optim(c(res1$pars,res1$coefs[1,]), posterior, method="BFGS", hessian=TRUE, control=list(parscale=sqrt(diag(Sig2_global)), maxit=100))
#   Estimated posterior covariance using the original optimizer to ensure the correct target posterior
					
# The center and variance estimates for teh Gaussian approximations based on the output from the optimizers.
					center_all = c(center_all, optimizer$par, optimizer2$par)
					sigma_all[[2*(i-1)+1]] = solve(optimizer$hessian)
					sigma_all[[2*(i)]] = solve(optimizer2$hessian)
					
# print("sigma_all ")	
# print(sigma_all)	
# print("center_all " )	
# print(center_all)
					X_k = c(X_k, rnorm(B, center_all[2*(i-1)+1], sqrt(sigma_all[[2*(i-1)+1]])) )		# Draw new samples based on old optimizer
					X_k = c(X_k, rnorm(B, center_all[2*i], sqrt(sigma_all[[2*(i)]])) )					# Draw new samples based on incremental based optimizer
					
				}# end of vector operations
				if (is.matrix(X_all)){	
					print("optimization starting from a value of ")
					print(X_imp)
					print("here3--matrix--")

#1. old optimizer
					optimizer = optim(X_imp, posterior,  hessian=TRUE)#,#method="BFGS",control=list(parscale=sqrt(diag(Sig2_global)), maxit=1000))
#			ptm.use = (proc.time() - ptm.opt)[3]
					print(paste("maximum posterior=", round(-optimizer$value,2),
								", likelihood=", round(likelihood(optimizer$par,log=T),2), 
								", prior=", round(log(prior(optimizer$par[1:4])),2), 
								", priorx0=", round(neglogpriorx0(optimizer$par[5:6]),2), 
								", time used=", round(ptm.use/60,2), 
								"minutes, convergence=", optimizer$convergence))
					center_all = rbind(center_all, optimizer$par)						# the center of new samples
					

					
#2. new optimizer		(the new smoothing based optimizer)
					
					
# first prep the parameters by naming them				
					pars=X_imp				
					names(pars) = proc$more$parnames
					res0 = inneropt(coefs, times=times, data=datamatrix, lik=lik, proc=proc, pars=pars, in.meth="nlminb", control.in=control.in)
					coefs = res0$coefs
# now optimize parameters
					res1 =  outeropt(data=datamatrix, times=times, pars=pars[1:4], coefs=res0$coefs, lik=lik, proc=proc, in.meth="nlminb", out.meth="nlminb",
										 control.in=control.in, control.out=control.out)
					                    
if(plots!=0){
#par(mfrow=c(2,2))
plot(times,data[,1]);
lines(times,lik$bvals%*%res1$coefs[,1],col=2);lines(times,lik$bvals%*%coefs[,1],col=3)
lines(times,lik$bvals%*%solve(t(lik$bvals)%*%lik$bvals)%*%t(lik$bvals)%*%data[,1],col=4)

plot(times,data[,2]);

lines(times,lik$bvals%*%res1$coefs[,2],col=2);lines(times,lik$bvals%*%coefs[,2],col=3)
lines(times,lik$bvals%*%solve(t(lik$bvals)%*%lik$bvals)%*%t(lik$bvals)%*%data[,2],col=4)

legend(2,1,lty=c(1,1,1),col=c(2,3,4),list("1","c","d"))

}
                    
					optimizer2 = optim(c(res1$pars,res1$coefs[1,]), posterior, method="BFGS", hessian=TRUE,
									  control=list(parscale=sqrt(diag(Sig2_global)), maxit=100))

					
					
					center_all = rbind(center_all, optimizer2$par)						# the center of new samples
					print("standard optimizer " )	
					print(optimizer$par)
					print(posterior(optimizer$par))
					print("hybrid optimizer step1 ")
					print(res1$par)
					print(posterior(c(res1$pars,res1$coefs[1,])))
					print("hybrid optimizer step2 ")
					print(optimizer2$par)
					print(posterior(optimizer2$par))
					
					
					
					if (min(eigen(optimizer$hessian)$values)>0)    	sigma_all[[2*(i-1)+1]] = solve(optimizer$hessian)						# the covariance of new samples
					if (min(eigen(optimizer2$hessian)$values)>0)    sigma_all[[2*(i)]] = solve(optimizer2$hessian)												# the covariance of new samples using the smoothing based optimizer
					
					if (min(eigen(optimizer$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
						eig = eigen(optimizer$hessian)
						eigen.values = eig$values
						eigen.values[which(eigen.values<0)] = 0
						hessian = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
						sigma_all[[2*(i-1)+1]] = solve(hessian + diag(1/diag(Sig2_global)) )
					}
					
					
					if (min(eigen(optimizer2$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
						eig = eigen(optimizer2$hessian)
						eigen.values = eig$values
						eigen.values[which(eigen.values<0)] = 0
						hessian = eig$vectors %*% diag(eigen.values) %*% t(eig$vectors)
						sigma_all[[2*i]] = solve(hessian + diag(1/diag(Sig2_global)) )
					}
					
					
############# NOTE original code used								X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[1]]) )			
############ meaning that the variance term was fixed even though it should have been changing.  I believe that this is an error in the original code				
					
					X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[2*(i-1)+1]]) )			# Draw new samples 
# print("X_k from optimizer")
# print(X_k)
############# END OF NOTE
					
					X_k = rbind(X_k, rmvnorm(B, optimizer2$par, sigma_all[[2*i]]) )			# Draw new samples using the smoothing based optimizer.
# print("X_k from smoothing based")
# print(X_k)
					
					
					print("here4----")
				}
				
				
				if(i==1){opt123=c(optimizer$par,res1$pars,res1$coefs[1,],optimizer2$par)
				}else{
					opt123=rbind(opt123,optimizer$par,res1$pars,res1$coefs[1,],optimizer2$par)
					
				}
				
#		print("which_remain 1")
#		print(which_remain)
# exclude the neighborhood of the local optima 
# first do this with the original optimizer
				if (is.matrix(X_all))	distance_remain = apply(abs(X_all[which_remain,]-optimizer$par),1,sum)
				if (is.vector(X_all))	distance_remain = abs(X_all[which_remain]-optimizer$par)
				label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
				which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
				which_remain = setdiff(which_remain, which_exclude)

				
				if(length(which_remain)>0){
# print("which_remain>0 length")
# print(which_remain)
# next do this with the new optimizer
					if (is.matrix(X_all))	distance_remain_alt = apply(abs(X_all[which_remain,]-optimizer2$par),1,sum)
					if (is.vector(X_all))	distance_remain_alt = abs(X_all[which_remain]-optimizer2$par)				
					label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
					which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
					which_remain = setdiff(which_remain, which_exclude)
				}
				
				
 print("which_remain3")
# print(which_remain)
				
				i=i+1
# print(paste("i",i))
			}  # end of "i" loop
 print("Here 3")
			if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
			if (is.vector(X_all))	X_all = c(X_all, X_k)		
		}
#print(k)
		
		
######## 
		if (k>1 | option.opt==0){
			print(paste("k>1, k=",k))
			important = which(Weights == max(Weights))
			if (length(important)>1)	important = important[1]
			if (is.matrix(X_all))	X_imp = X_all[important,]				# X_imp is the maximum weight input
			if (is.vector(X_all))	X_imp = X_all[important]
# print("X_imp")
# print(X_imp)
			if (is.matrix(X_all))	center_all = rbind(center_all, X_imp)
			if (is.vector(X_all))	center_all = c(center_all, X_imp)
# print("center_all")
# print(center_all)
			if (is.matrix(X_all))	distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )
			if (is.vector(X_all))	distance_all = abs(X_all-X_imp)			# Calculate the distances to X_imp
# print(paste("important",important))
			
			label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)		# Sort the distances
			which_var = label_nr$ix[1:B]								# Pick B inputs for covariance calculation
			if (is.matrix(X_all))	Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
			if (is.vector(X_all)){
				Weights_var = Weights[which_var]+1/length(X_all)
				Weights_var = Weights_var/sum(Weights_var)
				Sig2 = (X_all[which_var]-X_imp)^2 %*% Weights_var
			}
# print(paste("sigma_all",sigma_all))	
# print(paste("Sig2",Sig2))	
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
		
		
# print(paste("X_all",X_all))
 print("here it is")
		
		if (k>1){
			if (is.matrix(X_all)){
				gaussian_new = matrix(0, length(sigma_all), dim(X_all)[1] )
				gaussian_new[1:(dim(gaussian_all)[1]), 1:(dim(gaussian_all)[2])] = gaussian_all
				gaussian_new[length(sigma_all), ] = dmvnorm(X_all, X_imp, sigma_all[[D+k-1]])
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
	}  # end of k
	
	nonzero = which(Weights>0)
	which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
	if (is.matrix(X_all))	resample_X = X_all[which_X,]
	if (is.vector(X_all))	resample_X = X_all[which_X]
	options("warn"=0)
	return(list(stat=t(stat_all), resample=resample_X, center=center_all, optimres=opt123))
	
	
} # end of IMIS
