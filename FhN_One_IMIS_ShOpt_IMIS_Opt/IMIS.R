#########################################################################################################
#IMIS function to run the IMIS-Opt algorithm
#input:  
#            B                      - number of incremental particles
#            B.re                   - number of re-sampled particles
#            number_k               - number of iterations
#            D                      - number of initial points for the optimizaer to start
#            logging                - boolean whether to take log
#            data                   - data
#output:     a list of 
#            stat                   - statistics like number of unique particles, marginal likelihood..
#            X_all                  - the importance sampling distribution 
#            Weights                - weigths for each point in the importance sampling ditribution
#            sigma_all              - the list of matrices for the variances at the modes (optimization stage)
#                                     or highest weigths points (in the importance stage)
#            resample               - resampled points
#            center                 - modes (optimization stage) and highest weigths points (in the importance stage)
#            data                   - data
#########################################################################################################
IMIS <- function(B, B.re, number_k, D,logging=FALSE,data){
# Modifications from the original code.
# Enabled (optional) direct log likelihood/prior calculation
# Now re-scale the likelihood and prior shifting all values to avoing cases where there is enough data to ensure that the (not logged) prior and/or likelihood of the sampled values is nothing but zeros.
# Re-set all NaN weights to values of 0.  Therse arised when the prior(theta)==0 and we attempt to divide by this value (to get the envelop_all)
# Set up the BFGX optimizer instead of the Nelder Mead	
	B0 = B*10
	X_all = X_k = sample.prior(B0)				# Draw initial samples from the prior distribution
	if (is.vector(X_all))	Sig2_global = var(X_all)	# the prior covariance
	if (is.matrix(X_all))	Sig2_global = cov(X_all)	# the prior covariance
	stat_all = matrix(NA, 6, number_k)				# 6 diagnostic statistics at each iteration
	center_all = prior_all = like_all = NULL			# centers of Gaussian components, prior densities, and likelihoods
	sigma_all = list()						# covariance matrices of Gaussian components
	if (D>=1)	option.opt = 1					# use optimizer
	if (D==0){	option.opt = 0; D=1	}			# NOT use optimizer
	
	for (k in 1:number_k ){
		print(k)
		ptm.like = proc.time()
		prior_all = c(prior_all, prior(X_k,logs=logging))		# Calculate the prior densities
#		like_all = c(like_all, likelihood(X_k,data))		# Calculate the likelihoods
		like_all = c(like_all, likelihood(X_k,logs=logging,data=data))		# Calculate the likelihoods
		
		ptm.use = (proc.time() - ptm.like)[3]
		if (k==1)	{
			
			print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
			#if (logging){
			#	scalingfactor=max(like_all)
			#	like_all=like_all-scalingfactor
			#	prior_all=prior_all-max(prior_all)
			#}
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
		Weights[is.na(Weights)]=0
		stat_all[1,k] = log(mean(Weights))			# the raw marginal likelihood
		Weights = Weights / sum(Weights)			
		stat_all[2,k] = sum(1-(1-Weights)^B.re)		# the expected number of unique points
		stat_all[3,k] = max(Weights)				# the maximum weight
		stat_all[4,k] = 1/sum(Weights^2)			# the effictive sample size
		stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))	# the entropy relative to uniform
		stat_all[6,k] = var(Weights/mean(Weights))	# the variance of scaled weights
		if (k==1)	print("Stage   MargLike   UniquePoint   MaxWeight   ESS")
		print(c(k, round(stat_all[1:4,k], 3)))
		
		if (k==1 & option.opt==1){
			print("Here 1")
			if (is.matrix(X_all))	Sig2_global = cov(X_all[which(like_all>min(like_all)),])
			X_k = which_exclude = NULL					# exclude the neighborhood of the local optima 
			label_weight = sort(Weights, decreasing = TRUE, index=TRUE)
			which_remain = which(Weights>label_weight$x[B0]) 	# the candidate inputs for the starting points
			size_remain = length(which_remain)
			for (i in 1:D){
			  
				important = which_remain[which(Weights[which_remain]==max(Weights[which_remain]))]
				if (length(important)>1)	important = sample(important,1)	
				if (is.vector(X_all))	X_imp = X_all[important]
				if (is.matrix(X_all))	X_imp = X_all[important,]
        # Remove the selected input from candidates
				print(paste('important:',important))
				which_exclude = union( which_exclude, important )
				which_remain = setdiff(which_remain, which_exclude)
        #	posterior = function(theta){	-log(prior(theta))-log(likelihood(theta,data)) } 
				if (logging){
					posterior = function(theta){	-prior(theta,logs=logging)-likelihood(theta,logs=logging,data=data) }
				}else {
					posterior = function(theta){	-log(prior(theta))-log(likelihood(theta,data=data)) } 
				}
				print("here 1.1")
				if (is.vector(X_all)){	
					print("Here 2")
					print(sqrt(Sig2_global)/10)
					print(paste('X_imp:',X_imp))
					optimizer = optim(X_imp, posterior, method="BFGS", hessian=TRUE, 
									    control=list(parscale=sqrt(Sig2_global)/10,maxit=5000))
					center_all = c(center_all, optimizer$par)
					sigma_all[[i]] = solve(optimizer$hessian)
					X_k = c(X_k, rnorm(B, optimizer$par, sqrt(sigma_all[[i]])) )			# Draw new samples
				}
				if (is.matrix(X_all)){	
        # The rough optimizer uses the Nelder-Mead algorithm.
        # ptm.opt = proc.time()
         #optimizer = optim(X_imp, posterior, method="Nelder-Mead", control=list(maxit=10, parscale=sqrt(diag(Sig2_global))) )
					theta.NM = X_imp # ####################optimizer$par
					
# The more efficient optimizer uses the BFGS algorithm 
					optimizer = optim(theta.NM, posterior, method="BFGS", hessian=TRUE,
									  control=list(parscale=sqrt(diag(Sig2_global)), maxit=1000))
					ptm.use = (proc.time() - ptm.opt)[3]
#				print(paste("maximum posterior=", round(-optimizer$value,2), ", likelihood=", round(log(likelihood(optimizer$par,data)),2), 
					if (logging){
						print(paste("maximum posterior=", round(-optimizer$value,2), ", likelihood=", round(likelihood(optimizer$par,logs=logging,data=data),2), 
									", prior=", round(prior(optimizer$par,logs=logging),2), ", time used=", round(ptm.use/60,2), "minutes, convergence=", optimizer$convergence))
					}else{
						print(paste("maximum posterior=", round(-optimizer$value,2), ", likelihood=", round(log(likelihood(optimizer$par,data=data)),2), 
									", prior=", round(log(prior(optimizer$par)),2), ", time used=", round(ptm.use/60,2), "minutes, convergence=", optimizer$convergence))
						center_all = rbind(center_all, optimizer$par)						# the center of new samples
					}
					center_all = rbind(center_all, optimizer$par)						# the center of new samples
					if (min(eigen(optimizer$hessian)$values)>0)
					sigma_all[[i]] = solve(optimizer$hessian)						# the covariance of new samples
					if (min(eigen(optimizer$hessian)$values)<=0){						# If the hessian matrix is not positive definite, we define the covariance as following
						eigen.values = eigen(optimizer$hessian)$values
						eigen.values[which(eigen.values<0)] = 0
						hessian = eigen(optimizer$hessian)$vectors %*% diag(eigen.values) %*% t(eigen(optimizer$hessian)$vectors)
						sigma_all[[i]] = solve(hessian + diag(1/diag(Sig2_global)) )
					}
					print("here 2.1")
					X_k = rbind(X_k, rmvnorm(B, optimizer$par, sigma_all[[1]]) )			# Draw new samples
				}
# exclude the neighborhood of the local optima 
				distance_remain = abs(X_all[which_remain]-optimizer$par)
				label_dist = sort(distance_remain, decreasing = FALSE, index=TRUE)
				which_exclude = union( which_exclude, which_remain[label_dist$ix[1:floor(size_remain/D)]])
				which_remain = setdiff(which_remain, which_exclude)
			}
			print("Here 3")
			if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
			if (is.vector(X_all))	X_all = c(X_all, X_k)
		}
		print("here 3.1")
		if (k>1 | option.opt==0){
			
			important = which(Weights == max(Weights))
			if (length(important)>1)	important = important[1]
			if (is.matrix(X_all))	X_imp = X_all[important,]				# X_imp is the maximum weight input
			if (is.vector(X_all))	X_imp = X_all[important]
			if (is.matrix(X_all))	center_all = rbind(center_all, X_imp)
			if (is.vector(X_all))	center_all = c(center_all, X_imp)
			if (is.matrix(X_all))	distance_all = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)) )
			if (is.vector(X_all))	distance_all = abs(X_all-X_imp)			# Calculate the distances to X_imp
			label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)		# Sort the distances
			which_var = label_nr$ix[1:B]								# Pick B inputs for covariance calculation
			if (is.matrix(X_all))	Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
			if (is.vector(X_all)){
			  Weights[which(is.nan(Weights) | is.na(Weights))]=0
				Weights_var = Weights[which_var]+1/length(X_all)
				Weights_var = Weights_var/sum(Weights_var)
				Sig2 = (X_all[which_var]-X_imp)^2 %*% Weights_var
			}
			sigma_all[[D+k-1]] = Sig2
			if (is.matrix(X_all))	X_k = rmvnorm(B, X_imp, Sig2)				# Draw new samples
			if (is.vector(X_all))	X_k = rnorm(B, X_imp, sqrt(Sig2))			# Draw new samples
			if (is.matrix(X_all))	X_all = rbind(X_all, X_k)
			if (is.vector(X_all))	X_all = c(X_all, X_k)
		}
		print("here 4")
		if (k==1){
			gaussian_all = matrix(NA, D, B0+D*B)
			for (i in 1:D){
				if (is.matrix(X_all))	gaussian_all[i,] = dmvnorm(X_all, center_all[i,], sigma_all[[i]])
				if (is.vector(X_all))	gaussian_all[i,] = dnorm(X_all, center_all[i], sqrt(sigma_all[[i]]))
			}
		}
		if (k>1){
			gaussian_new = matrix(0, D+k-1, length(X_all) )
			if (is.matrix(X_all)){
				gaussian_new[1:(D+k-2), 1:(dim(X_all)[1]-B)] = gaussian_all
				gaussian_new[D+k-1, ] = dmvnorm(X_all, X_imp, sigma_all[[D+k-1]])
				for (j in 1:(D+k-2))	gaussian_new[j, (dim(X_all)[1]-B+1):dim(X_all)[1] ] = dmvnorm(X_k, center_all[j,], sigma_all[[j]])
			}
			if (is.vector(X_all)){
				gaussian_new[1:(D+k-2), 1:(length(X_all)-B)] = gaussian_all
				gaussian_new[D+k-1, ] = dnorm(X_all, X_imp, sqrt(sigma_all[[D+k-1]]))
				for (j in 1:(D+k-2))	gaussian_new[j, (length(X_all)-B+1):length(X_all) ] = dnorm(X_k, center_all[j], sqrt(sigma_all[[j]]))
			}
			gaussian_all = gaussian_new
		}
		print("here 5")
		if (stat_all[2,k] > (1-exp(-1))*B.re)	break
	} # end of k
	print("here 5.1")
	nonzero = which(Weights>0)
	which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
	if (is.matrix(X_all))	resample_X = X_all[which_X,]
	if (is.vector(X_all))	resample_X = X_all[which_X]
	
	return(list(stat=t(stat_all),  X_all=X_all,Weights=Weights,sigma_all=sigma_all,resample=resample_X, center=center_all,data=data))
}
# end of IMIS


