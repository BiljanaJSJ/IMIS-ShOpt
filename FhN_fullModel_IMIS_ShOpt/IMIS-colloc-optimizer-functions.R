




# This script is all about the CollocInfer functions with mods to let them hold initial states fixed.




#####################
outeropt.x0<-function (data, times, pars, coefs, lik, proc, in.meth = "nlminb", out.meth = "nlminb", control.in = list(), control.out = list(), active = 1:length(pars)) {
#coefs[1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))]=pars[lik$x0index]
	print("one")
#check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (file.exists("curcoefs.tmp")) {
        file.remove("curcoefs.tmp")
    }
    if (file.exists("optcoefs.tmp")) {
        file.remove("optcoefs.tmp")
    }
    if (file.exists("counter.tmp")) {
        file.remove("counter.tmp")
    }
    Ires = inneropt.x0(data, times, pars, coefs, lik, proc, in.meth, 
					control.in)
    ncoefs = matrix(Ires$coefs, dim(coefs))
    write.table(ncoefs, file = "optcoefs.tmp", col.names = FALSE, 
				row.names = FALSE)
    write.table(ncoefs, file = "curcoefs.tmp", col.names = FALSE, 
				row.names = FALSE)
		print("two")
    if (out.meth == "optim") {
        if (is.null(control.out$trace)) {
            control.out$trace = 6
        }
        if (is.null(control.out$maxit)) {
            control.out$maxit = 100
        }
        if (is.null(control.out$reltol)) {
            control.out$reltol = 1e-08
        }
        if (is.null(control.out$meth)) {
            control.out$meth = "BFGS"
        }
        ometh = control.out$meth
        control.out$meth = NULL
        res = optim(pars[active], ProfileErr.x0, allpars = pars, 
					times = times, data = data, coef = coefs, lik = lik, 
					proc = proc, active = active, hessian = T, in.meth = in.meth, 
					control.in = control.in, control = control.out, gr = ProfileDP.x0, 
					method = ometh)
        npar = res$par
    }
    else if (out.meth == "nlminb") {
        if (is.null(control.out$trace)) {
            control.out$trace = 10
        }
        if (is.null(control.out$eval.max)) {
            control.out$eval.max = 200
        }
        if (is.null(control.out$iter.max)) {
            control.out$iter.max = 100
        }
        if (is.null(control.out$rel.tol)) {
            control.out$rel.tol = 1e-08
        }
        res = nlminb(pars[active], ProfileErr.x0, allpars = pars, 
					 times = times, data = data, coef = coefs, lik = lik, 
					 proc = proc, in.meth = in.meth, control.in = control.in, 
					 control = control.out, gr = ProfileDP.x0, active = active)
        npar = res$par
    }
    else if (out.meth == "maxNR") {
        if (is.null(control.out$print.level)) {
            control.out$print.level = 2
        }
        if (is.null(control.out$iterlim)) {
            control.out$iterlim = 100
        }
        if (is.null(control.out$reltol)) {
            control.out$reltol = 1e-08
        }
        res = maxNR(ProfileErr1.x0, start = pars[active], allpars = pars, 
					times = times, data = data, coef = coefs, lik = lik, 
					proc = proc, in.meth = in.meth, control.in = control.in, 
					sgn = -1, active1 = active, grad = ProfileDP1.x0, print.level = control.out$print.level, 
					iterlim = control.out$iterlim)
        npar = res$estimate
    }
    else if (out.meth == "subplex") {
        if (is.null(control.out$maxit)) {
            control.out$maxit = 100
        }
        if (is.null(control.out$reltol)) {
            control.out$reltol = 1e-08
        }
        res = subplex(pars[active], ProfileErr.x0, control = control.out, 
					  hessian = FALSE, allpars = pars, times = times, data = data, 
					  coef = coefs, lik = lik, proc = proc, in.meth = in.meth, 
					  control.in = control.in, active = active)
        npar = res$par
    }
    else if (out.meth == "ProfileGN") {
        if (is.null(control.out$reltol)) {
            control.out$reltol = 1e-12
        }
        if (is.null(control.out$maxit)) {
            control.out$maxit = 50
        }
        if (is.null(control.out$maxtry)) {
            control.out$maxtry = 15
        }
        if (is.null(control.out$trace)) {
            control.out$trace = 1
        }
        res = Profile.GausNewt.x0(pars = pars, times = times, data = data, 
							   coefs = ncoefs, lik = lik, proc = proc, in.meth = in.meth, 
							   control.in = control.in, active = active, control = control.out,sgn=1)
        npar = res$pars[active]
        ncoefs = res$in.res$coefs
        g = res$in.res$df
        resid = res$in.res$f
    }
    else if (out.meth == "nls") {
        if (is.null(control.out$trace)) {
            control.out$trace = TRUE
        }
        if (is.null(control.out$maxiter)) {
            control.out$maxiter = 100
        }
        if (is.null(control.out$tol)) {
            control.out$tol = 1e-08
        }
        if (is.null(control.out$printEval)) {
            control.out$printEval = TRUE
        }
        if (is.null(control.out$warnOnly)) {
            control.out$warnOnly = TRUE
        }
        res = nls(~ProfileSSE.x0(pars, allpars, times, data, coefs, 
							  lik, proc, in.meth, control.in, active), data = list(allpars = pars, 
																		times = times, data = data, coefs = ncoefs, lik = lik, 
																		proc = proc, in.meth = in.meth, control.in = control.in, 
																		active = active), start = list(pars = pars[active]), 
																		trace = control.out$trace, control = control.out)
        npar = res$m$getPars()
        g = res$m$gradient()
        resid = res$m$resid()
    }
    else {
        stop("Unrecognized optimization method.")
    }
    if (file.exists("curcoefs.tmp")) {
        ncoefs = as.matrix(read.table("curcoefs.tmp"))
    }
    else {
        ncoefs = c()
    }
    if (file.exists("counter.tmp")) {
        counter = as.matrix(read.table("counter.tmp"))
    }
    else {
        counter = c()
    }
    if (file.exists("curcoefs.tmp")) {
        file.remove("curcoefs.tmp")
    }
    if (file.exists("optcoefs.tmp")) {
        file.remove("optcoefs.tmp")
    }
    if (file.exists("counter.tmp")) {
        file.remove("counter.tmp")
    }
    pars[active] = npar
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
    return(list(pars = pars, coefs = ncoefs, res = res, counter = counter))
}




##############
ProfileDP1.x0<-function (pars, allpars, times, data, coefs, lik, proc, in.meth = "house",control.in = NULL, sgn , sumlik = 1, active1 = 1:length(allpars)) {
print("ProfileDP1.x0")
    ProfileDP.x0(pars, allpars, times, data, coefs, lik, proc, in.meth = in.meth, 
			  control.in = control.in, sgn = sgn, sumlik = sumlik, active = active1)
}


##############fail
ProfileSSE.x0 <- function (pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb", control.in = NULL, active = 1:length(pars), dcdp = NULL, oldpars = NULL, use.nls = TRUE, sgn ){
    allpars[active] = pars
print("ProfileSSE.x0 ")
    f = ProfileErr.AllPar.x0(pars = allpars, times = times, data = data, 
						  coefs = coefs, lik = lik, proc = proc, in.meth = in.meth, 
						  control.in = control.in,sgn = sgn)
    if (use.nls) {
    }
    else {
        f$gradient = f$gradient[, active, drop = FALSE]
        f$dcdp = f$dcdp[, active, drop = FALSE]
    }
    return(f)
}



##############
ProfileErr1.x0<-function (pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb", 
control.in = NULL, sgn , active1 = 1:length(allpars)){
print("ProfileErr1.x0")  
	print(sgn)
	ProfileErr.x0(pars, allpars, times, data, coefs, lik, proc, 
			   in.meth = in.meth, control.in = control.in, sgn = sgn, active = active1)
}


#############
ProfileErr.x0<-function (pars, allpars, times, data, coefs, lik, proc, in.meth = "nlminb", 
control.in = NULL, sgn , active = 1:length(allpars)){
print("ProfileErr.x0")
		print(sgn)
#	coefs2=rep(0,length=dim(lik$bvals)[2]*(length(pars[lik$x0index])))
#	coefs2[1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))]=pars[lik$x0index]
#	for(lp in 1:length(pars[lik$x0index])){
#		coefs2[(2+(dim(lik$bvals)[2]*(lp-1))):(dim(lik$bvals)[2]*lp)]=coefs[(1+((dim(lik$bvals)[2]-1)*(lp-1))):((dim(lik$bvals)[2]-1)*lp)]		
#	}	
	allpars[active] = pars
    f = ProfileErr.AllPar.x0(allpars, times, data, coefs, lik, proc,  in.meth, control.in, sgn)
#print(f)
#print("ProfileErr.x0-out")
    return(f)
}







############### required function if using X0 as a parameter
ProfileErr.AllPar.x0 <- function (pars, times, data, coefs, lik, proc, in.meth = "nlminb", control.in = NULL, sgn ) {
	print("ProfileErr.AllPar.x0")
	print(sgn)
#	print(coefs)
	print(pars)
    if (file.exists("curcoefs.tmp")) {
        altcoefs = as.matrix(read.table("curcoefs.tmp"))
        if (!(length(altcoefs) == length(coefs))) {
            stop(paste("Variables in curcoefs.tmp do not conform;", 
					   " file exists from previous experiments?"))
        }
    } else {
        altcoefs = coefs
    }

	
	if (file.exists("counter.tmp")) {
        counter = read.table("counter.tmp")
        niter = counter[nrow(counter), 1]
    }else {
        counter = matrix(c(1, 0, pars), 1, length(pars) + 2)
        niter = 0
    }
	altdevals = as.matrix(lik$bvals %*% rbind(pars[lik$x0index],matrix(altcoefs,  ncol=(length(lik$x0index)+length(altcoefs))/ncol(lik$bvals))))
	
    colnames(altdevals) = proc$more$names
    f1 = SplineCoefsErr.x0(coefs, times, data, lik, proc, pars,sgn)
    f2 = SplineCoefsErr.x0(altcoefs, times, data, lik, proc, pars,sgn)
	
    if (f2 < f1) {
        coefs = altcoefs
    }
    Ires = inneropt.x0(data, times, pars, coefs, lik, proc, in.meth, 
					   control.in)
    ncoefs = Ires$coefs
    print("here in profileerr.allpar.x0===============================3")
	
	devals = as.matrix(lik$bvals %*% rbind(pars[lik$x0index],ncoefs))
    
	colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more))
    f2 = sum(lik$fn(data, times, altdevals, pars, lik$more))
    
	if (f <= f2) {
        write.table(ncoefs, file = "curcoefs.tmp", col.names = FALSE, 
					row.names = FALSE)
    }
    write.table(ncoefs, file = "optcoefs.tmp", col.names = FALSE, 
				row.names = FALSE)
    if (niter == 0) {
        counter[1, 2] = f
        write.table(counter, file = "counter.tmp", col.names = FALSE, 
					row.names = FALSE)
    }
    if (niter >= 1) {
        if (f < counter[niter, 2]) {
            counter = rbind(counter, c(niter + 1, f, pars))
            write.table(counter, file = "counter.tmp", col.names = FALSE, 
						row.names = FALSE)
        }
    }
    if (!is.null(lik$report)) {
		print(sgn)
    }
    f = sgn * f
#print(f)
#print(pars)
	print("ProfileErr.AllPar.x0 out")
    return(f)
#    return(list(value=f,gradient =ProfileDP.AllPar.x0(pars, times, data, ncoefs, lik, proc, in.meth = in.meth, control.in = control.in, sgn = sgn, sumlik = 1), npar = par))
}


ProfileDP.x0<-function (pars, allpars, times, data, coefs, lik, proc, in.meth = "house", control.in = NULL, sgn , sumlik = 1, active = 1:length(allpars)) {
	
	allpars[active] = pars
#print("ProfileDP.x0")
    g = ProfileDP.AllPar.x0(allpars, times, data, coefs, lik, proc, 
							in.meth, control.in, sgn, sumlik)
    if (sumlik) {
        names(g) = names(allpars)
        g = g[active]
    }
    else {
        colnames(g) = names(allpars)
        g = g[, active, drop = FALSE]
    }
#print(g)
#print("ProfileDP.x0-out")
    return(g)
}

#################
ProfileDP.AllPar.x0<-function (pars, times, data, coefs, lik, proc, in.meth = NULL, control.in = NULL, sgn , sumlik = 1){
	print("ProfileDP.AllPar.x0")
    if (file.exists("optcoefs.tmp")) {
        altcoefs = as.matrix(read.table("optcoefs.tmp"))
        if (!(length(altcoefs) == length(coefs))) {
            stop(paste("Variables in curcoefs.tmp do not conform;", 
					   "file exists from previous experiments? ---ProfileDP.AllPar.x0"))
        }else {
            coefs = altcoefs
        }
    }
    coefs = as.matrix(coefs)
    devals = as.matrix(lik$bvals %*% rbind(pars[lik$x0index],coefs))
    colnames(devals) = proc$more$names
    d2Hdc2 = SplineCoefsDC2sparse.x0(coefs, times, data, lik, proc, 
									 pars,sgn)
    d2Hdcdp = SplineCoefsDCDP.x0(coefs, times, data, lik, proc, 
								 pars,sgn)
    if (is.matrix(d2Hdc2)) {
        dcdp = -ginv(d2Hdc2) %*% d2Hdcdp
    }else {
        dcdp = -as.matrix(solve(d2Hdc2, d2Hdcdp))
    }
    if (sumlik) {
        dlikdc = as.vector(t(lik$bvals) %*% lik$dfdx(data, times, 
													 devals, pars, lik$more))
		
        df = t(dcdp) %*% (dlikdc[-c((0:(length(lik$x0index)-1)*dim(lik$bvals)[2])+1)]) + apply(lik$dfdp(data, times, 
																										devals, pars, lik$more), 2, sum)
        df = as.vector(sgn * df)
#print(max(abs(df)))
#print("ProfileDP.AllPar.x0 - out")
##print((df))
        return(df)
    }
    else {
        dlikdx = lik$dfdx(data, times, devals, pars, lik$more)
        dlikdp = lik$dfdp(data, times, devals, pars, lik$more)
        dlikdc = c()
        for (i in 1:dim(dlikdx)[2]) {
            dlikdc = cbind(dlikdc, as.matrix(diag(dlikdx[, i]) %*% 
											 lik$bvals))
        }
        df = dlikdc[,-c((0:(length(lik$x0index)-1)*dim(lik$bvals)[2])+1)] %*% dcdp + dlikdp
#print(apply(abs(df),1,max))
		
		print("ProfileDP.AllPar.x0 - out")
##print((df))
        return(as.matrix(df))
    }
}



#################
SplineCoefsErr.x0<-function(coefs, times, data, lik, proc, pars, sgn ){
#print("SplineCoefsErr.x0")
#	print(pars)
	coefs2 = rbind(pars[lik$x0index],matrix(coefs, ncol(lik$bvals)-1))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    f = sum(lik$fn(data, times, devals, pars, lik$more)) + proc$fn(coefs2, 
																   proc$bvals, pars, proc$more)
    if (!is.null(proc$report)) {
#print(f)
    }
	
#print(sgn * f)
	
#print("SplineCoefsErr.x0 now")
    return(sgn * f)
}




#################

SplineCoefsDC.x0<-function (coefs, times, data, lik, proc, pars, sgn ){
#print("SplineCoefsDC.x0")
	coefs2 = rbind(pars[lik$x0index],matrix(coefs, ncol(lik$bvals)-1))
#print(dim(t(lik$bvals)))
#	print("lik$bvals")
#	print(coefs)
#	print("coefs")
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    g = as.matrix(t(lik$bvals) %*% lik$dfdx(data, times, devals, pars, lik$more)) + 
	proc$dfdc(coefs2, proc$bvals, pars, proc$more)
    g = as.vector(g)
#	print(g)
#print("SplineCoefsDC.x0-out")
    return(sgn * g[-(1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index]))))])
# need to trim out the X0 coefs
}




#################

SplineCoefsDC2.x0<-function (coefs, times, data, lik, proc, pars, sgn ){
#print("SplineCoefsDC2.x0")
    result = as.matrix(SplineCoefsDC2sparse.x0(coefs, times, data, lik, proc, pars, sgn))
    return(result)
}


#################

SplineCoefsDCDP.x0<-function (coefs, times, data, lik, proc, pars, sgn ){
	print("SplineCoefsDCDP.x0")
	coefs2 = rbind(pars[lik$x0index],matrix(coefs, ncol(lik$bvals)-1))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    d2lik = lik$d2fdxdp(data, times, devals, pars, lik$more)
    H = c()
    for (i in 1:length(pars)) {
        H = cbind(H, as.vector(as.matrix(t(lik$bvals) %*% d2lik[, 
										 , i])))
    }
    H = H + proc$d2fdcdp(coefs2, proc$bvals, pars, proc$more)
#print(apply(H[-c((0:(length(lik$x0index)-1)*dim(lik$bvals)[2])+1),],2,max))
	print("SplineCoefsDCDP.x0-out")
    return(as.matrix(sgn * H[-c((0:(length(lik$x0index)-1)*dim(lik$bvals)[2])+1),]))
}





######################  
Profile.GausNewt.x0  <-function (pars, times, data, coefs, lik, proc, in.meth = "nlminb", control.in = NULL, active = 1:length(pars), 
control = list(reltol = 1e-06, maxit = 50, maxtry = 15, trace = 1),sgn=1){
print("Profile.GausNewt.x0 ")
    if (is.null(control)) {
        control = list()
    }
    if (is.null(control$reltol)) {
        control$reltol = 1e-12
    }
    if (is.null(control$maxit)) {
        control$maxit = 50
    }
    if (is.null(control$maxtry)) {
        control$maxtry = 15
    }
    if (is.null(control$trace)) {
        control$trace = 1
    }
#coefs[1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))]=pars[lik$x0index]
#   check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    res0 = ProfileSSE.AllPar.x0(pars, times, data, coefs, lik, proc, in.meth, control.in,sgn=sgn)
    ind = names(pars)
    if (is.null(ind)) {
        ind = 1:length(pars)
    }
    ind = ind[active]
    pars0 = pars
    pars1 = pars
	
	print("pars1---------------------------------------------")
    
    print(attributes(res0))
	F0 = sum(res0$value^2)
    F1 = F0
    iter = 0
    fundif = 1
	
	print("pars1.1---------------------------------------------")
    
    Jacobian = res0$gradient[, ind]
    residual = res0$value
    gradnorm0 = mean(abs(crossprod(Jacobian, residual)))
    gradnorm1 = 1
    dcdp = c()

	while (gradnorm1 > control$reltol & fundif > control$reltol & 
		   iter < control$maxit) {
        iter = iter + 1
#Dpars = lsfit(Jacobian, residual, int = FALSE)
        Dpars = solve(crossprod(Jacobian), crossprod(Jacobian,  residual))
        ntry = 0
        while (F1 >= F0 & t(Dpars) %*% Dpars > control$reltol & 
			   ntry < control$maxtry) {
            pars1[ind] = pars0[ind] - 0.5 * Dpars
print("pars1")
print(pars1)
            tres = ProfileSSE.AllPar.x0(pars1, times, data, res0$coefs,  lik, proc, in.meth, control.in,sgn=sgn)
            F1 = sum(tres$value^2)
            Dpars = Dpars/2
            ntry = ntry + 1
        }
		tres$coefs[1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))]=pars[lik$x0index]
        res0 = tres
		dcdp = tres$dcdp
        gradnorm1 = abs(mean(t(res0$gradient[, ind]) %*% res0$value))
        fundif = (F0 - F1)/abs(F0)
        pars0 = pars1
        res0 = tres
        F0 = F1
        gradnorm0 = gradnorm1
        if (control$trace > 0) {
print(c(iter, ntry, F0, pars0))
        }
    }
    newpars = pars0
print("Profile.GausNewt.x0 -out")
    return(list(pars = newpars, in.res = res0, value = F0))
}






#################
inneropt.x0<-function (data, times, pars, coefs, lik, proc, in.meth = "nlminb",   control.in = list()) {
#print("inneropt.x0")
# check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    if (in.meth == "SplineEst") {
        if (is.null(control.in$reltol)) {
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$maxtry)) {
            control.in$maxtry = 10
        }
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        res = SplineEst.NewtRaph.x0(coefs, times, data, lik, proc,
            pars, control.in, sgn=1)
        ncoefs = matrix(res$coefs, ncol= length(res$coefs)/(ncol(lik$bvals)-1))
    }
    else if (in.meth == "optim") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$maxit)) {
            control.in$maxit = 1000
        }
        if (is.null(control.in$reltol)) {
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$meth)) {
            control.in$meth = "BFGS"
        }
        if (is.null(control.in$reportHessian)) {
            control.in$reportHessian = TRUE
        }
        imeth = control.in$meth
        control.in$meth = NULL
        res = optim(coefs, SplineCoefsErr.x0, gr = SplineCoefsDC.x0, 
            hessian = control.in$reportHessian, control = control.in, 
            times = times, data = data, lik = lik, proc = proc, 
            pars = pars, sgn=1,method = in.meth)
        ncoefs = matrix(res$par, ncol=length(res$par)/(ncol(lik$bvals)-1))
    }
    else if (in.meth == "nlminb") {
        if (is.null(control.in$trace)) {
            control.in$trace = 0
        }
        if (is.null(control.in$eval.max)) {
            control.in$eval.max = 2000
        }
        if (is.null(control.in$iter.max)) {
            control.in$iter.max = 1000
        }
        if (is.null(control.in$rel.tol)) {
            control.in$rel.tol = 1e-12
        }
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2.x0
        }else {
            Hessian = NULL
        }
#print(coefs)
        res = nlminb(coefs, SplineCoefsErr.x0, gradient = SplineCoefsDC.x0, 
            hessian = Hessian, control = control.in, times = times, 
            data = data, lik = lik, proc = proc, pars = pars,sgn=1)
        ncoefs = matrix(res$par, ncol = length(res$par)/(ncol(lik$bvals)-1))
    }
    else if (in.meth == "maxNR") {
        if (is.null(control.in$print.level)) {
            control.in$print.level = 1
        }
        if (is.null(control.in$iterlim)) {
            control.in$iterlim = 1000
        }
        if (is.null(control.in$reltol)) {
            control.in$reltol = 1e-12
        }
        if (is.null(control.in$useHessian)) {
            Hessian = SplineCoefsDC2.x0
        }
        else {
            Hessian = NULL
        }
        res = maxNR(SplineCoefsErr.x0, coefs, times = times, data = data, 
            lik = lik, proc = proc, pars = pars, sgn = -1, grad = SplineCoefsDC.x0, 
            hess = Hessian, print.level = control.in$print.level, 
            iterlim = control.in$iterlim)
        ncoefs = matrix(res$estimate, ncol = length(res$estimate)/(ncol(lik$bvals)-1))
    }
    else if (in.meth == "trust") {
        if (is.null(control.in$rinit)) {
            control.in$rinit = 1
        }
        if (is.null(control.in$rmax)) {
            control.in$rmax = 100
        }
        if (is.null(control.in$parscale)) {
            control.in$parscale = abs(coefs)
        }
        if (is.null(control.in$iterlim)) {
            control.in$iterlim = 100
        }
        res = trust(SplineCoefsList.x0, coefs, rinit = control.in$rinit, 
            rmax = control.in$rmax, parscale = control.in$parscale, 
            iterlim = control.in$iterlim, times = time, data = data, 
            lik = lik, proc = proc, pars = pars, sgn = 1)
        ncoefs = matrix(res$argument, dim(coefs))
    }
    else {
        stop("Unknown optimizer specified")
    }
    if (!is.null(proc$more$names)) {
        colnames(ncoefs) = proc$more$names
    }
#print("inneropt.x0 - out")
    return(list(coefs = ncoefs, res = res))
}

#################

SplineCoefsList.x0<-function (coefs, times, data, lik, proc, pars, sgn){
#print("SplineCoefsList.x0")
    value = SplineCoefsErr.x0(coefs, times, data, lik, proc, pars, sgn)
    gradient = SplineCoefsDC.x0(coefs, times, data, lik, proc, pars, sgn)
    hessian = SplineCoefsDC2.x0(coefs, times, data, lik, proc, pars, sgn)
    return(list(value = value, gradient = gradient, hessian = hessian))
}



###########
SplineEst.NewtRaph.x0 <-function (coefs, times, data, lik, proc, pars, control = list(reltol = 1e-12, maxit = 1000, maxtry = 10, trace = 0,sgn)){
print("SplineEst.NewtRaph.x0")
	
    if (is.null(control)) {
        control = list()
    }
    if (is.null(control$reltol)) {
        control$reltol = 1e-12
    }
    if (is.null(control$maxit)) {
        control$maxit = 1000
    }
    if (is.null(control$maxtry)) {
        control$maxtry = 10
    }
    if (is.null(control$trace)) {
        control$trace = 0
    }
#    check.lik.proc.data.coefs(lik, proc, data, times, coefs)
    f0 = SplineCoefsErr.x0(coefs, times, data, lik, proc, pars,sgn=sgn)
	
    g = SplineCoefsDC.x0(coefs, times, data, lik, proc, pars,sgn=sgn)
	
    H = SplineCoefsDC2sparse.x0(coefs, times, data, lik, proc, pars,sgn=sgn)

    if (is.matrix(H)) {
        DC = -ginv(H) %*% g
    }
    else {
        DC = -as.matrix(solve(H, g))
    }
    gradnorm1 = 1
    fundif = 1
    iter = 0
    f1 = f0
    coefs0 = as.vector(coefs)

    while (gradnorm1 > control$reltol & fundif > 0 & iter < control$maxit) {
        iter = iter + 1
        if (is.matrix(H)) {
            DC = -ginv(H) %*% g
        }
        else {
            DC = -as.matrix(solve(H, g))
        }
        ntry = 0
        coefs1 = coefs0
        if (control$trace > 0) {
            print(c(f0, mean(abs(DC)), 0))
        }
        while ((f1 >= f0) & (t(DC) %*% DC > control$reltol) & 
			   (ntry < control$maxtry)) {
            coefs1 = coefs0 + DC
            f1 = SplineCoefsErr.x0(coefs1, times, data, lik, proc, 
								pars,sgn=sgn)
            DC = DC/2
            ntry = ntry + 1
        }
        coefs0 = coefs1
        g = SplineCoefsDC.x0(coefs0, times, data, lik, proc, pars)

        H = SplineCoefsDC2sparse.x0(coefs0, times, data, lik, proc, 
								 pars)
		
		
        gradnorm1 = mean(abs(DC))
        fundif = (f0 - f1)/abs(f0)
        f0 = f1
        if (control$trace > 1) {
            print(c(iter, ntry, f0, gradnorm1, fundif))
        }
    }
    if (control$trace > 0) {
        print(c(f0, gradnorm1, iter))
    }

    return(list(coefs = coefs0, g = g, value = f0, H = H))
}

##################

##################
SplineCoefsDC2sparse.x0<-function (coefs, times, data, lik, proc, pars, sgn ){
#print("SplineCoefsDC2sparse.x0")
	print(sgn)
	coefs2 = rbind(pars[lik$x0index],matrix(coefs, ncol(lik$bvals)-1))
    devals = as.matrix(lik$bvals %*% coefs2)
    colnames(devals) = proc$more$names
    d2lik = lik$d2fdx2(data, times, devals, pars, lik$more)
    H = list(len = ncol(lik$bvals))
    for (i in 1:dim(d2lik)[2]) {
        H[[i]] = list(len = ncol(devals))
        for (j in 1:dim(d2lik)[3]) {
            H[[i]][[j]] = t(lik$bvals) %*% diag(d2lik[, i, j]) %*% 
			lik$bvals
        }
    }
    H = blocks2mat(H)
    H = H + proc$d2fdc2(coefs2, proc$bvals, pars, proc$more)
##print(H[-(1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))),-(1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index]))))])	
#print("SplineCoefsDC2sparse.x0 - out")
    return(sgn * H[-(1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index])))),-(1+dim(lik$bvals)[2]*c(0:(-1+length(pars[lik$x0index]))))])
}





###################
ProfileSSE.AllPar.x0<-function(pars, times, data, coefs, lik, proc, in.meth = "nlminb", 
control.in = NULL, dcdp = NULL, oldpars = NULL, use.nls = FALSE, sgn=1){
    f1 = SplineCoefsErr.x0(coefs, times, data, lik, proc, pars,sgn)
    if (use.nls) {
        if (file.exists("curcoefs.tmp")) {
            altcoefs = as.matrix(read.table("curcoefs.tmp"))
            if (!(length(altcoefs) == length(coefs))) {
                stop(paste("Variables in curcoefs.tmp do not conform;", 
						   "file exists from previous experiments?"))
            }
        }else {
            altcoefs = coefs
        }
        if (file.exists("counter.tmp")) {
            counter = read.table("counter.tmp")
            niter = counter[nrow(counter), 1]
        }else {
            counter = matrix(c(1, 0, pars), 1, length(pars) +  2)
            niter = 0
        }
        f2 = SplineCoefsErr.x0(altcoefs, times, data, lik, proc, pars,sgn)
        if (f2 < f1) {
            coefs = altcoefs
            f1 = f2
        }
        altdevals = as.matrix(lik$bvals %*% rbind(pars[lik$x0index],matrix(altcoefs,ncol(lik$bvals)-1)))
        colnames(altdevals) = proc$more$names
    }
    if (!is.null(dcdp)) {
        tcoefs = as.vector(coefs) + dcdp %*% (pars - oldpars)
        f2 = SplineCoefsErr.x0(tcoefs, times, data, lik, proc, pars)
        if (f2 < f1) {
            coefs = tcoefs
            f1 = f2
        }
							  }
    Ires = inneropt.x0(data, times, pars, coefs, lik, proc, in.meth,control.in)
    ncoefs = Ires$coefs
    devals = as.matrix(lik$bvals %*%rbind(pars[lik$x0index], ncoefs))
    colnames(devals) = proc$more$names
    weights = checkweights(lik$more$weights, lik$more$whichobs,  data)
    f = as.vector(as.matrix(data - lik$more$fn(times, devals,  pars, lik$more$more)) * sqrt(weights))
    isnaf = is.na(f)
    f[isnaf] = 0
    dlikdp = lik$more$dfdp(times, devals, pars, lik$more$more)
    dlikdp = matrix(dlikdp, dim(dlikdp)[1] * dim(dlikdp)[2], dim(dlikdp)[3])
    dlikdx = lik$more$dfdx(times, devals, pars, lik$more$more)
    dlikdc = c()
    for (i in 1:dim(dlikdx)[2]) {
        tH = c()
        for (j in 1:dim(dlikdx)[3]) {
            tH = cbind(tH, as.matrix(diag(dlikdx[, i, j]) %*%  lik$bvals))
        }
        dlikdc = rbind(dlikdc, tH)
    }
	dlikdc=dlikdc[,-c(dim(lik$bvals)[2]*c(0:(dim(dlikdx)[2]-1))+1)]
    d2Hdc2 = SplineCoefsDC2sparse.x0(ncoefs, times, data, lik, proc, pars,sgn)
    d2Hdcdp = SplineCoefsDCDP.x0(ncoefs, times, data, lik, proc, pars,sgn)
    if (is.matrix(d2Hdc2)) {
        dcdp = ginv(d2Hdc2) %*% d2Hdcdp
    }else {
        dcdp = as.matrix(solve(d2Hdc2, d2Hdcdp))
    }
    df = dlikdc %*% dcdp + dlikdp
    df[isnaf, ] = 0
    colnames(df) = proc$more$parnames
    if (!is.null(lik$report)) {
        print(f)
    }
    f = sgn * f
    df = sgn * df
#if (use.nls) {
        tf = sum(lik$fn(data, times, devals, pars, lik$more))
        tf2 = sum(lik$fn(data, times, altdevals, pars, lik$more))
        if (tf <= tf2) {
            write.table(ncoefs, file = "curcoefs.tmp", row.names = FALSE, 
						col.names = FALSE)
            niter = counter[nrow(counter), 1]
        }
        if (niter == 0) {
            counter[1, 2] = tf
            write.table(counter, file = "counter.tmp", col.names = FALSE, 
						row.names = FALSE)
        }
        if (niter > 1) {
            if (tf < counter[niter, 2]) {
                counter = rbind(counter, c(niter + 1, tf, pars))
                write.table(counter, file = "counter.tmp", col.names = FALSE, 
							row.names = FALSE)
            }
        }
        write.table(ncoefs, file = "optcoefs.tmp", col.names = FALSE, 
					row.names = FALSE)
#       attr(f, "gradient") = df
#       return(f)
#   }
#   else {
        return(list(value = f, gradient = df, coefs = ncoefs, 
					dcdp = dcdp))
#   }
}
