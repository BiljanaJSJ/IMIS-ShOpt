######################################################################################
#make.FHN     - function that set up the FhN model
#  defines functions for evaluation of
#                  the ODE model(fhn.fun),
#                  the first derivative with respect ot x (fhn.dfdx)
#                  the first derivative with respect ot p (fhn.dfdp)
#                  the second derivative with respect ot x (fhn.d2fdx2)
#                  the second derivative with respect ot x and p (fhn.d2fdxdp)
######################################################################################

make.FHN<-function(){
    fhn.fun <- function(times, y, p, more) {

    names(p)=c("a","b","c","sdV","sdR")		
        r = y
        r[, "V"] = p["c"] * (y[, "V"] - y[, "V"]^3/3 + y[, "R"])
        r[, "R"] = -(y[, "V"] - p["a"] + p["b"] * y[, "R"])/p["c"]
		
        return(r)
    }
    fhn.fun.ode <- function(times, y, p) {

		names(p)=		names(p)=c("a","b","c","sdV","sdR")		
        r = y		
        dimnames(r) = dimnames(y)
        r["V"] = p["c"] * (y["V"] - y["V"]^3/3 + y["R"])
        r["R"] = -(y["V"] - p["a"] + p["b"] * y["R"])/p["c"]
	
        return(list(r))
    }
    fhn.dfdx <- function(times, y, p, more) {
        #print("fhn.dfdx")
        #print(p)
		names(p)=c("a","b","c","sdV","sdR")		
		r = array(0, c(dim(y), 2))
        dimnames(r) = list(NULL, c("V","R"), c("V","R"))
        r[, "V", "V"] = p["c"] - p["c"] * y[, "V"]^2
        r[, "V", "R"] = p["c"]
		
        r[, "R", "V"] = (-1/p["c"])
        r[, "R", "R"] = (-p["b"]/p["c"])
		
		
        return(r)
    }
	
    fhn.dfdp <- function(times, y, p, more) {
        #print("fhn.dfdp")
        #print(p)
		names(p)=c("a","b","c","sdV","sdR")
		r = array(0, c(dim(y), length(p)))
        dimnames(r) = list(NULL, c("V","R"), c("a","b","c","sdV","sdR"))
        r[, "V", "c"] = (y[, "V"] - y[, "V"]^3/3 + y[, "R"])
        r[, "R", "a"] = 1/p["c"]
        r[, "R", "b"] = (-y[, "R"]/p["c"])
        r[, "R", "c"] = ((y[, "V"] - p["a"] + p["b"] * y[, "R"])/(p["c"]^2))
		
		return(r)
    }
    fhn.d2fdx2 <- function(times, y, p, more) {
        #print("dfdx2")
        #print(p)
		names(p)=c("a","b","c","sdV","sdR")
		r = array(0, c(dim(y), 2, 2))
        dimnames(r) = list(NULL, c("V","R"), c("V","R"), c("V","R"))
        r[, "V", "V", "V"] = -2 * p["c"] * y[, "V"]
        return(r)
    }

	
	fhn.d2fdxdp <- function(times, y, p, more) {
        #print("d2fdxdp")
        #print(p)
		names(p)=c("a","b","c","sdV","sdR")
        r = array(0, c(dim(y), 2, length(p)))
        dimnames(r) = list(NULL, c("V","R"), c("V","R"), c("a","b","c","sdV","sdR"))
        r[, "V", "V", "c"] = 1 - y[, "V"]^2
        r[, "V", "R", "c"] = 1
        r[, "R", "V", "c"] = 1/p["c"]^2
        r[, "R", "R", "b"] = -1/p["c"]
        r[, "R", "R", "c"] = p["b"]/p["c"]^2
		
#		r[, , , "sd"] = 0
#		r[, , , "sd"] = 0
		#print(r)
		return(r)
    }

	
    return(list(fn = fhn.fun, fn.ode = fhn.fun.ode, dfdx = fhn.dfdx,dfdp = fhn.dfdp, d2fdx2 = fhn.d2fdx2, d2fdxdp = fhn.d2fdxdp))
}





###############
### SET UP THE PENALTY
# NOTE: there are two things of interest here,  The DAE that I'm using has 6 differential and 4 algebraic components and 
# that's why I have to alter the values of ddevals (the derivative of the spline basis fit) used to compare with the model output) to just devals
# (the evaluated spline fit) which then coincides with the model output.

make.SSEproc.FHN <- function(){
  SSEproc <- function(coefs, bvals, pars, more) {
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    f = make.SSElik()$fn(ddevals, more$qpts, devals, pars,  more)
    return(sum(f))
  }
  
  dSSEproc.dc <- function(coefs, bvals, pars, more) {
    
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    g1 = make.SSElik()$dfdx(ddevals, more$qpts, devals, pars, more)
    weights = checkweights(more$weights, more$whichobs, g1)
    g2 = weights * (ddevals - more$fn(more$qpts, devals,   pars, more$more))
    temp = as.matrix(t(bvals$bvals) %*% g1 + 2 * t(bvals$dbvals) %*% g2)
    
    
    g = as.vector(temp)
    
    return(g)
  }
  
  dSSEproc.dp <- function(coefs, bvals, pars, more) {
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    g = make.SSElik()$dfdp(ddevals, more$qpts, devals, pars, more)
    g = apply(g, 2, sum)
    return(g)
    
  }
  
  d2SSEproc.dc2 <- function(coefs, bvals, pars, more) {
    
    
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    H1 = make.SSElik()$d2fdx2(ddevals, more$qpts, devals,  pars, more)
    H2 = more$dfdx(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[,  , 1, drop = FALSE])
    H = list(len = dim(bvals$bvals)[2])
    
    
    for (i in 1:dim(devals)[2]) {
      H[[i]] = list(len = dim(devals))
      for (j in 1:dim(devals)[2]) {
        H[[i]][[j]] = t(bvals$bvals) %*% diag(H1[, i, j]) %*% bvals$bvals - 2 * t(bvals$dbvals) %*% diag(H2[, i, j] * weights) %*% bvals$bvals - 	2 * t(bvals$bvals) %*% diag(H2[, j, i] * weights) %*% bvals$dbvals
        
        
      }
      H[[i]][[i]] = H[[i]][[i]] + 2 * t(bvals$dbvals) %*% diag(weights) %*% bvals$dbvals
    }
    
    H = blocks2mat(H)
    return(H)
    
  }
  
  d2SSEproc.dcdp <- function(coefs, bvals, pars, more) {
    
    devals = as.matrix(bvals$bvals %*% coefs)
    ddevals = as.matrix(bvals$dbvals %*% coefs)
    colnames(devals) = more$names
    colnames(ddevals) = more$names
    names(pars) = more$parnames
    H1 = make.SSElik()$d2fdxdp(ddevals, more$qpts, devals, pars, more)
    H2 = 2 * more$dfdp(more$qpts, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, H1[, , 1, drop = FALSE])
    H = c()
    for (i in 1:length(pars)) {
      H = cbind(H, as.vector(as.matrix(t(bvals$bvals) %*%  H1[, , i] - t(bvals$dbvals) %*% (weights * H2[,, i]))))
    }
    return(H)
    
  }
  
  
  return(list(fn = SSEproc, dfdc = dSSEproc.dc, dfdp = dSSEproc.dp, 
              d2fdc2 = d2SSEproc.dc2, d2fdcdp = d2SSEproc.dcdp))
}

##################
### Set up the likelihood with sd as a model parameter

make.SSElik.with.prior.and.sd   <-function(){
    SSE <- function(data, times, devals, pars, more) {
		if(pars[3]<=0| pars[4]<=0| pars[5]<=0){return(Inf)
        }else{
            fdevals = more$fn(times, devals, pars, more$more)
            difs = data - fdevals
            difs[is.na(difs)] = 0
            weights = checkweights(.5/pars[4:5], more$whichobs, difs)
            f = apply(weights * difs^2+matrix(.5*log(pars[4:5]),ncol=dim(difs)[2],nrow=dim(difs)[1],byrow=T), 1, sum)

            return(f+neglogprior(pars)/length(f))
        }
    }
    dSSE.dx <- function(data, times, devals, pars, more) {
# print("likSSE.dx")
		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0        
        weights = checkweights(.5/pars[4:5], more$whichobs, difs)
        difs = weights * difs
        dfdx = more$dfdx(times, devals, pars, more$more)
        g = c()
        for (i in 1:dim(dfdx)[3]) {
            g = cbind(g, apply(difs * dfdx[, , i], 1, sum))
        }
        return(-2 * g)
    }
    dSSE.dy <- function(data, times, devals, pars, more) {
#print("likSSE.dy")
		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
        weights = checkweights(.5/pars[4:5], more$whichobs, difs)
		difs = weights * difs
        return(2 * difs)
    }
    dSSE.dp <- function(data, times, devals, pars, more) {
#print("likSSE.dp")

		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
		
        weights = checkweights(.5/pars[4:5], more$whichobs, difs)
        wsse = difs*weights*difs
		difs = weights * difs
		dfdp = more$dfdp(times, devals, pars, more$more)
        g = c()
        for (i in 1:dim(dfdp)[3]) {
            g = cbind(g, apply(difs * dfdp[, , i], 1, sum))
        }
	
		g[,4] = g[,4]+apply(-wsse*weights+.5/pars[4:5], 1, sum)
		
        return(g+matrix(dneglogpriordpar(pars)/dim(g)[1],dim(g)[1],dim(g)[2],byrow=T))
        
        
    }
    d2SSE.dx2 <- function(data, times, devals, pars, more) {
#	print("likSSE.dx2")
		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        dfdx = more$dfdx(times, devals, pars, more$more)
        d2fdx2 = more$d2fdx2(times, devals, pars, more$more)
        difs[is.na(difs)] = 0
		
        weights = checkweights(.5/pars[4:5], more$whichobs, difs)
		difs = weights * difs
        H = array(0, c(dim(devals), dim(devals)[2]))
        for (i in 1:dim(d2fdx2)[3]) {
            for (j in 1:dim(d2fdx2)[4]) {
                H[, i, j] = apply(-difs * d2fdx2[, , i, j] +   weights * dfdx[, , j] * dfdx[, , i], 1, sum)
            }
        }
        return(2 * H)
    }
    d2SSE.dxdy <- function(data, times, devals, pars, more) {
#        print("likSSE.dxdy")
		dfdx = more$dfdx(times, devals, pars, more$more)
		
        weights = checkweights(.5/pars[4:5], more$whichobs, dfdx[,, 1])
		weights[is.na(data)] = 0
        for (i in 1:dim(dfdx)[3]) {
            dfdx[, , i] = weights * dfdx[, , i]
        }
        return(-2 * aperm(dfdx, c(1, 3, 2)))
    }
    d2SSE.dy2 <- function(data, times, devals, pars, more) {
#		print("likSSE.dy2")
        r = array(0, c(dim(data), dim(data)[2]))
        ind = cbind(rep(1:dim(data)[1], dim(data)[2]), kronecker(cbind(1:dim(data)[2],  1:dim(data)[2]), rep(1, dim(data)[1]))) 
		weights = checkweights(.5/pars[4:5], more$whichobs, data)
        r[ind] = weights
        return(2 * r)
    }
    d2SSE.dxdp <- function(data, times, devals, pars, more) {
#		print("likSSE.dxdp")
#   print(pars)
		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        difs[is.na(difs)] = 0
		
        weights = checkweights(.5/pars[4:5], more$whichobs, difs)
        difs = weights * difs
        dfdx = more$dfdx(times, devals, pars, more$more)
        dfdp = more$dfdp(times, devals, pars, more$more)
        d2fdxdp = more$d2fdxdp(times, devals, pars, more$more)
        H = array(0, c(dim(devals), length(pars)))
        for (i in 1:dim(d2fdxdp)[3]) {
            for (j in 1:dim(d2fdxdp)[4]) {
                H[, i, j] = apply(-difs * d2fdxdp[, , i, j] +   weights * dfdx[, , i] * dfdp[, , j], 1, sum)
            }
            H[,i,4 ] = H[,i , 4] - apply(difs * ((dfdx[, , i])/pars[4]),1,sum)
            H[,i,5 ] = H[,i , 5] - apply(difs * ((dfdx[, , i])/pars[5]),1,sum)
        }
        
        
    

#print("dhereH")
        return(2 * H)
    }
    d2SSE.dydp <- function(data, times, devals, pars, more) {
		#print("likSSE.dydp")
        dfdp = more$dfdp(times, devals, pars, more$more)
		
		weights = checkweights(.5/pars[4:5], more$whichobs, dfdp[, , 1])
		fdevals = more$fn(times, devals, pars, more$more)
        difs = data - fdevals
        weights[is.na(data)] = 0
        difs[is.na(difs)] = 0
		dfdx = more$dfdx(times, devals, pars, more$more)
		g=dfdp*0
		for (i in 1:dim(dfdp)[3]) {
            g[, , i] = weights * dfdp[, , i]   
        }
		g[, , 4] = g[, , 4] + weights/pars[4] *difs
		g[, , 5] = g[, , 5] + weights/pars[5] *difs        
        return(-2 * g)
    }
    return(list(fn = SSE, dfdx = dSSE.dx, dfdy = dSSE.dy, dfdp = dSSE.dp, 
				d2fdx2 = d2SSE.dx2, d2fdxdy = d2SSE.dxdy, d2fdy2 = d2SSE.dy2, 
				d2fdxdp = d2SSE.dxdp, d2fdydp = d2SSE.dydp))
}


