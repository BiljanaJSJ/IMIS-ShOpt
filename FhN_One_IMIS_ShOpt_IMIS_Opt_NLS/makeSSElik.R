


########################################################################
#Set the CollocInfer functions
########################################################################



make.SSElik <- function () 
{
  SSE <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    #weights = mat(weights)
    f = apply(weights * difs^2, 1, sum)
    return(f)
  }
  dSSE.dx <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    #weights = mat(weights)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    g = c()
    for (i in 1:dim(dfdx)[3]) {
      g = cbind(g, apply(difs * dfdx[, , i], 1, sum))
    }
    return(-2 * g)
  }
  dSSE.dy <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    #weights = mat(weights)
    difs = weights * difs
    return(2 * difs)
  }
  dSSE.dp <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    weights = mat(weights)
    difs = weights * difs
    dfdp = more$dfdp(times, devals, pars, more$more)
    g = c()
    for (i in 1:dim(dfdp)[3]) {
      g = cbind(g, apply(difs * dfdp[, , i], 1, sum))
    }
    return(-2 * g)
  }
  d2SSE.dx2 <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    dfdx = more$dfdx(times, devals, pars, more$more)
    d2fdx2 = more$d2fdx2(times, devals, pars, more$more)
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    #weights = mat(weights)
    difs = weights * difs
    H = array(0, c(dim(devals), dim(devals)[2]))
    for (i in 1:dim(d2fdx2)[3]) {
      for (j in 1:dim(d2fdx2)[4]) {
        H[, i, j] = apply(-difs * d2fdx2[, , i, j] + 
                            weights * dfdx[, , j] * dfdx[, , i], 1, sum)
      }
    }
    return(2 * H)
  }
  d2SSE.dxdy <- function(data, times, devals, pars, more) {
    dfdx = more$dfdx(times, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, mat(dfdx[, 
                                                                 , 1]))
    #weights = mat(weights)
    weights[is.na(data)] = 0
    for (i in 1:dim(dfdx)[3]) {
      dfdx[, , i] = weights * dfdx[, , i]
    }
    return(-2 * aperm(dfdx, c(1, 3, 2)))
  }
  d2SSE.dy2 <- function(data, times, devals, pars, more) {
    r = array(0, c(dim(data), dim(data)[2]))
    ind = cbind(rep(1:dim(data)[1], dim(data)[2]), kronecker(cbind(1:dim(data)[2], 
                                                                   1:dim(data)[2]), rep(1, dim(data)[1])))
    weights = checkweights(more$weights, more$whichobs, mat(data))
    #weights = mat(weights)
    r[ind] = weights
    return(2 * r)
  }
  d2SSE.dxdp <- function(data, times, devals, pars, more) {
    fdevals = more$fn(times, devals, pars, more$more)
    difs = data - fdevals
    difs[is.na(difs)] = 0
    weights = checkweights(more$weights, more$whichobs, mat(difs))
    #weights = mat(weights)
    difs = weights * difs
    dfdx = more$dfdx(times, devals, pars, more$more)
    dfdp = more$dfdp(times, devals, pars, more$more)
    d2fdxdp = more$d2fdxdp(times, devals, pars, more$more)
    H = array(0, c(dim(devals), length(pars)))
    for (i in 1:dim(d2fdxdp)[3]) {
      for (j in 1:dim(d2fdxdp)[4]) {
        H[, i, j] = apply(-difs * d2fdxdp[, , i, j] + 
                            weights * dfdx[, , i] * dfdp[, , j], 1, sum)
      }
    }
    return(2 * H)
  }
  d2SSE.dydp <- function(data, times, devals, pars, more) {
    dfdp = more$dfdp(times, devals, pars, more$more)
    weights = checkweights(more$weights, more$whichobs, mat(dfdp[, 
                                                                 , 1]))
    #weights = mat(weights)
    weights[is.na(data)] = 0
    for (i in 1:dim(dfdp)[3]) {
      dfdp[, , i] = weights * dfdp[, , i]
    }
    return(-2 * dfdp)
  }
  return(list(fn = SSE, dfdx = dSSE.dx, dfdy = dSSE.dy, dfdp = dSSE.dp, 
              d2fdx2 = d2SSE.dx2, d2fdxdy = d2SSE.dxdy, d2fdy2 = d2SSE.dy2, 
              d2fdxdp = d2SSE.dxdp, d2fdydp = d2SSE.dydp))
}