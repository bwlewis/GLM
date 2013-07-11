# This is a super-minimalist stripped-down GLM routine. We don't compute
# any statistics other than the model coefficients to keep the code as simple
# as possible. The idea is to highlight key aspects of the one standard
# GLM computational approach.

minimalist_glm.fit =
function(y, x, family=gaussian(), weights, maxit=25, tol=1e-08)
{
  if (missing(y)) stop("Argument y is missing")
  if (missing(x)) stop("Argument x is missing") 
  nobs  = length(y)
  nvar  = ncol(x) 
  if (missing(weights)) weights = rep(1, nobs)   


# Initialization...(I think the way family works is truly weird.
# Its initialize function expects certain variables in the local environment,
# and places a new local variable called mustart. Pretty strange.)
  etastart = start = mustart = NULL
  eval(family$initialize)
  eta          = family$linkfun(mustart)
  mu           = mustart
  beta         = rep(0,nvar)
  iter         = 0
  dev          = sum(family$dev.resids(y, mu, weights)) 
  rtol         = tol + 1
 
  while((rtol>tol) && (iter<maxit))
  {
    iter       = iter + 1
    dev0       = dev
    varmu      = family$variance(mu)
    mu.eta.val = family$mu.eta(eta)    
    z          = eta + (y - mu)/mu.eta.val    
    W          = (weights*mu.eta.val*mu.eta.val)/varmu 

    XTX = crossprod(x,(W * x))
    XTz = t(crossprod(W*z,x))

    beta  = solve(XTX, XTz, tol=2*.Machine$double.eps)
    eta   = drop(x %*% beta)
    mu    = family$linkinv(eta)
    dev   = sum(family$dev.resids(y, mu, weights)) 
    rtol  = max(abs(dev0-dev)/(abs(dev)+0.1))
  }

  list(coefficients=drop(beta), iter=iter)
}
