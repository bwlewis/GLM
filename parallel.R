parallel_glm.fit =
function(y, x, np=1, family=gaussian(), weights, maxit=25, tol=1e-08)
{
  if (missing(y)) stop("Argument y is missing")
  if (missing(x)) stop("Argument x is missing") 
  nobs  = length(y)
  nvar  = ncol(x) 
  if(np>nobs) stop("np too big")
  if (missing(weights)) weights = rep(1, nobs)   

# Initialization...(I think the way family works is truly weird)
  etastart = start = mustart = NULL
  eval(family$initialize)
  eta          = family$linkfun(mustart)
  mu           = mustart
  beta         = rep(0,nvar)
  iter         = 0
  dev          = sum(family$dev.resids(y, mu, weights)) 
  rtol         = tol + 1

# Row partitions used for parallel computation of the the products
  P = lapply(1:np,function(j) {h=floor(nrow(x)/np);((j-1)*h+1):max((j*h),(j==np)*nrow(x))})
# A function that combines results from the foreach workers
  combiner = function(x,y)
  {
    list(XTWX=x$XTWX + y$XTWX,
         XTWz=x$XTWz + y$XTWz)
  }
# Let's distribute the big model matrix data partitions to our workers.
# NOTE: WE ASSUME THAT THE NUMBER OF WORKERS EXACTLY MATCHES np here.
# Better approaches are available, but this one is simple.
# We'll do this just once.
  it = function(j)
  {
    if(j<=np) return(x[P[[j]],,drop=FALSE])
    stop("StopIteration")
  }
  foreach(X=iter(it)) %dopar%
  {
    rm(list="X",envir=globalenv())
    assign("X",X, envir=globalenv())
  }
  
  while((rtol>tol) && (iter<maxit))
  {
    iter       = iter + 1
    dev0       = dev
    varmu      = family$variance(mu)
    mu.eta.val = family$mu.eta(eta)    
    z          = eta + (y - mu)/mu.eta.val    
    W          = (weights*mu.eta.val*mu.eta.val)/varmu 

    PX = foreach(j=P,
                  .combine=combiner,
                  .inorder=FALSE) %dopar%
         {
          list(XTWX=crossprod(X, (W[j] * X)),
               XTWz=t(crossprod(W[j]*z[j],X)))
         }

    beta  = solve(PX$XTWX, PX$XTWz, tol=2*.Machine$double.eps)
    eta   = drop(x %*% beta)
    mu    = family$linkinv(eta)
    dev   = sum(family$dev.resids(y, mu, weights)) 
    rtol  = max(abs(dev0-dev)/(abs(dev)+0.1))
  }

  list(coefficients=drop(beta), iter=iter)
}
