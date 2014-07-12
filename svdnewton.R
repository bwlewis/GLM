svdsubsel <- function(A,k=ncol(A))
{
  S <- svd(scale(A,center=FALSE,scale=TRUE), k)
  n <- which(svd(A)$d < 2*.Machine$double.eps)[1]
  if(!is.na(n) && k>=n)
  {
    k <- n - 1
    warning("k was reduced to match the rank of A")
  }
  Q <- qr( t(S$v[,1:k]) ,LAPACK=TRUE)
  sort(Q$pivot[1:k],decreasing=FALSE)
}

irls_svdnewton =
function(A, b, family=binomial, maxit=25, tol=1e-08, weights=rep(1,nrow(A)),
         rank_deficiency=c("select columns","minimum norm","error"))
{
  rank_deficiency = match.arg(rank_deficiency)
  m = nrow(A)
  n = ncol(A)
  select = 1:n
  if(any(weights==0)) A[weights,]=0
  S = svd(A)
  tiny_singular_values = S$d/S$d[1] < tol
  k = sum(tiny_singular_values)
  if(k>0) 
  {
    if(rank_deficiency=="select columns")
    {
      warning("Numerically rank-deficient model matrix")
# NB This is a different selection method than R's default glm.fit uses.
# See https://bwlewis.github.io/GLM and https://bwlewis/github.io/GLM/svdss.html
      select = svdsubsel(A,n-k)
      S = svd(A[,select,drop=FALSE])
    } else if(rank_deficiency=="error")
    {
      stop("Near rank-deficient model matrix")
    }
  }
  t = rep(0,m)
  s = rep(0,length(select))
  good = weights > 0
  for(j in 1:maxit)
  {
    g       = family()$linkinv(t[good])
    varg    = family()$variance(g)
    if(any(is.na(varg))) stop("NAs in variance of the inverse link function")
    if(any(varg==0)) stop("Zero value in variance of the inverse link function")
    gprime  = family()$mu.eta(t[good])
    if(any(is.na(gprime))) stop("NAs in the inverse link function derivative")
    z       = rep(0,m)
    W       = rep(1,m)
    z[good] = t[good] + (b[good] - g) / gprime
    W[good] = sqrt(weights[good] * as.vector(gprime^2 / varg))

    good   = W > .Machine$double.eps*2
    if(sum(good)<m) warning("Tiny weights encountered")
    s_old   = s
    C   = chol(crossprod(S$u[good,,drop=FALSE], W[good]*S$u[good,,drop=FALSE]))
    s   = forwardsolve(t(C), crossprod(S$u[good,,drop=FALSE],W[good]*z[good]))
    s   = backsolve(C,s)
    t   = rep(0,m)
    t[good] = S$u[good,,drop=FALSE] %*% s
    if(sqrt(crossprod(s - s_old)) < tol) break
  }
  x = rep(NA, n)
  if(rank_deficiency=="minimum norm") S$d[tiny_singular_values] = Inf
  x[select] = S$v %*% ((1/S$d) * crossprod(S$u,t))
  list(coefficients=x,iterations=j)
}
