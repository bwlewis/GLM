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
function(A, b, family=binomial, maxit=25, tol=1e-08,
         rank_deficiency=c("select columns","minimum norm","error"))
{
  rank_deficiency = match.arg(rank_deficiency)
  m = nrow(A)
  n = ncol(A)
  select = 1:n
  S = svd(A)
  tiny_singular_values = S$d/S$d[1] < tol
  k = sum(tiny_singular_values)
  if(k>0) 
  {
    if(rank_deficiency=="select columns")
    {
      warning("Near rank-deficient model matrix")
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
  good = rep(TRUE,m)
  for(j in 1:maxit)
  {
    g       = family()$linkinv(t[good])
    gprime  = family()$mu.eta(t[good])
    z       = rep(0,m)
    W       = rep(0,m)
    z[good] = t[good] + (b[good] - g) / gprime
    W[good] = as.vector(gprime^2 / family()$variance(g))
    good   = W > .Machine$double.eps*2 & abs(W) < Inf
XXX Follow R's convetion here for error detection and good selection
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
