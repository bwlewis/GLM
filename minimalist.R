minimalist_glm =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  x = rep(0,ncol(A))
  for(j in 1:maxit)
  {
    eta    = A %*% x
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    ATWA   = crossprod(A,(W * A))
    ATWz   = t(crossprod(W*z,A))
    xold   = x
    x      = solve(ATWA, ATWz, tol=2*.Machine$double.eps)
    if(sqrt(crossprod(x-xold))<tol) break
  }
  list(coefficients=x,iterations=j)
}
