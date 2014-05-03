irls =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  x = rep(0,ncol(A))
  for(j in 1:maxit)
  {
    eta   = A %*% x
    g     = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z     = eta + (b - g) / gprime
    W     = as.vector(gprime^2 / family()$variance(g))
    xold  = x
    x     = solve(crossprod(A,W*A), crossprod(A,W*z), tol=2*.Machine$double.eps)
    if(sqrt(crossprod(x-xold)) < tol) break
  }
  list(coefficients=x,iterations=j)
}

irls_qrnewton =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  r  = b
  QR = qr(A)
  Q  = qr.Q(QR)
  R  = qr.R(QR)
  for(j in 1:maxit)
  {
    eta    = b - r
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wmin   = min(W)
    if(wmin < sqrt(.Machine$double.eps)) warning("Tiny weights encountered, likely numerical problems!")
    rold   = r
    C   = chol(crossprod(Q, W*Q))
    s   = forwardsolve(t(C), crossprod(Q,W*z))
    s   = backsolve(C,s)
    r      = b - Q %*% s
    if(sqrt(crossprod(r-rold)) < tol) break
  }
  x = backsolve(R, crossprod(Q,b-r))
  list(coefficients=x,iterations=j)
}

irls_svdnewton =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  r  = b
  S  = svd(A)
  for(j in 1:maxit)
  {
    eta    = b - r
    g      = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z      = eta + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wmin   = min(W)
    if(wmin < sqrt(.Machine$double.eps)) warning("Tiny weights encountered, likely numerical problems!")
    rold   = r
    C   = chol(crossprod(S$u, W*S$u))
    s   = forwardsolve(t(C), crossprod(S$u,W*z))
    s   = backsolve(C,s)
    r      = b - S$u %*% s
    if(sqrt(crossprod(r-rold)) < tol) break
  }
  x = S$v %*% ((1/S$d) * crossprod(S$u,b-r))
  list(coefficients=x,iterations=j)
}
