irls_basic =
function(X, y, family=binomial(), maxit=25, tol=1e-08)
{
  devold = 0
  eta  = rep(0, length(y))
  w = rep(1, nrow(X))
  mu  = family$linkinv(eta)
  dev = sum(family$dev.resids(y, mu, w))
  for(j in 1:maxit)
  {
    dmu_deta = family$mu.eta(eta)
    z        = eta + (y - mu) / dmu_deta
    W        = sqrt(family$variance(mu))
    beta = qr.solve(W * X, W * z)
    eta = drop(X %*% beta)
    mu  = family$linkinv(eta)
    dev = sum(family$dev.resids(y, mu, w))
    if(abs(dev - devold) / (0.1 + abs(dev)) < tol) break
    devold = dev
  }
  list(coefficients=beta,iterations=j)
}


irls_qrnewton =
function(X, y, family=binomial(), maxit=25, tol=1e-08)
{
  devold = 0
  eta  = rep(0, length(y))
  w = rep(1, nrow(X))
  mu = family$linkinv(eta)
  dev = sum(family$dev.resids(y, mu, w))
  QR = qr(X)
  Q  = qr.Q(QR)
  R  = qr.R(QR)
  for(j in 1:maxit)
  {
    dmu_deta = family$mu.eta(eta)
    z      = eta + (y - mu) / dmu_deta
    W      = drop(family$variance(mu))
    C   = chol(crossprod(Q, W*Q))
    eta = Q %*% backsolve(C, forwardsolve(t(C), crossprod(Q,W*z)))
    mu = family$linkinv(eta)
    dev = sum(family$dev.resids(y, mu, w))
    if(abs(dev - devold) / (0.1 + abs(dev)) < tol) break
    devold = dev
  }
  beta = drop(backsolve(R, crossprod(Q,eta)))
  list(coefficients=beta,iterations=j)
}


x = c(0.5,0.75,1,1.25,1.5,1.75,1.75,2,2.25,2.5,2.75,3,3.25,3.5,4,4.25,4.5,4.75,5,5.5)
y = c(0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1)
logistic = glm.fit(cbind(1, x), y, family=binomial())
irlsqr = irls_qrnewton(cbind(1, x), y, family=binomial())
irlsbasic = irls_basic(cbind(1, x), y, family=binomial())
rbind(coef(logistic), coef(irlsqr), coef(irlsbasic))

set.seed(1)
x = matrix(rnorm(10000*800),10000)
y = sample(c(0,1),10000, replace=TRUE)
print(system.time(glm.fit(cbind(1, x), y, family=binomial())))
print(system.time(irls_qrnewton(cbind(1, x), y, family=binomial())))


