# Generalized Linear Model Experiments
#
# Compare a logistic regression computed using the usual maximum likelihood
# estimate with an alternative model computed directly by a nonlinear least
# squares estimate.


# We use the 'Contraception' data set from the mlmRev package in favor of a
# made-up example. You might need to install the mlmRev package to obtain
# these data. We mostly follow the very concise introduction to generalized
# linear models by Doug Bates available here:
# http://www.stat.wisc.edu/courses/st849-bates/lectures/GLMH.pdf
data("Contraception",package="mlmRev")
nobs = nrow(Contraception)

# Compute the usual maximum likeliehood estimate. We use a formula suggested
# by Doug Bates in the note cited above, worth reading! Note that we also
# specify x=TRUE, which builds and returns a model matrix for us. We'll use
# that later.
MLE = glm(formula = use ~ age + I(age^2) + urban + livch,
          family = binomial, x=TRUE, data=Contraception)

# Here is the standard model summary provided by GLM:
print(summary(MLE))

# Let's compare the maximum likelihood estimate with a generic nonlinear
# least squares solution...
# Extract the response y and model matrix x from the MLE model object.
x = M$x
y = M$y
logistic = function(x) 1/(1 + exp(-x))
f = function(beta)
{
   w  = x %*% beta
   lw = logistic(w)
   crossprod(y - lw)[1]
}
NNLS = optim(rep(0,ncol(x)), fn=f, method="BFGS")
names(NNLS$par) = colnames(x)

# We can compare the coefficients estimated by each approach. They're
# somewhat different.
print(data.frame(MLE=coef(MLE), NNLS=NNLS$par))




# Let's compare bootstrapped estimates to get a sense of how
# these two approaches vary.

t1 = proc.time()
mleboot = replicate(500,
{
  i = sample(nobs,nobs,replace=TRUE)
  glm.fit(x=x[i,],y=y[i], family=binomial())$coef
})
print(proc.time()-t1)

fboot = function(i)
{
  function(beta)
  {
    w  = x[i,] %*% beta
    lw = logistic(w)
    crossprod(y[i] - lw)[1]
  }
}
t1 = proc.time()
nnlsboot = replicate(500,
{
  i  = sample(nobs,nobs,replace=TRUE)
  f = fboot(i)
  optim(rep(0,ncol(x)), fn=f, method="BFGS")$par
})
print(proc.time()-t1)

dev.off()
split.screen(c(2,1))
xrange = c(min(min(mleboot[3,]),min(nnlsboot[3,])),
           max(max(mleboot[3,]),max(nnlsboot[3,])))
screen(1)
hist(mleboot[3,],main="Maximum likelihood estimate",col=4,xlim=xrange,breaks=25)
screen(2)
hist(nnlsboot[3,],main="Direct nonlinear optimization",col=4,xlim=xrange,breaks=25)
