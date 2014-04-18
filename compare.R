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
x = MLE$x
y = MLE$y
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

# Note that the residual norm will always be smaller for the NNLS estimate
# because that's what that approach minimizes. Also note that the MLE approach
# arrives at an estimate *much* faster.


# Let's compare bootstrapped estimates to get better a sense of how these two
# estimation approaches vary...First, let's obtain bootstrapped estimates of
# the standard MLE solution:
set.seed(1)
mleboot = replicate(1000,
{
  i = sample(nobs,nobs,replace=TRUE)
  glm.fit(x=x[i,],y=y[i], family=binomial())$coef
})

# And now with the nonlinear least squares solution:
fboot = function(i)
{
  function(beta)
  {
    w  = x[i,] %*% beta
    lw = logistic(w)
    crossprod(y[i] - lw)[1]
  }
}
nnlsboot = replicate(1000,
{
  i  = sample(nobs,nobs,replace=TRUE)
  f = fboot(i)
  optim(rep(0,ncol(x)), fn=f, method="BFGS")$par
})

# Let's compare histograms of one of the more significant model variables in
# this problem, the square of age:
jpeg(file="compare.jpg",quality=100,width=800,height=800)
split.screen(c(2,1))
xrange = c(min(min(mleboot[3,]),min(nnlsboot[3,])),
           max(max(mleboot[3,]),max(nnlsboot[3,])))
screen(1)
hist(mleboot[3,],main="age^2 maximum likelihood estimate",col=4,xlim=xrange,breaks=25)
screen(2)
hist(nnlsboot[3,],main="age^2 direct nonlinear optimization",col=4,xlim=xrange,breaks=25)
dev.off()
