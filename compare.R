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
nnls_glm = function(A,b,g)
{
  f = function(x) crossprod(b - g(A %*% x))[]
  ans = optim(rep(0,ncol(A)), fn=f, method="BFGS")$par
  names(ans) = colnames(A)
  ans
}

# Extract the response y and model matrix x from the MLE model object.
A = MLE$x    # model matrix
b = MLE$y    # response vector
NNLS = nnls_glm(A, b, binomial()$linkinv)

# We can compare the coefficients estimated by each approach. They're
# somewhat different.
print(data.frame(MLE=coef(MLE), NNLS=NNLS))

# Note that the residual norm will always be smaller for the NNLS estimate
# because that's what that approach minimizes. Also note that the MLE approach
# arrives at an estimate much faster.

# Let's compare bootstrapped estimates to get better a sense of how these two
# estimation approaches vary...
set.seed(1)
compare_bootstrap = replicate(1000,
{
  i = sample(nobs,nobs,replace=TRUE)
  c(mle=glm.fit(x=A[i,],y=b[i], family=binomial())$coef,
       nnls=nnls_glm(A[i,], b[i], binomial()$linkinv))
})
# (The bootstrap results are returned in a matrix where each column represents
# a resampled set of estimates. The first 7 rows contain the MLE estimates and
# the next 7 the NNLS estimates.)

# Let's compare histograms of one of the more significant model variables in
# this problem, the square of age (in the 3rd and 10th rows):
jpeg(file="compare.jpg",quality=100,width=800,height=800)
split.screen(c(2,1))
xrange = c(min(min(compare_bootstrap[3,]),min(compare_bootstrap[10,])),
           max(max(compare_bootstrap[3,]),max(compare_bootstrap[10,])))
screen(1)
hist(compare_bootstrap[3,],main="age^2 maximum likelihood estimate",col=4,xlim=xrange,breaks=25)
screen(2)
hist(compare_bootstrap[10,],main="age^2 direct nonlinear optimization",col=4,xlim=xrange,breaks=25)
dev.off()
