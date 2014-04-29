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
Contraception$use = as.numeric(Contraception$use) - 1
nobs = nrow(Contraception)

# Compute the usual maximum likeliehood estimate. We use a formula suggested
# by Doug Bates in the note cited above, worth reading! Note that we also
# specify x=TRUE, which builds and returns a model matrix for us. We'll use
# that later.
BIN = glm(formula = use ~ age + I(age^2) + urban + livch,
          family = binomial(), x=TRUE, data=Contraception)

GAU = glm(formula = use ~ age + I(age^2) + urban + livch, start=coef(l),
          family = gaussian(link="logit"), x=TRUE, data=Contraception)

# Here is the standard model summary provided by GLM:
print(summary(BIN))

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
A = BIN$x    # model matrix
b = BIN$y    # response vector
NNLS = nnls_glm(A, b, binomial()$linkinv)

# We can compare the coefficients estimated by each approach. They're
# somewhat different.
print(data.frame(BIN=coef(BIN),GAU=coef(GAU),NNLS=NNLS))
