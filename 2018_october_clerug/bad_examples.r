source("glm.svd.r")
options(warn=-1)

n = 1000
p = 30
reg.method = "column projection"
LAPACK = FALSE  # FALSE -> use RRQR modified SVD subset selection
                # TRUE -> use SVD subset selection
tol = 1e-12

set.seed(1)
X = matrix(rnorm(n*p), n)
X[sample(n*p, n*p*0.5)] = 0
X = X %*% diag(exp(-1.3*seq_len(p)))
X[,15] = rowSums(X[,-15]) / 30
beta = rep(0, p)
beta[seq_len(3)] = 1
y = X %*% beta
g = glm.fit(as.matrix(X), as.vector(y), family=gaussian(), intercept=FALSE, singular.ok=TRUE, control=list(trace=TRUE))
print(coef(g))
s = glm.svd(as.matrix(X), y, family=gaussian, reg.method=reg.method, tol=tol, LAPACK=LAPACK)
print(coef(s))
cat("glm rel error ",drop(sqrt(crossprod(na.omit(coef(g) - beta)))), "\n")
cat("irls rel error ",drop(sqrt(crossprod(na.omit(coef(s) - beta)))), "\n")


# the next two illustrate the assumption in the dqrdc2 method
# Here, X_3 = X_4 = X_5 and X_3 is a component of the solution
set.seed(1)
X = matrix(rnorm(n*p), n)
X[sample(n*p, n*p*0.5)] = 0
X[,4] = X[,3]
X[,5] = X[,3]
beta = rep(0, p)
beta[seq_len(3)] = 1
y = X %*% beta
g = glm.fit(as.matrix(X), as.vector(y), family=gaussian(), intercept=FALSE, singular.ok=TRUE, control=list(trace=TRUE))
cat("glm rel error ",drop(sqrt(crossprod(na.omit(coef(g) - beta)))), "\n")

# Here, X_3 = X_4 = X_5 and X_5 is a component of the solution, whoops!
# Any one of these is as good as the other. R picks in the order that
# the columns appear.
set.seed(1)
X = matrix(rnorm(n*p), n)
X[sample(n*p, n*p*0.5)] = 0
X[,4] = X[,3]
X[,5] = X[,3]
beta = rep(0, p)
beta[1] = 1
beta[2] = 1
beta[5] = 1
y = X %*% beta
g = glm.fit(as.matrix(X), as.vector(y), family=gaussian(), intercept=FALSE, singular.ok=TRUE, control=list(trace=TRUE))
cat("glm rel error ",drop(sqrt(crossprod(na.omit(coef(g) - beta)))), "\n")
