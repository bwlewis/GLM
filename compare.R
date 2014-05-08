source("implementations.R")
data("Contraception",package="mlmRev")
# Model estimated with R's glm function, returning model matrix and response
# in $x and $y, respectively:
R_GLM = glm(formula = use ~ age + I(age^2) + urban + livch, family = binomial, x=TRUE, data=Contraception)
# Model estimated with our radically stripped-down minimalist implementation:
mini = irls(R_GLM$x, R_GLM$y, family=binomial)
print(data.frame(R_GLM=coef(R_GLM), minimalist=coef(mini)))

iqrn = irls_qrnewton(R_GLM$x, R_GLM$y, family=binomial)
print(data.frame(R_GLM=coef(R_GLM), qr_newton=coef(iqrn)))

isvdn = irls_svdnewton(R_GLM$x, R_GLM$y, family=binomial)
print(data.frame(R_GLM=coef(R_GLM), svd_newton=coef(isvdn)))

# Let's test the sparse-aware IRLS example. But we need some data prep for it
# first. The 1st three columns of our model matrix are dense:
library("Matrix")
A_dense = Matrix(R_GLM$x[,1:3], sparse=FALSE)
# The next four columns are sparse:
A_sparse = Matrix(R_GLM$x[,4:7], sparse=TRUE)
isparse = irls_sparse(A_dense, A_sparse, R_GLM$y, family=binomial)
print(data.frame(R_GLM=coef(R_GLM), irls_sparse=coef(isparse)))

# Let's test the incremental implementation...
# Write out the model matrix to a data file for the incremental example.
write.table(R_GLM$x, file="data.csv", sep=",", col.names=FALSE, row.names=FALSE)
inc = irls_incremental("data.csv", 500, R_GLM$y, family=binomial)
print(data.frame(R_GLM=coef(R_GLM), incremental=coef(inc)))
