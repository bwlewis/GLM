# Investigate matrix products involving mixed sparse and dense matrices.
# We compute 'crossprod(A)' two ways:
# 1. Using the usual dense linear algebra code path
# 2. Computing sparse and dense parts separately and then combining
#
# The perofrmance difference you will see is highly sensitive to your
# CPU architecture, the BLAS library that you're using, the size of
# the problem, and the ratio of sparse to dense portions.

# Set up an example:
m  = 100000  # Number of rows
ns = 900     # Number of sparse columns
nd = 100     # Number of dense columns
n  = nd + ns

# The 1st nd columns of A will be dense, the rest sparse (about 1% fill-in).
set.seed(1)
A = cbind(
      matrix(rnorm(m*nd),nrow=m),
      matrix(sample(0:1,prob=c(0.99,0.01),size=ns*m,replace=TRUE), nrow=m))

# Time forming the matrix cross product
t1  = proc.time()
ATA = crossprod(A)
t1  = proc.time() - t1

# Split A into dense and sparse parts
library("Matrix")
A_dense = Matrix(A[,1:nd],sparse=FALSE)
A_sparse = Matrix(A[,(nd+1):n],sparse=TRUE)

# Time split version
t2  = proc.time()
X11 = crossprod(A_dense)
X12 = crossprod(A_dense,A_sparse)
X21 = crossprod(A_sparse,A_dense)
X22 = crossprod(A_sparse)
X   = matrix(0,nrow=n,ncol=n)
X[1:nd,1:nd] =     as.matrix(X11)
X[1:nd,(nd+1):n] = as.matrix(X12)
X[(nd+1):n,1:nd] = as.matrix(X21)
X[(nd+1):n,(nd+1):n] = as.matrix(X22)
t2 = proc.time() - t2

# Print our timings:
cat("Time to compute the usual crossprod(A):\n")
print(t1)
cat("Time to compute crossprod(A) by splitting into sparse and dense parts:\n")
print(t2)

# Verify that we compute the same thing, within expected numerical accuracy:
cat("Frobenius norm of the difference of the two methods output:\n")
print(sqrt(norm(ATA - X,"F")))
