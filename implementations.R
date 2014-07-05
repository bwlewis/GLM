# Example iteratively re-weighted least squares (IRLS) implementations
# Mike Kane & Bryan Lewis, 2013-2014.
#
# The implementations generally follow the same input/output pattern.  They
# take as inputs a model matrix A, a response vector b whose length is the
# number of rows of A, an R 'family' function that defines the error
# distribution family and link function, a maximum number of iterations, and an
# iteration convergence tolerance. The methods produce a list with two
# elements, the model coefficients and the number of iterations.

# The most basic IRLS method, and the shortest implementation we could come
# up with. This method solves the normal equations associated with a weighted
# least squares problem in each iteration.
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

# A method discussed by O'Leary that uses a QR factorization of the model
# matrix. This method should be much more numerically stable in the face of
# ill-conditioned model matrices than the simple method defined above.  If the
# QR method used uses Givens rotations, this method is numerically stable for
# stiff problems too.
irls_qrnewton =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  t  = 0
  QR = qr(A)
  Q  = qr.Q(QR)
  R  = qr.R(QR)
  for(j in 1:maxit)
  {
    g      = family()$linkinv(t)
    gprime = family()$mu.eta(t)
    z      = t + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wmin   = min(W)
    if(wmin < sqrt(.Machine$double.eps))
      warning("Tiny weights encountered")
    told   = t
    C   = chol(crossprod(Q, W*Q))
    s   = forwardsolve(t(C), crossprod(Q,W*z))
    s   = backsolve(C,s)
    t      = Q %*% s
    if(sqrt(crossprod(t - told)) < tol) break
  }
  x = backsolve(R, crossprod(Q,t))
  list(coefficients=x,iterations=j)
}

# The next method is a minor variation on the QR Newton method defined above
# that uses the SVD instead. It exhibits similar numerical stability and can
# definitively check model matrix rank deficiency, at the cost of computing
# the SVD instead of the QR factorization up front.
irls_svdnewton =
function(A, b, family=binomial, maxit=25, tol=1e-08)
{
  t = 0
  S  = svd(A)
  if(min(S$d)/max(S$d)<tol) warn("Near rank-deficient model matrix")
  for(j in 1:maxit)
  {
    g      = family()$linkinv(t)
    gprime = family()$mu.eta(t)
    z      = t + (b - g) / gprime
    W      = as.vector(gprime^2 / family()$variance(g))
    wmin   = min(W)
    if(wmin < sqrt(.Machine$double.eps))
      warning("Tiny weights encountered")
    told   = t
    C   = chol(crossprod(S$u, W*S$u))
    s   = forwardsolve(t(C), crossprod(S$u,W*z))
    s   = backsolve(C,s)
    t      = S$u %*% s
    if(sqrt(crossprod(t - told)) < tol) break
  }
  x = S$v %*% ((1/S$d) * crossprod(S$u,t))
  list(coefficients=x,iterations=j)
}


# Sparse weighted cross product helper function
# Input: Dense Matrix A_dense, sparse Matrix A_sparse, weights W,
# where A = [A_dense, A_sparse] and length W=ncol(A).
# Output: Dense representation of crossprod(A,W*A)
sp_wt_cross = function(A_dense, A_sparse, W)
{
  nd = ncol(A_dense)
  ns = ncol(A_sparse)
  n  = nd + ns
  ATWA = matrix(0, nrow=n, ncol=n)
  WA_dense  = W*A_dense
  WA_sparse = W*A_sparse
  ATWA[1:nd,1:nd] = as.matrix(crossprod(A_dense,WA_dense))
  ATWA[1:nd,(nd+1):n] = as.matrix(crossprod(A_dense,WA_sparse))
  ATWA[(nd+1):n,1:nd] = as.matrix(crossprod(A_sparse,WA_dense))
  ATWA[(nd+1):n,(nd+1):n] = as.matrix(crossprod(A_sparse,WA_sparse))
  ATWA
}

# Example IRLS implementation that can take advantage of sparse model matrices.
# Here we assume that the model matrix A is already permuted and partitioned
# into A = [A_dense, A_sparse] dense and sparse columns.  The response vector b
# must already be permuted on input to correspond to the matrix splitting.
irls_sparse =
function(A_dense, A_sparse, b, family=binomial, maxit=25, tol=1e-08)
{
  nd = ncol(A_dense)
  ns = ncol(A_sparse)
  n  = nd + ns
  x = rep(0, n)
  for(j in 1:maxit)
  {
    eta   = as.vector(A_dense %*% x[1:nd] + A_sparse %*% x[(nd+1):n])
    g     = family()$linkinv(eta)
    gprime = family()$mu.eta(eta)
    z     = eta + (b - g) / gprime
    W     = as.vector(gprime^2 / family()$variance(g))
    xold  = x
    ATWA  = sp_wt_cross(A_dense,A_sparse,W)
    wz    = W*z
    ATWz  = c(as.vector(crossprod(A_dense, wz)), as.vector(crossprod(A_sparse, wz)))

    C   = chol(ATWA, pivot=TRUE)
    if(attr(C,"rank")<ncol(C)) stop("Rank-deficiency detected.")
    p   = attr(C, "pivot")
    s   = forwardsolve(t(C), ATWz[p])
    x   = backsolve(C,s)[p]

    if(sqrt(crossprod(x-xold)) < tol) break
  }
  list(coefficients=x,iterations=j)
}



# The following simple iterator function is required by irls_incremental below.
# This function iterates by rows through a delimited text file nrows at a time,
# returning NULL at the end. For more sophisticated iterators, see the
# iterators or lazy.frame packages. Use init=TRUE argument to initialize/reset
# iterator. Example:
# chunk = iterator("mydata.csv")
# chunk(init=TRUE)
# chunk() ... until it returns NULL
iterator = function(filename, nrows=100, sep=",")
{
  function(init=FALSE)
  {
    if(init)
    {
      f <<- file(filename)
      open(f)
      return(NULL)
    }
    tryCatch(
      as.matrix(read.table(f, sep=sep, nrows=nrows)),
      error=function(e)
      {
        close(f)
        NULL
      })
  }
}

irls_incremental =
function(filename, chunksize, b, family=binomial, maxit=25, tol=1e-08)
{
  x     = NULL
  chunk = iterator(filename, nrows=chunksize) # a basic data file iterator
  for(j in 1:maxit)
  {
    k = 1                                     # Track the rows
    chunk(init=TRUE)                          # initialize the iterator
    A     = chunk()                           # get first chunk of model matrix
# Initialize first time through (after ascertaining ncol(A)):
    if(is.null(x)) x = rep(0,ncol(A))
    ATWA = matrix(0,ncol(A),ncol(A))
    ATWz = rep(0,ncol(A))
    while(!is.null(A))                        # iterate
    {
      eta    = A %*% x
      g      = family()$linkinv(eta)
      gprime = family()$mu.eta(eta)
      z      = eta + (b[k:(k+nrow(A)-1)] - g) / gprime
      k      = k + nrow(A)
      W      = as.vector(gprime^2 / family()$variance(g))
      ATWz   = ATWz + crossprod(A,W*z)
      ATWA   = ATWA + crossprod(A,W*A)
      A      = chunk()    # Next chunk
    }
    xold  = x
    C     = chol(ATWA, pivot=TRUE)
    if(attr(C,"rank")<ncol(C)) stop("Rank-deficiency detected.")
    p     = attr(C, "pivot")
    x     = backsolve(C,forwardsolve(t(C),ATWz[p]))[p]
    if(sqrt(crossprod(x-xold)) < tol) break
  }
  list(coefficients=x,iterations=j)
}
