# Let's split the matrix products up into a parallel setting...
require("doSNOW")

# Computers to participate in cluster
nodes = c("localhost","localhost")
# Create a SNOW cluster and register it with foreach
cl = makeCluster(nodes, type="SOCK")
registerDoSNOW(cl)

source("parallel.R")
n = 17
np= length(nodes)

X = matrix(rnorm(n*5),n)
y = rep(1,n)

ans = parallel_glm.fit(y=y,x=X,np=np)

stopCluster(cl)
