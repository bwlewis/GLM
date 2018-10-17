tls = function(X, y)
{
  p = ncol(X)
  s = svd(cbind(X, y))
  -s$v[1:p, p+1]  / s$v[p+1, p+1]
}

set.seed(3)
n = 10
x = matrix(rnorm(n), ncol=1)
y = 0.2 * x + runif(n, -1)
b = coef(lm.fit(x, y)) # ols fit
bt = tls(x, y)

#pdf("~/Desktop/ch10_tls.pdf", width = 6, height = 6)
opar = par(mar = c(0,0,0,0))#), mfrow=c(1,2))

plot(x, y, asp=1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", pch=19)
#legend("bottomright", legend=c("Ordinary least squares fit"), fill=c(1), bty="n", cex=2)
abline(0, b, lwd=2, col=1, lty=2)
abline(0, bt, lwd=3, col=1, lty=1)
py = b * x
for(i in 1:n) lines(x=c(x[i], x[i]), y=c(y[i], py[i]), col=1, lwd=2, lty=3)

#plot(x, y, asp=1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", pch=19)
#legend("bottomright", legend=c("Total least squares fit"), fill=c(4), bty="n", cex=2)
#abline(0, b, lwd=1, lty=2)
#abline(0, bt, lwd=3, col=1, lty=1)

px = sqrt(1/(1+bt^2))
py = bt * px
theta = atan(py/px)
c = cos(theta)
s = sin(theta)
P = matrix(c(c, s, -s, c), 2, byrow=TRUE)
G = matrix(c(c, -s, s, c), 2, byrow=TRUE)
z = cbind(x, y) %*% t(P)
z = cbind(z[, 1], 0) %*% t(G)
for(i in 1:n) lines(c(x[i], z[i,1]), c(y[i], z[i,2]), col=1, lwd=2, lty=1)
#dev.off()
par(opar)
