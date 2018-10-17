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

jpeg(file="linear.jpg", quality=100, width=500, height=500)
plot(x, y, asp=1, xlab="x", ylab="y", bty="n", pch=19)
grid()
abline(0, b, lwd=2, col=4)
dev.off()

jpeg(file="not_linear.jpg", quality=100, width=500, height=500)
plot(x, y, asp=1, xlab="x", ylab="y", bty="n", pch=19)
grid()
xx=seq(from=-2, to=max(x), length.out=100)
e=coef(lm.fit(cbind(1,exp(x)), y))
lines(xx, e[1] + e[2]*exp(xx), type='l', col=2, lwd=2)
dev.off()

jpeg(file="whichlinear.jpg", quality=100, width=500, height=500)
plot(x, y, asp=1, xlab="x", ylab="y", bty="n", pch=19)
grid()
abline(0, b, lwd=2, col=4)
abline(-0.5, 0.5*b, lwd=2, col=3)
abline(0, -0.1*b, lwd=2, col=2)
dev.off()


jpeg(file="ols.jpg", quality=100, width=500, height=500)
plot(x, y, asp=1,  xlab="x", ylab="y", bty="n", pch=19)
grid()
#legend("bottomright", legend=c("Ordinary least squares fit"), fill=c(1), bty="n", cex=2)
abline(0, b, lwd=2, col=4, lty=1)
py = b * x
for(i in 1:n) lines(x=c(x[i], x[i]), y=c(y[i], py[i]), col=1, lwd=2, lty=2)
dev.off()
