x=c(0.5,0.75,1,1.25,1.5,1.75,1.75,2,2.25,2.5,2.75,3,3.25,3.5,4,4.25,4.5,4.75,5,5.5)
y=c(0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1)

jpeg(file="logistic.jpg", quality=100, width=500, height=500)
plot(x,y,xlab="Hours of study", ylab="Probability of passing", xlim=c(0,6), ylim=c(-0.25,1.25))
grid()
abline(coef(lm.fit(cbind(1, x), y)), col=4, lwd=2, lty=2)

l = coef(glm.fit(cbind(1, x), y, family=binomial()))
xx = seq(0, 6, length.out=100)
p = 1/(1 + exp(-(l[2]*xx + l[1])))

lines(xx, p, col=2, lwd=2)
dev.off()
