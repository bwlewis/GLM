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
