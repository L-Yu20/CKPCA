## plot
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(MASS)
library(moments)
library("wbs")
n = 200
px = 30 
tau = 0.5
an = floor(n^0.5)
af = 3  

## Balanced dataset
n1 = n/2
n2 = n

c2 = n2 - n1
c1 = n1

set.seed(1111)

## data setting
## data setting
# rho = 0.5
mu1 = rep(0, px)
sigma1 = diag(1,px)
for (i in 1:4) {
  sigma1[i, i] = 100
}
for (i in 5:px) {
  sigma1[i, i] = 0.01
}


x = matrix(nrow = n, ncol = px)
x[(1):(n1), ] = mvrnorm(n = c1, mu1, sigma1)
x[(n1 + 1):(n2), ] = mvrnorm(n = c2, mu1, sigma1)

## set change point
k = px
u1 = matrix(rep(0, k * c2), nrow = c2, ncol = k)
u1[, 5] = 10
x[1:(n1), 1:k] = x[1:(n1), 1:k]
x[(n1 + 1):(n2), 1:k] = x[(n1 + 1):(n2), 1:k] + u1

## cpca dimension reduction
Mx = cov(x) / n * (n - 1)
sigma = matrix(rep(0, px * px), nrow = px, ncol = px)
zs = floor(n / an)
for (j in 1:(zs - 1)) {
  sigma = sigma + cov(x[((j - 1) * an + 1):(j * an), ])
}
sigma = sigma + cov(x[(((zs - 1) * an) + 1):n, ])
sigma = sigma / zs
Dx = Mx - sigma
spec = eigen(Dx)
evalue = rep(0, px + 1)
va = abs(spec$values)
evalue[1:px] = va
evector = spec$vectors
ccn = log(log(n)) * sqrt(px / n) / 5
evalueplus = evalue + ccn
lambda1 = evalueplus[2:(px + 1)] / evalueplus[1:px]
lambda2 = which(lambda1[1:(px - 1)] < 0.5)
aaa = min(lambda1[1:(px - 1)])
if (aaa > 0.5) {
  hatq = 1
} else {
  hatq = max(lambda2)
}
B = evector[, 1:(hatq)]
cc1 = t(B) %*% t(x)

### pca method
Dx = Mx
spec = eigen(Dx)
evalue = rep(0, px + 1)
va = abs(spec$values)
evalue[1:px] = va
evector = spec$vectors
total = sum(evalue)
to = 0.95 * total
cto = matrix(rep(0, px), nrow = px, ncol = 1)
for (i in 1:px) {
  for (j in 1:i) {
    cto[i] = cto[i] + evalue[j]
  }
}
lo = which(cto > to)
hatq = min(lo)
B = evector[, 1:(hatq)]
cc2 = t(B) %*% t(x)

## without dimension reduction
cc4 = t(x)

kk1 = c(n1)

## plot the figure
plot(x[, 1], cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in mean", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')

plot(cc1[1, 1:n], ylab = "f(X)", cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in mean (CPCA)", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')

plot(cc2[1, 1:n], ylab = "f(X)", cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in mean (PCA)", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')


