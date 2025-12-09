## plot
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(MASS)
library(moments)
library("wbs")
n = 800
px = 200 ## px = 100 or 200
tau = 0.5
an = floor(n^0.5)
af = 3  

## Balanced dataset
n1 = 100
n2 = 200
n3 = 300
n4 = 400
n5 = 500
n6 = 600
n7 = 700
n8 = 800

c2 = n2 - n1
c4 = n4 - n3
c6 = n6 - n5
c8 = n8 - n7
c1 = n1
c3 = n3 - n2
c5 = n5 - n4
c7 = n7 - n6

set.seed(1111)

## data setting
rho = 0.5
mu1 = rep(0, px)
sigma1 = matrix(nrow = px, ncol = px)
for (i in 1:px) {
  for (j in 1:px) {
    sigma1[i, j] = rho^(abs(i - j))  ## case2
  }
}
for (i in 1:px) {
  sigma1[i, i] = 2.5
}
x = matrix(nrow = n, ncol = px)
tt1 = px * c2
tt2 = px * c4
tt3 = px * c6
tt4 = px * c8
t1 = runif(n = tt1, min = -3, max = 3)
t2 = runif(n = tt2, min = -3, max = 3)
t3 = runif(n = tt3, min = -3, max = 3)
t4 = runif(n = tt4, min = -3, max = 3)
x[(1):(n1), ] = mvrnorm(n = c1, mu1, sigma1)
x[(n1 + 1):(n2), ] = matrix(t1, nrow = c2, ncol = px)
x[(n2 + 1):(n3), ] = mvrnorm(n = c3, mu1, sigma1)
x[(n3 + 1):(n4), ] = matrix(t2, nrow = c4, ncol = px)
x[(n4 + 1):(n5), ] = mvrnorm(n = c5, mu1, sigma1)
x[(n5 + 1):(n6), ] = matrix(t3, nrow = c6, ncol = px)
x[(n6 + 1):(n7), ] = mvrnorm(n = c7, mu1, sigma1)
x[(n7 + 1):(n8), ] = matrix(t4, nrow = c8, ncol = px)

## using thumb method to choose bandwidth
vv = var(x)
v0 = 0
for (i in 1:px) {
  v0 = v0 + vv[i, i]
}
v0 = v0 / px
m = 0.8
sig = sqrt(m * px * v0)
sig0 = sig

## ckpca
k = matrix(rep(0, n * n), nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    aa = x[i, ] - x[j, ]
    bb = sum(aa^2)
    k[i, j] = exp(-bb / (2 * sig^2))
  }
}
op = matrix(rep(1, n * n), nrow = n, ncol = n)
zs = floor(n / an)
a0 = array(dim = c(n, an, zs))
for (i in 1:zs) {
  b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
  a1 = diag(rep(1, an))
  b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
  b0 = rbind(b1, a1, b2)
  a0[, , i] = b0
}
c0 = array(dim = c(n, an, zs))
for (i in 1:zs) {
  b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
  a1 = matrix(rep(1, an * an), nrow = an, ncol = an)
  b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
  b0 = rbind(b1, a1, b2)
  c0[, , i] = b0
}
d0 = array(dim = c(n, an, zs))
for (i in 1:zs) {
  d0[, , i] = a0[, , i] - c0[, , i] / an
  d0[, , i] = d0[, , i]
}
U = matrix(rep(0, n * n), nrow = n, ncol = n)
for (i in 1:zs) {
  U = U + d0[, , i] %*% t(d0[, , i]) / (an - 1)
}
bn = n - zs * an
b1 = matrix(rep(0, (n - bn) * an), nrow = n - bn, ncol = bn)
a1 = diag(rep(1, bn))
bb0 = rbind(b1, a1)
b1 = matrix(rep(0, (n - bn) * bn), nrow = n - bn, ncol = bn)
a1 = matrix(rep(1, bn * bn), nrow = bn, ncol = bn)
bb1 = rbind(b1, a1)
dd0 = (bb0 - bb1 / bn)
U = U + dd0 %*% t(dd0) / (bn - 1)
U = U / zs
l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
l2 = diag(1, n)
L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
M = (L - U)
kk0 = (L - U) %*% k
spec = eigen(kk0)
evalue = rep(0, n + 1)
va = abs(spec$values)
evalue[1:n] = va
evector = spec$vectors
ccn = log(log(n)) * sqrt(1 / n) / 5
evalueplus = evalue + ccn
lambda1 = evalueplus[2:(px + 1)] / evalueplus[1:px]
lambda2 = which(lambda1[1:(px - 1)] < 0.5)
aaa = min(lambda1[1:(px - 1)])
if (aaa > 0.5) {
  hatq = 1
} else {
  hatq = max(lambda2[which(lambda2 < px)])
}
B = evector[, 1:(hatq)]
ab = sqrt(sum(B^2))
cc1 = k %*% B
cc1 = t(cc1)

## kernel pca
kcc = matrix(rep(0, n * n), nrow = n, ncol = n)
for (i in 1:n) {
  for (j in 1:n) {
    aa = x[i, ] - x[j, ]
    bb = sum(aa^2)
    kcc[i, j] = exp(-bb / (2 * sig0^2))
  }
}
ke0 = (L) %*% kcc
spec1 = eigen(ke0)
evalue1 = rep(0, n + 1)
va1 = abs(spec1$values)
evalue1[1:n] = va1
evector1 = spec1$vectors
total = sum(evalue1)
to = 0.95 * total
cto = matrix(rep(0, n), nrow = n, ncol = 1)
for (i in 1:n) {
  for (j in 1:i) {
    cto[i] = cto[i] + evalue1[j]
  }
}
lo = which(cto > to)
hatq = min(lo)
B1 = evector1[, 1:(hatq)]
cc2 = kcc %*% B1
cc2 = t(cc2)

## without dimension reduction
cc4 = t(x)


kk1 = c(n1, n2, n3, n4, n5, n6, n7)

## plot the figure
plot(x[, 1], cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in distribution", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')

plot(cc1[1, 1:n], ylab = "f(X)", cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in distribution (CKPCA)", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')

plot(cc2[1, 1:n], ylab = "f(X)", cex.axis = 1.5, cex.lab = 1.5)
title(main = "Changes in distribution (KPCA)", cex.main = 1.5)
abline(v = kk1, lty = 2, col = 'red')


