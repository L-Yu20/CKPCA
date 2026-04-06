##genetics data
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(MASS)
library(moments)
library("wbs")

## Setting seed
set.seed(1111)

data(ACGH, package = "ecp")
x = ACGH$data
dimx = dim(x)
n = dimx[1]
px = dimx[2]
tau = 0.5
an = floor(n^0.5)

##bandwidth
vv = var(x)
v0 = 0
for (i in 1:px) {
  v0 = v0 + vv[i, i]
}
v0 = v0 / px
m = 0.8
sig = sqrt(m * px * v0)

## Kernel matrix
dist_mat = as.matrix(dist(x)) ^ 2
k = exp(-dist_mat / (2 * sig ^ 2))

##dimension reduction
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
  hatq = max(lambda2)
}
B = evector[, 1:(hatq)]
ab = sqrt(sum(B^2))
cc1 = k %*% B

##the lower dimensional data
cc1 = t(cc1)

##ecp30
output1 = e.divisive(t(cc1), R = 499, alpha = 1, min.size = 30)
kk1 = output1$estimates

##plot the figure
kk = length(kk1)
kk1 = kk1[2:(kk - 1)]

par(mar = c(5, 4.7, 2, 2) + 0.1)

plot(cc1[1, 1:2215], main = "The E-Divisive method after dimension reduction",
     xlab = "Index",
     ylab = expression(f[1*n](X[i])) )
abline(v = kk1, lty = 2, col = 'red')
