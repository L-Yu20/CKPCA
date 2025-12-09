##macroeconomic data
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")

## Setting seed
set.seed(1111)

dataset1 = read.csv(file = 'E:\\ckpca_code\\dataset\\macro\\1.csv', header = F)
dataset4 = read.csv(file = 'E:\\ckpca_code\\dataset\\macro\\4.csv', header = F)
dataset5 = read.csv(file = 'E:\\ckpca_code\\dataset\\macro\\5.csv', header = F)
dataset6 = read.csv(file = 'E:\\ckpca_code\\dataset\\macro\\6.csv', header = F)

##Data processing
x1 = dataset1[, 2:4]
x4 = dataset4[, 2:15]
x5 = dataset5[, 2:50]
x6 = dataset6[, 2:27]
x1 = as.matrix(x1)
x4 = as.matrix(x4)
x5 = as.matrix(x5)
x6 = as.matrix(x6)
x4 = log(x4)
x5 = log(x5)
x6 = log(x6)
n = 364

xx5 = matrix(nrow = n, ncol = 49)
for (i in 1:49) {
  for (j in 1:n) {
    xx5[j, i] = x5[(j + 1), i] - x5[j, i]
  }
}

xx6 = matrix(nrow = n + 1, ncol = 26)
for (i in 1:26) {
  for (j in 1:(n + 1)) {
    xx6[j, i] = x6[(j + 1), i] - x6[j, i]
  }
}

xxx6 = matrix(nrow = n, ncol = 26)
for (i in 1:26) {
  for (j in 1:n) {
    xxx6[j, i] = xx6[(j + 1), i] - xx6[j, i]
  }
}

x1 = x1[1:n, ]
x4 = x4[1:n, ]
x5 = xx5[1:n, ]
x6 = xx6[1:n, ]
x = cbind(x1, x4, x5, x6)

dimx = dim(x)
n = dimx[1]
px = dimx[2]
x = as.matrix(x)

##using thumb method to choose bandwidth
vv = var(x)
v0 = 0
for (i in 1:px) {
  v0 = v0 + vv[i, i]
}
v0 = v0 / px
m = 2
sig = sqrt(m * px * v0)
an = floor(sqrt(n))

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
output1 = e.divisive(t(cc1), R = 499, alpha = 1)
kk1 = output1$estimates

##plot the figure
cc1 = ts(t(cc1), frequency = 12, start = 1992 + 2 / 12)
plot(cc1[, 1], ann = F, xaxt = "n")
axis(1, at = c(1995, 2000, 2005, 2010, 2015, 2020),
     labels = c("1995", "2000", "2005", "2010", "2015", "2020"))
title(xlab = "Time")
title("The U.S. Macroeconomic data")
kk2 = kk1
kk1 = kk1 / 12 + 1992 + 2 / 12
kk = length(kk1)
kk1 = kk1[2:(kk - 1)]
abline(v = kk1, lty = 2, col = 'red')
