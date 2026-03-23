###corrected kernel pca clustering for kmeans
library(MASS)
library(cluster)
library('mvtnorm')
library(tictoc)
times = 1000
v = 0
px = 10
n = 600

##Balanced dataset
nn1 = 200
nn2 = 400
nn3 = 600

##Imbalanced dataset
# nn1 = 300
# nn2 = 500
# nn3 = 600

d1 = nn1
d2 = nn2 - nn1
d3 = nn3 - nn2
n0 = n / 3
repeat1 = 30
kinds = 3
an = floor(n ^ 0.5)
ri0 = matrix(ncol = 1, nrow = times)
ri1 = matrix(ncol = 1, nrow = times)
ri2 = matrix(ncol = 1, nrow = times)
ri3 = matrix(ncol = 1, nrow = times)
cluster1 = matrix(nrow = n, ncol = repeat1)

tic()

for (time in 1:times){
  cat("Running iteration:", time, "\n")
  set.seed(1111 + time)
  ##data setting
  mu1 = rep(0, px)
  sigma1 = diag(rep(1, px))
  x = mvrnorm(n, mu1, sigma1)
  t1 = runif(n = d1, min = 0, max = 1)
  t2 = runif(n = d2, min = 2, max = 3)
  t3 = runif(n = d3, min = 4, max = 5)
  t = c(t1, t2, t3)
  for (i in 1:n) {
    x[i, ] = x[i, ] / sqrt(sum(x[i, ] ^ 2))
  }
  for (i in 1:n) {
    x[i, ] = t[i] * x[i, ]
  }
  
  ##bandwidth
  vv = var(x)
  v0 = 0
  for (i in 1:px) {
    v0 = v0 + vv[i, i]
  }
  v0 = v0 / px
  m = 0.8
  sig = sqrt(m * px * v0)
  
  ##calculate k
  k = matrix(rep(0, n * n), nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      aa = x[i, ] - x[j, ]
      bb = sum(aa ^ 2)
      k[i, j] = exp(-bb / (2 * sig ^ 2))
    }
  }
  op = matrix(rep(1, n * n), nrow = n, ncol = n)
  l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
  l2 = diag(1, n)
  L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
  ke0 = (L) %*% k
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
  cc = k %*% B1
  cc = t(cc)
  ###kmeans
  cluster2 = kmeans(t(cc), kinds)$cluster
  
  ###
  rest1 = kmeans(x, kinds)
  cluster1[, 1] = rest1$cluster
  
  ##iteration
  for (order in 2:repeat1){
    sigma = 0
    b1 = which(cluster1[, order - 1] == 1)  
    b2 = which(cluster1[, order - 1] == 2)
    b3 = which(cluster1[, order - 1] == 3)
    cc1 = length(b1)
    cc2 = length(b2)
    cc3 = length(b3)
    ccc = c(cc1, cc2, cc3)
    g1 = matrix(rep(0, n * ccc[1]), nrow = n, ncol = ccc[1])
    g2 = matrix(rep(0, n * ccc[2]), nrow = n, ncol = ccc[2])
    g3 = matrix(rep(0, n * ccc[3]), nrow = n, ncol = ccc[3])
    h1 = matrix(rep(0, n * ccc[1]), nrow = n, ncol = ccc[1])
    h2 = matrix(rep(0, n * ccc[2]), nrow = n, ncol = ccc[2])
    h3 = matrix(rep(0, n * ccc[3]), nrow = n, ncol = ccc[3])
    for (j in 1:ccc[1]){
      g1[b1[j], j] = 1
      h1[b1[j], ] = 1
    }
    for (j in 1:ccc[2]){
      g2[b2[j], j] = 1
      h2[b2[j], ] = 1
    }
    for (j in 1:ccc[3]){
      g3[b3[j], j] = 1
      h3[b3[j], ] = 1
    }
    n1 = cc1
    n2 = cc2
    n3 = cc3
    U = matrix(rep(0, n * n), nrow = n, ncol = n)
    U = U + (g1 - h1 / n1) %*% t(g1 - h1 / n1) / (n - 3)
    U = U + (g2 - h2 / n2) %*% t(g2 - h2 / n2) / (n - 3)
    U = U + (g3 - h3 / n3) %*% t(g3 - h3 / n3) / (n - 3)
    kk0 = (L - U) %*% k
    ke0 = L %*% k
    M = L - U
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
    if (aaa > 0.5){
      hatq = 1
    } else{
      hatq = max(lambda2[which(lambda2 < px)])
    }
    B = evector[, 1:(hatq)]
    cc = k %*% B
    cc = t(cc)  
    rest1 = kmeans(t(cc), kinds)
    cluster1[, order] = rest1$cluster
  }
  
  
  
  #rand index
  rri = matrix(nrow = 1, ncol = repeat1)
  for (kk in 2:repeat1){
    indexx = matrix(nrow = 1, ncol = n)
    indexx0 = matrix(nrow = 1, ncol = n)
    indexx = cluster1[, kk]
    indexx0 = cluster1[, kk - 1]
    all0 = n * (n - 1) / 2
    tp = 0
    tn = 0
    for (jj in 2:n){
      for (hh in 1:(jj - 1)){
        if (indexx[jj] == indexx[hh] & indexx0[jj] == indexx0[hh])
          tp = tp + 1 else if (indexx[jj] != indexx[hh] & indexx0[jj] != indexx0[hh])
            tn = tn + 1
      }
    }
    rri[kk] = (tp + tn) / all0
  }
  for (i in 2:(repeat1 - 1)){
    if (rri[i] > 0.999 & rri[i + 1] > 0.999){
      break
    }
  }
  st = i
  
  ##our method
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0[1:nn1] = 1
  index0[(nn1 + 1):nn2] = 2
  index0[(nn2 + 1):nn3] = 3
  randindex = matrix(nrow = 1, ncol = repeat1)
  all0 = n * (n - 1) / 2
  h = st
  index = cluster1[, h]
  tp = 0
  tn = 0
  for (jj in 2:n){
    for (hh in 1:(jj - 1)){
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri0[time] = (tp + tn) / all0
  
  ##dim-1
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0[1:nn1] = 1
  index0[(nn1 + 1):nn2] = 2
  index0[(nn2 + 1):nn3] = 3
  index = cluster2
  tp = 0
  tn = 0
  for (jj in 2:n){
    for (hh in 1:(jj - 1)){
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri1[time] = (tp + tn) / all0
  
  ##kmeans(x)
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0[1:nn1] = 1
  index0[(nn1 + 1):nn2] = 2
  index0[(nn2 + 1):nn3] = 3
  randindex = matrix(nrow = 1, ncol = repeat1)
  all0 = n * (n - 1) / 2
  index = cluster1[, 1]
  tp = 0
  tn = 0
  for (jj in 2:n){
    for (hh in 1:(jj - 1)){
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri2[time] = (tp + tn) / all0
  
}


##our method
mri0 = mean(ri0)
sd0 = sd(ri0)
lower0 = quantile(ri0, 0.025)
upper0 = quantile(ri0, 0.975)
##dim-1
mri1 = mean(ri1)
sd1 = sd(ri1)
lower1 = quantile(ri1, 0.025)
upper1 = quantile(ri1, 0.975)
##kmeans(x)
mri2 = mean(ri2)
sd2 = sd(ri2)
lower2 = quantile(ri2, 0.025)
upper2 = quantile(ri2, 0.975)

##
result = rbind(
  c(mri0, sd0, lower0, upper0),
  c(mri1, sd1, lower1, upper1),
  c(mri2, sd2, lower2, upper2)
)

result
write.csv(result, file = paste0("E:/kmeans_px", px, ".csv"))



toc()
