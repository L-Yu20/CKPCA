##corrected kernel pca change point for example2
library(tictoc)
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(MASS)
library(moments)
library("wbs")
library(RcppCNPy)
remove(list = ls())

tic()

#### Basic parameter settings
n = 800        ## number of samples  
px = 100       ## experiment dimension  
tau = 0.5      ## tau for TRR
case = 1
times = 1000   ## number of experiments  
an = floor(n ^ 0.5)  ## betan  
type = "ba"    # "ba" for balanced, "im" for imbalanced  

## Dynamic path and file name
save_path = "E:/"
file_tag = paste0("e1", "c", case, type, px)     # e.g., "e1c1ba100"
## Loading path
npy_file = paste0(save_path, "ckpca_code/change/simu_data/", file_tag, ".npy")
## Saving path
file_name = paste0(save_path, "simu", file_tag, ".csv")

##Data loading
flat_vec = npyLoad(npy_file)
big_array1 = array(flat_vec, dim = c(n, px, times))


##balanced or imbalanced
if (type == "ba") {
  n1 = 100; n2 = 200; n3 = 300; n4 = 400
  n5 = 500; n6 = 600; n7 = 700; n8 = 800
} else if (type == "im") {
  n1 = 30;  n2 = 170; n3 = 350; n4 = 440
  n5 = 520; n6 = 630; n7 = 710; n8 = 800
} else {
  stop("type must be either 'ba' or 'im'")
}

##True change point
true_cps = c(n1, n2, n3, n4, n5, n6, n7, n8)

##Function to calculate the Rand Index
source("E:/ckpca_code/randi.R")

###Save the change point locations and the number 
###of change points for various methods
number1 = matrix(nrow = 1, ncol = times)
position1 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number2 = matrix(nrow = 1, ncol = times)
position2 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number3 = matrix(nrow = 1, ncol = times)
position3 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number4 = matrix(nrow = 1, ncol = times)
position4 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number5 = matrix(nrow = 1, ncol = times)
position5 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number6 = matrix(nrow = 1, ncol = times)
position6 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number7 = matrix(nrow = 1, ncol = times)
position7 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)

for (kkk in 1:times) {
  cat("Running iteration:", kkk, "\n")
  set.seed(1111 + kkk)
  ##Get the data from the kkk-th experiment
  x1 = big_array1[, , kkk]
  x = x1
  
  ##using thumb method to choose bandwidth
  vv = var(x)
  v0 = 0
  for (i in 1:px) {
    v0 = v0 + vv[i, i]
  }
  v0 = v0 / px
  m = 0.8
  sig = sqrt(m * px * v0)
  
  
  ###ckpca
  ## Kernel matrix
  dist_mat = as.matrix(dist(x)) ^ 2
  k = exp(-dist_mat / (2 * sig ^ 2))
  ## Compute matrices L and U
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
  if (bn > 1) {
    b1 = matrix(rep(0, (n - bn) * an), nrow = n - bn, ncol = bn)
    a1 = diag(rep(1, bn))
    bb0 = rbind(b1, a1)
    b1 = matrix(rep(0, (n - bn) * bn), nrow = n - bn, ncol = bn)
    a1 = matrix(rep(1, bn * bn), nrow = bn, ncol = bn)
    bb1 = rbind(b1, a1)
    dd0 = (bb0 - bb1 / bn)
    U = U + dd0 %*% t(dd0) / (bn - 1)
  } else {
    # Skip last segment when it's too short to compute variance (bn = 1)
    warning("Last segment size is 1; skipping to avoid division by zero.")
  }
  ### U
  U = U / zs
  l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
  l2 = diag(1, n)
  ### L
  L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
  M = (L - U)
  ### kn = (L - U) K
  kk0 = (L - U) %*% k
  ## Eigen decomposition
  spec = eigen(kk0)
  evalue = rep(0, n + 1)
  va = abs(spec$values)
  evalue[1:n] = va
  evector = spec$vectors
  ## TRR
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
  ab = sqrt(sum(B ^ 2))
  ## Projection
  cc1 = k %*% B
  cc1 = t(cc1)
  ## kernel PCA
  ke0 = (L) %*% k
  spec1 = eigen(ke0)
  evalue1 = rep(0, n + 1)
  va1 = abs(spec$values)
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
  cc2 = k %*% B1
  cc2 = t(cc2)
  ### Mahalanobis matrix
  Tx = matrix(nrow = n, ncol = n * px)
  a = matrix(nrow = n * n, ncol = px)
  Mx = matrix(nrow = px, ncol = px)
  for (i in 1:px) {
    tt = x[, i] %*% matrix(rep(1, n), nrow = 1, ncol = n) - matrix(rep(1, n), nrow = n, ncol = 1) %*% t(x[, i])
    Tx = matrix(tt, nrow = 1)
    a[, i] = Tx
  }
  Mx = t(a) %*% a / (n * (n - 1))
  sigma = matrix(rep(0, px * px), nrow = px, ncol = px)
  zs = floor(n / an)
  for (j in 1:(zs - 1)) {
    sigma = sigma + cov(x[((j - 1) * an + 1):(j * an), ])
  }
  sigma = sigma + cov(x[(((zs - 1) * an) + 1):n, ])
  sigma = sigma / zs
  Dx = Mx - 2 * sigma
  spec = eigen(Dx)
  evalue = rep(0, px + 1)
  va = abs(spec$values)
  evalue[1:px] = va
  evector = spec$vectors
  ccn = log(log(n)) * sqrt(px / n) / 2
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
  cc3 = t(B) %*% t(x)
  ### without dimension reduction
  cc4 = t(x)
  
  ## E-Divisive after dimension reduction
  
  ## ecp_c
  output1 = e.divisive(t(cc1), R = 499, alpha = 1)
  kk1 = output1$estimates
  number1[kkk] = length(kk1)
  position1[kkk, 1:(number1[kkk])] = kk1
  
  ## ecp_p
  output2 = e.divisive(t(cc2), R = 499, alpha = 1)
  kk2 = output2$estimates
  number2[kkk] = length(kk2)
  position2[kkk, 1:(number2[kkk])] = kk2
  
  ## ecp_m
  output3 = e.divisive(t(cc3), R = 499, alpha = 1)
  kk3 = output3$estimates
  number3[kkk] = length(kk3)
  position3[kkk, 1:(number3[kkk])] = kk3
  
  ## ecp
  output4 = e.divisive(t(cc4), R = 499, alpha = 1)
  kk4 = output4$estimates
  number4[kkk] = length(kk4)
  position4[kkk, 1:(number4[kkk])] = kk4
  
  ## sbs_c
  hatq = dim(cc1)[1]
  if (hatq > 1) {
    kk5 = sbs.alg(cc1)$ecp
    number5[kkk] = length(kk5)
    if (number5[kkk] > 0) {
      position5[kkk, 1:number5[kkk]] = kk5
    }
  } else {
    xx = t(cc1)
    w = wbs(xx)
    w.cpt = changepoints(w, penalty = "bic.penalty")
    kk5 = sort(w.cpt$cpt.ic$bic.penalty)
    number5[kkk] = length(kk5)
    if (number5[kkk] > 0) {
      position5[kkk, 1:number5[kkk]] = kk5
    }
  }
  
  ## sbs_p
  hatq = dim(cc2)[1]
  if (hatq > 1) {
    kk6 = sbs.alg(cc2)$ecp
    number6[kkk] = length(kk6)
    if (number6[kkk] > 0) {
      position6[kkk, 1:number6[kkk]] = kk6
    }
  } else {
    xx = t(cc2)
    w = wbs(xx)
    w.cpt = changepoints(w, penalty = "bic.penalty")
    kk6 = sort(w.cpt$cpt.ic$bic.penalty)
    number6[kkk] = length(kk6)
    if (number6[kkk] > 0) {
      position6[kkk, 1:number6[kkk]] = kk6
    }
  }
  
  ## sbs_m
  hatq = dim(cc3)[1]
  if (hatq > 1) {
    kk7 = sbs.alg(cc3)$ecp
    number7[kkk] = length(kk7)
    if (number7[kkk] > 0) {
      position7[kkk, 1:number7[kkk]] = kk7
    }
  } else {
    xx = t(cc3)
    w = wbs(xx)
    w.cpt = changepoints(w, penalty = "bic.penalty")
    kk7 = sort(w.cpt$cpt.ic$bic.penalty)
    number7[kkk] = length(kk7)
    if (number7[kkk] > 0) {
      position7[kkk, 1:number7[kkk]] = kk7
    }
  }
  
  
}





## Save the estimated number of change points (hatk), MSE and Rand Index
mk = rep(0, 7)     # estimated number of change points
me = rep(0, 7)     # root mean square error
ri = rep(0, 7)     # mean Rand index
ci1 = rep(0, 7)    # 2.5% quantile of Rand index
ci2 = rep(0, 7)    # 97.5% quantile of Rand index

### ----- Method: ecp_c -----
numberr1 = number1
positionn1 = position1
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):he1] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[1] = mean(numberr1)
me[1] = sqrt(mean((numberr1 - 7)^2))
res1 = randi(positionn1, numberr1, n, times, true_cps)
ri[1] = res1$mean
ci1[1] = res1$CI_lower
ci2[1] = res1$CI_upper

### ----- Method: ecp_p -----
numberr1 = number2
positionn1 = position2
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):he1] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[2] = mean(numberr1)
me[2] = sqrt(mean((numberr1 - 7)^2))
res2 = randi(positionn1, numberr1, n, times, true_cps)
ri[2] = res2$mean
ci1[2] = res2$CI_lower
ci2[2] = res2$CI_upper

### ----- Method: ecp_m -----
numberr1 = number3
positionn1 = position3
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):he1] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[3] = mean(numberr1)
me[3] = sqrt(mean((numberr1 - 7)^2))
res3 = randi(positionn1, numberr1, n, times, true_cps)
ri[3] = res3$mean
ci1[3] = res3$CI_lower
ci2[3] = res3$CI_upper

### ----- Method: ecp -----
numberr1 = number4
positionn1 = position4
for (k in 1:times) {
  kk = positionn1[k, ]
  he1 = length(which(kk > 0))
  if (he1 == 2) {
    positionn1[k, ] = 0
  } else {
    positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
    positionn1[k, (he1 - 1):he1] = c(0, 0)
  }
}
numberr1 = numberr1 - 2
mk[4] = mean(numberr1)
me[4] = sqrt(mean((numberr1 - 7)^2))
res4 = randi(positionn1, numberr1, n, times, true_cps)
ri[4] = res4$mean
ci1[4] = res4$CI_lower
ci2[4] = res4$CI_upper

### ----- Method: sbs_c -----
numberr1 = number5
positionn1 = position5
mk[5] = mean(numberr1)
me[5] = sqrt(mean((numberr1 - 7)^2))
res5 = randi(positionn1, numberr1, n, times, true_cps)
ri[5] = res5$mean
ci1[5] = res5$CI_lower
ci2[5] = res5$CI_upper

### ----- Method: sbs_p -----
numberr1 = number6
positionn1 = position6
mk[6] = mean(numberr1)
me[6] = sqrt(mean((numberr1 - 7)^2))
res6 = randi(positionn1, numberr1, n, times, true_cps)
ri[6] = res6$mean
ci1[6] = res6$CI_lower
ci2[6] = res6$CI_upper

### ----- Method: sbs_m -----
numberr1 = number7
positionn1 = position7
mk[7] = mean(numberr1)
me[7] = sqrt(mean((numberr1 - 7)^2))
res7 = randi(positionn1, numberr1, n, times, true_cps)
ri[7] = res7$mean
ci1[7] = res7$CI_lower
ci2[7] = res7$CI_upper

### Save all results
res = cbind(mk, me, ri, ci1, ci2)
write.csv(res, file_name)
res



toc()
