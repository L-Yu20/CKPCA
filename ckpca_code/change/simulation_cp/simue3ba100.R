## corrected kernel pca change point for example 3
library(tictoc)
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(MASS)
library(moments)
library("wbs")
library(changepoint.geo)
library(RcppCNPy)
remove(list = ls())

tic()

#### Basic parameter settings
n = 800   ## number of samples
px = 100  ## experiment dimension
tau = 0.5  ## tau for TRR
times = 1000  ## number of experiments
an = floor(n ^ 0.5) ## betan
type = "ba"  # "ba" for balanced, "im" for imbalanced

## Paths
save_path = "E:/"

## Auto-generate file tags and paths
file_tag = paste0("e3", type, px)                    # e.g., "e3ba200"
## Loading path
npy_file = paste0(save_path, "ckpca_code/change/simu_data/", file_tag, ".npy")
## Saving path
file_name = paste0(save_path, "simu", file_tag, ".csv")

##Data loading
flat_vec = npyLoad(npy_file)
big_array1 = array(flat_vec, dim = c(n, px, times))

## balanced or imbalanced
if (type == "ba") {
  n1 = 100; n2 = 200; n3 = 300; n4 = 400
  n5 = 500; n6 = 600; n7 = 700; n8 = 800
} else if (type == "im") {
  n1 = 30;  n2 = 170; n3 = 350; n4 = 440
  n5 = 520; n6 = 630; n7 = 710; n8 = 800
} else {
  stop("type must be either 'ba' or 'im'")
}

## True change point
true_cps = c(n1, n2, n3, n4, n5, n6, n7, n8)

## Function to calculate the Rand Index
source("E:/ckpca_code/randi.R")

### Save the change point locations and the number 
### of change points for various methods
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
number8 = matrix(nrow = 1, ncol = times)
position8 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number9 = matrix(nrow = 1, ncol = times)
position9 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number10 = matrix(nrow = 1, ncol = times)
position10 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)
number11 = matrix(nrow = 1, ncol = times)
position11 = matrix(rep(0, 50 * times), nrow = times, ncol = 50)

for (kkk in 1:times) {
  cat("Running iteration:", kkk, "\n")
  set.seed(1111 + kkk)
  
  ## Get the data from the kkk-th experiment
  x1 = big_array1[ , , kkk]
  x = x1
  
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
  
  ## Mahalanobis matrix
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
  
  ### change point detection
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
  
  ## sbs
  hatq = dim(cc4)[1]
  if (hatq > 1) {
    kk8 = sbs.alg(cc4)$ecp
    number8[kkk] = length(kk8)
    if (number8[kkk] > 0) {
      position8[kkk, 1:(number8[kkk])] = kk8
    }
  } else {
    xx = t(cc4)
    w = wbs(xx)
    w.cpt = changepoints(w, penalty = "bic.penalty")
    kk8 = sort(w.cpt$cpt.ic$bic.penalty)
    number8[kkk] = length(kk8)
    if (number8[kkk] > 0) {
      position8[kkk, 1:(number8[kkk])] = kk8
    }
  }
  
  ## inspect method
  kk9 = inspect(t(x))$changepoints[, 1]
  number9[kkk] = length(kk9)
  if (number9[kkk] > 0) {
    position9[kkk, 1:(number9[kkk])] = kk9
  }
  
  ### geomcp
  ans2 = geomcp(x)
  kk10 = ans2@dist.cpts
  number10[kkk] = length(kk10)
  if (number10[kkk] > 0) {
    position10[kkk, 1:(number10[kkk])] = kk10
  }
  
  ### dcbs
  kk11 = dcbs.alg(t(x))$ecp
  number11[kkk] = length(kk11)
  if (number11[kkk] > 0) {
    position11[kkk, 1:(number11[kkk])] = kk11
  }
}

## Save the estimated number of change points (hatk), MSE and Rand Index
mk = rep(0, 11)     # estimated number of change points
me = rep(0, 11)     # root mean square error
ri = rep(0, 11)     # mean Rand index
ci1 = rep(0, 11)    # 2.5% quantile of Rand index
ci2 = rep(0, 11)    # 97.5% quantile of Rand index

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
me[1] = sqrt(mean((numberr1 - 7) ^ 2))
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
me[2] = sqrt(mean((numberr1 - 7) ^ 2))
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
me[3] = sqrt(mean((numberr1 - 7) ^ 2))
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
me[4] = sqrt(mean((numberr1 - 7) ^ 2))
res4 = randi(positionn1, numberr1, n, times, true_cps)
ri[4] = res4$mean
ci1[4] = res4$CI_lower
ci2[4] = res4$CI_upper

### ----- Method: sbs_c -----
numberr1 = number5
positionn1 = position5
mk[5] = mean(numberr1)
me[5] = sqrt(mean((numberr1 - 7) ^ 2))
res5 = randi(positionn1, numberr1, n, times, true_cps)
ri[5] = res5$mean
ci1[5] = res5$CI_lower
ci2[5] = res5$CI_upper

### ----- Method: sbs_p -----
numberr1 = number6
positionn1 = position6
mk[6] = mean(numberr1)
me[6] = sqrt(mean((numberr1 - 7) ^ 2))
res6 = randi(positionn1, numberr1, n, times, true_cps)
ri[6] = res6$mean
ci1[6] = res6$CI_lower
ci2[6] = res6$CI_upper

### ----- Method: sbs_m -----
numberr1 = number7
positionn1 = position7
mk[7] = mean(numberr1)
me[7] = sqrt(mean((numberr1 - 7) ^ 2))
res7 = randi(positionn1, numberr1, n, times, true_cps)
ri[7] = res7$mean
ci1[7] = res7$CI_lower
ci2[7] = res7$CI_upper

### ----- Method: sbs -----
numberr1 = number8
positionn1 = position8
mk[8] = mean(numberr1)
me[8] = sqrt(mean((numberr1 - 7) ^ 2))
res8 = randi(positionn1, numberr1, n, times, true_cps)
ri[8] = res8$mean
ci1[8] = res8$CI_lower
ci2[8] = res8$CI_upper

## inspect
numberr1 = number9
positionn1 = position9
mk[9] = mean(numberr1)
me[9] = sqrt(mean((numberr1 - 7) ^ 2))
res9 = randi(positionn1, numberr1, n, times, true_cps)
ri[9] = res9$mean
ci1[9] = res9$CI_lower
ci2[9] = res9$CI_upper

### ----- Method: Geomcp -----
numberr1 = number10
positionn1 = position10
mk[10] = mean(numberr1)
me[10] = sqrt(mean((numberr1 - 7) ^ 2))
res10 = randi(positionn1, numberr1, n, times, true_cps)
ri[10] = res10$mean
ci1[10] = res10$CI_lower
ci2[10] = res10$CI_upper

## DCBS
numberr1 = number11
positionn1 = position11
mk[11] = mean(numberr1)
me[11] = sqrt(mean((numberr1 - 7) ^ 2))
res11 = randi(positionn1, numberr1, n, times, true_cps)
ri[11] = res11$mean
ci1[11] = res11$CI_lower
ci2[11] = res11$CI_upper

### Save all results
res = cbind(mk, me, ri, ci1, ci2)
write.csv(res, file_name)
res

toc()
