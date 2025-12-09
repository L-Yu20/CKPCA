###kernel pca change point
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
n = 800       ## number of samples
px = 200      ## experiment dimension
tau = 0.5     ## tau for TRR
times = 1000  ## number of experiments
an = floor(n^0.5)  ## betan
af = 4        ## a = 4 or 6
hh = c(0.4, 0.8, 1.2, 1.6, 2)  ## candidate set for m
hh1 = length(hh)
file_path = "E:/ckpca_code/change/simu_data/e2a4ba200.npy"

## Save the hatk, MSE, and Rand index with 95% CI
mk = rep(0, 5)     ## hatk (estimated number of change points)
me = rep(0, 5)     ## RMSE of hatk
ri = rep(0, 5)     ## mean Rand index
ci1 = rep(0, 5)    ## lower bound of 95% CI for Rand index
ci2 = rep(0, 5)    ## upper bound of 95% CI for Rand index

##Balanced dataset
n1 = 100; n2 = 200; n3 = 300; n4 = 400
n5 = 500; n6 = 600; n7 = 700; n8 = 800

##True change point
true_cps = c(n1, n2, n3, n4, n5, n6, n7, n8)

##Loading data
flat_vec = npyLoad(file_path)
big_array1 = array(flat_vec, dim = c(n, px, times))

##Function to calculate the Rand Index
source("E:/ckpca_code/randi.R")

for (hh0 in 1:hh1){
  
  ###simulation
  number1 = matrix(nrow = 1, ncol = 1000)
  position1 = matrix(rep(0, 10000), nrow = 1000, ncol = 50)
  
  for (kkk in 1:times){
    cat("Running iteration:", kkk, "\n")
    cat("m:", hh[hh0], "\n")
    set.seed(1111 + kkk)
    
    ##Get the data from the kkk-th experiment
    x1 = big_array1[ , , kkk]
    x = x1
    
    ##using thumb method to choose bandwidth
    vv = var(x)
    v0 = 0
    for (i in 1:px){
      v0 = v0 + vv[i, i]
    }
    v0 = v0 / px
    m = hh[hh0]
    sig = sqrt(m * px * v0)  ##
    sig0 = sig
    sig
    
    ###ckpca
    ## Kernel matrix
    dist_mat = as.matrix(dist(x))^2
    k = exp(-dist_mat / (2 * sig^2))
    ##Compute matrices L and U
    op = matrix(rep(1, n * n), nrow = n, ncol = n)
    zs = floor(n / an)
    a0 = array(dim = c(n, an, zs))
    for (i in 1:zs){
      b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
      a1 = diag(rep(1, an))
      b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
      b0 = rbind(b1, a1, b2)
      a0[,,i] = b0
    }
    c0 = array(dim = c(n, an, zs))
    for (i in 1:zs){
      b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
      a1 = matrix(rep(1, an * an), nrow = an, ncol = an)
      b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
      b0 = rbind(b1, a1, b2)
      c0[,,i] = b0
    }
    d0 = array(dim = c(n, an, zs))
    for (i in 1:zs){
      d0[,,i] = a0[,,i] - c0[,,i] / an
      d0[,,i] = d0[,,i]
    }
    U = matrix(rep(0, n * n), nrow = n, ncol = n)
    for (i in 1:zs){
      U = U + d0[,,i] %*% t(d0[,,i]) / (an - 1)
    }
    bn = n - zs * an
    if (bn > 1) {
      b1  =  matrix(rep(0, (n - bn) * an), nrow = n - bn, ncol = bn)
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
    ###U
    U = U / zs
    l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
    l2 = diag(1, n)
    ###L
    L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
    M = (L - U)
    ###kn = (L-U)K
    kk0 = (L - U) %*% k
    ## Eigen decomposition
    spec = eigen(kk0)
    evalue = rep(0, n + 1)
    va = abs(spec$values)
    evalue[1:n] = va
    evector = spec$vectors
    ##TRR
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
    ab = sqrt(sum(B^2))
    ## Projection
    cc1 = k %*% B
    cc1 = t(cc1)
    
    # # # ###E-Divisive after dimension reduction
    ##ecp_c
    output1 = e.divisive(t(cc1), R = 499, alpha = 1)
    kk1 = output1$estimates
    number1[kkk] = length(kk1)
    position1[kkk, 1:(number1[kkk])] = kk1
    
  }
  
  ### For E-Divisive after dimension reduction
  ### ecp_c
  numberr1 = number1
  positionn1 = position1
  
  for (k in 1:times) {
    kk = positionn1[k, ]
    he1 = length(which(kk > 0))
    if (he1 == 2) {
      positionn1[k, ] = 0
    } else {
      positionn1[k, 1:(he1 - 2)] = kk[2:(he1 - 1)]
      positionn1[k, (he1 - 1):(he1)] = c(0, 0)
    }
  }
  
  numberr1 = numberr1 - 2
  
  # Average number of detected change points
  mk[hh0] = mean(numberr1[1:times])
  
  # Compute RMSE for change point count
  mse1 = 0
  for (j in 1:times) {
    a1 = (numberr1[j] - 7)^2
    mse1 = mse1 + a1
  }
  mse1 = sqrt(mse1 / times)
  me[hh0] = mse1
  
  # Rand index and 95% CI
  posi = positionn1[1:times, ]
  num = numberr1[1:times]
  ri_result = randi(posi, num, n, times, true_cps)
  ri[hh0] = ri_result$mean
  ci1[hh0] = ri_result$CI_lower
  ci2[hh0] = ri_result$CI_upper
  
}

### Save all results
res = cbind(mk, me, ri, ci1, ci2)
write.csv(res, "E:/choosem.csv")
res

toc()
