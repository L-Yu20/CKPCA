### Choose cn (ckpca + ecp) evaluation

library(tictoc)
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(moments)
library("wbs")
library(RcppCNPy)

remove(list = ls())
tic()

#### Basic settings
n = 800        # number of samples
px = 100       # dimension
tau = 0.5      # tau for TRR
times = 2   # number of experiments
type = "ba"    # "ba" for balanced, "im" for imbalanced
case = 1

## Dynamic path and file name
save_path = "E:/"
file_tag = paste0("e1", "c", case, type, px)     # e.g., "e1c1ba100"
## Loading path
npy_file = paste0(save_path, "ckpca_code/change/simu_data/", file_tag, ".npy")
## Saving path
file_name = paste0(save_path, "choose_cn", file_tag, ".csv")

## Function to calculate the Rand Index
source("E:/ckpca_code/randi.R")

##Data loading
flat_vec = npyLoad(npy_file)
big_array1 = array(flat_vec, dim = c(n, px, times))

## Candidate block sizes (beta_n)
cn0 = c(
  0.02,0.2,2 
)
lan = length(cn0)

## balanced or imbalanced
if (type == "ba") {
  n1 = 100; n2 = 200; n3 = 300; n4 = 400
  n5 = 500; n6 = 600; n7 = 700; n8 = 800
} else if (type == "im") {
  n1 = 30; n2 = 170; n3 = 350; n4 = 440
  n5 = 520; n6 = 630; n7 = 710; n8 = 800
} else {
  stop("type must be either 'ba' or 'im'")
}

## True change points
true_cps = c(n1, n2, n3, n4, n5, n6, n7, n8)

## Store results
number1 = matrix(0, nrow = lan, ncol = times)
position1 = array(0, dim = c(lan, times, 50))

## Save the hatk, MSE, and Rand index with 95% CI
mk = rep(0, lan)   # estimated number of change points
me = rep(0, lan)   # RMSE
ri = rep(0, lan)   # Rand index
ci1 = rep(0, lan)  # lower bound of 95% CI for Rand index
ci2 = rep(0, lan)  # upper bound of 95% CI for Rand index

for (kkk in 1:times) {
  
  ## Get the data from the kkk-th experiment
  x1 = big_array1[ , , kkk]
  x = x1
  
  ## using thumb method to choose bandwidth
  vv = var(x)
  v0 = 0
  for (i in 1:px) {
    v0 = v0 + vv[i, i]
  }
  v0 = v0 / px
  m = 0.8
  sig = sqrt(m * px * v0)
  
  ### ckpca
  ## Kernel matrix
  dist_mat = as.matrix(dist(x))^2
  k = exp(-dist_mat / (2 * sig^2))
  
  for (ajk in 1:lan) {
    set.seed(1111 + kkk)
    cat("Running iteration:", kkk, "\n")
    cat("cn:", cn0[ajk], "\n")
    
    ## get betan
    # an = an0[ajk]
    an = floor(n^0.5)  ## betan
    
    ## Compute matrices L and U
    op = matrix(rep(1, n * n), nrow = n, ncol = n)
    zs = floor(n / an)
    a0 = array(dim = c(n, an, zs))
    for (i in 1:zs) {
      b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
      a1 = diag(rep(1, an))
      b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
      b0 = rbind(b1, a1, b2)
      a0[,,i] = b0
    }
    c0 = array(dim = c(n, an, zs))
    for (i in 1:zs) {
      b1 = matrix(rep(0, (i - 1) * an * an), nrow = (i - 1) * an, ncol = an)
      a1 = matrix(rep(1, an * an), nrow = an, ncol = an)
      b2 = matrix(rep(0, (n - i * an) * an), nrow = (n - i * an), ncol = an)
      b0 = rbind(b1, a1, b2)
      c0[,,i] = b0
    }
    d0 = array(dim = c(n, an, zs))
    for (i in 1:zs) {
      d0[,,i] = a0[,,i] - c0[,,i] / an
      d0[,,i] = d0[,,i]
    }
    U = matrix(rep(0, n * n), nrow = n, ncol = n)
    for (i in 1:zs) {
      U = U + d0[,,i] %*% t(d0[,,i]) / (an - 1)
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
      warning("Last segment size is 1; skipping to avoid division by zero.")
    }
    ### U
    U = U / zs
    l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
    l2 = diag(1, n)
    ### L
    L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
    M = (L - U)
    ### get kn = (L - U)K
    kk0 = (L - U) %*% k
    ## Eigen decomposition
    spec = eigen(kk0)
    evalue = rep(0, n + 1)
    va = abs(spec$values)
    evalue[1:n] = va
    evector = spec$vectors
    ## TRR
    ####
    ####change this cn
    ccn = cn0[ajk]*log(log(n)) * sqrt(1 / n)
    ####
    ####
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
    ## Projection
    cc1 = k %*% B
    cc1 = t(cc1)
    
    ## E-divisive
    output1 = e.divisive(t(cc1), R = 499, alpha = 1)
    est_cp = output1$estimates
    number1[ajk, kkk] = length(est_cp)
    position1[ajk, kkk, 1:length(est_cp)] = est_cp
  }
}

## Evaluation
for (ajk in 1:lan) {
  numberr1 = number1[ajk, ]
  positionn1 = position1[ajk, , ]
  
  for (k in 1:times) {
    pos = positionn1[k, ]
    he1 = length(which(pos > 0))
    if (he1 == 2) {
      positionn1[k, ] = 0
    } else {
      positionn1[k, 1:(he1 - 2)] = pos[2:(he1 - 1)]
      positionn1[k, (he1 - 1):(he1)] = c(0, 0)
    }
  }
  
  numberr1 = numberr1 - 2
  mk[ajk] = mean(numberr1)
  me[ajk] = sqrt(mean((numberr1 - 7)^2))
  ri_result = randi(positionn1, numberr1, n, times, true_cps)
  ri[ajk] = ri_result$mean
  ci1[ajk] = ri_result$CI_lower
  ci2[ajk] = ri_result$CI_upper
}

### Save all results
res = cbind(mk, me, ri, ci1, ci2)
write.csv(res, file_name)
res

toc()

