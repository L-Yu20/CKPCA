remove(list = ls())
##mnist data for sc
library(MASS)
library(cluster)
library('mvtnorm')
library("mclust")
library("fpc")
library("R.matlab")
times = 50
set.seed(22)
ri2 = matrix(ncol = 1, nrow = times)
ri1 = matrix(ncol = 1, nrow = times)
ri0 = matrix(ncol = 1, nrow = times)
dataset1 = read.csv(file = 'E:\\ckpca_code\\dataset\\mnist\\mnist_test.csv', header = F)
dataset2 = read.csv(file = 'E:\\ckpca_code\\dataset\\mnist\\mnist_train.csv', header = F)
dataset = rbind(dataset1, dataset2)
y0 = dataset[, 1]
w6 = which(y0 == 6)
w8 = which(y0 == 8)
w9 = which(y0 == 9)
nl1 = 300
nl2 = 300
nl3 = 300
cl1 = nl1
cl2 = nl1 + nl2
cl3 = nl1 + nl2 + nl3
n = nl1 + nl2 + nl3
y1 = rep(0, n)
y1[1:cl1] = 1
y1[(cl1 + 1):cl2] = 2
y1[(cl2 + 1):cl3] = 3

SpectralClustering <- function(data, num_clusters) {
  #' Conducts spectral clustering algorithm.
  #'
  #' @description This function implements spectral clustering algorithm.
  #'
  #' It conducts k-means on top k left singular vectors of the data matrix.
  #'
  #' @param data matrix. Input data matrix for clustering.
  #' @param num_clusters int. Number of clusters.
  #' @usage SpectralClustering(data, num_clusters)
  #' @return A vector containing the cluster assignment after iterations.
  #' @references T. Liu, Y. Lu, B. Zhu, H. Zhao (2021). High-dimensional Clustering via Feature Selection with
  #' Applications to Single Cell RNA-seq Data.
  #' @examples
  #' synthetic_data <- GenerateSyntheticData(n=10, p=10, s=5, k=2, signal_strength=1, noise_type="gaussian")
  #' label.est <- SpectralClustering(synthetic_data$data, 2)
  #' @export
  data.svd = svd(data)
  r = min(dim(data)[1], dim(data)[2], num_clusters)
  kmeans.result = kmeans(
    data.svd$u[, 1:r],
    centers = num_clusters,
    iter.max = 200,
    nstart = 100
  )
  return(kmeans.result$cluster)
}


for (time in 1:times){
  ## Random sampling
  rj1 = sample(1:length(w6), nl1, replace = FALSE)
  rj2 = sample(1:length(w8), nl2, replace = FALSE)
  rj3 = sample(1:length(w9), nl3, replace = FALSE)
  
  ww6 = w6[rj1]
  ww8 = w8[rj2]
  ww9 = w9[rj3]
  
  x1 = dataset[ww6, 2:785]
  x2 = dataset[ww8, 2:785]
  x3 = dataset[ww9, 2:785]
  
  x = rbind(x1, x2, x3)
  x = as.matrix(x)
  
  dimx = dim(x)
  px = dimx[2]
  n = dimx[1]
  
  for (i in 1:n) {
    if (sum(x[i, ]) == 0) {
      x[i, ] = x[i, ]
    } else {
      x[i, ] = scale(x[i, ])
    }
  }
  
  repeat1 = 30
  kinds = 3
  an = floor(n^0.5)
  
  ## Bandwidth
  vv = var(x)
  v0 = 0
  for (i in 1:px) {
    v0 = v0 + vv[i, i]
  }
  v0 = v0 / px
  m = 0.8
  sig = sqrt(m * px * v0)
  
  
  ## Initial value
  cluster1 = matrix(nrow = n, ncol = repeat1)
  k = matrix(rep(0, n * n), nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      aa = x[i, ] - x[j, ]
      bb = sum(aa^2)
      k[i, j] = exp(-bb / (2 * sig^2))
    }
  }
  
  op = matrix(rep(1, n * n), nrow = n, ncol = n)
  l1 = matrix(rep(1, n * n), nrow = n, ncol = n)
  l2 = diag(1, n)
  L = (l2 - l1 / n) %*% (l2 - l1 / n) / n
  
  ke0 = L %*% k
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
  B1 = evector1[, 1:hatq]
  
  cc = k %*% B1
  cc = t(cc)
  cc = Re(cc)
  
  cluster1[, 1] = SpectralClustering(x, kinds)
  
  ## spe
  cluster2 = SpectralClustering(t(cc), kinds)
  
  
  ## iteration
  for (order in 2:repeat1) {
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
    
    if (cc1 == 0) {
      g1 = g1
      h1 = h1
    } else {
      for (j in 1:ccc[1]) {
        g1[b1[j], j] = 1
        h1[b1[j], ] = 1
      }
    }
    
    if (cc2 == 0) {
      g2 = g2
      h2 = h2
    } else {
      for (j in 1:ccc[2]) {
        g2[b2[j], j] = 1
        h2[b2[j], ] = 1
      }
    }
    
    if (cc3 == 0) {
      g3 = g3
      h3 = h3
    } else {
      for (j in 1:ccc[3]) {
        g3[b3[j], j] = 1
        h3[b3[j], ] = 1
      }
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
    lambda1 = evalueplus[2:(n + 1)] / evalueplus[1:n]
    lambda2 = which(lambda1[1:(n - 1)] < 0.5)
    aaa = min(lambda1[1:(n - 1)])
    
    if (aaa > 0.5) {
      hatq = 1
    } else {
      hatq = max(lambda2[which(lambda2 < px)])
    }
    
    B = evector[, 1:hatq]
    cc = k %*% B
    cc = t(cc)
    cc = Re(cc)
    
    cluster1[, order] = SpectralClustering(t(cc), kinds)
  }
  
  
  
  ## rand index
  rri = matrix(nrow = 1, ncol = repeat1)
  for (kk in 2:repeat1) {
    indexx = matrix(nrow = 1, ncol = n)
    indexx0 = matrix(nrow = 1, ncol = n)
    indexx = cluster1[, kk]
    indexx0 = cluster1[, kk - 1]
    
    all0 = n * (n - 1) / 2
    tp = 0
    tn = 0
    
    for (jj in 2:n) {
      for (hh in 1:(jj - 1)) {
        if (indexx[jj] == indexx[hh] & indexx0[jj] == indexx0[hh])
          tp = tp + 1 else if (indexx[jj] != indexx[hh] & indexx0[jj] != indexx0[hh])
            tn = tn + 1
      }
    }
    
    rri[kk] = (tp + tn) / all0
  }
  
  for (i in 2:repeat1) {
    if (rri[i] > 0.999) {
      break
    }
  }
  st = i
  
  ## our method
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0 = y1
  randindex = matrix(nrow = 1, ncol = repeat1)
  
  all0 = n * (n - 1) / 2
  h = st
  index = cluster1[, h]
  tp = 0
  tn = 0
  
  for (jj in 2:n) {
    for (hh in 1:(jj - 1)) {
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri0[time] = (tp + tn) / all0
  
  ## dim-1
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0 = y1
  index = cluster2
  tp = 0
  tn = 0
  
  for (jj in 2:n) {
    for (hh in 1:(jj - 1)) {
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri1[time] = (tp + tn) / all0
  
  ## sc(x)
  index = matrix(nrow = 1, ncol = n)
  index0 = matrix(nrow = 1, ncol = n)
  index0 = y1
  index = cluster1[, 1]
  tp = 0
  tn = 0
  
  for (jj in 2:n) {
    for (hh in 1:(jj - 1)) {
      if (index[jj] == index[hh] & index0[jj] == index0[hh])
        tp = tp + 1 else if (index[jj] != index[hh] & index0[jj] != index0[hh])
          tn = tn + 1
    }
  }
  ri2[time] = (tp + tn) / all0
  
 
}





## our method
mri0 = mean(ri0)
sd0 = sd(ri0)

## dim-1
mri1 = mean(ri1)
sd1 = sd(ri1)

## sc(x)
mri2 = mean(ri2)
sd2 = sd(ri2)


result = rbind(
  c(mri0, sd0),
  c(mri1, sd1),
  c(mri2, sd2)
)

write.csv(result, file = "E:/mnsc.csv")


