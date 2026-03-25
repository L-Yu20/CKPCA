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
library(abind)

## Basic parameter settings
n = 800
times = 1000
u = 0.5

## choose model
## Choose model
px = 200             ## px = 100 or 200
type = "ba"          ## "ba" for balanced, "im" for imbalanced
save_path = "E:/ckpca_code/change/simu_data/"

## Auto-generate filename
file_tag = paste0("e3", type, px)        # e.g., "e3ba200"
file_name = paste0(save_path, file_tag, ".npy")

# Filename explanation:
# e3  stands for Example 3
# ba  stands for balanced data
# 100 stands for px = 100

## "ba" for balanced, "im" for imbalanced
if (type == "ba") {
  n1 = 100; n2 = 200; n3 = 300; n4 = 400
  n5 = 500; n6 = 600; n7 = 700; n8 = 800
} else if (type == "im") {
  n1 = 30;  n2 = 170; n3 = 350; n4 = 440
  n5 = 520; n6 = 630; n7 = 710; n8 = 800
} else {
  stop("type must be either 'ba' or 'im'")
}

# Distances between adjacent change points
c1 = n1;       c2 = n2 - n1; c3 = n3 - n2; c4 = n4 - n3
c5 = n5 - n4;  c6 = n6 - n5; c7 = n7 - n6; c8 = n8 - n7

## Generate data
mat_list = list()
for (kkk in 1:times) {
  set.seed(1111 + kkk)
  ## data setting
  rho = 0.5
  mu1 = rep(0, px)
  sigma1 = matrix(nrow = px, ncol = px)
  for (i in 1:px) {
    for (j in 1:px) {
      sigma1[i, j] = rho ^ (abs(i - j))
    }
  }
  for (i in 1:px) {
    sigma1[i, i] = 1
  }
  x = matrix(nrow = n, ncol = px)
  x[(1):(n1), ] = mvrnorm(n = c1, mu1, sigma1)
  x[(n1 + 1):(n2), ] = mvrnorm(n = c2, mu1, sigma1)
  x[(n2 + 1):(n3), ] = mvrnorm(n = c3, mu1, sigma1)
  x[(n3 + 1):(n4), ] = mvrnorm(n = c4, mu1, sigma1)
  x[(n4 + 1):(n5), ] = mvrnorm(n = c5, mu1, sigma1)
  x[(n5 + 1):(n6), ] = mvrnorm(n = c6, mu1, sigma1)
  x[(n6 + 1):(n7), ] = mvrnorm(n = c7, mu1, sigma1)
  x[(n7 + 1):(n8), ] = mvrnorm(n = c8, mu1, sigma1)
  
  ## set change point
  k = floor(sqrt(px))
  u1 = matrix(rep(0, k * c2), nrow = c2, ncol = k)
  u2 = matrix(rep(0, k * c4), nrow = c4, ncol = k)
  u3 = matrix(rep(0, k * c6), nrow = c6, ncol = k)
  u4 = matrix(rep(0, k * c8), nrow = c8, ncol = k)
  u1[, 1:5] = u
  u1[, 6:10] = u / 2
  u2[, 1:5] = u
  u2[, 6:10] = u / 3
  u3[, 1:5] = u / 2
  u3[, 6:10] = u
  u4[, 1:5] = u / 3
  u4[, 6:10] = u 
  x[1:(n1), 1:k] = x[1:(n1), 1:k]
  x[(n1 + 1):(n2), 1:k] = x[(n1 + 1):(n2), 1:k] + u1
  x[(n2 + 1):(n3), 1:k] = x[(n2 + 1):(n3), 1:k]
  x[(n3 + 1):(n4), 1:k] = x[(n3 + 1):(n4), 1:k] + u2
  x[(n4 + 1):(n5), 1:k] = x[(n4 + 1):(n5), 1:k]
  x[(n5 + 1):(n6), 1:k] = x[(n5 + 1):(n6), 1:k] + u3
  x[(n6 + 1):(n7), 1:k] = x[(n6 + 1):(n7), 1:k]
  x[(n7 + 1):(n8), 1:k] = x[(n7 + 1):(n8), 1:k] + u4
  
  ## Save the data from each loop iteration
  mat_list[[kkk]] = x
}

## Save the data
big_array = abind(mat_list, along = 3)
print(dim(big_array))
npySave(file_name, big_array)


