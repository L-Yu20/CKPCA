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

## Choose model
px = 200             ## px = 100 or 200
af = 4               ## degrees of freedom, af = 4 or 6
type = "ba"          ## "ba" for balanced, "im" for imbalanced
save_path = "E:/ckpca_code/change/simu_data/"  ## base path

## Auto-generate filename
file_tag = paste0("e2", "a", af, type, px)      # e.g., "e2a4ba200"
file_name = paste0(save_path, file_tag, ".npy") # full path
# Filename explanation:
# e2  stands for Example 2
# a4   stands for af = 4
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
  ## distribution
  x = matrix(nrow = n, ncol = px)
  x[(1) : (n1), ]       = mvrnorm(n = c1, mu1, sigma1)
  x[(n1 + 1) : (n2), ]  = rmvt(n = c2, sigma1, df = af)
  x[(n2 + 1) : (n3), ]  = mvrnorm(n = c3, mu1, sigma1)
  x[(n3 + 1) : (n4), ]  = rmvt(n = c4, sigma1, df = af)
  x[(n4 + 1) : (n5), ]  = mvrnorm(n = c5, mu1, sigma1)
  x[(n5 + 1) : (n6), ]  = rmvt(n = c6, sigma1, df = af)
  x[(n6 + 1) : (n7), ]  = mvrnorm(n = c7, mu1, sigma1)
  x[(n7 + 1) : (n8), ]  = rmvt(n = c8, sigma1, df = af)
  
  ## Save the data from each loop iteration
  mat_list[[kkk]] = x
}

## Save the data
big_array = abind(mat_list, along = 3)
print(dim(big_array))
npySave(file_name, big_array)


