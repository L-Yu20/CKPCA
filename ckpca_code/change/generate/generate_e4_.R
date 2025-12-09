library(tictoc)
library(MASS)
library('ecp')
library('mvtnorm')
library("InspectChangepoint")
library("hdbinseg")
library(moments)
library("wbs")
library(RcppCNPy)
library(abind)

## Basic parameter settings
n = 800
times = 2

## Choose model
px = 100             ## px = 100 or 200
case = 4
type = "ba"          ## "ba" for balanced, "im" for imbalanced
save_path = "E:/ckpca_code/change/simu_data/"  ## base path
## Auto-generate file name
file_tag = paste0("e1", "c", case, type, px)      # e.g., "e1c1ba100"
file_name = paste0(save_path, file_tag, ".npy")        # full path to .npy
# Filename explanation:
# e1  stands for experiment 1
# c1  stands for case 1
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

## Distances between adjacent change points
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
      if (case == 1) {
        sigma1[i, j] = rho  # case 1: constant correlation rho
      } else if (case == 2) {
        sigma1[i, j] = rho ^ abs(i - j)  # case 2: AR(1)-type correlation
      } else {
        stop("Invalid case value")
      }
    }
  }
  
  for (i in 1:px) {
    sigma1[i, i] = 2.5
  }
  
  x = matrix(nrow = n, ncol = px)

  x[1:n1, ] = mvrnorm(n = c1, mu1, sigma1)
  x[(n1 + 1):n2, ] = mvrnorm(n = c1, mu1, sigma1)
  x[(n2 + 1):n3, ] = mvrnorm(n = c3, mu1, sigma1)
  x[(n3 + 1):n4, ] = matrix(t2, nrow = c4, ncol = px)
  x[(n4 + 1):n5, ] = mvrnorm(n = c5, mu1, sigma1)
  x[(n5 + 1):n6, ] = matrix(t3, nrow = c6, ncol = px)
  x[(n6 + 1):n7, ] = mvrnorm(n = c7, mu1, sigma1)
  x[(n7 + 1):n8, ] = matrix(t4, nrow = c8, ncol = px)
  
  ## Save the data from each loop iteration
  mat_list[[kkk]] = x
}

### Save the data
big_array = abind(mat_list, along = 3)
print(dim(big_array)) 
npySave(file_name, big_array)
