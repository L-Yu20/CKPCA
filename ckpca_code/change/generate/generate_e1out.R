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
times = 1000

## Choose model
px = 100             ## px = 100 or 200
case = 1
save_path = "E:/ckpca_code/change/simu_data/"  ## base path
## Auto-generate output file name
file_tag_out = paste0("e1", "c", case, "out", px)  # e.g., "e1c1out100"
file_name = paste0(save_path, file_tag_out, ".npy")     # full path to .npy file
# Filename explanation:
# e1  stands for Example 1
# c1  stands for case 1
# out  stands for data with outliers
# 100 stands for px = 100

## imbalanced
n1 = 30;  n2 = 170; n3 = 350; n4 = 440
n5 = 520; n6 = 630; n7 = 710; n8 = 800

## Distances between adjacent change points
c1 = n1;       c2 = n2 - n1; c3 = n3 - n2; c4 = n4 - n3
c5 = n5 - n4;  c6 = n6 - n5; c7 = n7 - n6; c8 = n8 - n7

mat_list = list()

for (kkk in 1:times) {
  set.seed(1111 + kkk)
  
  rho = 0.5
  mu1 = rep(0, px)
  sigma1 = matrix(nrow = px, ncol = px)
  
  for (i in 1:px) {
    for (j in 1:px) {
      if (case == 1) {
        sigma1[i, j] = rho
      } else if (case == 2) {
        sigma1[i, j] = rho ^ abs(i - j)
      } else {
        stop("Invalid case value")
      }
    }
  }
  
  for (i in 1:px) {
    sigma1[i, i] = 2.5
  }
  
  x = matrix(nrow = n, ncol = px)
  tt1 = px * c2
  tt2 = px * c4
  tt3 = px * c6
  tt4 = px * c8
  
  t1 = runif(n = tt1, min = -3, max = 3)
  t2 = runif(n = tt2, min = -3, max = 3)
  t3 = runif(n = tt3, min = -3, max = 3)
  t4 = runif(n = tt4, min = -3, max = 3)
  
  x[1:n1, ] = mvrnorm(n = c1, mu1, sigma1)
  x[(n1 + 1):n2, ] = matrix(t1, nrow = c2, ncol = px)
  x[(n2 + 1):n3, ] = mvrnorm(n = c3, mu1, sigma1)
  x[(n3 + 1):n4, ] = matrix(t2, nrow = c4, ncol = px)
  x[(n4 + 1):n5, ] = mvrnorm(n = c5, mu1, sigma1)
  x[(n5 + 1):n6, ] = matrix(t3, nrow = c6, ncol = px)
  x[(n6 + 1):n7, ] = mvrnorm(n = c7, mu1, sigma1)
  x[(n7 + 1):n8, ] = matrix(t4, nrow = c8, ncol = px)
  
  ### Add outliers (5% of rows in each segment, 5% of columns)
  ac1 = floor(0.05 * n1);       ac2 = floor(0.05 * (n2 - n1))
  ac3 = floor(0.05 * (n3 - n2)); ac4 = floor(0.05 * (n4 - n3))
  ac5 = floor(0.05 * (n5 - n4)); ac6 = floor(0.05 * (n6 - n5))
  ac7 = floor(0.05 * (n7 - n6)); ac8 = floor(0.05 * (n8 - n7))
  
  aa1 = sample(1:n1, ac1, replace = FALSE)
  aa2 = sample((n1 + 1):n2, ac2, replace = FALSE)
  aa3 = sample((n2 + 1):n3, ac3, replace = FALSE)
  aa4 = sample((n3 + 1):n4, ac4, replace = FALSE)
  aa5 = sample((n4 + 1):n5, ac5, replace = FALSE)
  aa6 = sample((n5 + 1):n6, ac6, replace = FALSE)
  aa7 = sample((n6 + 1):n7, ac7, replace = FALSE)
  aa8 = sample((n7 + 1):n8, ac8, replace = FALSE)
  
  ac0 = floor(0.05 * px)
  
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa1, aa00] = x[aa1, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa2, aa00] = x[aa2, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa3, aa00] = x[aa3, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa4, aa00] = x[aa4, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa5, aa00] = x[aa5, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa6, aa00] = x[aa6, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa7, aa00] = x[aa7, aa00] + 5
  aa00 = sample(1:px, ac0, replace = FALSE); x[aa8, aa00] = x[aa8, aa00] + 5
  
  mat_list[[kkk]] = x
}

big_array = abind(mat_list, along = 3)
print(dim(big_array))
npySave(file_name, big_array)


