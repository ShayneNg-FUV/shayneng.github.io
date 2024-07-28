---
layout: page
title: First Project
description: Written in R
img: assets/img/12.jpg
importance: 1
category: work
related_publications: false
---

Here's the code:

{% raw %}
```r
---
title: |
  | Jackknife and Regularized Jackknife Estimation 
author: |
  | Dung (Shayne) Nguyen
date: "April 19, 2024"
output: html_notebook
---
```
# 2SLSE, JIVE & RJIVE

## 0. Setting up the Environment
```r
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, results = 'hide')

options(scipen = 10)

library(pacman)

p_load(MASS, Rmisc)

rm(list=ls())
```

## 1. Data Generating Process (D.G.P.)

### 1-1. Constructing the Data Simulation Function via D.G.P.
```r
generate_data<- function(n, K, Pi) {
  # Z: instrumental variable
  # [n*K]
  Z <- matrix(sample(c(0.5, -0.5), n*K, replace = TRUE, prob = c(0.5, 0.5)),
                     nrow = n, ncol = K)
  
  # Z_i: the i^th row of Z
  # [1*K]
  # E((Z_i)'(Z_i)) = 0.25 * I_K
  # [K*K, diagonal]
  exp_ZiTZi <- diag(0.25, nrow = K)
  print(exp_ZiTZi)
  # Pi (T.B.D.)
  # [K*1]
  
  sigma_e_sq <- 0.2
  sigma_e <- sqrt(sigma_e_sq)
  
  mu_sq <- 30
  sigma_U_sq <- n * t(Pi) %*% exp_ZiTZi %*% Pi / mu_sq
  sigma_U <- sqrt(sigma_U_sq)
  
  
  rho <- 0.6 # correlation of e & U
  sigma_e_U <- rho * sigma_e * sigma_U
  
  # Covariance of e & U
  # [2*2]
  cov <- matrix(c(sigma_e_sq, sigma_e_U, sigma_e_U, sigma_U_sq),
                nrow = 2, ncol = 2)
  
  # Error
  # [n*2]
  err <- mvrnorm(n = n, mu = c(0, 0), Sigma = cov)
  # U: error of stage 1
  # [n*1]
  U <- matrix(err[, 1])
  # e: error of stage 2
  # [n*1]
  e <- matrix(err[, 2])
  
  # X: independent variable
  # [n*1]
  X <- Z %*% Pi + U
  
  delta_0 <- 1
  # Y: dependent variable
  # [n*1]
  Y <- delta_0 * X + e
  
  return(list(Z = Z, X = X, Y = Y,
              U = U, e = e))
}
```

### 1-2. Implementing the Data Simulation Function to Generate a new dataset
```r
n <- 100
K <- 95

# pi = 1 in row 1 ~ 0.4K (i.e. 38)
# pi = 0 otherwise
Pi_0 <- matrix(0, nrow = 95, ncol = 1)
Pi_0[1:(K*0.4),] <- 1

generate_data(n = n, K = K, Pi = Pi_0)
```

## 2. Defining the Two-Stage Least Squares Estimator (2SLSE)
```r
tslse <- function(data) {
  Z <- data$Z
  X <- data$X
  Y <- data$Y
  
  # Stage 1
  # X\hat = Z (Z' Z)^(-1) Z' X
  # [n*1]
  X_pred <- Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X
  
  # Stage 2
  # delta\hat = (X\hat' X\hat)^(-1) X\hat' Y
  # [1*1]
  beta_hat <- solve(t(X_pred) %*% X_pred) %*% t(X_pred) %*% Y
  
  return(beta_hat)
}
```


## 3. Defining the Jackknife Instrumental Variables Estimator (JIVE)
```r
jive <- function(data, alpha_hat) {
  Z <- data$Z
  X <- data$X
  Y <- data$Y
  
  # P = Z (Z' Z)^(-1) Z'
  # [n*n]
  P <- Z %*% solve(t(Z) %*% Z) %*% t(Z)         
  
  # Note: G = 1
  
  # V\bar = (P - alpha\hat I_n) X
  # [n*1]
  V <- (P - alpha_hat * diag(1, nrow = n)) %*% X
  
  # delta\tilde = (V\bar' X)^(-1) (V\bar' Y)
  # [1*1]
  delta_pred <- solve(t(V) %*% X) %*% (t(V) %*% Y)
  
  return(delta_pred)
}
```



## 4. Defining the Regularized Jackknife Instrumental Variables Estimator (RJIVE)
```r
rjive <- function(data) {
  Z <- data$Z
  X <- data$X
  Y <- data$Y
  
  n <- nrow(Z)
  K <- ncol(Z)
  
  # C: constant of proportionality (sample standard deviation of X_i)
  C <- sd(X)
  # gamma: scalar penalty parameter
  # gamma^(1/2) = C * K^(1/2)
  gamma_sqrt <- C * sqrt(K)
  # Lambda = gamma^(1/2) I_K
  # [K*K, diagonal]
  Lambda <- gamma_sqrt * diag(1, nrow = K)
  
  # Pi\hat_aggr: Aggregation of Pi\hat_Lambda_-i
  # [K*n]
  Pi_aggr <- matrix(0, nrow = K, ncol = n)
  for (i in 1:n) {
    # Z_-i: sub-Z (without the i^th row)
    # [(n-1)*K]
    Z_ni <- Z[-i,]
    # X_-i: sub-X (without the i^th row)
    # [(n-1)*1]
    X_ni <- X[-i,]
    # Pi\hat_Lambda_-i = (Z_-i' Z_-i + Lambda' Lambda)^(-1) (Z_-i' X_-i)
    # [K*1]
    Pi_aggr[, i] <-
      solve(t(Z_ni) %*% Z_ni + t(Lambda) %*% Lambda) %*% (t(Z_ni) %*% X_ni)
  }
  
  # sum_1 = sum[(Pi_-i)' (Z_i)' (X_i)']
  # [1*1]
  sum_1 <- 0
  for (i in 1:n) {
    sum_1 <- sum_1 + t(Pi_aggr[, i]) %*% Z[i,] %*% t(X[i,])
  }
  
  # sum_2 = sum[(Pi_-i)' (Z_i)' Y_i]
  # [1*1]
  sum_2 <- 0
  for (i in 1:n) {
    sum_2 <- sum_2 + t(Pi_aggr[, i]) %*% Z[i,] %*% Y[i,]
  }
  
  # delta\tilde = (sum_1)^(-1) sum_2
  # [1*1]
  delta_pred <- solve(sum_1) %*% sum_2
  
  return(delta_pred)
}
```

## 5. Monte Carlo (MC) Simulation

### 5-1. Constructing the MC Simulation Function

```r 
#message=TRUE, warning=TRUE, include=TRUE
mc_simulation <- function(R, n, K, Pi) {
  results <- list()

  for (r in 1:R){
    
    data <- generate_data(n, K, Pi)
    
    # Deriving estimators of each method
    delta_hat_jive_0 <- jive(data, alpha_hat =  0)
    delta_hat_jive_01 <- jive(data, alpha_hat = 0.1)
    delta_hat_jive_05 <- jive(data, alpha_hat = 0.5)
    delta_hat_jive_1 <- jive(data, alpha_hat = 1)
    delta_hat_rjive <- rjive(data)
    delta_hat_tsls <- tslse(data)
    
    #Storing results of each method
    results[[r]] <- list(
      delta_hat_jive_0 = delta_hat_jive_0,
      delta_hat_jive_01 = delta_hat_jive_01,
      delta_hat_jive_05 = delta_hat_jive_05,
      delta_hat_jive_1 = delta_hat_jive_1,
      delta_hat_rjive = delta_hat_rjive,
      delta_hat_tsls = delta_hat_tsls)
  }
  
  #Extract estimates from results for each estimator using sapply to simplify the result to a vector
  estimates_jive_0 <- sapply(results, function(x) x$delta_hat_jive_0)
  estimates_jive_01 <- sapply(results, function(x) x$delta_hat_jive_01)
  estimates_jive_05 <- sapply(results, function(x) x$delta_hat_jive_05)
  estimates_jive_1 <- sapply(results, function(x) x$delta_hat_jive_1)
  estimates_rjive <- sapply(results, function(x) x$delta_hat_rjive)
  estimates_tsls <- sapply(results, function(x) x$delta_hat_tsls)
  
  # Calculate the bias for each estimator (difference from the true value)
  median_bias_jive_0 <- sapply(estimates_jive_0, function(x) (x- 1))
  median_bias_jive_01 <- sapply(estimates_jive_01, function(x) (x- 1))
  median_bias_jive_05 <- sapply(estimates_jive_05, function(x) (x- 1))
  median_bias_jive_1 <- sapply(estimates_jive_1, function(x) (x- 1))
  median_bias_rjive <- sapply(estimates_rjive, function(x) (x- 1))
  median_bias_tsls <- sapply(estimates_tsls, function(x) (x- 1))
  
  # Calculate the Mean Absolute Difference for each estimator
  mad_jive_0 <- sapply(estimates_jive_0, function(x) abs(x- 1))
  mad_jive_01 <- sapply(estimates_jive_01, function(x) abs(x- 1))
  mad_jive_05 <- sapply(estimates_jive_05, function(x) abs(x- 1))
  mad_jive_1 <- sapply(estimates_jive_1, function(x) abs(x-1))
  mad_rjive <- sapply(estimates_rjive, function(x) abs(x- 1))
  mad_tsls <- sapply(estimates_tsls, function(x) abs(x - 1))
  

  # Create a data frame summarizing the median bias and MAD for each estimator
  summary_stats <- data.frame(
    Estimator = c("JIVE alpha_hat = 0","JIVE alpha_hat = 0.1","JIVE alpha_hat = 0.5","JIVE alpha_hat = 1", "RJIVE", "TSLS"),
    Median_Bias = c(median(median_bias_jive_0),median(median_bias_jive_01),median(median_bias_jive_05),median(median_bias_jive_1),median(median_bias_rjive) ,median(median_bias_tsls)),
    MAD = c(median(mad_jive_0),median(mad_jive_01),median(mad_jive_05),median(mad_jive_1), median(mad_rjive), median(mad_tsls))
  )
  return(summary_stats)
}

```


### 5-2. Defining *R*, *n*, *K* & *$\Pi$*, and Implementing MC Simulation
```r
set.seed(123)

R <- 1500
n <- 100
K <- 95

# pi = 1 in row 1 ~ 0.4K (i.e. 38)
# pi = 0 otherwise
Pi_1 <- matrix(0, nrow = K, ncol = 1)
Pi_1[1:(K*0.4),] <- 1

mc_1 <- mc_simulation(R = R, n = n, K = K, Pi = Pi_1)
mc_1
```

## 6. Repeating Q5 with Different *$\Pi$*
```r
set.seed(123)

R <- 1500
n <- 100
K <- 95

# pi = 1 in row 1 ~ 5
# pi = 0 otherwise
Pi_2 <- matrix(0, nrow = K, ncol = 1)
Pi_2[1:5,] <- 1

mc_2 <- mc_simulation(R = R, n = n, K = K, Pi = Pi_2)
mc_2
```

{% endraw %}
