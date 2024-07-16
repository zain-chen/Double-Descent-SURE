## load the needed package
library(MASS)
library(nloptr)
library(parallel)

## set the simulation parameters
set.seed(1002)
beta = runif(100)                               # real coefficients
beta = beta/sqrt(sum(beta^2))                    # convert to a unit vector                                      # number of simulations
# number of samples
N = seq(50, 150, 1)

## conduct the simulation
simulate <- function(N){
  M = 10
  test_MSE = matrix(nrow = length(N), ncol = M)
  test_SURE = matrix(nrow = length(N), ncol = M)
  for (i in seq_along(N)){
    n = N[i]
    print(n)
    for (m in 1:M){
      print(m)
      # generate training data
      X = replicate(100, rnorm(n))
      e = rnorm(n, sd = 1)
      y = X %*% beta + e
      
      # estimate the beta_hat
      sure <- function(beta, X, y, n){ 
        # sigma = 1, H is hat matrix
        #if(n < 1000){H = X %*% ginv(X)}
        #else{H = X %*% solve(t(X) %*% X) %*% t(X)}
        H = X %*% ginv(X)
        sure = t(y - X %*% beta) %*% (y - X %*% beta) - n + 2*sum(diag(H)) 
        return(as.numeric(sure))
      }
      
      beta_hat_init = ginv(X) %*% y # pseudo inverse
      opts <- list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1.0e-8)
      objective_function <- function(beta) sure(beta, X, y, n)
      result <- nloptr(x0 = beta_hat_init, 
                       eval_f = objective_function, 
                       opts = opts)
      
      # generate test set
      X_test = replicate(100, rnorm(1000))
      e_test = rnorm(100, sd = 1)
      y_test = X_test %*% beta + e_test
      
      # measure model accuracy
      beta_hat = result$solution
      preds_test = X_test %*% beta_hat
      test_MSE[i, m] = sqrt(mean((y_test - preds_test)^2))
      test_SURE[i, m] = sure(beta_hat, X_test, y_test, n)
    }
  }
  result <- list(test_MSE, test_SURE)
  return(result)
}

result <- simulate(N)
test_MSE <- result[[1]]
test_SURE <- result[[2]]

plot(100/N,rowMeans(test_MSE), ylab = "Test MSE", xlab = "p / n", 
     pch = 4)
lines(100/N, rowMeans(test_MSE), col = "red", lwd = 2)

plot(100/N,rowMeans(test_SURE), ylab = "Test SURE Risk", xlab = "p / n", pch = 4)
lines(100/N, rowMeans(test_SURE), col = "red", lwd = 2)

plot(N,rowMeans(test_MSE), ylab = "Test MSE", xlab = "n", 
     pch = 4)
lines(N, rowMeans(test_MSE), col = "red", lwd = 2)

plot(N,rowMeans(test_SURE), ylab = "Test SURE Risk", xlab = "n", pch = 4)
lines(N, rowMeans(test_SURE), col = "red", lwd = 2)