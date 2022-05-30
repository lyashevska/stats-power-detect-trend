#####
## working with lambda
get_mean_ziP <-
  function(lambda, theta) {
    mu <-
      (1 - exp(-exp(theta[1] + exp(theta[2]) * log(lambda)))) * lambda * exp(lambda) / (exp(lambda) - 1)
    return(mu)
  }
get_diff <- function(lambda_new, lambda_orig, theta, p_change) {
  ## p_change = -0.1 is a 10% decrease for example
  mu_orig <- get_mean_ziP(lambda_orig, theta)
  mu_new <- (1 + p_change) * mu_orig
  dmu <- mu_new - get_mean_ziP(lambda_new, theta)
  return(dmu)
}

## dmu^2
get_diff2 <- function(lambda_new, lambda_orig, theta, p_change) {
  ## p_change = -0.1 is a 10% decrease for example
  mu_orig <- get_mean_ziP(lambda_orig, theta)
  mu_new <- (1 + p_change) * mu_orig
  dmu <- mu_new - get_mean_ziP(lambda_new, theta)
  return(dmu ^ 2)
}

rzipch <- function(gamma, p_change) {
  y <- gamma
  
  lambda <- exp(gamma)
  # theta from fit
  theta <- mod$family$getTheta()
  n <- length(y)
  # start if p_change is not 1
  if (p_change != 0) {
    lambda_new <- sapply(lambda, function(z) {
      root <- try(uniroot(
        get_diff,
        interval = c(1e-14, z + 0.5),
        lambda_orig = z,
        theta = theta,
        p_change = p_change
      )$root,
      silent = TRUE)
      ## FIX 3
      ## assum error is due to very small numbers
      if (class(root) == "try-error") {
        1e-14
      } else{
        root
      }
    })
  }
  
  gamma <- log(lambda)
  eta <- theta[1] + exp(theta[2]) * gamma
  p <- 1 - exp(-exp(eta))
  ind <- p > runif(n)
  y[!ind] <- 0
  np <- sum(ind)
  y[ind] <- qpois(runif(np, dpois(0, lambda[ind]), 1), lambda[ind])
  ## FIX 2
  # very small number to 0
  y[!is.finite(y)] <- 0
  y
}


rzipch2 <- function(gamma, p_change) {
  y <- gamma
  
  lambda <- exp(gamma)
  # theta from fit
  theta <- mod$family$getTheta()
  n <- length(y)
  # start if p_change is not 1
  if (p_change != 0) {
    lambda_new <- sapply(lambda, function(z) {
      fit <- optim(
        par = (1 + p_change) * z,
        ## starting value
        fn = get_diff2,
        lower = 0,
        upper = 500,
        lambda_orig = z,
        theta = theta,
        p_change = p_change,
        method = "Brent"
      )
      ## FIX 3 assume error is due to very small numbers
      if (fit$convergence == 0) {
        fit$par
      } else{
        NA
      }
    })
  } else{
    lambda_new <- lambda ## no change in mean
  }
  gamma <- log(lambda_new)
  eta <- theta[1] + exp(theta[2]) * gamma
  p <- 1 - exp(-exp(eta))
  ind <- p > runif(n)
  y[!ind] <- 0
  np <- sum(ind)
  y[ind] <-
    qpois(runif(np, dpois(0, lambda_new[ind]), 1), lambda_new[ind])
  ## FIX 2 hack!
  y[!is.finite(y)] <- 0
  y
}
