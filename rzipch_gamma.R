## working with gamma
get_mean_ziP <- function(gamma, theta){
    mu <- (1 - exp(-exp(theta[1] + exp(theta[2]) * gamma))) * exp(gamma) * exp(exp(gamma)) / (exp(exp(gamma)) - 1)
    return(mu)
}
## now a function to return the difference in means for a percentage change
get_diff <- function(gamma_new, gamma_orig, theta, p_change){
    ## p_change = -0.1 is a 10% decrease for example
    mu_orig <- get_mean_ziP(gamma_orig, theta)
    mu_new <- (1 + p_change) * mu_orig
    dmu <- mu_new - get_mean_ziP(gamma_new, theta)
    return(dmu)
}
rzipch <- function(gamma, p_change) {
	y <- gamma; 
	# theta from fit
	theta <- mod$family$getTheta()
	n <- length(y)
	# start if p_change is not 1
	if (p_change !=0) {
	# update gamma and lambda
#	gamma <- sapply(gamma, function(z){if(abs(z) < 0.2){uniroot(get_diff, interval = c(-6, 6), 
#					gamma_orig = z, theta = theta, p_change = p_change)$root}else{log(exp(z) * (1 + p_change)) } })
#
	# poisson
	gamma <- sapply(gamma, function(z){log(exp(z) * (1 + p_change))})
	} 
	lambda <- exp(gamma)
	eta <- theta[1] + exp(theta[2])*gamma
	p <- 1- exp(-exp(eta))
	ind <- p > runif(n)
	y[!ind] <- 0
	np <- sum(ind)
	y[ind] <- qpois(runif(np,dpois(0,lambda[ind]),1),lambda[ind])
	y
} 


