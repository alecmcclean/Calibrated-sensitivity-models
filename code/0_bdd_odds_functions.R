####################################################
### Functions for bounded odds ratio analysis
####################################################

expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(p) log(p / (1 - p))

aplus <- function(a) pmax(a, 0)
aminus <- function(a) -1 * pmin(a, 0)


###################################
### Calculate thetas and nus

# Loss function
theta_optim <- function(coef_theta, y, X, Gamma) {
  
  preds <- as.matrix(X) %*% as.matrix(coef_theta) 
  loss_vals <- 0.5 * aplus(y - preds)^2 + Gamma * aminus(y - preds)^2
  return(mean(loss_vals))
  
}

nu_optim <- function(coef_nu, y, X, theta, Gamma) {
  
  preds <- as.matrix(X) %*% as.matrix(coef_nu)
  loss_vals <- 0.5 * (1 + (Gamma - 1) * ifelse(y > theta, 1, 0) - preds)^2
  return(mean(loss_vals))
  
}

CalculateIFlow <- function(A, Y, theta1, theta0, nu1, nu0, pi, Gamma) {
  psi1 <- aplus(Y - theta1) - Gamma * aminus(Y - theta1)
  psi0 <- Gamma * aplus(Y - theta0) - aminus(Y - theta0)
  treated <- A * Y + (1 - A) * theta1 + A * (psi1 * (1 - pi) / (nu1 * pi))
  control <- (1 - A) * Y + A * theta0 + (1 - A) * (psi0 * pi / (nu0 * (1 - pi)))
  return(treated - control)  
}

CalculateIFhigh <- function(A, Y, theta1, theta0, nu1, nu0, pi, Gamma) {
  psi1 <- Gamma * aplus(Y - theta1) - aminus(Y - theta1)
  psi0 <- aplus(Y - theta0) - Gamma * aminus(Y - theta0)
  treated <- A * Y + (1 - A) * theta1 + A * (psi1 * (1 - pi) / (nu1 * pi))
  control <- (1 - A) * Y + A * theta0 + (1 - A) * (psi0 * pi / (nu0 * (1 - pi)))
  return(treated - control)  
}
