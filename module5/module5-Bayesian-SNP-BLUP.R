
# Generate some example data
n <- 100 # Number of individuals
m <- 1000 # Number of SNPs
Z <- matrix(rnorm(n * m), nrow = n) # SNP matrix (design matrix)
mu = 1
alpha_true <- rnorm(m) # True SNP effects
y <- mu + Z %*% alpha_true + rnorm(n) # Phenotype data

# Chr 1 of QTL-MAS 2012 data
phen = read.csv("phen.csv")
y = phen$milk
geno = read.table("chr1.raw",head=T)
Z = as.matrix(geno[,-(1:6)])
n = nrow(Z)
m = ncol(Z)

# Gibbs sampling
# Load required libraries
require(invgamma)

n_iterations <- 1000

# Initialize parameters
alpha <- matrix(0, nrow=m, ncol=n_iterations) # Initialize SNP effects
var_alpha <- rep(1e-2, n_iterations) # Initialize SNP effect variance
var_e <- rep(1, n_iterations) # Initialize residual variance
mu <- rep(1, n_iterations) # Initialize population mean

# Compute err
err = y - mu[1] - Z %*% alpha[,1]

for (iter in 1:(n_iterations-1)) {
  # Update alpha given others
  for (j in 1:m) {
    err = err + Z[,j] * alpha[j, iter]
    dj = sum(Z[,j] * Z[,j]) + var_e[iter] / var_alpha[iter]
    alpha[j, iter+1] <- rnorm(1, mean = sum(Z[,j] * err)/dj, sd = sqrt(var_e[iter]/dj))
    err = err - Z[,j] * alpha[j, iter+1]
  }
  
  # Update mu given others
  err = err + mu[iter]
  mu[iter+1] = rnorm(1, mean = mean(err), sd = sqrt(var_e[iter]/n))
  err = err - mu[iter+1]

  # Update var_e given others
  var_e[iter+1] <- rinvgamma(1, shape = n / 2, rate = sum(err^2) / 2)

  # Update var_alpha given others
  var_alpha[iter+1] <- rinvgamma(1, shape = m / 2, rate = sum(alpha[,iter+1]^2) / 2)
}

# MCMC convergence diagnostics
library(coda)

var_alpha = mcmc(var_alpha)
summary(var_alpha)
plot(var_alpha)

first_few_alpha = mcmc(alpha[,1:3])
summary(first_few_alpha)
plot(first_few_alpha)

