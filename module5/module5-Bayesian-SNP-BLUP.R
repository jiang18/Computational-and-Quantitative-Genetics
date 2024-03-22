
# PLINK raw genotype data
geno = read.table("./wheat/wheat.raw",head=T)
# Phenotype data
phen = read.csv("./wheat/wheat.csv")

y = phen$T1
Z = as.matrix(geno[,-(1:6)])
Z = scale(Z, scale=F)       # # Center genotypes. Can help convergence. 
n = nrow(Z)
m = ncol(Z)

# frequency and 2pq
freq = apply(geno[,-(1:6)], 2, mean)/2
het = 2*(1-freq) * freq

# Gibbs sampling
# Load required libraries
require(invgamma)

burnin <- 1000       # number of initial iterations to be disregarded
n_iterations <- 5000 # total number of iterations
thinning <- 10       # use every 10-th sample when calculating mean values

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

  print(paste("Completed MCMC iteration", iter))
}

# MCMC convergence diagnostics
library(coda)

var_alpha = mcmc(var_alpha)
summary(var_alpha)
plot(var_alpha)

var_e = mcmc(var_e)
summary(var_e)
plot(var_e)

first_few_alpha = mcmc(t(alpha[1:3,]))
summary(first_few_alpha)
plot(first_few_alpha)

# Calculate SNP effects and variance components
sample_idx = seq(burnin, n_iterations, thinning)
est.alpha = apply(alpha[, sample_idx], 1, mean)
est.var_e = mean(var_e[sample_idx])
est.var_g = sum(het) * mean(var_alpha[sample_idx])  # genetic variance = sum(2pq) * var_alpha
est.hsq = est.var_g / (est.var_g + est.var_e)

# Calculate GEBVs
gebv = as.matrix(geno[,-(1:6)]) %*% est.alpha
