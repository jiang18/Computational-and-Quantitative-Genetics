
# PLINK raw genotype data
geno = read.table("./wheat/wheat.raw",head=T)
# Phenotype data
phen = read.csv("./wheat/wheat.csv")

y = phen$T1
Z = as.matrix(geno[,-(1:6)])
Z = scale(Z, scale=F)       # Center genotypes. Can help convergence. 
n = nrow(Z)
m = ncol(Z)

# frequency and 2pq
freq = apply(geno[,-(1:6)], 2, mean)/2
het = 2*(1-freq) * freq

blup = gibbs_blup(y, Z)    # The function of gibbs_blup is defined below.

# MCMC convergence diagnostics
library(coda)
var_e = mcmc(blup$var_e)
summary(var_e)
plot(var_e)
first_few_alpha = mcmc(t(blup$alpha[1:3,]))
summary(first_few_alpha)
plot(first_few_alpha)

# Calculate h2 estimate
est.var_g = sum(het) * blup$estimate$var_alpha  # genetic variance = sum(2pq) * var_alpha
est.hsq = est.var_g / (est.var_g + blup$estimate$var_e)

# Calculate GEBVs
gebv = as.matrix(geno[,-(1:6)]) %*% blup$estimate$alpha

##################################################################################
# Gibbs sampling for the following model: y = u + Z alpha + e
# Input: y and Z
# Options:
#  burnin is the number of initial iterations to be disregarded.
#  n_iterations is total number of iterations.
#  thinning=10 means every 10-th sample is used for calculating mean values.
# Output: MCMC samples for alpha, var_e, var_alpha, and their estimates
#  Return a list.
##################################################################################
gibbs_blup <- function (y, Z, n_iterations=5000, burnin=1000, thinning=10) {
  # Load required libraries
  require(invgamma)
  
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
  sample_idx = seq(burnin, n_iterations, thinning)
  est.alpha = apply(alpha[, sample_idx], 1, mean)
  est.var_e = mean(var_e[sample_idx])
  est.var_alpha = mean(var_alpha[sample_idx])

  return(list(alpha=alpha, var_alpha=var_alpha, var_e=var_e, estimate=list(alpha=est.alpha, var_e=est.var_e, var_alpha=est.var_alpha)))
}

