
# PLINK raw genotype data
geno = read.table("./wheat/wheat.raw",head=T)
# Phenotype data
phen = read.csv("./training.csv")
colnames(phen)[1] = "IID"

data = merge(phen, geno, by="IID")

y = data$T1
Z = as.matrix(data[,-(1:10)])
Z = scale(Z, scale=F)   # Center genotypes. Can help convergence.

# what we have known
hsq = 0.3
freq = apply(geno[,-(1:6)], 2, mean)/2
het = 2*(1-freq) * freq

# hyperparameters
df0 = 4
scale0 = var(y)*hsq / sum(het) * (df0 - 2) / df0

BayesA = gibbs_bayesa(y, Z, df0, scale0)    # The function of gibbs_bayesa is defined below.


# MCMC convergence diagnostics
library(coda)

var_e = mcmc(BayesA$var_e)
summary(var_e)
plot(var_e)
first_few_alpha = mcmc(t(BayesA$alpha[1:3,]))
summary(first_few_alpha)
plot(first_few_alpha)

# Calculate GEBVs
gebv = as.matrix(geno[,-(1:6)]) %*% BayesA$estimate$alpha
gebv = cbind(geno$IID, gebv)
colnames(gebv) = c("ID", "GEBV")

val = read.csv("./validation.csv")
merged = merge(val, gebv)
cor(merged$T1, merged$GEBV)


##################################################################################
# Gibbs sampling of BayesA for the following model: y = u + Z alpha + e
# Input: y and Z
# Options:
#  hyperparamters (df0 and scale0) in scaled-inv-chisquare prior distribution.
#  burnin is the number of initial iterations to be disregarded.
#  n_iterations is total number of iterations.
#  thinning=10 means every 10-th sample is used for calculating mean values.
# Output: MCMC samples and the estimates
#  Return a list including alpha, var_e, var_alpha, and a nested list of their estimates
##################################################################################
gibbs_bayesa <- function (y, Z, df0=4, scale0=1e-3, n_iterations=5000, burnin=1000, thinning=10) {
  # Gibbs sampling
  # Load required libraries
  require(invgamma)

  n = nrow(Z)
  m = ncol(Z)

  # Initialize parameters
  alpha <- matrix(0, nrow=m, ncol=n_iterations) # Initialize SNP effects
  var_alpha <- matrix(1e-3, nrow=m, ncol=n_iterations) # Initialize SNP effect variance
  var_e <- rep(var(y)*0.7, n_iterations) # Initialize residual variance
  mu <- rep(mean(y), n_iterations) # Initialize population mean

  # Compute err
  err = y - mu[1] - Z %*% alpha[,1]

  for (iter in 1:(n_iterations-1)) {
    # Update alpha and var_alpha given others
    for (j in 1:m) {
      # alpha
      err = err + Z[,j] * alpha[j, iter]
      dj = sum(Z[,j] * Z[,j]) + var_e[iter] / var_alpha[j, iter]
      alpha[j, iter+1] <- rnorm(1, mean = sum(Z[,j] * err)/dj, sd = sqrt(var_e[iter]/dj))
      err = err - Z[,j] * alpha[j, iter+1]

      # var_alpha
      var_alpha[j, iter+1] <- rinvgamma(1, shape = (1 + df0)/ 2, rate = (alpha[j,iter+1]^2 + scale0*df0)/ 2)
    }
    
    # Update mu given others
    err = err + mu[iter]
    mu[iter+1] = rnorm(1, mean = mean(err), sd = sqrt(var_e[iter]/n))
    err = err - mu[iter+1]

    # Update var_e given others
    var_e[iter+1] <- rinvgamma(1, shape = n / 2, rate = sum(err^2) / 2)

    print(paste("Completed MCMC iteration", iter))
  }
  sample_idx = seq(burnin, n_iterations, thinning)
  est.alpha = apply(alpha[, sample_idx], 1, mean)
  est.var_e = mean(var_e[sample_idx])
  est.var_alpha = apply(var_alpha[, sample_idx], 1, mean)

  return(list(alpha=alpha, var_alpha=var_alpha, var_e=var_e, estimate=list(alpha=est.alpha, var_e=est.var_e, var_alpha=est.var_alpha)))
}


