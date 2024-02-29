# Simulation
n = 1000
x = runif(n, min=-1, max=1)
x = scale(x, scale=T)
e = rnorm(n)
b = 2
y = x * b + e

# Linear regression
fit = lm(y ~ x - 1)
summary(fit)

# MCMC
require(invgamma)

niter = 10000
burnin = 1000
step = 10
mcmc = data.frame(b=rep(NA, niter), vare=rep(NA, niter))
# initial values
mcmc$b[1] = 0.1
mcmc$vare[1] = 0.1

b_mean = sum(x*y) / sum(x*x)
x_sumsq = sum(x*x)
for(i in 2:niter) {
    mcmc$b[i] = rnorm(n=1, mean=b_mean, sd=sqrt(mcmc$vare[i-1]/x_sumsq))

    error = y - x *  mcmc$b[i]
    mcmc$vare[i] = rinvgamma(n=1, shape=n/2-1, rate=sum(error*error)/2)
}

# Check MCMC samples
plot(xlim=c(1,niter), ylim=c(min(mcmc$b),max(mcmc$b)), type = "n", x=0, y= 0, xlab="Iteration", ylab="Sample")
points(1:niter, mcmc$b)

plot(xlim=c(1,niter), ylim=c(min(mcmc$vare),max(mcmc$vare)), type = "n", x=0, y= 0, xlab="Iteration", ylab="Sample")
points(1:niter, mcmc$vare)

hist(mcmc$vare[seq(burnin, niter, step)])
