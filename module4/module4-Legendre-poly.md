library(orthopolynom)
poly = legendre.polynomials(n=3, normalized=F)

x = runif(100000, min=-1, max=1)
p = matrix(unlist(polynomial.values(poly, x)), nrow=length(x))

# correlations are close to 0 because of orthogonality
cor(p)

