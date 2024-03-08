library(orthopolynom)
poly = legendre.polynomials(n=1, normalized=F)

env = read.table("rnm.par")
x = env[-c(1:2),1]
# x = scaleX(x, u=-1, v=1)
p = matrix(unlist(polynomial.values(poly, x)), nrow=length(x))

source("../demo/mph_functs.R")
G = read_grm("example1-3-2")

for(i in 1:2) {
    for(j in i:2) {
        if(i == 1 && j == 1) {
            next
        } else if (i == j) {
           newG = ( p[,i] %*% t(p[,j]) ) * G
        } else {
           newG = ( p[,i] %*% t(p[,j]) + p[,j] %*% t(p[,i])) * G
        }
        write_grm(paste(i,j,sep="x"), colnames(newG), newG)
        print(c(i, j))
    }
}


