# 1. Avoid unnecessary computations
# Slow
A <- matrix(1:9, nrow = 3)
B <- t(A) %*% A
C <- solve(A) %*% B

# Fast
A <- matrix(1:9, nrow = 3)
B <- crossprod(A)  # Equivalent to t(A) %*% A
C <- solve(A, B)   # Solves A * C = B directly


# 2. Use parallel computation
library(parallel)

# Create example matrices
A <- matrix(rnorm(1000000), nrow = 1000)
B <- matrix(rnorm(1000000), nrow = 1000)

# Define a function for matrix multiplication
multiply_row <- function(i) {
  A[i, ] %*% B
}

# Create a cluster and export matrices to worker nodes
cl <- makeCluster(4)  # Create a cluster with 4 cores
clusterExport(cl, c("A", "B"))  # Export matrices to worker nodes

# Perform parallel computation
result <- parLapply(cl, seq_len(nrow(A)), multiply_row)

# Stop the cluster
stopCluster(cl)

# Convert the result to a matrix
C <- do.call(rbind, result)

