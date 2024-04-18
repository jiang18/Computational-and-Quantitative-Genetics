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
A <- matrix(rnorm(10000 * 10000), nrow = 10000)
B <- matrix(rnorm(10000 * 10000), nrow = 10000)

# Define a function for matrix multiplication
multiply_block <- function(block_indices) {
  start_row <- block_indices[1]
  end_row <- block_indices[2]
  A[start_row:end_row, ] %*% B
}

# Parallel computation
start_time <- Sys.time()

# Create a cluster and export matrices to worker nodes
cl <- makeCluster(8)  # Create a cluster with 8 cores
clusterExport(cl, c("A", "B"))  # Export matrices to worker nodes

# Define the block size and create block indices
block_size <- 1000
num_blocks <- ceiling(nrow(A) / block_size)
block_indices <- matrix(c(seq(1, nrow(A), block_size),
                          pmin(seq(block_size, nrow(A), block_size), nrow(A))),
                        ncol = 2, byrow = TRUE)

# Perform parallel computation
result <- parLapply(cl, split(block_indices, 1:nrow(block_indices)), multiply_block)

# Stop the cluster
stopCluster(cl)

# Convert the result to a matrix
C_parallel <- do.call(rbind, result)

end_time <- Sys.time()
parallel_time <- end_time - start_time

# Serial computation
start_time <- Sys.time()
C_serial <- A %*% B
end_time <- Sys.time()
serial_time <- end_time - start_time

# Compare the results
identical(C_parallel, C_serial)

cat("Parallel time:", parallel_time, "\n")
cat("Serial time:", serial_time, "\n")

# 3. Use RcppEigen
library(SMUT)

# Create example matrices
A <- matrix(rnorm(10000 * 10000), nrow = 10000)
B <- matrix(rnorm(10000 * 10000), nrow = 10000)

# Serial computation using RcppEigen
start_time <- Sys.time()
C_serial <- A %*% B
end_time <- Sys.time()
serial_time <- end_time - start_time

# Parallel computation using RcppEigen
start_time <- Sys.time()
C_parallel <- eigenMapMatMult(A, B)
end_time <- Sys.time()
parallel_time <- end_time - start_time

# Compare the results
identical(C_parallel, C_serial)

cat("Serial time (RcppEigen):", serial_time, "\n")
cat("Parallel time (RcppEigen):", parallel_time, "\n")
