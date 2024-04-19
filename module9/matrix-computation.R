########################################################################
# 1. Avoid unnecessary computations
########################################################################
# Slow
A <- matrix(1:9, nrow = 3)
B <- t(A) %*% A
C <- solve(A) %*% B

# Fast
A <- matrix(1:9, nrow = 3)
B <- crossprod(A)  # Equivalent to t(A) %*% A
C <- solve(A, B)   # Solves A * C = B directly

########################################################################
# 2. Use parallel computation
########################################################################
library(parallel)

n = 10000

A <- matrix(rnorm(n*n), n, n)
B <- matrix(rnorm(n*n), n, n)

# Define a function for matrix multiplication
multiply_block <- function(block_indices) {
  start_row <- block_indices[1]
  end_row <- block_indices[2]
  A[start_row:end_row, ] %*% B
}

# File paths for matrices A and B
file_A <- "A.rds"
file_B <- "B.rds"

# Save matrices A and B to RDS files in the current working directory
saveRDS(A, file = file_A)
saveRDS(B, file = file_B)

# Parallel computation
start_time <- Sys.time()

# Set the number of cores
n_cores <- 10

# Create a cluster and export file paths to worker nodes
cl <- makeCluster(n_cores)
clusterExport(cl, c("file_A", "file_B"))

# Define the block size and create block indices
block_size <- round(n / n_cores)
num_blocks <- ceiling(n / block_size)
block_indices <- matrix(c(seq(1, nrow(A), block_size), pmin(seq(block_size, nrow(A), block_size), nrow(A))), ncol = 2, byrow = FALSE)

# Perform parallel computation
result <- parLapplyLB(cl, split(block_indices, 1:nrow(block_indices)), function(idx) {
  A <- readRDS(file_A)
  B <- readRDS(file_B)
  A[idx[1]:idx[2], ] %*% B
})

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

########################################################################
# 3. Use RcppEigen
########################################################################
library(Rcpp)
library(RcppEigen)

n = 10000

A <- matrix(rnorm(n*n), n, n)
B <- matrix(rnorm(n*n), n, n)

library(microbenchmark)
sourceCpp("matmul.cpp")
microbenchmark(A %*% B, matmult(A, B, n_cores=14), times=1)
