# Load the parallel package
library(parallel)

# Function to perform a heavy computation
heavy_computation <- function(x) {
  result <- 0
  for (i in 1:1000000) {
    result <- result + sqrt(x) * i
  }
  return(result)
}

# Create a vector of input values
input_values <- 1:100

# Using lapply (sequential computation)
start_time <- Sys.time()
result_lapply <- lapply(input_values, heavy_computation)
end_time <- Sys.time()
lapply_time <- end_time - start_time

# Using mclapply (parallel computation)
start_time <- Sys.time()
result_mclapply <- mclapply(input_values, heavy_computation, mc.cores = 4)
end_time <- Sys.time()
mclapply_time <- end_time - start_time

# Print the results and execution times
# print("Result using lapply:")
# print(result_lapply)
print(paste("Execution time with lapply:", lapply_time))

# print("Result using mclapply:")
# print(result_mclapply)
print(paste("Execution time with mclapply:", mclapply_time))
