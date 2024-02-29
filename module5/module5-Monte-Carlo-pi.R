# Function to estimate pi
estimate_pi <- function(num_points) {
  # Generate random points within a square of side 2
  points.x <- runif(n = num_points, min = -1, max = 1)
  points.y <- runif(n = num_points, min = -1, max = 1)

  # Calculate distance from the origin for each point
  distances <- sqrt(points.x^2 + points.y^2)
  
  # Count points within the circle (radius 1)
  in_circle <- sum(distances <= 1)
  
  # Estimate pi using the ratio of points inside the circle and total points
  pi_estimate <- 4 * in_circle / num_points
  return(pi_estimate)
}

# Set number of points for simulation
num_points <- 100e6

# Estimate pi
pi_est <- estimate_pi(num_points)

# Print the estimated value
cat("Estimated value of pi:", pi_est, "\n")
