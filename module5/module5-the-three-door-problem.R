# Function to simulate a single game
simulate_game <- function(switch = TRUE) {
  # Place car randomly behind one of the doors
  car_door <- sample(1:3, 1)
  
  # Player chooses a door randomly
  chosen_door <- sample(1:3, 1)
  
  # Host reveals a goat behind one of the remaining doors
  revealed_door <- ifelse(chosen_door == car_door, sample(setdiff(1:3, chosen_door), 1), setdiff(1:3, c(chosen_door, car_door)))

  # Player either stays or switches doors
  if (switch) {
    chosen_door <- setdiff(1:3, c(chosen_door, revealed_door))[1]
  }
  
  # Return TRUE if player wins, FALSE otherwise
  return(chosen_door == car_door)
}

# Set number of simulations
n_simulations <- 100000

# Simulate games with switching and without switching
wins_switch <- sum(sapply(rep(TRUE,n_simulations), simulate_game))
wins_stay <- sum(sapply(rep(FALSE,n_simulations), simulate_game))

# Calculate win percentages
win_rate_switch <- wins_switch / n_simulations
win_rate_stay <- wins_stay / n_simulations

# Print results
cat("Win rate with switching:", win_rate_switch, "\n")
cat("Win rate with staying:", win_rate_stay, "\n")
