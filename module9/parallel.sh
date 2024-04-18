#!/bin/bash

# Function to perform a simple computation
compute() {
  local index=$1
  local result=$((index * index))
  echo "Result for index $index: $result"
}

# Number of parallel tasks
num_tasks=5

# Perform parallel computation
for ((i=1; i<=num_tasks; i++)); do
  compute $i & 
done

# Wait for all tasks to complete
wait

echo "All tasks completed."
