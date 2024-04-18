#!/bin/bash

# The directory containing the files
source_dir="/lustre19/expphy/volatile/halla/sbs/seeds/simc/p657sf_sbs9_sbs70p_simc/simcout"

# The directory to move the files to
target_dir="/lustre19/expphy/volatile/halla/sbs/seeds/simc/p657sf_sbs9_sbs70p_simc/trash_run"

# The number threshold, provided as the first command-line argument
number_threshold=$1

# Check if the number threshold was provided
if [ -z "$number_threshold" ]; then
  echo "Please provide a number threshold as an argument."
  exit 1
fi

echo "Files to be moved:"

# Array to hold files that match the criteria
declare -a files_to_move

# Iterate over files in the source directory
for file in "$source_dir"/*; do
  # Extract the number followed by a period from the file name
  # Adjusted the regular expression to match a number followed directly by a period
  number=$(echo "$file" | grep -oE '[0-9]+\.' | grep -oE '[0-9]+' | head -1)
  
  # Check if the number is greater than the provided threshold
  if [[ "$number" -gt "$number_threshold" ]]; then
    # Add file to the list
    files_to_move+=("$file")
    echo "$file"
  fi
done

# Check if there are any files to move
if [ ${#files_to_move[@]} -eq 0 ]; then
  echo "No files to move."
  exit 0
fi

# Show target directory
echo "Target directory: $target_dir"

# Ask for user confirmation
read -p "Are you sure you want to move the above files? (y/n) " answer

# Check user's answer
case $answer in
  [Yy]* )
    for file in "${files_to_move[@]}"; do
      mv "$file" "$target_dir"
    done
    echo "Files have been moved."
    ;;
  [Nn]* )
    echo "Operation canceled."
    exit 0
    ;;
  * )
    echo "Please answer yes or no."
    ;;
esac
