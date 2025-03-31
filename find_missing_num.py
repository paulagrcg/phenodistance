import os
import re

# Define the folder path and filename pattern
data_folder = "/home/pg520/phenodistance/data"
pattern = re.compile(r"to_analyse(\d+)")

# Get all numbers from filenames matching the pattern
numbers_in_files = []
for filename in os.listdir(data_folder):
    match = pattern.match(filename)
    if match:
        numbers_in_files.append(int(match.group(1)))

# Define the full range of numbers
full_range = set(range(0, 12459))

# Find missing numbers
numbers_in_files = set(numbers_in_files)
missing_numbers = sorted(full_range - numbers_in_files)

# Output the missing numbers
print("Missing numbers:", missing_numbers, len(missing_numbers))
# Save the missing numbers to a npy
import numpy as np
np.save("missing_numbers.npy", missing_numbers)
