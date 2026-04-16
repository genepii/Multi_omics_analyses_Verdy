#!/bin/bash

# Run the Python script
python3.9 concatenate_diamond_table.py

# Run the first shell script
bash correct_diamond.sh

# Run the second shell script
bash lca_viral.sh

echo "All scripts have been executed."

