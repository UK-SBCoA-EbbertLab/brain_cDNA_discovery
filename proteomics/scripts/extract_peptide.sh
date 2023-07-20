#!/bin/bash

# Check if input and output file arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: $0 <input_file.fasta> <output_file.txt>"
  exit 1
fi

input_file=$1
output_file=$2

# Check if input file exists
if [ ! -f "$input_file" ]; then
  echo "Input file '$input_file' not found."
  exit 1
fi

# Extract peptide sequences and write them to the output file
grep -v '^>' "$input_file" | awk 'BEGIN{RS=">"; FS="\n"} NR>1{print $1}' > "$output_file"

echo
