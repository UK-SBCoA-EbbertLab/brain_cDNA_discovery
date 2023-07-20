#!/bin/bash
# Check if input and output file arguments are provided
if [ $# -ne 3 ]; then
	echo "Usage: $0 <input_file.txt> <input_file2.txt> <output_file.fa>"
  exit 1
fi
# Open the files for reading
exec 3<$1
exec 4<$2

while read -u 3 t && read -u 4 p; do
    echo "$t" >> $3
    echo "$p" >> $3
done

# Close the file descriptors
exec 3<&-
exec 4<&-

