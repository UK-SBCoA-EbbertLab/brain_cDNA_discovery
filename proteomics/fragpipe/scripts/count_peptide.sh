#!/bin/bash
# Check if input and output file arguments are provided
if [ $# -ne 3 ]; then
        echo "Usage: $0 <peptide_input.txt> <search_file.fa> <output_file.txt>"
  exit 1
fi

exec 4<"$1"
echo Start
while read -u4 p ; do
    echo "$p" ":" $(grep -c "$p" $2)>> $3
done
