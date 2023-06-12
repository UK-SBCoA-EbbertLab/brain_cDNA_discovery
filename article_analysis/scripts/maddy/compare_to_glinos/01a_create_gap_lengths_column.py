import sys
import os

# the merged bedfile is passed
merged_bedfile_path = sys.argv[1]
outfile = os.path.splitext(merged_bedfile_path)[0] + ".txt"

with open(outfile, "w") as outf:
    with open(merged_bedfile_path) as bedfile:
        previous_end = -1
        current_chrom = -1
        # for each line in the bedfile, we want to calculate the distance from the previous feature to the current feature
        for line in bedfile:
            line = line.strip().split()
            # if the current chromosome doesn't match the previous one, we know that we have moved on and we don't calculate a gap (since there isn't a previous feature to compare to)
            if current_chrom != line[0]:
                current_chrom = line[0]
                previous_end = int(line[2])
                st_end_diff = "NA"
            else:
                st_end_diff = str(int(line[1]) - previous_end)
                previous_end = int(line[2])
            line.append(st_end_diff)
            outf.write("\t".join(line) + "\n")



