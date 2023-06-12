ensembl_file="../gtfs/Homo_sapiens.GRCh38.105.gtf"
outfile="../gtfs/Homo_sapiens.GRCh38.105.simplified.gtf"

with open(outfile, 'w') as out_ensembl:
    with open(ensembl_file) as ensembl:
        for line in ensembl:
            split_line = line.strip().split('\t')
            if line.startswith("#"):
                continue
            elif split_line[2] == "transcript" or split_line[2] == "exon":
                gene_id = ""
                transcript_id = ""
                exon_number = ""
                attributes = split_line[8].split(";")
                for a in attributes:
                    a = a.strip()
                    if a.startswith("gene_id"):
                        gene_id = a + ";"
                    elif a.startswith("transcript_id"):
                        transcript_id = a + ";"
                    elif a.startswith("exon_number"):
                        exon_number = a + ";"
                    else: 
                        continue

                if exon_number != "":
                    split_line[8] = " ".join([gene_id, transcript_id, exon_number])
                else:
                    split_line[8] = " ".join([gene_id, transcript_id])

                out_ensembl.write('\t'.join(split_line) + '\n')



