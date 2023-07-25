import sys

d = {}
g = {}


# create a dictionary of all the header lines
with open(sys.argv[1]) as df:
    for line in df:
        print(line)
        (gene, transcript) = line.strip().split('\t')
        if gene not in g.keys():
            g[gene] = gene
        if transcript.startswith('BambuTx') and gene in g.keys():
            gene = g[gene]
        d[transcript] = (gene, transcript, gene)

# grab the original header lines and replace them with the formatted uniprot-esqu header lines
with open(sys.argv[3], 'w') as ou:
    with open(sys.argv[2]) as f:
        for line in f:
            if line.startswith(">"):
                line=line.strip(">").strip('\n')
                if line.startswith("var"):
                    print(line)
                    line_list=line.split()
                    rf_id=line_list[0].split("_")
                    print(line_list)
                    d[line_list[1]] = (d[line_list[1]][0], f'{rf_id[1]}_{rf_id[2]}', d[line_list[1]][2])
                    line=line_list[1]
                gene_tuple=d[line]
                print(gene_tuple)
                ou.write(f'>db|{gene_tuple[1]}|{gene_tuple[0]} {gene_tuple[2]} OS=Homo sapien OX=9606 GN={gene_tuple[2]}')
                ou.write('\n')
            else:
                ou.write(line)
                    

