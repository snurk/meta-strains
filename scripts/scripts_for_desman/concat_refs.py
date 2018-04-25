import sys
import os
from Bio import SeqIO
from collections import Counter


output = sys.argv[1]
inputs = sys.argv[2:]

with open(output, "w") as out_file:
    for cur_input in inputs:
        cur_strain = os.path.splitext(os.path.basename(cur_input))[0]
        fasta_sequences = SeqIO.parse(open(cur_input), 'fasta')

        counter = 1
        for fasta in fasta_sequences:
            fasta.id = cur_strain+'_'+str(counter)
            fasta.description = ""
            #fasta.seq = str(fasta.seq).replace('K', 'G')
            print(Counter(str(fasta.seq)))
            SeqIO.write(fasta, out_file, "fasta")
            counter += 1 
