#!/usr/bin/env python
from __future__ import print_function

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

if len(sys.argv) < 5: 
    print("Usage: %s <number of Ns to insert> <contigs_file> <joined seq name> <output>" % sys.argv[0])
    sys.exit(1) 

f_n = sys.argv[2]
input_seq_iterator = SeqIO.parse(open(f_n, "r"), "fasta")

delim="N" * int(sys.argv[1])

joined_str=delim.join(str(record.seq) for record in input_seq_iterator)
 
output_handle = open(sys.argv[4], "w")

joined=SeqIO.SeqRecord(Seq(joined_str, generic_dna), sys.argv[3], '', '')
SeqIO.write(joined, output_handle, "fasta")
output_handle.close()
