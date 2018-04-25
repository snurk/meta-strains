import sys
from Bio import SeqIO
 

def cp_long_seq(input_file, output_file):
    input_handle = open(input_file, 'rU')
    output_handle = open(output_file, "w")
    
    for record in SeqIO.parse(input_handle, "fasta"):
        if int(record.id) % 2 == 0:
            output_handle.write(record.id+'\t'+str(len(record.seq))+'\n')                        
    
    input_handle.close()
    output_handle.close()


if __name__=="__main__":
    cp_long_seq(sys.argv[1], sys.argv[2])
