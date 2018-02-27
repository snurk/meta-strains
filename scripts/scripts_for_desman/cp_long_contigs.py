import sys
from Bio import SeqIO
 

def cp_long_seq(input_file, output_file, min_len):
    long_sequences = []
    input_handle = open(input_file, 'rU')
    output_handle = open(output_file, "w")
    
    total_sum = 0

    for record in SeqIO.parse(input_handle, "fasta"):
        if len(record.seq) >= min_len:
            long_sequences.append(record)
            total_sum += len(record.seq)
                                                                  
    print("Found %i long sequences with %i bp" % (len(long_sequences), total_sum))

    SeqIO.write(long_sequences, output_handle, "fasta")
    input_handle.close()
    output_handle.close()


if __name__=="__main__":
    cp_long_seq(sys.argv[1], sys.argv[2], int(sys.argv[3]))
