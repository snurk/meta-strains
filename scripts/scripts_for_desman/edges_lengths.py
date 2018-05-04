import sys

def cp_long_seq(input_file, output_file):
    input_handle = open(input_file, 'rU')
    output_handle = open(output_file, "w")
    

    for line in input_handle:
        line = line.strip().split()
        line_type = line[0]
        if line_type != 'S':
            break
        e_id = line[1]
        e_seq = line[2]

        output_handle.write(e_id+'\t'+str(len(e_seq))+'\n')                        
    
    input_handle.close()
    output_handle.close()


if __name__=="__main__":
    cp_long_seq(sys.argv[1], sys.argv[2])
