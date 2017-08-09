import sys


with open(sys.argv[1], "r") as file_varscan, \
        open(sys.argv[2], "w") as file_R_matrix, open(sys.argv[3], "w") as file_X_matrix:

    num_of_sites = 0
    num_of_bad_sites = 0

    for line in file_varscan.readlines()[1:]:
        num_of_sites += 1

        columns = line.strip().split()
        samples = columns[10:]
        cur_r, cur_x = [], []
        flag_2_alleles = True

        for sample in samples:
            info = sample.split(':')

            if info[1] == '-' or info[2] == '-' or info[3] == '-':
                flag_2_alleles = False
                break

            if int(info[1]) - (int(info[2]) + int(info[3])) == 0:  # <= 3:
                cur_r.append(info[1])
                cur_x.append(info[3])
            else:
                flag_2_alleles = False
                break

        if flag_2_alleles:
            file_R_matrix.write(" ".join(cur_r) + '\n')
            file_X_matrix.write(" ".join(cur_x) + '\n')
        else:
            num_of_bad_sites += 1

print("Number of sites:", num_of_sites)
print("Number of bad sites:",  num_of_bad_sites)
print("Number of biallelic sites:",  num_of_sites - num_of_bad_sites)
