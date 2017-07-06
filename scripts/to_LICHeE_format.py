import sys
import numpy as np


r_matrix = np.loadtxt(sys.argv[1], dtype=int, delimiter=' ')
x_matrix = np.loadtxt(sys.argv[2], dtype=int, delimiter=' ')

#r_matrix = np.delete(r_matrix, [2, 5, 6, 10], axis=1)
#x_matrix = np.delete(x_matrix, [2, 5, 6, 10], axis=1)

num_of_samples = len(r_matrix[0])

with open(sys.argv[3], "w") as f:
    f.write("#chr\tposition\tdescription\tNormal\t"
            + "\t".join(["S" + str(i) for i in range(1, num_of_samples+1)])
            + "\n")
    for i in range(len(r_matrix)):
        s = "0" + '\t' + str(i) + '\t' + "None" + '\t' + "0.0" + '\t'
        s += "\t".join(map(str, (x_matrix[i] / r_matrix[i]).tolist())) + '\n'
        f.write(s)


# #chr    position    description    Normal    S1    S2    S3    S4
# 17      123456      A/T DUSP19     0.0       0.1   0.2   0.25  0.15
# 11      341567      C/G MUC16      0.0       0.4   0.09  0.38  0.24
# 9       787834      A/C OR2A14     0.0       0.35  0.14  0.17  0.48