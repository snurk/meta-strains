import sys
import numpy as np


r_matrix = np.loadtxt(sys.argv[1], dtype=int, delimiter=' ')
x_matrix = np.loadtxt(sys.argv[2], dtype=int, delimiter=' ')

# ssm_data
with open(sys.argv[3], "w") as f:
    f.write("id\tgene\ta\td\tmu_r\tmu_v\n")
    for i in range(len(r_matrix)):
        s = "s" + str(i) + '\t' + "g" + str(i) + '\t'
        s += ",".join(map(str, x_matrix[i].tolist())) + '\t'
        s += ",".join(map(str, r_matrix[i].tolist())) + '\t'
        s += "0.999\t0.999\n"
        f.write(s)

# cnv_data
f = open(sys.argv[4], "w")
f.close()