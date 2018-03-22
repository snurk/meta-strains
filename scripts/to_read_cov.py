
# coding: utf-8

import sys

with open(sys.argv[1], "r") as f_r, open(sys.argv[2], "w") as f_w:
    f_w.write("Gene," + ','.join(["sample"+str(i) for i in range(1,12)]) + '\n')
    for line in f_r:
        line = line.strip().split()
        f_w.write('e'+line[0]+',')
        line = line[1:]
        cov = [str(float(k) * 100 / 46) for k in line]
        f_w.write(','.join(cov) + '\n')

