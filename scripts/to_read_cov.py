
# coding: utf-8

import argparse
import sys

parser = argparse.ArgumentParser(description="Convert k-mer coverage to read coverage")
parser.add_argument('input',  type=str, help="Input file")
parser.add_argument('output',  type=str, help="Output file") 
parser.add_argument('sample_names',  type=str, help="List of sample names separeted by commas")
parser.add_argument('-r', type=int, help="Read length (default 100)", default=100)
parser.add_argument('-k', type=int, help="Size of K (default 77)", default=77)
args = parser.parse_args()


with open(sys.argv[1], "r") as f_r, open(sys.argv[2], "w") as f_w:
    first_line = False
    for line in f_r:
        line = line.strip().split()
        if not first_line:
            f_w.write("Gene," + args.sample_names + '\n')
            first_line = True
        f_w.write('e'+line[0]+',')
        line = line[1:]
        cov_r = [str(float(cov_k) * args.r / (args.r - args.k + 1)) for cov_k in line]
        f_w.write(','.join(cov_r) + '\n')

