import sys
import numpy as np


def filter_by_coverage(row):
    return (min_coverage < row).all() and (row < max_coverage).all()


r_matrix = np.loadtxt(sys.argv[1], dtype=int, delimiter=' ')
x_matrix = np.loadtxt(sys.argv[2], dtype=int, delimiter=' ')

min_coverage = np.percentile(r_matrix, 25, axis=0)
max_coverage = np.percentile(r_matrix, 75, axis=0)

bool_arr = np.array([filter_by_coverage(row) for row in r_matrix])

r_matrix = r_matrix[bool_arr]
x_matrix = x_matrix[bool_arr]

np.savetxt(sys.argv[3], r_matrix, fmt='%s')
np.savetxt(sys.argv[4], x_matrix, fmt='%s')
