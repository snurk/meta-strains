import sys
import numpy as np


def filter_by_coverage(row):
    return (min_coverage < row).all() and (row < max_coverage).all()


r_matrix = np.loadtxt(sys.argv[1], dtype=int, delimiter=' ')
x_matrix = np.loadtxt(sys.argv[2], dtype=int, delimiter=' ')

print("Sites before filtering:", len(r_matrix))

min_coverage = np.percentile(r_matrix, 25, axis=0)
max_coverage = np.percentile(r_matrix, 75, axis=0)

bool_arr = np.array([filter_by_coverage(row) for row in r_matrix])
r_matrix = r_matrix[bool_arr]
x_matrix = x_matrix[bool_arr]
print("Sites after coverage filtering:", len(r_matrix))

ind = [i for i in range(len(x_matrix)) if not ((np.abs(r_matrix[i, :] - x_matrix[i, :]) <= 3).all())]
r_matrix = r_matrix[ind]
x_matrix = x_matrix[ind]
print("Sites after filtering non-sense:", len(r_matrix))

ind = [i for i in range(len(x_matrix)) if not ((x_matrix[i, :] <= 3).all())]
r_matrix = r_matrix[ind]
x_matrix = x_matrix[ind]
print("Sites after filtering low coverage:", len(r_matrix))

# chose about 200 random sites
# if len(r_matrix) > 200:
#     sample = np.random.choice(len(r_matrix), size=200)
#     sample = np.unique(sample)
#     r_matrix = r_matrix[sample]
#     x_matrix = x_matrix[sample]

np.savetxt(sys.argv[3], r_matrix, fmt='%s')
np.savetxt(sys.argv[4], x_matrix, fmt='%s')
