import sys

with open(sys.argv[1]) as f, open(sys.argv[1], "w") as f2:
    i = -1
    for line in f:
        i += 1
        if i % 4 != 0:
            f2.write(line)
            continue
        line = line.split()
        line[0] = line[0][:-2]
        line = " ".join(line)
        f2.write(line + '\n')