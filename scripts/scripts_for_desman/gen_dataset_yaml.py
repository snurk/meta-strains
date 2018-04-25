import sys

n = int(sys.argv[1])
left_reads = sys.argv[2:n+2]
right_reads = sys.argv[n+2:2*n+2]

for i in range(n):
    print('- "left reads":')
    print('  - "%s"' % left_reads[i])
    print('  "orientation": "fr"')
    print('  "right reads":')
    print('  - "%s"' % right_reads[i])
    print('  "type": "paired-end"')

