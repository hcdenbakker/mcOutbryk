#!/usr/bin/env python3

from sys import argv


header = []
variants =[]
SNVsites = []
with open(argv[1], 'r') as gz:
    for line in gz:
        line = line.rstrip('\n')
        if line[1] == '#':
            header.append(line.rstrip('\n'))
        elif line[1] =='C':
                columns = line.rstrip('\n')
        else:
            # line = line.rstrip('\n')
            variant = line.split('\t')
            variants.append(variant[:10])
for h in header:
    print(h)
print('\t'.join(columns.split('\t')[:10]))
for v in variants:
    print('\t'.join(v))
