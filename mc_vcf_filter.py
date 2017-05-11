#!/usr/bin/env python3

from sys import argv


header = []
variants =[]
SNVsites = []
with open(argv[1], 'r') as gz:
    for line in gz:
        line = line.rstrip('\n')
        if line[0] == '#':
            header.append(line.rstrip('\n'))
        else:
            # line = line.rstrip('\n')
            variant = line.split('\t')
            # print(line.decode('utf8').rstrip('\n').split('\t'))
            try:
                if int(variant[-1].split(':')[1]) == 0:
                    continue
                else:
                    if (len(variant[4]) == 1 and len(variant[3]) == 1) and (
                            int(variant[-1].split(':')[0]) / int(variant[-1].split(':')[1]) < 0.1):
                        variants.append(variant)
                        SNVsites.append(int(variant[1]))
            except ValueError:
                continue
for h in header:
    print(h)
for v in variants:
    print('\t'.join(v))
