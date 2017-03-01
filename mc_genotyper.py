#!/usr/bin/env python3

from sys import argv




consensus = ''
with open(argv[1], 'r') as gz:
    for line in gz:
        line = line.rstrip('\n')
        if line[0] == '#':
            continue
        else:
            # line = line.rstrip('\n')
            variant = line.split('\t')
            # print(line.decode('utf8').rstrip('\n').split('\t'))
            if len(variant[4]) > 1:
                continue
            else:
                if variant[-1] == '.:.':
                    consensus += '?'
                elif variant[-1] == '0:0':
                    consensus += '-'
                else:
                    cov = int(variant[-1].split(':')[0]) / (
                        int(variant[-1].split(':')[0]) + int(variant[-1].split(':')[1]))
                    if cov > 0.9:
                        consensus += variant[3]
                    elif cov < 0.1:
                        consensus += variant[4]
                    else:
                        consensus += 'N'

print(consensus)
