#!/usr/bin/env python3

import subprocess
import os
from sys import argv

def remove_dummy(vcf_in, vcf_out):
    header = []
    variants =[]
    SNVsites = []
    with open(vcf_in, 'r') as gz:
        for line in gz:
            line = line.rstrip('\n')
            if line[1] == '#':
                header.append(line.rstrip('\n'))
            elif line[1] =='C':
                    columns = line.rstrip('\n')
            else:
                # line = line.rstrip('\n')
                variant = line.split('\t')
                variants.append(variant)
    with open(vcf_out, 'w') as out:
        for h in header:
            out.write(h)
        out.write('\t'.join(columns.split('\t')[:9]) + '\t' + columns.split('\t')[10])
        for v in variants:
            out.write('\t'.join(v[:9]) + '\t' + v[10])



def geno_cov(mccortex, reference, calls_vcf, raw_ctx, outdir):
    ref = os.path.basename(reference)
    sample = os.path.basename(raw_ctx).rstrip('.ctx')
    subprocess.call(
        [mccortex + " vcfcov -q -n 50M -m 1G -r " + outdir + "/ref/" + ref + " -f -o " + outdir + "/" + sample + ".cov.vcf " +
         calls_vcf + " " + outdir + "/raw/" + raw_ctx + ".ctx" ], stdout=subprocess.PIPE, shell=True)

def main():
    mccortex = argv[1]
    reference = argv[2]
    calls_vcf = argv[3]
    raw_ctx = argv[4]
    outdir = argv[5]
    geno_cov(mccortex, reference, calls_vcf, raw_ctx, outdir)

if __name__ == '__main__':
    main()
