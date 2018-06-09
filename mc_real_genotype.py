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
            out.write(h + '\n')
        out.write('\t'.join(columns.split('\t')[:9]) + '\t' + columns.split('\t')[10] + '\n')
        for v in variants:
            out.write('\t'.join(v[:9]) + '\t' + v[10] + '\n')


def average_kcov(infile):
    consensus = ''
    with open(infile, 'r') as gz:
        kcovs=[]
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
                        continue
                    elif variant[-1] == '0:0':
                        continue
                    else:
                        kcovs.append(int(variant[-1].split(':')[0]) + int(variant[-1].split(':')[1]))

    return sum(kcovs)/len(kcovs)

def geno_cov(mccortex, reference, calls_vcf, raw_ctx, outdir):
    ref = os.path.basename(reference)
    sample = os.path.basename(raw_ctx).rstrip('.ctx')
    subprocess.call(
        [mccortex + " vcfcov -q -n 50M -m 1G -r " + outdir + "/ref/" + ref + " -f -o " + outdir + "/" + sample + ".cov.vcf " +
         calls_vcf + " " + outdir + "/raw/" + raw_ctx + ".ctx" ], stdout=subprocess.PIPE, shell=True)
    kcov = average_kcov(outdir + "/" + sample + ".cov.vcf")
    remove_dummy(outdir + "/" + sample + ".cov.vcf", outdir + "/" + sample + ".cov1.vcf")
    subprocess.call(
        [mccortex + " vcfgeno -q -f --out " + outdir + "/" + sample + ".geno.vcf --kcov " + str(kcov) + " --ploidy 2 " +
         outdir + "/" + sample + ".cov1.vcf"], stdout=subprocess.PIPE, shell=True)

def main():
    mccortex = argv[1]
    reference = argv[2]
    calls_vcf = argv[3]
    raw_ctx = argv[4]
    outdir = argv[5]
    geno_cov(mccortex, reference, calls_vcf, raw_ctx, outdir)

if __name__ == '__main__':
    main()
