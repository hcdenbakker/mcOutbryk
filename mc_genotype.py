#!/usr/bin/env python3

import subprocess
import os
from sys import argv

def geno_cov(reference, calls_vcf, raw_ctx, outdir):
    ref = os.path.basename(reference)
    sample = os.path.basename(raw_ctx).rstrip('.ctx')
    subprocess.call(
        ["mccortex63 vcfcov -m 12G -n 500M -r " + outdir + "/ref/" + ref + " -f -o " + outdir + "/" + sample + ".cov.vcf " +
         calls_vcf + " " + outdir + "/raw/" + raw_ctx + ".ctx" ], stdout=subprocess.PIPE, shell=True)

def main():
    reference = argv[1]
    calls_vcf = argv[2]
    raw_ctx = argv[3]
    outdir = argv[4]
    geno_cov(reference, calls_vcf, raw_ctx, outdir)

if __name__ == '__main__':
    main()
