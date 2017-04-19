#!/usr/bin/env python3


import os
import subprocess
from sys import argv


def bubble_call(clean_binary, reference_binary, reference, outdir):
    #bubble calls cleaned binaries in parallel and does filtering of resulting vcfs
    bin = os.path.basename(clean_binary)
    ref_bwa = os.path.basename(reference)
    subprocess.call(["cd " + outdir + "; mccortex63 join -m 8G -o " + clean_binary + "_plus_" + reference_binary +" ref/" + reference_binary +
                     " clean/" + clean_binary], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; mccortex63 bubbles -m 8G -o " + clean_binary + "_bubbles.txt.gz --haploid 0 -S -f " +
                     clean_binary + "_plus_" + reference_binary], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; bash mccortex_print_flanks.sh " + clean_binary + "_bubbles.txt.gz " +
                     "> " + clean_binary + ".flanks" ], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; bwa mem ref/" + ref_bwa + " "+ clean_binary + ".flanks > " +
                    clean_binary +".sam"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; mccortex63 calls2vcf -F" + clean_binary +".sam -o " + clean_binary +
                     "_bub.vcf " + clean_binary + "_bubbles.txt.gz ref/" + ref_bwa], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; vcf-sort " + clean_binary + "_bub.vcf | vcfuniq >" + clean_binary + "_bub.uniq.vcf"],
                    stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; mccortex63 vcfcov -m 8G --out " + clean_binary + "_bub.cov.vcf --ref ref/" +
                     ref_bwa + " " + clean_binary + "_bub.uniq.vcf clean/" + clean_binary], stdout=subprocess.PIPE, shell=True)

    subprocess.call(["cd " + outdir + "; mc_vcf_filter.py " + clean_binary + "_bub.cov.vcf > " + clean_binary +
                     ".final.vcf"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cd " + outdir + "; rm clean/" + clean_binary + " " + clean_binary + "_plus_" + reference_binary + " " +
                     clean_binary + "_bubbles.txt.gz " + clean_binary + ".flanks " + clean_binary + ".sam " +
                     clean_binary + "_bub.uniq.vcf " + clean_binary + "_bub.vcf " + clean_binary + "_bub.cov.vcf "],
                     stdout=subprocess.PIPE, shell=True)

def main():
    sample_ctx = argv[1]
    reference_ctx = argv[2]
    reference = argv[3]
    outdir = argv[4]
    bubble_call(sample_ctx, reference_ctx, reference, outdir)

if __name__ == '__main__':
    main()
