#!/usr/bin/env python3

# todo: create log
# todo: flag, log, exclude highly divergent samples (> 10,000 SNPs)


import os
import re
import subprocess
from argparse import (ArgumentParser)


def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='MCOutbryk - de novo SNP and ambiguous site Caller')
    # inputs
    parser.add_argument(
        '--reference', nargs='+', type=str, required=True,
        help='Reference genome in fasta format, include path to file if file is not in same directory ')
    parser.add_argument(
        '--sample_list', nargs='+', type=str, required=True,
        help='Text file, one isolate/sample per line ')
    parser.add_argument(
        '--results_dir', nargs='+', type=str, required=True,
        help='Directory which will contain results after the analysis ')
    parser.add_argument(
        '--processes', nargs='+', type=str, required=True,
        help='number of processes to run in parallel')
    return parser.parse_args()


def perform_actions_list(sample_list, outdir):
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    subprocess.call(["rm " + outdir + "/to_be_processed.txt"], stdout=subprocess.PIPE, shell=True)
    to_be_processed = open(outdir + '/to_be_processed.txt', 'w')
    for sample in sample_list:
        regex1 = re.compile(sample + '_S._L001_R1_001.fastq.gz')
        regex2 = re.compile(sample + '_S.._L001_R1_001.fastq.gz')
        matches1 = [string for string in files if re.match(regex1, string)]
        matches2 = [string for string in files if re.match(regex2, string)]
        if len(matches1) == 1:
            to_be_processed.write(sample + '\t' + matches1[0] + '\t' + matches1[0].rstrip('_R1_001.fastq.gz') + '001_R2_001.fastq.gz\n')
        if len(matches2) == 1:
            to_be_processed.write(sample + '\t' + matches2[0] + '\t' + matches2[0].rstrip('_R1_001.fastq.gz') + '001_R2_001.fastq.gz\n')
        if str(sample + '_1.fastq.gz') in files:
            to_be_processed.write(sample + '\t' + sample + '_1.fastq.gz' + '\t' + sample + '_2.fastq.gz\n')

    to_be_processed.close()


def create_ref_binary(outdir, reference):
    ref = os.path.basename(reference).rstrip('.fasta')
    ref_bwa = os.path.basename(reference)
    subprocess.call(
        ["mccortex63 build -s " + ref + " -k 33 -m 8GB -f -t 32 -1 " + reference + " " + outdir + "/ref/" + ref
         + ".ctx"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cp " + reference + " " + outdir + "/ref; bwa index " + outdir + "/ref/" + ref_bwa],
                    stdout=subprocess.PIPE, shell=True)


def get_list(infile):
    samples = []
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        samples.append(line.strip('\n'))
    return samples


def create_consensus(infile):
    consensus = ''
    with open(infile, 'r') as gz:
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

    return consensus


def main():
    args = parse_args()
    reference = args.reference[0]
    samples = get_list(args.sample_list[0])
    sample_list = args.sample_list[0]
    outdir = args.results_dir[0]
    num_procs = args.processes[0]
    ref = os.path.basename(reference).rstrip('.fasta')
    ref_ctx = str(ref + ".ctx")
    create_ref_binary(outdir, reference)
    perform_actions_list(samples, outdir)
    subprocess.call(["cat " + outdir + "/to_be_processed.txt | parallel -j " +
                     num_procs + " --colsep '\t' mccortex63 build -s {1} -k 33 -Q 15 -p -n 500M -m 16GB  -2 {2}:{3} " + outdir + "/raw/{1}.ctx"],
                    stdout=subprocess.PIPE, shell=True)
    subprocess.call([
        "cat " + outdir + "/to_be_processed.txt | parallel -j " + num_procs + " --colsep '\t' mccortex63 clean -m 8G -o " +
        outdir + "/clean/{1}.ctx " + outdir + "/raw/{1}.ctx"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["find " + outdir +"/clean -size 0 -delete"], stdout=subprocess.PIPE, shell=True)
    clean_files = [f.rstrip('.ctx') for f in os.listdir('./'+ outdir +'/clean')]
    not_cleaned = set(samples) - set(clean_files)
    uncleaned = open(outdir + '/low_coverage_samples.txt', 'w')
    [uncleaned.write(uc + '\n') for uc in not_cleaned]
    uncleaned.close()
    cleaned = open(outdir + '/cleaned.txt', 'w')
    [cleaned.write(cl + '\n') for cl in clean_files]
    cleaned.close()
    subprocess.call([
        "cat " + outdir + "/cleaned.txt | parallel -j " + num_procs + " mcBubble.py {}.ctx " + ref_ctx + " " + reference + " " + outdir]
        , stdout=subprocess.PIPE, shell=True)
    subprocess.call(["for f in " + outdir + "/*.final.vcf; do bgzip $f; done"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["for f in " + outdir + "/*.final.vcf.gz; do bcftools index $f; done"], stdout=subprocess.PIPE,
                    shell=True)
    subprocess.call(["bcftools merge " + outdir + "/*.final.vcf.gz > " + outdir + "/final.merged.vcf"],
                    stdout=subprocess.PIPE, shell=True)
    subprocess.call(["merged_stripper.py " + outdir + "/final.merged.vcf > " + outdir + "/final.stripped.vcf"],
                    stdout=subprocess.PIPE, shell=True)
    subprocess.call([
        "cat " + sample_list + "| parallel mc_genotype.py " + reference + " " + outdir + "/final.stripped.vcf {} " + outdir]
        , stdout=subprocess.PIPE, shell=True)
    fasta_out = open(outdir + '/MC_consensus.fasta', 'w')
    for sample in samples:
        fasta_out.write('>' + sample + '\n')
        fasta_out.write(create_consensus(outdir + "/" + sample + '.cov.vcf') + '\n')


if __name__ == '__main__':
    main()

