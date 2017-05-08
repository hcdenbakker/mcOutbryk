#!/usr/bin/env python3

# todo: write samples with empty cov.vcf files to log
# todo: start from different points, inputs (fastq, fastq.gz, bam, (raw) ctx)
# todo: Create a quicker, but dirtier way to create final SNP-matrix: merge initial VCFs, and change to multifasta

import os
import re
import subprocess
from argparse import (ArgumentParser)
from time import localtime, strftime


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
    parser.add_argument('--delete_highly_divergent',
                        choices=['yes', 'no'], default= 'yes', help='delete highly divergent isolates automatically (default: yes)')
    parser.add_argument('--div_cutoff', nargs='?', const=5000, type=int, default = 5000,
                        help='number of SNPs used as cutoff to define highly divergent isolates (default: 5000)')

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
        ["mccortex63 build -q -s " + ref + " -k 33 -m 8GB -f -t 32 -1 " + reference + " " + outdir + "/ref/" + ref
         + ".ctx"], stdout=subprocess.PIPE, shell=True)
    subprocess.call(["cp " + reference + " " + outdir + "/ref; bwa index " + outdir + "/ref/" + ref_bwa],
                    stdout=subprocess.PIPE, shell=True)


def get_list(infile):
    samples = []
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        samples.append(line.strip('\n'))
    return samples


def summary(vcf):
    name = os.path.basename(vcf)
    count = 0
    with open(vcf, 'r') as gz:
        for line in gz:
            line = line.rstrip('\n')
            if line[0] == '#': continue
            else: count += 1
    return name, count


def create_master_vcf(list, outdir):
    #get minimum information to create new vcf from list[0]
    with open(outdir +'/'+ list[0], 'r') as m:
        for line in m:
            if line.startswith('##contig='):
                contig = line
                chrom = line.split(',')[0].split('=')[-1]
                print(chrom)

    #create dict to keep track of position, ref, alt
    master_dict = {}
    for vcf in list:
        with open(outdir +'/'+ vcf, 'r') as v:
            for line in v:
                if line[0] == '#':
                    continue
                else:
                    variant = line.split('\t')
                    if len(variant[4]) > 1:
                        continue
                    else:
                        if int(variant[1]) in master_dict:
                            if variant[4] in master_dict[int(variant[1])][1]:
                                continue
                            else:
                                master_dict[int(variant[1])][1].append(variant[4])
                        else:
                            master_dict[int(variant[1])] = (variant[3], [variant[4]])
    #Make the list and build the vcf
    sites = [a for a in master_dict]
    with open(outdir +'/simple_merge.vcf', 'w') as merged_vcf:
        merged_vcf.write('##fileformat=VCFv4.2\n')
        merged_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' +
                         '##FORMAT=<ID=K33R,Number=A,Type=Integer,Description="Coverage on ref (k=33): sum(kmer_covs) / exp_num_kmers">\n' +
                         '##FORMAT=<ID=K33A,Number=A,Type=Integer,Description="Coverage on alt (k=33): sum(kmer_covs) / exp_num_kmers">\n')
        merged_vcf.write(contig)
        merged_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDummy\n')
        for s in sorted(sites):
            #CP002002.1	1415	.	C	T	.	PASS	BUBBLE=357;K33	K33R:K33A	.:.
            merged_vcf.write(chrom + '\t' + str(s) + '\t' + '.\t' + str(master_dict[s][0]) + '\t' +
                             ','.join(master_dict[s][1]) + '\t' + '.' + '\t' + 'PASS' +
                             '\t' + '.' + '\t' + 'K33R:K33A' + '\t' + '.:.\n')
            #print(s, master_dict[s][0], master_dict[s][1])






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

def create_consensus_from_geno(infile):
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
                    v = variant[-1].split(':')
                    if v[0] == './.':
                        consensus += '?'
                    elif variant[1] == '0' == variant[2]:
                        consensus += '-'
                    elif v[0] == '0/0' and v[-1] == '0':
                        consensus += '?'
                    elif v[0] == '0/0' and int(v[-1]) > 0:
                        consensus += variant[3]
                    elif v[0] == '1/1' and v[-1] == '0':
                        consensus += '?'
                    elif v[0] == '1/1' and int(v[-1]) > 0:
                        consensus += variant[4]
                    elif v[0] == '0/1' and v[-1] == '0':
                        consensus += '?'
                    elif v[0] == '0/1' and int(v[-1]) > 0:
                        consensus += 'N'

    return consensus


def main():
    args = parse_args()
    reference = args.reference[0]
    samples = get_list(args.sample_list[0])
    sample_list = args.sample_list[0]
    outdir = args.results_dir[0]
    num_procs = args.processes[0]
    high_filter = args.delete_highly_divergent
    cutoff = args.div_cutoff
    print(cutoff)
    subprocess.call(['mkdir ' + outdir], stdout=subprocess.PIPE, shell=True)
    log = open(outdir + '/mcOutbryk.log', 'w')
    log.write(strftime("%Y-%m-%d %H:%M:%S", localtime()) + ': Start\n')
    log.close()
    ref = os.path.basename(reference).rstrip('.fasta')
    ref_ctx = str(ref + ".ctx")
    create_ref_binary(outdir, reference)
    perform_actions_list(samples, outdir)
    subprocess.call(["cat " + outdir + "/to_be_processed.txt | parallel --gnu -j " +
                     num_procs + " --colsep '\t' mccortex63 build -q -s {1} -k 33 -Q 15 -p -n 500M -m 16GB  -2 {2}:{3} " + outdir + "/raw/{1}.ctx"],
                    stdout=subprocess.PIPE, shell=True)
    log = open(outdir + '/mcOutbryk.log', 'a')
    log.write(strftime("%Y-%m-%d %H:%M:%S", localtime()) + ': Raw graphs created \n')
    log.close()
    subprocess.call([
        "cat " + outdir + "/to_be_processed.txt | parallel --gnu -j " + num_procs + " --colsep '\t' mccortex63 clean -q -m 8G -o " +
        outdir + "/clean/{1}.ctx " + outdir + "/raw/{1}.ctx"], stdout=subprocess.PIPE, shell=True)
    log = open(outdir + '/mcOutbryk.log', 'a')
    log.write(strftime("%Y-%m-%d %H:%M:%S", localtime()) + ': Samples cleaned' '\n')
    subprocess.call(["find " + outdir +"/clean -size 0 -delete"], stdout=subprocess.PIPE, shell=True)
    clean_files = [f[:-4] for f in os.listdir('./'+ outdir +'/clean')]
    not_cleaned = set(samples) - set(clean_files)
    uncleaned = open(outdir + '/low_coverage_samples.txt', 'w')
    [uncleaned.write(uc + '\n') for uc in not_cleaned]
    uncleaned.close()
    log.write('Samples not cleaned:' + str(len(not_cleaned)) + '\n')
    if len(not_cleaned) > 0:
        [log.write(uc + '\n') for uc in not_cleaned]
    log.close()
    cleaned = open(outdir + '/cleaned.txt', 'w')
    [cleaned.write(cl + '\n') for cl in clean_files]
    cleaned.close()
    subprocess.call([
        "cat " + outdir + "/cleaned.txt | parallel --gnu -j " + num_procs + " mcBubble.py {}.ctx " + ref_ctx + " " + reference + " " + outdir]
        , stdout=subprocess.PIPE, shell=True)
    individual_vcfs = [f for f in os.listdir('./' + outdir ) if f.endswith('.final.vcf')]
    log = open(outdir + '/mcOutbryk.log', 'a')
    log.write('Called SNPs per sample:'+ '\n')
    SNV_dict ={}
    for vcf in individual_vcfs:
        vcf_sum = summary(outdir + '/' + vcf)
        log.write(str(vcf_sum[0]) +': ' + str(vcf_sum[1]) + '\n')
        SNV_dict[vcf_sum[0][:-14]] = int(vcf_sum[1])
    highly_divergent = [f for f in SNV_dict if SNV_dict[f] > cutoff]
    if len(highly_divergent) > 0:
        log.write('Highly divergent isolates detected!')
        print('highly divergent isolates detected!')
        for f in highly_divergent:
            print(f)
            log.write(f + '\n')
        if high_filter == 'yes':
            log.write('Highly divergent isolates will not be included in variant list!'+ '\n')
            print('Highly divergent isolates will not be included in variant list!')
            for f in highly_divergent:
                subprocess.call(["rm " + outdir + '/' + f + '.ctx.final.vcf'], stdout=subprocess.PIPE, shell=True)
            samples = list(set(samples) - set(highly_divergent))
            print(samples)
        else:
            log.write('You have chosen to include highly divergent isolates in construction variant list!\n')
            print('highly divergent isolates will be included in further analyses!')
    create_master_vcf([f for f in os.listdir('./' + outdir) if f.endswith('.final.vcf')], outdir)
    subprocess.call(['rm '+outdir + "/*.final.vcf"],stdout=subprocess.PIPE, shell=True)
    subprocess.call([
        "cat " + sample_list + "| parallel --gnu mc_genotype.py " + reference + " " + outdir + "/simple_merge.vcf {} " + outdir]
        , stdout=subprocess.PIPE, shell=True)
    fasta_out = open(outdir + '/MC_consensus.fasta', 'w')
    cov_vcfs = [f for f in os.listdir('./' + outdir) if f.endswith('.cov.vcf')]
    for vcf in cov_vcfs:
        if os.stat(outdir + "/" + vcf).st_size == 0:
            print(vcf + ' is an empty file!')
        else:
            if vcf[:-8] in highly_divergent:
                fasta_out.write('>' + vcf[:-8]  +'_highly_divergent' + '\n')
                fasta_out.write(create_consensus(outdir + "/" + vcf) + '\n')
            else:
                fasta_out.write('>' + vcf[:-8]  + '\n')
                fasta_out.write(create_consensus(outdir + "/" + vcf) + '\n')
    subprocess.call(['rm '+outdir + "/*.final.vcf.gz; rm " + outdir + "/*.final.vcf.gz.csi"]
                    , stdout=subprocess.PIPE, shell=True)
    subprocess.call([
        "cat " + sample_list + "| parallel --gnu mc_real_genotype.py " + reference + " " + outdir + "/simple_merge.vcf {} " + outdir]
        , stdout=subprocess.PIPE, shell=True)
    fasta_out_geno = open(outdir + '/MC_consensus_geno.fasta', 'w')
    geno_vcfs = [f for f in os.listdir('./' + outdir) if f.endswith('.geno.vcf')]
    for vcf in geno_vcfs:
        if os.stat(outdir + "/" + vcf).st_size == 0:
            print(vcf + ' is an empty file!')
        else:
            if vcf[:-9] in highly_divergent:
                fasta_out_geno.write('>' + vcf[:-9] + '_highly_divergent' + '\n')
                fasta_out_geno.write(create_consensus_from_geno(outdir + "/" + vcf) + '\n')
            else:
                fasta_out_geno.write('>' + vcf[:-9] + '\n')
                fasta_out_geno.write(create_consensus_from_geno(outdir + "/" + vcf) + '\n')
    subprocess.call(["rm " + outdir + "/*.cov1.vcf"], stdout=subprocess.PIPE, shell=True)

    log.write(strftime("%Y-%m-%d %H:%M:%S", localtime()) + ': End ')


if __name__ == '__main__':
    main()
