#call consensus alleles given a list of sites/reference genome

from sys import argv
from argparse import (ArgumentParser)
import subprocess
import os

'''
assemblies = get_list(argv[1])
    reference = argv[2]
    variant_vcf = argv[3]
'''
def parse_args():
    "Parse the input arguments, use '-h' for help."
    parser = ArgumentParser(description='mcAssembly - calling known variable sites in assemblies')
    # inputs
    parser.add_argument(
        '--reference', nargs='+', type=str, required=True,
        help='Reference genome in fasta format, include path to file if file is not in same directory ')
    parser.add_argument(
        '--assembly_list', nargs='+', type=str, required=True,
        help='Text file, one assembly name per line ')
    parser.add_argument(
        '--SNP_vcf', nargs='+', type=str, required=True,
        help='VCF with sites to be called')

    return parser.parse_args()


def get_list(file):
    assembly_list =[]
    with open(file) as f:
        for line in f:
            assembly_list.append(line.strip('\n'))
    return assembly_list

def create_binaries(assemblies):
    for file in assemblies:
        base_name = '.'.join(file.split('.')[:-1])
        subprocess.call(
            ["mccortex63 build -q -s " + base_name + " -k 33 -m 8GB -f -t 32 -1 " + file + " assemblies/" + base_name
             + ".ctx"], stdout=subprocess.PIPE, shell=True)


def create_consensus(infile):
    consensus = ''
    with open(infile, 'r') as gz:
        for line in gz:
            line = line.rstrip('\n')
            if line[0] == '#':
                continue
            else:
                variant = line.split('\t')
                if len(variant[4]) > 1:
                    continue
                else:
                    if variant[-1] == '.:.':
                        consensus += '?'
                    elif variant[-1] == '0:0':
                        consensus += '-'
                    elif variant[-1] == '1:0':
                        consensus += variant[3]
                    elif variant[-1] == '0:1':
                        consensus += variant[4]
                    else:
                        consensus += 'N'

    return consensus


def call_consensus(reference, variant_vcf):
    #create list of binaries
    binaries = [f for f in os.listdir('./assemblies') if f.endswith('.ctx')]
    #create 'coverage' vcfs for binaries
    for biny in binaries:
        base_name = '.'.join(biny.split('.')[:-1])
        subprocess.call(
            [
                "mccortex63 vcfcov -q -m 1G -r " + reference + " -f -o " + base_name + ".cov.vcf " +
                variant_vcf + " assemblies/"+ biny], stdout=subprocess.PIPE, shell=True)
    #call consensus based on 'coverage' vcfs
    multi_fasta = open('MC_assemblies.fasta', 'w')
    for biny in binaries:
        multi_fasta.write('>' + '.'.join(biny.split('.')[:-1]) + '\n')
        multi_fasta.write(create_consensus('.'.join(biny.split('.')[:-1]) + ".cov.vcf") + '\n')




def main():
    args = parse_args()
    assemblies = get_list(args.assembly_list[0])
    reference = args.reference[0]
    variant_vcf = args.SNP_vcf[0]
    create_binaries(assemblies)
    call_consensus(reference, variant_vcf)

if __name__ == '__main__':
    main()