# mcOutbryk
SNP calling pipeline using mccortex

Inspired by/follow up of  https://github.com/iqbal-lab/outbryk.git

### requirements:

vcftools

vcflib

bwa

python3

mccortex

### Installation:

Clone this repository:
```
git clone https://github.com/hcdenbakker/mcOutbryk
```
Add the mcOutbryk folder to your path

Clone mccortex:
```
git clone --recursive https://github.com/mcveanlab/mccortex.git
```
Compile the mccortex63 binary:
```
cd mccortex

make MAXK=63 all
make MAXK=31 all
```

Add the mccortex63 binary and the mccortex/scripts folder to your path 

### Execution
Execute the pipeline in the folder with illumina fastq.gz files of the isolates you are interested in.
Create a 1 column list with the sample names and  decide on an appropriate reference sequence.

For help run:
`MCOutbryk. py -h`

```
MCOutbryk - de novo SNP and ambiguous site Caller

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE [REFERENCE ...]
                        Reference genome in fasta format, include path to file
                        if file is not in same directory
  --sample_list SAMPLE_LIST [SAMPLE_LIST ...]
                        Text file, one isolate/sample per line
  --results_dir RESULTS_DIR [RESULTS_DIR ...]
                        Directory which will contain results after the
                        analysis
  --processes PROCESSES [PROCESSES ...]
                        number of processes to run in parallel
  --delete_highly_divergent {yes,no}
                        delete highly divergent isolates automatically
                        (default: yes)
  --div_cutoff [DIV_CUTOFF]
                        number of SNPs used as cutoff to define highly
                        divergent isolates (default: 5000)
  --K [K]               K-mer size (default: 33)
  --SNP_vcf SNP_VCF [SNP_VCF ...]
                        VCF with sites to be called; if this is given, the
                        steps to create a SNP sites vcf de novo from the
                        isolates in the sample list will be skipped.
```

To run mcOutbryk:
```
MCOutbryk.py --reference reference.fasta --sample_list your_sample_list.txt --results_dir directory_for_results --processes parallel_processes
```
Warning: while this pipeline is generally heavy on the RAM usage, buidling the uncleaned graphs takes about 10 Gb of memory per process,
it is recommended to perform 1 process per 10 Gb of memory (I run 10 parallel processes on a 128 Gb machine) 

The pipeline does the following:

1. Checks if raw data for samples in the list are available in the execution directory

2. Create raw and cleaned graphs

3. Bubble call each individual isolate vs reference (and delete cleaned graphs and some intermediary files)

4. Create vcfs with trusted SNPs and create a 'pseudo-vcf' of all trusted SNPs

5. Use 'pseudo-vcf' to genotype (= call SNPs) of all isolates using their raw graph
  
6. Create SNP sites multi-fasta file (stored in resultsdir/MC_consensus.fasta; ambiguous calls: N, gaps: -, impossible to call: ?)); this file can be used as direct input for FastTree

