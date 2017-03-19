# mcOutbryk
SNP calling pipeline using mccortex

Inspired by/follow up of  https://github.com/iqbal-lab/outbryk.git

### requirements:

bcftools

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
```

Add the mccortex63 binary and the mccortex/scripts folder to your path 

### Execution
Execute the pipeline in the folder with illumina fastq.gz files of the isolates you are interested in.
Create a 1 column list with the sample names and  decide on an appropriate reference sequence.

To run mcOutbryk:
```
MCOutbryk.py --reference reference.fasta --sample_list your_sample_list.txt --results_dir directory_for_results --processes parallel_processes
```
Warning: while this pipeline is generally heavy on the RAM ussage, buidling the uncleaned graphs takes about 10 Gb of memory per process,
it is recommended to perform 1 process per 10 Gb of memory (I run 10 parallel processes on a 128 Gb machine) 

The pipeline does the following:

1. Checks if raw data for samples in the list are available in the execution directory

2. Create raw and cleaned graphs

3. Bubble call each individual isolate vs reference (and delete cleaned graphs and some intermediary files)

4. Create vcfs with trusted SNPs and create a vcf of all trusted SNPs with bcftools merge

5. Use combined lists to genotype all isolates with raw graphsresults stored in .cov.vcf files
  
6. Create SNP sites multi-fasta file (stored in resultsdir/MC_consensus.fasta; ambiguous calls: N, gaps: -, impossible to call: ?)); this file can be used as direct input for FastTree

