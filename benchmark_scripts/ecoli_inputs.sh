#!/bin/bash

# This script houses various paths to files related to this organism.  When supplied as the first argument upon calling the feature quantification scripts, these variables will be sourced into the calling script and used by it

out_dir=/###OUTDIR###/benchmarks/ecoli

nuc_genome_fna=${out_dir}/E2348_full.fa
express_fna=$nuc_genome_fna
hisat_fna=$nuc_genome_fna
nuc_cds_fna=${out_dir}/E2348_full.fa.cds.fna
gff3=/###DATADIR###/fadu/ecoli/E2348_full.gff3
fastq1=/###DATADIR###/fadu/ecoli/SRR2601722_1.fastq
fastq2=/###DATADIR###/fadu/ecoli/SRR2601722_2.fastq
bam=/###DATADIR###/fadu/ecoli/SRR2601722.sortedbyposition.bam
cds_bam=/###DATADIR###/fadu/ecoli/SRR2601722.cds.sortedbyname.bam
SRR=SRR2601722
maxins=865

threads=1
feat_type=CDS
attr_id=ID

express_stranded=''
fadu_stranded=no
featurecounts_stranded=0
htseq_stranded=no
kallisto_stranded=
salmon_stranded=IU
