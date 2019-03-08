#!/bin/bash

# This script houses various paths to files related to this organism.  When supplied as the first argument upon calling the feature quantification scripts, these variables will be sourced into the calling script and used by it

out_dir=/###OUTDIR###/benchmarks/echaffeensis

nuc_genome_fna=${out_dir}/Canis_familiaris.CanFam3.1.dna.toplevel.geark.combined.fna
express_fna=$nuc_genome_fna
hisat_fna=$nuc_genome_fna
nuc_cds_fna=${out_dir}/geark.fna.cds.combined.fna
gff3=/###DATADIR###/fadu/echaffeensis/geark.gff3
fastq1=/###DATADIR###/fadu/echaffeensis/SRR1188323_1.fastq
fastq2=/###DATADIR###/fadu/echaffeensis/SRR1188323_2.fastq
bam=/###DATADIR###/fadu/echaffeensis/SRR1188323.sortedbyposition.bam
cds_bam=/###DATADIR###/fadu/echaffeensis/SRR1188323.cds.sortedbyname.bam
SRR=SRR1188323
maxins=459

threads=1
feat_type=CDS
attr_id=ID

express_stranded=''
fadu_stranded=no
featurecounts_stranded=0
htseq_stranded=no
kallisto_stranded=
salmon_stranded=IU
