#!/bin/bash

# This script houses various paths to files related to this organism.  When supplied as the first argument upon calling the feature quantification scripts, these variables will be sourced into the calling script and used by it

out_dir=/###OUTDIR###/benchmarks/wbm

hisat_fna=${out_dir}/brugia_malayi.PRJNA10729.WBPS9.genomic.fa
nuc_genome_fna=${out_dir}/GCF_000008385.1_ASM838v1_genomic.fna
nuc_cds_fna=${out_dir}/GCF_000008385.1_ASM838v1_genomic.fna.cds.combined.fna
express_fna=$nuc_cds_fna
gff3=/###DATADIR###/wBm/GCF_000008385.1_ASM838v1_genomic.gff
fastq1=/###DATADIR###/fadu/wbm/SRR5192555_1.fastq
fastq2=/###DATADIR###/fadu/wbm/SRR5192555_2.fastq
bam=/###DATADIR###/fadu/wbm/SRR5192555.sortedbyposition.bam
cds_bam=/###DATADIR###/fadu/wbm/SRR5192555.cds.sortedbyposition.bam
SRR=SRR5192555
maxins=1003

threads=1
feat_type=CDS
attr_id=ID

express_stranded=--rf-stranded
fadu_stranded=reverse
featurecounts_stranded=2
htseq_stranded=reverse
kallisto_stranded=--rf-stranded
salmon_stranded=ISR
