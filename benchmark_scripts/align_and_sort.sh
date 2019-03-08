#!/bin/bash

# align_and_sort.sh - This script aligns the FASTQ files with HiSAT2, then sorts the resulting BAM with Samtools

script_dir=`dirname "$0"`
org_inputs=$1
source $script_dir/$org_inputs

working_dir=$2

/usr/local/packages/hisat2-2.1.0/hisat2 -p 1 -x "$hisat_fna" -1 "$fastq1" -2 "$fastq2" --no-spliced-alignment -X "$maxins" | /usr/local/packages/samtools-1.3.1/bin/samtools view -bhSo "${working_dir}/${SRR}.bam" -

/usr/local/packages/samtools-1.9/bin/samtools sort "${working_dir}/${SRR}.bam" -o "${working_dir}/${SRR}.sortedbyposition.bam"

exit 0
