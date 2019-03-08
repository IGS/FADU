#!/bin/bash

# express_align_and_sort.sh - This script aligns the FASTQ files with HiSAT2, then sorts the resulting BAM with Samtools.  This version of the script is used to obtain the correct inputs for eXpress
script_dir=`dirname $0`
org_inputs=$1
source $script_dir/$org_inputs

working_dir=$2

/usr/local/packages/hisat2-2.1.0/hisat2 -p 1 -x "$express_fna" -1 "$fastq1" -2 "$fastq2" --no-spliced-alignment -X "$maxins" | /usr/local/packages/samtools-1.3.1/bin/samtools view -bhSo "${working_dir}/${SRR}.cds.bam" -

/usr/local/packages/samtools-1.9/bin/samtools sort "${working_dir}/${SRR}.cds.bam" -o "${working_dir}/${SRR}.cds.sortedbyposition.bam"

exit 0
