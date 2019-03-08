#!/bin/bash

# kallisto.sh - Script to execute 'kallisto index' and 'kallisto quant'
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

/usr/local/packages/kallisto-0.44.0/kallisto index -i "${nuc_cds_fna}.kallisto.index" --make-unique "$nuc_cds_fna"

/usr/local/packages/kallisto-0.44.0/kallisto quant $kallisto_stranded -t "$threads" -i "${nuc_cds_fna}.kallisto.index" -o "${working_dir}" "$fastq1" "$fastq2"

exit 0
