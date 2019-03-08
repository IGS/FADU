#!/bin/bash

# salmon.sh - Script to execute 'salmon index' and 'salmon quant'
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

/usr/local/packages/salmon-0.7.2/bin/salmon index -t "$nuc_cds_fna" -i "${nuc_cds_fna}.salmon.index"

/usr/local/packages/salmon-0.7.2/bin/salmon quant -i "${nuc_cds_fna}.salmon.index"  -l "$salmon_stranded" -1 "$fastq1" -2 "$fastq2" -o "${working_dir}/salmon.counts" -p "$threads"

exit 0
