#!/bin/bash

# express.sh - Performs both an alignment step and calls eXpress using the resulting BAM output
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

bash ${script_dir}/express_align_and_sort.sh "$org_inputs" "$working_dir" 

/usr/local/packages/express-1.5.1/express -o ${working_dir} "$express_stranded" "$nuc_cds_fna" "$cds_bam" 

exit 0
