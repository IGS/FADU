#!/bin/bash

# featurecounts_overlap.sh - Performs both an alignment step and calls featureCounts using the resulting BAM output.  This call of featureCounts uses the default algorithm
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

bash ${script_dir}/align_and_sort.sh "$org_inputs" "$working_dir"

/usr/local/packages/subread-1.6.3-source/bin/featureCounts -p --primary -a "$gff3" -o "${working_dir}/featurecounts_default.counts" -t "$feat_type" -g "$attr_id" -s "$featurecounts_stranded" "$bam"

exit 0
