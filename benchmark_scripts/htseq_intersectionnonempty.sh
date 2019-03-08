#!/bin/bash
# htseq_intersectionnonempty.sh - Performs both an alignment step and calls HTSeq using the resulting BAM output.  This call of HTSeq uses the "intersection-nonempty" algorithm
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

bash ${script_dir}/align_and_sort.sh "$org_inputs" "$working_dir"

#Uses python3
python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$feat_type" -i "$attr_id" -s "$htseq_stranded" "$bam" "$gff3" > "${working_dir}/htseq_intersection_nonempty.counts"

exit 0
