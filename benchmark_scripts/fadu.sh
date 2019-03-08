#!/bin/bash

# fadu.sh - Performs both an alignment step and calls FADU using the resulting BAM output
# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

bash ${script_dir}/align_and_sort.sh "$org_inputs" "$working_dir"

# Uncomment this line if Julia cannot find the installed modules, assuming they are in the user's home directory
#export JULIA_DEPOT_PATH="~/.julia/"
/usr/local/bin/julia /home/sadkins/devel/git/FADU/fadu.jl --remove_multimapped -g "$gff3" -b "$bam" -o "${working_dir}" -s "$fadu_stranded" -f "$feat_type" -a "$attr_id"
exit 0
