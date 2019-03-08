#!/bin/bash

hisat2_index.sh - Used to build HiSAT2 indexes with gene references

# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs

mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

/usr/local/packages/hisat2-2.1.0/hisat2-build $hisat_fna $hisat_fna
