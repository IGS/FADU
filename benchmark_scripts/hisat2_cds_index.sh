#!/bin/bash

# hisat_cds_index.sh - Builds the Hisat2 index for CDS references

# $1 is one of the three organism input files

org_inputs=$1
source $org_inputs
mkdir -p $out_dir
script_dir=`dirname $0`
working_dir=$(mktemp -d -p "$out_dir")

/usr/local/packages/hisat2-2.1.0/hisat2-build $express_fna $express_fna
