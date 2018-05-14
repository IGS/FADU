#!/bin/bash

gff3="./GCF_000008385.1_ASM838v1_genomic.gff"
tempdir="./FADU_tmp"
outputdir="./FADU"
bam="./test.bam"

python3 ../fadu.py --bam_file $bam --gff3 $gff3 --tmp_dir $tempdir --output_dir $outputdir --stranded yes --feature_type CDS --attribute_type ID --count_by fragment --keep_only_properly_paired
