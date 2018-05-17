#!/bin/bash

gff3="./brugia_malayi.PRJNA10729.WBPS9.canonical_geneset.gtf"
#gff3="./GCF_000008385.1_ASM838v1_genomic.gff"
tempdir="./FADU_tmp"
outputdir="./FADU"
bam="./agss_BmV18hpi_a.sortedbyposition.bam"
#bam="./test.bam"

#python3 ../fadu.py --bam_file $bam --gff3 $gff3 --tmp_dir $tempdir --output_dir $outputdir --stranded yes --feature_type CDS --attribute_type ID --count_by fragment --keep_only_properly_paired --debug="DEBUG"
python3 ../fadu.py --bam_file $bam --gtf $gff3 --tmp_dir $tempdir --output_dir $outputdir --stranded reverse --feature_type exon --attribute_type gene_id --count_by fragment --rm_multimapped_reads --debug="DEBUG"
