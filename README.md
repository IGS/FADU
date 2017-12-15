# FADU
Feature Aggregate Depth Utility

## Authors:
* Shaun Adkins (sadkins@som.umaryland.edu)
* Matthew Chung (mattchung@umaryland.edu)

## Description
Generate counts of reads that map to non-overlapping portions of genes

## Requirements
Python 3, Pysam version 0.12.0.1

## Input
* One of the following:
** A single BAM file of reads
** A file listing the paths of BAM files
* One of the following:
** A GFF3-formatted annotation file
** A GTF-formatted annotation file

## Output
* Tab-delimited file containing gene statistics ending in 'uniquebp.stats.tsv'
** These are the following fields:
*** Name of the contig
*** Strand
*** Name of the feature
*** length of the feature
*** Number of bases in the feature the do not overlap with another feature of the same feature type
*** Percent of non-overlapping bases in feature compared to feature length
* Separate tab-delimited file for each BAM input ending in '.counts'
** THese are the following fields:
*** Name of feature
*** Readcounts (total depth of non-overlapping bases in feature / length of reads)

```
usage: fadu.py [-h]
               (--bam_file /path/to/alignment.bam | --bam_list /path/to/bam.list)
			   (--gff3 /path/to/annotation.gff3 | --gtf /path/to/annotation.gtf)
               --output_dir /path/to/output/dir [--tmp_dir /path/to/tmp/dir]
               [--stranded {yes,no,reverse}] [--feature_type FEATURE_TYPE]
               [--attribute_type ATTRIBUTE_TYPE] [--num_cores NUM_CORES]
               [--debug DEBUG/INFO/WARNING/ERROR/CRITICAL]


optional arguments:
  -h, --help            show this help message and exit
  --bam_file /path/to/alignment.bam, -b /path/to/alignment.bam
                        Path to BAM file. Choose between --bam_file or
                        --bam_list.
  --bam_list /path/to/bam.list, -B /path/to/bam.list
                        List file of BAM-formatted files (no SAM at this
                        time). Choose between --bam_file or --bam_list.
  --gff3 /path/to/annotation.gff3, -g /path/to/annotation.gff3
                        Path to GFF3-formatted annotation file. Choose between
                        --gff3 or --gtf.
  --gtf /path/to/annotation.gtf, -G /path/to/annotation.gtf
                        Path to GTF-formatted annotation file. Choose between
                        --gff3 or --gtf.
  --output_dir /path/to/output/dir, -o /path/to/output/dir
                        Required. Directory to store output.
  --tmp_dir /path/to/tmp/dir, -t /path/to/tmp/dir
                        Directory to store temporary output. If not provided,
                        the output directory will serve as the temporary
                        directory also
  --stranded {yes,no,reverse}, -s {yes,no,reverse}
                        Indicate if BAM reads are from a strand-specific assay
  --feature_type FEATURE_TYPE, -f FEATURE_TYPE
                        Which GFF3/GTF feature type (column 3) to obtain
                        readcount statistics for. Default is 'gene'. Case-
                        sensitive.
  --attribute_type ATTRIBUTE_TYPE, -a ATTRIBUTE_TYPE
                        Which GFF3/GTF attribute type (column 9) to obtain
                        readcount statistics for. Default is 'ID'. Case-
                        sensitive.
  --num_cores NUM_CORES, -n NUM_CORES
                        Number of cores to spread processes to when processing
                        BAM list.
  --debug DEBUG/INFO/WARNING/ERROR/CRITICAL, -d DEBUG/INFO/WARNING/ERROR/CRITICAL
                        Set the debug level
```
