---
layout: default
---

# Description
Generate counts of reads that map to non-overlapping portions of genes

# Authors
* Shaun Adkins (sadkins@som.umaryland.edu)
* Matthew Chung (mattchung@umaryland.edu)

# Requirements
Python 3, Pysam version 0.12.0.1

# Input
* One of the following:
  * A single BAM file of reads (sorted by position)
  * A file listing the paths of BAM files (each file sorted by position)
* One of the following:
  * A GFF3-formatted annotation file
  * A GTF-formatted annotation file

# Output
* Tab-delimited file containing gene statistics ending in 'uniquebp.stats.tsv'
  * These are the following fields:
    * Name of the contig
    * Strand
    * Name of the feature
    * length of the feature
    * Number of bases in the feature the do not overlap with another feature of the same feature type
    * Percent of non-overlapping bases in feature compared to feature length
* Separate tab-delimited file for each BAM input ending in '.counts'
  * These are the following fields:
    * Name of feature
    * Readcounts (total depth of non-overlapping bases in feature / length of reads)

```
usage: fadu.py [-h]
               (--bam_file /path/to/alignment.bam | --bam_list /path/to/bam.list)
               (--gff3 /path/to/annotation.gff3 | --gtf /path/to/annotation.gtf)
               --output_dir /path/to/output/dir [--tmp_dir /path/to/tmp/dir]
               [--stranded {yes,no,reverse}] [--feature_type FEATURE_TYPE]
               [--attribute_type ATTRIBUTE_TYPE] [--count_by {read,fragment}]
               [--keep_only_properly_paired] [--rm_multimapped_reads]
               [--num_cores NUM_CORES]
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
  --count_by {read,fragment}, -c {read,fragment}
                        How to count the reads when performing depth
                        calculations. Default is 'read'.
  --keep_only_properly_paired
                        Enable flag to remove any reads that are not properly
                        paired from the depth count statistics.
  --rm_multimapped_reads
                        Enable flag to remove any multimapped reads ('NH' tag > 1)
  --num_cores NUM_CORES, -n NUM_CORES
                        Number of cores to spread processes to when processing
                        BAM list.
  --debug DEBUG/INFO/WARNING/ERROR/CRITICAL, -d DEBUG/INFO/WARNING/ERROR/CRITICAL
                        Set the debug level
```

# Special Notes
## Notes about kept reads for depth and coverage calculations
### Default
By default, all reads that map to the reference are kept (SAM flag 0x4).  Unmapped reads would not factor into depth calculations, and so they are removed in advance so they do not factor into calculating the average read length
### With --keep\_only\_properly\_paired enabled
If this option is enabled, then only properly paired reads are kept (SAM flag 0x2)

## Note about optical or PCR duplicate reads
Because `samtools depth` excludes any reads with the 0x400 bit flag, FADU will also throw out these reads.  If you want duplicates to be included in the analysis, ensure no reads have this bit flag enabled.  This includes refraining from running tools to detect duplicates, such as PicardTools "MarkDuplicates" utility.

# Frequently Asked Questions
## Is FADU designed to work with both prokaryotic and eukaryotic samples?
FADU currently does not support eukaryotic samples, or at least those that come from Tophat or HISAT2 aligners.  

# Issues
If you have a question or suggestion about FADU, feel free to e-mail either of the authors above, or create an issue on this GitHub page.

