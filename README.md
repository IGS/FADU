# FADU
Feature Aggregate Depth Utility

## Github Pages site
https://igs.github.io/FADU/

# Description
Generate counts of reads that map to non-overlapping portions of genes

# Authors
* Shaun Adkins (sadkins@som.umaryland.edu)
* Matthew Chung (mattchung@umaryland.edu)

# Requirements
Julia - v0.7 or later (v1.0 or later preferred)
Modules
* ArgParse
* BioAlignments - v1.0.0
* GenomicFeatures - v1.0.0

# Input
* A single BAM file of reads (sorted by position)
* A GFF3-formatted annotation file

# Output
The output is a tab-delimited file whose name is the BAM input file's basename, ending in '.counts.txt'.  For example if the name of the BAM file was 'sample123.bam', the output file name is 'sample123.counts.txt'

The output file has the following fields:
* Name of feature
* Number of BAM fragments that aligned to this feature's non-overlapping bases
* Fractionalized fragment counts 
  * For each aligned fragment, the depth is (intersect of fragment coords to feature coords that do not overlap with other features) / length of the fragment

```
usage: fadu.jl -b /path/to/file.bam -g /path/to/annotation.gff3
               -o /path/to/output/dir [-s STRANDED] [-f FEATURE_TYPE]
               [-a ATTRIBUTE_TYPE]

optional arguments:
  -b, --bam_file /path/to/file.bam
                        Path to BAM file (SAM is not supported).
  -g, --gff3_file /path/to/annotation.gff3
                        Path to GFF3-formatted annotation file (GTF is not supported).
  -o, --output_dir /path/to/output/dir/
                        Directory to write the output.
  -s, --stranded STRANDED
                        Indicate if BAM reads are from a
                        strand-specific assay. Choose between
                        'yes', 'no', or 'reverse'. (default: "no")
  -f, --feature_type FEATURE_TYPE
                        Which GFF3/GTF feature type (column 3) to
                        obtain. readcount statistics for.
                        Case-sensitive. (default: "gene")
  -a, --attribute_type ATTRIBUTE_TYPE
                        Which GFF3/GTF feature type (column 9) to
                        obtain. readcount statistics for.
                        Case-sensitive. (default: "ID")
  --version             show version information and exit
  -h, --help            show this help message and exit
```

# Special Notes
## Notes about kept reads for calculations
When processing the input BAM file, any reads that are not properly paired or are not primary mappings will be thrown out.  This is determined by the SAM-format bit-flag.

# Frequently Asked Questions
## Is FADU designed to work with both prokaryotic and eukaryotic samples?
FADU currently does not support eukaryotic samples, or at least those that come from Tophat or HISAT2 aligners.

## When running FADU.jl, I get a "too many arguments" error.  How do I fix this?
A safe way to pass in values to the options is to enclose any string arguments (like file paths or strandedness) in quotations so that it ensures Julia's argument parser reads the parameter as an option value instead of a command-type argument.

## Can I use FADU to calculate counts for features using read depth instead of fragment depth?
Currently at this time, only fragment depth is supported.  Performing calculations using read depth may be implemented in the future though.

# Issues
If you have a question or suggestion about FADU, feel free to e-mail either of the authors above, or create an issue [here](https://github.com/IGS/FADU/issues)

