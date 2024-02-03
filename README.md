# FADU

Feature Aggregate Depth Utility

## Github Pages site

https://igs.github.io/FADU/

## Description

Most current available quantification tools for transcriptomics analyses have been designed using human data sets, in which full-length transcript sequences are available. In most prokaryotic systems, full-length transcript sequences are unavailable, causing prokaryotic transcriptomics analyses being performed using coding sequences instead. In contrast to eukaryotes, prokaryotes contain polycistronic transcripts and when genes are quantified based on coding sequences instead of transcript sequences, this leads to an increased abundance of ambiguous multi-gene fragments, especially between genes in operons. Here we describe FADU, a quantification tool designed specifically for prokaryotic RNA-Seq analyses with the purpose of quantifying operonic genes while minimizing the pitfalls associated with improperly assigning fragment counts from ambiguous transcripts.

## Authors

* Shaun Adkins (sadkins@som.umaryland.edu)
* Matthew Chung (mattchung@umaryland.edu)

## Current FADU release

v1.8.3

## Requirements

Julia - v1.4.2 or later

### Packages

NOTE: Packages installation instructions can be found at https://docs.julialang.org/en/v1.4/stdlib/Pkg/index.html#Getting-Started-1

* ArgParse
* BGZFStreams.jl - v0.3.0
* GenomicFeatures.jl - v2.0.0
* GFF3 - v0.2.0
* Indexes.jl - v0.1.1
* StructArrays.jl - v0.4.4
* XAM.jl - v0.2.6
* BED.jl - v0.3.0 (if you wish to exclude regions)

### NOTE on versions

FADU v1.7 and earlier will not work with Julia v1.4.2 and has different package requirements.  If using one of those versions, please use the recommended version of Julia and recommended packages listed on the README of that particular release.

### OS Requirements

FADU is supported for both the Linux and Mac OSX operating systems, and has been tested on the following operating system versions:

* Linux - CentOS 7
* Mac OSX - High Sierra (v10.13.6)
* MacOS M1 architecture- Monterrey (v12.4 - used Julia v1.7.3 for this test)

## Input

* A single BAM file of reads (sorted by position)
* A GFF3-formatted annotation file

## Output

The output is a tab-delimited file whose name is the BAM input file's basename, ending in '.counts.txt'.  For example if the name of the BAM file was 'sample123.bam', the output file name is 'sample123.counts.txt'

The output file has the following fields:

* Name of feature
* Number of bases that do not overlap with other features (uniq\_len)
* Number of BAM records that aligned to this feature's non-overlapping bases
  * Aligned reads not part of a fragment get 0.5 count, and fragments get 1 count
* Fractionalized feature counts
  * For each aligned fragment, the depth is (intersect of fragment coords to feature coords that do not overlap with other features) / length of the fragment
  * Aligned reads that are not part of a fragment get their depth contribution cut in half
* TPM count in scientific notation

## Installation instructions

First `cd` to the directory of choice to install FADU

```bash
git clone https://github.com/IGS/FADU.git
cd FADU
```

## Example command

For this example we will create a directory called "fadu\_output" to write output into.  This example uses a BAM alignment file of subsampled Wolbachia reads.  This example should take about 60 seconds to run.

```bash
mkdir ./fadu_output
julia ./fadu.jl -M -g "./test_data/GCF_000008385.1_ASM838v1_genomic.gff" -b "./test_data/SRR5192555.100000x.sortedbyposition.bam" -o "./fadu_output" -s "reverse" -f "CDS" -a "ID"
```
The resulting output file should be called "SRR5192555.100000x.sortedbyposition.counts.txt" and located in the "fadu\_output" directory.

Here are the counts for the first 10 features:

```bash
featureID   uniq_len    num_alignments  counts  tpm
cds0    1017    4.5 2.81    1744.38
cds1    1194    3.5 2.48    1310.50
cds10   1161    0.0 0.00    0.00
cds100  591 0.0 0.00    0.00
cds1000 741 8.0 7.46    6358.81
cds1001 850 5.0 4.45    3310.38
cds1002 829 1.0 0.48    363.86
cds1003 1167    0.0 0.00    0.00
cds1004 1164    0.0 0.00    0.00
cds1005 816 0.0 0.00    0.00
```

## Usage

```bash
usage: fadu.jl -b /path/to/file.bam -g /path/to/annotation.gff3
               -o /path/to/output/dir/ [--no_output_header]
               [-x /path/to/regions.bed] [-s STRANDED]
               [-f FEATURE_TYPE] [-a ATTRIBUTE_TYPE] [-p]
               [-m MAX_FRAGMENT_SIZE] [-M] [-e EM_ITER] [--version]
               [-h]

Generate counts of reads that map to non-overlapping portions of genes

optional arguments:
  -b, --bam_file /path/to/file.bam
                        Path to BAM file (SAM is not supported).
  -g, --gff3_file /path/to/annotation.gff3
                        Path to GFF3-formatted annotation file (GTF is
                        not supported).
  -o, --output_dir /path/to/output/dir/
                        Directory to write the output.
  --no_output_header    If enabled, do not write the header line to
                        the output file.
  -x, --exclude_regions_file /path/to/regions.bed
                        Path to BED file containing regions to exclude
                        from analysis.  Any alignments that overlap
                        these regions will be ignored.
  -s, --stranded STRANDED
                        Indicate if BAM reads are from a
                        strand-specific assay. Choose between 'yes',
                        'no', or 'reverse'. (default: "no")
  -f, --feature_type FEATURE_TYPE
                        Which GFF3/GTF feature type (column 3) to
                        obtain. readcount statistics for.
                        Case-sensitive. (default: "gene")
  -a, --attribute_type ATTRIBUTE_TYPE
                        Which GFF3/GTF feature type (column 9) to
                        obtain. readcount statistics for.
                        Case-sensitive. (default: "ID")
  -p, --keep_only_proper_pairs
                        If enabled, keep only properly paired reads
                        when performing calculations.
  -m, --max_fragment_size MAX_FRAGMENT_SIZE
                        If the fragment size of properly-paired reads
                        exceeds this value, process pair as single
                        reads instead of as a fragment. Setting this
                        value to 0 will make every fragment pair be
                        processed as two individual reads. If
                        --keep_only_proper_pairs is enabled, then any
                        fragment exceeding this value will be
                        discarded. (type: Int64, default: 1000)
  -M, --remove_multimapped
                        If enabled, remove any reads or fragments that
                        are mapped to multiple regions of the genome,
                        indiated by the 'NH' attribute being greater
                        than 1.  Otherwise, use EM algorithm to
                        quantify reads after all other reads are
                        counted.
  -e, --em_iterations EM_ITER
                        Number of iterations to perform EM algorithm
                        on multimapped read depth. Only applies if
                        --remove_multimapped flag is not passed (is
                        disabled) (type: Int64, default: 1)
  --version             show version information and exit
  -h, --help            show this help message and exit
```

## Frequently Asked Questions

### Is FADU designed to work with both prokaryotic and eukaryotic samples?

FADU currently does not support eukaryotic samples

### When running FADU.jl, I get a "too many arguments" error.  How do I fix this?

A safe way to pass in values to the options is to enclose any string arguments (like file paths or strandedness) in quotations so that it ensures Julia's argument parser reads the parameter as an option value instead of a command-type argument.

### Can I use FADU to calculate counts for features using read depth instead of fragment depth?

Currently at this time, only fragment depth is supported.  Performing calculations using read depth may be implemented in the future though.

### How can I use FADU output as input for DESeq2?

Perhaps the easiest possibility would be to extract the 'featureID' and 'counts' columns and structure them like the output from htseq-count.  From there you can follow the corresponding instructions from the DESeq2 manual on using htseq-count output as input.

### My dataset consists of only single-end reads, and the assigned counts only appear to be half of what was expected.  What gives?

In the case of singletons, FADU gives them half the weight compared to a paired-end read in the final counts.  One reason FADU does this is because we noticed that if we treated single-end and properly-paired-end reads with the same weight, there was some overcounting done with respect to certain genes or features.

### Is there a Docker variation of FADU that I can run instead of installing Julia?

Yes! Currently, there is no Docker image for FADU on Docker Hub as their free tier removes any images that have been inactive for 6+ months. Instead, you can build the Docker image as follows (this assumes you are in the FADU root directory):

```bash
docker build -t fadu .
```

This Docker image uses Julia 1.7 and installs the various packages required to run FADU. Once the image has been built you can create and run a FADU container by running the following:

```bash
docker run -v $PWD:/opt/FADU fadu <fadu_options>
```

This command assumes that your input files are in your current working directory and that your output will be written to a space relative to your current working directory as well.

## Issues

If you have a question or suggestion about FADU, feel free to e-mail either of the authors above, or create an issue [here](https://github.com/IGS/FADU/issues)

