# CHANGELOG

## v1.6

* FEATURE - Added i"--em\_iterations" option to use the Expectation-Maximization algorithm to re-add and allocate multimapped alignment records based on the proportion of the feature counts derived from singly-mapped reads.
* FEATURE -  I recently learned that the BioAlignments package takes advantage of multithreading with respect to reading the BAM file.  So therefore, FADU supports multithreading.  This can be achieved by running "export JULIA_NUM_THREADS=#" before running FADU where # is the number of threads you want to use.
* DELETION - Due to the rewriting, the "--chunk_size" option has been removed.
* Changed nested dictionary information for each feature (num_alignments, feat_counts, coords_set) into a mutable struct, which shaved off some slight runtime
* Rewriting how FADU processes overlaps.  Instead of reading by chunks of BAM alignments and processing all overlaps of those alignments to all features, FADU will process overlaps of all BAM alignment records for a feature-by-feature basis.  This results in faster performance speed.
* Created new "include" file alignment\_overlap.jl, which implements a custom version of the "eachoverlaps" function to handle both fragments and reads accordingly.
* Created new "include" file bam\_record.jl, which houses some functions that work on Record types.
* Renamed some variables and functions to adhere to Julia styling conventions
* Rewrote step to get nonoverlapping coordinates for each reference sequence in the GFF, which should speed up that step, particularly with larger annotation files

## v1.5

* FEATURE - If the --remove\_multimapped option is not enabled, multimapped reads/fragments will be added to feature quantification using the Expectation-Maximization algorithm
* DELETION - Removing --min\_mapping\_quality argument.  All multimapped reads will be 0 or 1, so this argument seemed unnecessary.
* Small optimization in creating the GenomicFeatures IntervalCollection, by parsing out just the feature types wanted from the get go instead of later on.
* Removed "uniq_coords" argument from the "process_overlaps" function since it was unused.
* Removed "feature" argument from the "compute_align_feat_set" function since it was unused.
* Renamed 'feat_depth' keys to 'feat_counts' and 'counter' keys to 'num_alignments' to better align to proper terminology

## v1.4

* BUGFIX - Corrected TPM counts
* BUGFIX - Corrected overcounting of the individual alignment-feature overlap counts
* FEATURE - Added --min\_mapping\_quality argument to only keep reads that surpass this min quality score
* More tweaks for speed
* Renamed a bunch of "fragment" variables to "alignment" since many of these variables also represented reads

## v1.3

* BUGFIX - Fixed bug with --max\_fragment\_size where reverse-stranded properly paired reads where still being used as a fragment since the template length reported is negative
* BUGFIX - Changed way fragment length is calculated.  Instead of using read1 of the proper pair, the read with the negative template length is used.  Realized that using template length to give fragment size can be incorrect in some cases, particularly when soft-clipping is involved.  So by using the rightmost read, the fragment can be calculated by taking the rightmost position of this read and the leftmost position (listed as 'nextposition' in the BAM.Record) of the leftmost read
* FEATURE - Adding --remove\_multimapped as an option.  Enabling it removes any reads mapping to multiple regions in the genome.
* FEATURE - Adding --chunk_size as an option.  Sets the number of validated reads/fragments that are read into memory before processing overlaps with features.
* Changed the custom BAM record IntervalCollection to simply read 'R' or 'F' as the metadata instead of "read" or "fragment".  This should hopefully be more efficient as Char types are used instead of String types.

## v1.2

* FEATURE - Added two new output columns.  1) Nonoverlapping feature length. 2) TPM (in scientific notation)
* FEATURE - Adding --max\_fragment\_size argument to handle extremely long (incorrect) fragments, and process them as reads instead.
* Specifying some type information for variables.
* Moving flag variables in assign_read_to_strand outside the function to become their own functions


## v1.1

* BUGFIX - Correctly verifies BAM filepath provided in arguments does exist.
* FEATURE - Added --keep\_only\_proper\_pairs option.  Before, only properly paired reads were kept.  Now any discordant reads and singletons will be kept unless this option is specified.  Discordant reads and singletons will carry half the weight of a fragment, since both reads are processed individually.

## v1.0

* Initial version of the Julia rewrite.  FADU was rewritten from Python to Julia in order to have it run faster and use less memory.
* Algorithm constructs IntervalTree data structures from both the annotation features and the BAM fragment records. After this, the algorithm determines the fractionalized fragment depth for each feature by iterating through the overlaps of both sets of interval trees.
