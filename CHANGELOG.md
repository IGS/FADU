# CHANGELOG

## v1.9.1

* BUGFIX - Updating BGZFStreams and XAM packages to address a BoundsError issue (https://github.com/IGS/FADU/issues/14)
* Updating various packages to reflect versioning I cited in (https://github.com/BioJulia/BGZFStreams.jl/issues/30)

## v1.9

* FEATURE - Adding --no_output_header option to not print the header line in output. Useful for when you want to pipe this output to another tool.
* FEATURE - Adding --exclude_regions option to exclude any regions found in a passed-in BED file.
* FEATURE - Adding --log_file and --log_level options to display or hide log messages of varying severity. Note that currently --log_level=DEBUG will print the line number, as I haven't figured out how to disable that for debug-level messages yet.
* Optimized some pieces of code for speed (thanks Github Copilot). This should increase the speed where a GFF feature has 100,000s of overlaps with the BAM alignments

## v1.8.3

* BUGFIX - Changed thrown FileNotFoundException to SystemError as the former does not exist in Julia
* BUGFIX - Correcting refindex check (see commit https://github.com/BioJulia/XAM.jl/pull/12/commits/8b7a2ecc2b80ca7717b65d1489c99f50a0bf0d10). This addresses issue #4.
* Added descriptors to the project.toml file

## v1.8.2

* BUGFIX - converting vectorized operation to array comprehension to fix an error

## v1.8.1

* BUGFIX - Fixed edge case where the alignment overlap iterator would attempt to read a stream block past the end of the BAM.Reader file

## v1.8

* Update packages to be compatible with Julia v1.4.2

## v1.7

* Changed the Set of nonoverlapping coords to a BitSet, which significantly sped up execution
* Rewrote some of the EM algorithm to use a dictionary of StructArrays, which speeds up execution
* Other various performance tweaks

## v1.6.1

* BUGFIX - Fixing issue with generating uniq coordinates dictionary using unstranded data if the "+" strand key for a sequence ID is not already present.
* Added a try/catch block for improperly-formatted GFF3 files.  If the file is improper formatting (such as non-URL-escaped '=' in a attribute tag key), the GFF3 Reader from BioJulia will throw an error, which is something out of the control of FADU.jl

## v1.6

* FEATURE - Added "--em\_iterations" option to use the Expectation-Maximization algorithm to re-add and allocate multimapped alignment records based on the proportion of the feature counts derived from singly-mapped reads.
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
