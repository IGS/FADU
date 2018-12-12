# CHANGELOG

## v1.4
* BUGFIX - Corrected TPM counts
* FEATURE - Added --min\_mapping\_quality argument to only keep reads that surpass this min quality score
* More tweaks for speed

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
