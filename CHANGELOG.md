# CHANGELOG

* More optimization of a couple functions
  * Fragment interval trees metadata will just be a character for read (R) or fragment (F) rather than the string of that word

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
