# CHANGELOG

## v1.1
* BUGFIX - Correctly verifies BAM filepath provided in arguments does exist.
* FEATURE - Added --keep\_only\_proper\_pairs option.  Before, only properly paired reads were kept.  Now any discordant reads and singletons will be kept unless this option is specified.  Discordant reads and singletons will carry half the weight of a fragment, since both reads are processed individually.

## v1.0
* Initial version of the Julia rewrite.  FADU was rewritten from Python to Julia in order to have it run faster and use less memory.
* Algorithm constructs IntervalTree data structures from both the annotation features and the BAM fragment records. After this, the algorithm determines the fractionalized fragment depth for each feature by iterating through the overlaps of both sets of interval trees.
