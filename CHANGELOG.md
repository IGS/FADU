# CHANGELOG

## v1.0
* Initial version of the Julia rewrite.  FADU was rewritten from Python to Julia in order to have it run faster and use less memory.
* Algorithm constructs IntervalTree data structures from both the annotation features and the BAM fragment records. After this, the algorithm determines the fractionalized fragment depth for each feature by iterating through the overlaps of both sets of interval trees.
