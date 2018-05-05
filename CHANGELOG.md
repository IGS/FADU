NOTE: I haven't really kept track of changes onto now, so this will generally start post-1.0 release

#v1.1
* Will now throw out optical or PCR duplicate reads.  `samtools depth` excludes them anyways, and keeping them in the analysis can cause errors later on
* Will now throw out multi-mapped reads.  Keeping these reads can overestimate the readcounts of some genes, as the read will now count for two alignments.
* Fixed small symlinking bug where the script will error if the symlinked source and destination are the same
* Will create temporary directory if it currently does not exist
* Improved error handling in some areas of the script (raising errors instead of using assertions, for example)
* Properly paired BAM files will also be indexed
* Fixed issue where the range of fragment depth was not including the ending coordinate
* Ran the code through a couple of Python linters to clean it up