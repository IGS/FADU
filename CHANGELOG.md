# CHANGELOG

## v1.2

* Rewrote a large portion of the BAM processing parts of the script
  * Far fewer intermediate files will be outputted.  Now, only the depth file, BAM index, and the symlinked BAM file will be in the temporary output directory.  The main output directory results are still the same.
  * FADU now reads through the BAM file only once.  In v1.1, when counting by fragments, FADU would read through a list of alignments twice, once before running "samtools depth" and again (properly paired reads only) for adjusting the depth.  This should result in faster execution times.
  * Changed the nesting order of keys and values in the "depth_dict" dictionary.  Made "strand" higher and "coord" nested within.
* Eliminated the use of "samtools depth" in the script, so coordinate depth is calculated and adjusted on-the-fly as alignment records are processed.
* If the user elects to count by fragments, then only the fragment depth will be returned to file, not read and fragment depth.

## v1.1

* Will now throw out optical or PCR duplicate reads.  `samtools depth` excludes them anyways, and keeping them in the analysis can cause errors later on
* Adding option to throw out multi-mapped reads.  Keeping these reads can overestimate the readcounts of some genes, as the read will now count for two alignments.
* Fixed small symlinking bug where the script will error if the symlinked source and destination are the same
* Will create temporary directory if it currently does not exist
* Improved error handling in some areas of the script (raising errors instead of using assertions, for example)
* Properly paired BAM files will also be indexed
* Fixed issue where the range of fragment depth was not including the ending coordinate
* Ran the code through a couple of Python linters to clean it up
* Found issues if the BAM file was being operated on as a relative path, so ensuring absolute paths are read.
* Fixing bug where KeyError was encountered if properly paired reads overlap and there is a gap in the overlap leading to 0-depth for some coordinates
