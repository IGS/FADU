# Notes about the benchmarking scripts

NOTE: These benchmarking scripts are outdated and for a previous version of FADU

## Omissions
Our benchmarking of these scripts took place on the Sun Grid Engine (now known as Oracle Grid Engine).  The scripts had extra code in place to schedule the alignment and quantification jobs onto the grid, and this code has been removed.  In addition, certain file paths were either generalized or had a placeholder like ###DATADIR### inserted to scrub out references to our filesystem structure.

## Organism scripts
The following 3 script files contain paths to various files and settings that were used when benchmarking for that organism:

echaffeensis\_inputs.sh
ecoli\_inputs.sh
wbm\_inputs.sh

When any of these scripts is passed as an argument to any of the "quantification" scripts that script will source the variables of the organism script and use them in the quantification program commands.  If you wish to benchmark against your own organism, you can create a script file that is patterned after these.

## Running a benchmark script

The following scripts represent the following quantification programs and algorithms (if applicable) that were benchmarked against:

express.sh
fadu.sh
featurecounts\_fractionoverlap.sh
featurecounts\_overlap.sh
featurecounts.sh
htseq\_intersectionnonempty.sh
htseq\_intersectionstrict.sh
htseq\_union\_nonunique.sh
htseq\_union.sh
kallisto.sh
salmon.sh

To call one of these scripts run:
```bash
bash <program.sh> <organism.sh>
```
