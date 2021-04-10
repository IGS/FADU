# FADU

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
	- [Software](#software)
	- [Directories](#directories)
	- [Create directories](#create-directories)
- [Assess RNA quantification tool performance using simulated prokaryote data](#assess-rna-quantification-tool-performance-using-simulated-prokaryote-data)
	- [Simulate E. coli K-12 substr. MG1655 RNA-Seq data with operon annotations](#simulate-e-coli-k-12-substr-mg1655-rna-seq-data-with-operon-annotations)
		- [Download and format reference files](#download-and-format-reference-files)
		- [Edit GTF file to contain operon annotations](#edit-gtf-file-to-contain-operon-annotations)
			- [Set R inputs](#set-r-inputs)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
			- [Parse operon information for E. coli K-12 substr MG1655 from OperonDB](#parse-operon-information-for-e-coli-k-12-substr-mg1655-from-operondb)
			- [Edit E. coli K-12 substr. MG1655 gtf to include operon annotations](#edit-e-coli-k-12-substr-mg1655-gtf-to-include-operon-annotations)
		- [Simulate RNA-Seq FASTA files using operon GTF](#simulate-rna-seq-fasta-files-using-operon-gtf)
			- [Set R inputs](#set-r-inputs-1)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
			- [Create gene FASTA reference](#create-gene-fasta-reference)
			- [Create transcript FASTA reference using OperonDB predictions](#create-transcript-fasta-reference-using-operondb-predictions)
			- [Simulate RNA-Seq experiment using transcript FASTA reference](#simulate-rna-seq-experiment-using-transcript-fasta-reference)
		- [Convert simulated FASTAs to FASTQs](#convert-simulated-fastas-to-fastqs)
		- [Quantify simulated FASTQs using alignment-based quantification tools](#quantify-simulated-fastqs-using-alignment-based-quantification-tools)
			- [Create HISAT2 indexes and align simulated reads](#create-hisat2-indexes-and-align-simulated-reads)
			- [Sort and index BAM files](#sort-and-index-bam-files)
			- [Process GFF to contain only gene feature](#process-gff-to-contain-only-gene-feature)
			- [Remove Windows end-of-line characters from R-generated FASTAs](#remove-windows-end-of-line-characters-from-r-generated-fastas)
			- [Quantify simulated FASTQ files using alignment-based tools \(eXpress, FADU, featureCounts, HTSeq\)](#quantify-simulated-fastq-files-using-alignment-based-tools-express-fadu-featurecounts-htseq)
		- [Quantify simulated FASTQs using alignment-free quantification tools](#quantify-simulated-fastqs-using-alignment-free-quantification-tools)
			- [Create reference indexes](#create-reference-indexes)
			- [Quantify simulated FASTQ files using alignment-free tools \(kallisto, Salmon\)](#quantify-simulated-fastq-files-using-alignment-free-tools-kallisto-salmon)
	- [Determine whether FADU better predicts operon counts compared to other tools in a simulated data set](#determine-whether-fadu-better-predicts-operon-counts-compared-to-other-tools-in-a-simulated-data-set)
		- [Set R inputs](#set-r-inputs-2)
		- [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
		- [Load R functions](#load-r-functions)
		- [Establish list of samples](#establish-list-of-samples)
		- [Store counts data from each tool in a list consisting of dataframes](#store-counts-data-from-each-tool-in-a-list-consisting-of-dataframes)
		- [Parse operon information for E. coli K-12 substr MG1655 from OperonDB](#parse-operon-information-for-e-coli-k-12-substr-mg1655-from-operondb-1)
		- [Read and prep GTF reference files](#read-and-prep-gtf-reference-files)
		- [Conduct hiearchical clustering analysis on counts for simulated datasets](#conduct-hiearchical-clustering-analysis-on-counts-for-simulated-datasets)
		- [Identify simulated differentially expressed genes from list of differentially expressed transcripts](#identify-simulated-differentially-expressed-genes-from-list-of-differentially-expressed-transcripts)
		- [Conduct differential expression analysis](#conduct-differential-expression-analysis)
			- [DESeq2](#deseq2)
			- [edgeR](#edger)
		- [Construct table of DE stats](#construct-table-of-de-stats)
			- [DESeq2](#deseq2-1)
			- [edgeR](#edger-1)
		- [Plot stacked bar charts from DE stats](#plot-stacked-bar-charts-from-de-stats)
			- [Plot legend](#plot-legend)
			- [DESeq2](#deseq2-2)
			- [edgeR](#edger-2)
			- [Combine DESeq2 and edgeR plots](#combine-deseq2-and-edger-plots)
		- [Compare counts for each tool to expected counts from simulation](#compare-counts-for-each-tool-to-expected-counts-from-simulation)
			- [Plot the distributions of the actual versus expected values for each tool](#plot-the-distributions-of-the-actual-versus-expected-values-for-each-tool)
- [Assess RNA quantification tool performance using actual RNA-Seq data](#assess-rna-quantification-tool-performance-using-actual-rna-seq-data)
	- [Download and construct reference files](#download-and-construct-reference-files)
		- [Download E. chaffeensis, E. coli, and wBm reference files](#download-e-chaffeensis-e-coli-and-wbm-reference-files)
		- [Create CDS nucleotide fasta files](#create-cds-nucleotide-fasta-files)
			- [Input Sets](#input-sets)
			- [Commands](#commands)
		- [Create combined reference files for B. malayi + wBm](#create-combined-reference-files-for-b-malayi--wbm)
		- [Remove extra wBm entry in combined B. malayi + wBm reference](#remove-extra-wbm-entry-in-combined-b-malayi--wbm-reference)
	- [Download FASTQ data from SRA](#download-fastq-data-from-sra)
	- [Quantify FASTQs using alignment-based quantification tools](#quantify-fastqs-using-alignment-based-quantification-tools)
		- [Create HISAT2 indexes and align simulated reads](#create-hisat2-indexes-and-align-simulated-reads-1)
		- [Sort and index BAM files](#sort-and-index-bam-files-1)
		- [Quantify simulated FASTQ files using alignment-based tools \(eXpress, FADU, featureCounts, HTSeq\)](#quantify-simulated-fastq-files-using-alignment-based-tools-express-fadu-featurecounts-htseq-1)
	- [Quantify FASTQs using alignment-free quantification tools](#quantify-fastqs-using-alignment-free-quantification-tools)
		- [Create reference indexes](#create-reference-indexes-1)
		- [Quantify simulated FASTQ files using alignment-free tools \(kallisto, Salmon\)](#quantify-simulated-fastq-files-using-alignment-free-tools-kallisto-salmon-1)
	- [Prepare wBm files for more specific analyses](#prepare-wbm-files-for-more-specific-analyses)
		- [Split wBm genome BAM file by strandedness](#split-wbm-genome-bam-file-by-strandedness)
		- [Index stranded wBm BAM files](#index-stranded-wbm-bam-files)
		- [Calculate depth for stranded wBm BAM files](#calculate-depth-for-stranded-wbm-bam-files)
	- [Determine whether FADU better quantifies compared to other tools in actual data sets](#determine-whether-fadu-better-quantifies-compared-to-other-tools-in-actual-data-sets)
		- [Set R inputs](#set-r-inputs-3)
		- [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo-1)
		- [Load R functions](#load-r-functions-1)
		- [Establish list of samples](#establish-list-of-samples-1)
		- [Store counts data from each tool in a list consisting of dataframes](#store-counts-data-from-each-tool-in-a-list-consisting-of-dataframes-1)
		- [Remove B. malayi genes from wBm analysis](#remove-b-malayi-genes-from-wbm-analysis)
		- [Plot FADU log2counts v other tools' log2counts](#plot-fadu-log2counts-v-other-tools-log2counts)
			- [Set plot order](#set-plot-order)
			- [Plot FADU default v other tools](#plot-fadu-default-v-other-tools)
			- [Plot FADU 10em v other tools](#plot-fadu-10em-v-other-tools)
			- [Plot FADU no_multimapped v other tools](#plot-fadu-no_multimapped-v-other-tools)
			- [Create figure legend](#create-figure-legend)
		- [Conduct MA analysis](#conduct-ma-analysis)
			- [Select quantification methods for MA analysis](#select-quantification-methods-for-ma-analysis)
			- [Identify the methods over-/under-counted relative to default FADU](#identify-the-methods-over-under-counted-relative-to-default-fadu)
			- [Identify the genes over-/under-counted relative to default FADU](#identify-the-genes-over-under-counted-relative-to-default-fadu)
			- [Plot MA plots](#plot-ma-plots)
		- [Examine counts for genes within a specific operon](#examine-counts-for-genes-within-a-specific-operon)
			- [Set plot order](#set-plot-order-1)
			- [Set operon genes and start and stop coordinates](#set-operon-genes-and-start-and-stop-coordinates)
			- [Read and subset depth over operon region](#read-and-subset-depth-over-operon-region)
			- [Calculate normalized relative count values](#calculate-normalized-relative-count-values)
			- [Plot depth over operon region](#plot-depth-over-operon-region)
			- [Plot normalized relative count values of operon region](#plot-normalized-relative-count-values-of-operon-region)
			- [Create figure legend](#create-figure-legend-1)
			- [Combine plots](#combine-plots)
		- [Examine counts for genes under-counted by FADU](#examine-counts-for-genes-under-counted-by-fadu)
			- [Set under-counted genes and start and stop coordinates](#set-under-counted-genes-and-start-and-stop-coordinates)
			- [Read and subset depth over under-counted gene region](#read-and-subset-depth-over-under-counted-gene-region)
			- [Plot depth over under-counted gene region](#plot-depth-over-under-counted-gene-region)
			- [Plot count values of operon region](#plot-count-values-of-operon-region)
- [Generate scripts for benchmarking](#generate-scripts-for-benchmarking)
	- [Create benchmarking scripts directory](#create-benchmarking-scripts-directory)
	- [Create benchmarking scripts](#create-benchmarking-scripts)
	- [Remove extra scripts that do not need to be benchmarked](#remove-extra-scripts-that-do-not-need-to-be-benchmarked)
	- [Download Shaun's benchmarking stats](#download-shauns-benchmarking-stats)
	- [Plot benchmarking data](#plot-benchmarking-data)
		- [Set R inputs](#set-r-inputs-4)
		- [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo-2)
		- [Load R functions](#load-r-functions-2)
		- [Convert wallclock time units to seconds](#convert-wallclock-time-units-to-seconds)
		- [Parse benchmarking data to create timing data frames for each dataset](#parse-benchmarking-data-to-create-timing-data-frames-for-each-dataset)
		- [Parse benchmarking data to create VMem data frames for each dataset](#parse-benchmarking-data-to-create-vmem-data-frames-for-each-dataset)
		- [Set plot order](#set-plot-order-2)
		- [Create timing plots](#create-timing-plots)
		- [Create VMem plots](#create-vmem-plots)
		- [Combine timing and VMem plots](#combine-timing-and-vmem-plots)

<!-- /MarkdownTOC -->

# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software

```{bash, eval = F}
SCRIPTS_BIN_DIR=/home/mattchung/scripts

JULIA_BIN_DIR=/usr/local/bin
JULIA_LIB_DIR=/home/mattchung/.julia
PYTHON_BIN_DIR=/usr/local/packages/python-3.5.2/bin
PYTHON_LIB_DIR=/usr/local/packages/python-3.5.2/lib

BBTOOLS_BIN_DIR=/usr/local/packages/bbtools-38.47
DOS2UNIX_BIN_DIR=/usr/bin
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SEQTK_BIN_DIR=/usr/local/packages/seqtk-1.2/bin

EXPRESS_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/express-1.5.1
FADU_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/FADU_v1.7
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
KALLISTO_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/kallisto_v0.46.1
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
SUBREAD_BIN_DIR=/usr/local/packages/subread-1.6.4/bin
```

## Directories

```{bash, eval = F}
REFERENCES_DIR=/local/projects-t3/EBMAL/mchung_dir/fadu/references
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/fadu/
OUTPUT_DIR=/local/projects-t3/EBMAL/mchung_dir/fadu/output
```

## Create directories

```{bash, eval = F}
mkdir "$WORKING_DIR"/simulation/
mkdir "$WORKING_DIR"/simulation/reads
```

# Assess RNA quantification tool performance using simulated prokaryote data

## Simulate E. coli K-12 substr. MG1655 RNA-Seq data with operon annotations

### Download and format reference files
```{bash, eval = F}
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz

gunzip "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff.gz
gunzip "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gtf.gz

cp "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna "$REFERENCES_DIR"/NC_000913.3.fa

grep -v "#" "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff > "$REFERENCES_DIR"/temp
mv "$REFERENCES_DIR"/temp "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff
```

### Edit GTF file to contain operon annotations

#### Set R inputs

R inputs from bash should be:

REFERENCES.DIR = "$REFERENCES_DIR"

```{R}
REFERENCES.DIR <- "Z:/EBMAL/mchung_dir/fadu/references"
```
#### Load packages and view sessionInfo

```{R}
library(stringr)
library(rlist)

sessionInfo()
```

```{R, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] stringr_1.4.0 rlist_0.4.6.1

loaded via a namespace (and not attached):
[1] compiler_3.5.1    magrittr_1.5      tools_3.5.1       yaml_2.2.0        stringi_1.4.3    
[6] data.table_1.12.2 knitr_1.22        xfun_0.6
```

#### Parse operon information for E. coli K-12 substr MG1655 from OperonDB
```{R}
operon_page <- readLines("http://operondb.cbcb.umd.edu/cgi-bin/operondb/pairs.cgi?genome_id=376",warn=F)
operon_delimit <- c(grep("coords:",operon_page),length(operon_page))

operon_list <- list()
while(length(operon_delimit) > 1){
	operon_page_subset <- operon_page[(operon_delimit[1]+1):operon_delimit[2]]
	operon_page_subset <- gsub(".*term=","",grep("entrez",operon_page_subset,value = T))
	operon_list <- list.append(operon_list,
								unique(sort(gsub(". .*","",operon_page_subset))))
	operon_delimit <- operon_delimit[-1]
}
```

#### Edit E. coli K-12 substr. MG1655 gtf to include operon annotations
```{R}
gtf <- read.delim(paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.gtf"), header = F, comment.char = "#")
gtf <- gtf[gtf[,3] == "gene",]
gtf[,9] <- gsub(";.*",";",gtf[,9])

for(i in 1:length(operon.list)){
	if(length(operon_list[[i]]) > 0){
		start <- min(as.numeric(gtf[gtf[,9] %in% paste0("gene_id ", operon_list[[i]],";"),4]))
		stop <- max(as.numeric(gtf[gtf[,9] %in% paste0("gene_id ", operon_list[[i]],";"),5]))
		strand <- as.character(unique(gtf[gtf[,9] %in% paste0("gene_id ",operon_list[[i]],";"),7]))
		
		gtf <- gtf[!(gtf[,9] %in% paste0("gene_id ",operon_list[[i]],";")),]
		
		gtf <- as.data.frame(rbind(gtf,
								  c("NC_000913.3",
									"RefSeq",
									"gene",
									start,
									stop,
									".",
									strand,
									".",
									paste0("gene_id operon",str_pad(i, 3, pad = "0"),";"))))
	}
}

gtf <- gtf[order(as.numeric(as.character(gtf[,4]))),]

write.table(gtf,
			paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.operon.gtf"),
			row.names = F,
			col.names = F,
			quote = F,
			sep = "\t")
```

### Simulate RNA-Seq FASTA files using operon GTF

#### Set R inputs

R inputs from bash should be:

REFERENCES.DIR = "$REFERENCES_DIR"
WORKING.DIR = "$WORKING_DIR"/simulation/reads

```{R}
REFERENCES.DIR <- "Z:/EBMAL/mchung_dir/fadu/references"
WORKING.DIR <- "Z:/EBMAL/mchung_dir/fadu/"
```

#### Load packages and view sessionInfo

```{R}
library(Biostrings)
library(polyester)

sessionInfo()
```

```{R, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Biostrings_2.50.2   XVector_0.22.0      IRanges_2.16.0      S4Vectors_0.20.1    BiocGenerics_0.28.0 polyester_1.9.7    

loaded via a namespace (and not attached):
[1] zlibbioc_1.28.0  compiler_3.5.1   limma_3.38.3     tools_3.5.1      logspline_2.1.13 yaml_2.2.0
```

#### Create gene FASTA reference
```{R}
gene_fasta <- seq_gtf(gtf = paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.gtf"), 
					  seqs = REFERENCES.DIR,
					  feature = "transcript",
					  exononly = F,
					  idfield = "gene_id", 
					  attrsep = "; ")
writeXStringSet(gene_fasta, paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.gene.fna"))
```

#### Create transcript FASTA reference using OperonDB predictions

```{R}
gene_fasta <- seq_gtf(gtf = paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_iseed=enomic.operon.gtf"), 
					  seqs = REFERENCES.DIR,
					  feature = "transcript", 
					  exononly = F,
					  idfield = "gene_id", 
					  attrsep = "; ")
writeXStringSet(gene_fasta, paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.transcript.fna"))
```

#### Simulate RNA-Seq experiment using transcript FASTA reference

RNA-Seq experiment contains 2 biological groups with 2 replicates each and 300 up- and 300 down-regulated genes.
```{R}
num_reps <- c(2,2)

fold_change_matrix <- matrix(1,
							 ncol = 1,
							 nrow = count_transcripts(paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.operon.gtf"),
													  fasta = F,identifier="gene_id"))

set.seed(5)
fold_change_matrix[sample(1:nrow(fold_change_matrix),600)[1:300],1] <- 2
fold_change_matrix[sample(1:nrow(fold_change_matrix),600)[301:600],1] <- 0.5

#fold_change_matrix <- fold_change_matrix[1:987,]
simulate_experiment(fasta=paste0(REFERENCES.DIR,"/GCF_000005845.2_ASM584v2_genomic.transcript.fna"),
					outdir = paste0(WORKING.DIR,"/simulation/reads"),
					num_reps = num_reps,
					fold_changes = fold_change_matrix,
					strand_specific = T,
					gzip = T,
					seed = 5)
```
### Convert simulated FASTAs to FASTQs

##### Inputs
```{bash, eval = F}
READS_DIR="$WORKING_DIR"/simulation/reads
```

##### Commands
```{bash, eval = F}
for SAMPLE in $(find "$READS_DIR" -name "*[.]fasta.gz" | sed "s/.*\\///g" | sort | sed "s/_[12].fasta.gz//g" | uniq)
do
	"$BBTOOLS_BIN_DIR"/reformat.sh in="$READS_DIR"/"$SAMPLE"_1.fasta.gz in2="$READS_DIR"/"$SAMPLE"_2.fasta.gz out="$READS_DIR"/"$SAMPLE"_1.fastq.gz out2="$READS_DIR"/"$SAMPLE"_2.fastq.gz qfake=35
done
```

### Quantify simulated FASTQs using alignment-based quantification tools

#### Create HISAT2 indexes and align simulated reads
##### Inputs
```{bash, eval = F}
READS_DIR="$WORKING_DIR"/simulation/reads
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
THREADS=16
```

##### Commands
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/bam

"$HISAT2_BIN_DIR"/hisat2-build --large-index "$NUC_GENOME_FNA" "$NUC_GENOME_FNA"
"$HISAT2_BIN_DIR"/hisat2-build --large-index "$NUC_GENE_FNA" "$NUC_GENE_FNA"

for SAMPLE in $(find "$READS_DIR" -name "*[.]fasta.gz" | sed "s/.*\\///g" | sort | sed "s/_[12].fasta.gz//g" | uniq)
do
	echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENOME_FNA" -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/bam/"$SAMPLE".genome.bam -" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N hisat2 -wd "$WORKING_DIR"/bam

	echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENE_FNA" -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/bam/"$SAMPLE".transcript.bam -" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N hisat2 -wd "$WORKING_DIR"/bam
done
```

#### Sort and index BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
THREADS=16
```

##### Commands
```{bash, eval = F}
for SAMPLE in $(find "$BAM_DIR" -name "*[.]bam" | sed "s/.*\\///g" | sort | sed "s/[.].*//g" | uniq)
do
  "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SAMPLE".genome.bam -@ "$THREADS" -o "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam
  "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SAMPLE".transcript.bam -@ "$THREADS" -n -o "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam
  "$SAMTOOLS_BIN_DIR"/samtools index "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -@ "$THREADS"
done 
```

#### Process GFF to contain only gene feature

##### Inputs
```{bash, eval = F}
GFF3=/local/projects-t3/EBMAL/mchung_dir/fadu/references/GCF_000005845.2_ASM584v2_genomic.gff
```

```{bash, eval = F}
awk '$3 == "gene" {print $0}' "$GFF3" > $(echo "$GFF3" | sed "s/[.]gff/.gene.gff/g")""
```

#### Remove Windows end-of-line characters from R-generated FASTAs

##### Inputs
```{bash, eval = F}
## Input Set 1
FNA=/local/projects-t3/EBMAL/mchung_dir/fadu/references/GCF_000005845.2_ASM584v2_genomic.gene.fna

## Input Set 2
FNA=/local/projects-t3/EBMAL/mchung_dir/fadu/references/GCF_000005845.2_ASM584v2_genomic.transcript.fna
```

```{bash, eval = F}
"$DOS2UNIX_BIN_DIR"/dos2unix "$FNA"
```

#### Quantify simulated FASTQ files using alignment-based tools (eXpress, FADU, featureCounts, HTSeq)

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.gff
FEAT_TYPE=gene
ATTR_ID=ID
THREADS=16
STRANDED=yes
```

##### Commands
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/quant

for SAMPLE in $(find "$BAM_DIR" -name "*[.]bam" | sed "s/.*\\///g" | sort | sed "s/[.].*//g" | uniq)
do
if [ "$STRANDED" == reverse ]
then

	echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

elif [ "$STRANDED" == yes ]
then

	echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

else

	echo -e ""$EXPRESS_BIN_DIR"/express -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

fi
done
```
### Quantify simulated FASTQs using alignment-free quantification tools

#### Create reference indexes
##### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
```

##### Commands
```{bash, eval = F}
"$KALLISTO_BIN_DIR"/kallisto index -i "$NUC_GENE_FNA".kallisto.index --make-unique "$NUC_GENE_FNA"
"$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$NUC_GENE_FNA" -i "$NUC_GENE_FNA".salmon.index 
```

#### Quantify simulated FASTQ files using alignment-free tools (kallisto, Salmon)

##### Inputs
```{bash, eval = F}
READS_DIR="$WORKING_DIR"/simulation/reads
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
THREADS=16
STRANDED=yes
```

##### Commands
```{bash, eval = F}
for SAMPLE in $(find "$READS_DIR" -name "*[.]fasta.gz" | sed "s/.*\\///g" | sort | sed "s/_[12].fasta.gz//g" | uniq)
do
	echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/quant/"$SAMPLE".salmon_default.counts --validateMappings" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/

	echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/quant/"$SAMPLE".salmon_optimized.counts --validateMappings  --allowDovetail" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/

	if [ "$STRANDED" == reverse ]
	then
		echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --rf-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/
	elif [ "$STRANDED" == yes ]
	then
		echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --fr-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/

	else
		echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" -o "$WORKING_DIR"/quant/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/
	fi
done
```

## Determine whether FADU better predicts operon counts compared to other tools in a simulated data set

### Set R inputs

```{R}
QUANT_DIR <- "Z:/EBMAL/mchung_dir/fadu/quant"
OUTPUT_DIR <- "Z:/EBMAL/mchung_dir/fadu/output"
REFERENCES_DIR <- "Z:/EBMAL/mchung_dir/fadu/references"

EXPECTED_COUNTS_PATH <- "Z:/EBMAL/mchung_dir/fadu/simulation/reads/sim_counts_matrix.rda"
EXPECTED_DE_PATH <- "Z:/EBMAL/mchung_dir/fadu/simulation/reads/sim_tx_info.txt"
```

### Load R packages and view sessionInfo

```{R}
library(dendextend)
library(ggdendro)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(pvclust)
library(reshape2)
library(rlist)
library(see)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] see_0.3.0                   rlist_0.4.6.1               reshape2_1.4.3             
 [4] pvclust_2.0-0               gridExtra_2.3               ggplot2_3.2.0              
 [7] ggdendro_0.1-20             DESeq2_1.22.2               SummarizedExperiment_1.12.0
[10] DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0         
[13] Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
[16] IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0        
[19] dendextend_1.12.0           cowplot_1.0.0               edgeR_3.24.3               
[22] limma_3.38.3               

loaded via a namespace (and not attached):
 [1] bitops_1.0-6           bit64_0.9-7            insight_0.8.0          RColorBrewer_1.1-2    
 [5] tools_3.5.0            backports_1.1.4        R6_2.4.0               rpart_4.1-13          
 [9] Hmisc_4.2-0            DBI_1.0.0              lazyeval_0.2.2         colorspace_1.4-1      
[13] nnet_7.3-12            withr_2.1.2            tidyselect_0.2.5       bit_1.1-14            
[17] compiler_3.5.0         htmlTable_1.13.1       labeling_0.3           bayestestR_0.4.0      
[21] scales_1.0.0           checkmate_1.9.4        genefilter_1.64.0      ggridges_0.5.2        
[25] stringr_1.4.0          digest_0.6.20          foreign_0.8-70         XVector_0.22.0        
[29] base64enc_0.1-3        pkgconfig_2.0.2        htmltools_0.3.6        htmlwidgets_1.3       
[33] rlang_0.4.0            rstudioapi_0.10        RSQLite_2.1.2          acepack_1.4.1         
[37] dplyr_0.8.3            RCurl_1.95-4.12        magrittr_1.5           GenomeInfoDbData_1.2.0
[41] Formula_1.2-3          parameters_0.3.0       Matrix_1.2-14          Rcpp_1.0.2            
[45] munsell_0.5.0          viridis_0.5.1          stringi_1.4.3          MASS_7.3-51.4         
[49] zlibbioc_1.28.0        plyr_1.8.4             grid_3.5.0             blob_1.2.0            
[53] crayon_1.3.4           lattice_0.20-35        splines_3.5.0          annotate_1.60.1       
[57] locfit_1.5-9.1         zeallot_0.1.0          knitr_1.23             pillar_1.4.2          
[61] effectsize_0.0.1       geneplotter_1.60.0     XML_3.98-1.20          glue_1.3.1            
[65] latticeExtra_0.6-28    data.table_1.12.2      vctrs_0.2.0            gtable_0.3.0          
[69] purrr_0.3.2            assertthat_0.2.1       xfun_0.8               xtable_1.8-4          
[73] survival_2.41-3        viridisLite_0.3.0      tibble_2.1.3           AnnotationDbi_1.44.0  
[77] memoise_1.1.0          cluster_2.0.7-1       
```

### Load R functions
```{R}
formalize_names <- function(vector){
	vector <- gsub("express_default","eXpress",vector)
	vector <- gsub("express_optimized","eXpress\n-B10 --no-bias-correct",vector)

	vector <- gsub("fadu_default","FADU",vector)
	vector <- gsub("fadu_10em","FADU\n--em_iterations 10",vector)
	vector <- gsub("fadu_nomm","FADU\n--remove_multimapped",vector)

	 vector <- gsub("htseq_unionnonuniqueall","HTSeq\n-m union --nonunique all",vector)
	vector <- gsub("htseq_union","HTSeq\n-m union",vector)
	vector <- gsub("htseq_intersectionstrict","HTSeq\n-m intersection-strict",vector)
	vector <- gsub("htseq_intersectionnonempty","HTSeq\n-m intersection-nonempty",vector)
	

	vector <- gsub("featurecounts_default","featureCounts",vector)
	vector <- gsub("featurecounts_overlap","featureCounts\n-O",vector)
	vector <- gsub("featurecounts_fractionaloverlap","featureCounts\n-O --fraction",vector)

	vector <- gsub("salmon_default","Salmon\n--validateMappings",vector)
	vector <- gsub("salmon_optimized","Salmon\n--validateMappings --allowDovetail",vector)

	vector <- gsub("kallisto_default","kallisto",vector)
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 
```

### Establish list of samples

```{R}
samples <- c("sample_01","sample_02","sample_03","sample_04")
```

### Store counts data from each tool in a list consisting of dataframes

```{R}
quant_tools <- c("fadu_default",
                 "fadu_10em",
                 "fadu_nomm",
                 "express_default",
                 "express_optimized",
                 "featurecounts_default",
                 "featurecounts_overlap",
                 "featurecounts_fractionaloverlap",
                 "htseq_union",
                 "htseq_intersectionstrict",
                 "htseq_intersectionnonempty",
                 "htseq_unionnonuniqueall",
                 "kallisto_default",
                 "salmon_default",
                 "salmon_optimized")

counts <- list()
for(i in 1:length(samples)){
  counts[[i]] <- as.data.frame(matrix(nrow = nrow(read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))),
                                 ncol = length(quant_tools)))
  rownames(counts[[i]]) <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,1]
  colnames(counts[[i]]) <- quant_tools
  
  counts[[i]]$fadu_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  counts[[i]]$fadu_10em <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_10em.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  counts[[i]]$fadu_nomm <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_nomm.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  
  counts[[i]]$express_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".express_default.counts/results.xprs"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".express_default.counts/results.xprs"))[,2]),8]
  counts[[i]]$express_optimized <- read.delim(paste0(QUANT_DIR,"/",samples[i],".express_optimized.counts/results.xprs"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".express_optimized.counts/results.xprs"))[,2]),8]
  
  counts[[i]]$featurecounts_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_default.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_default.counts"),comment.char="#")[,1]),7]
  counts[[i]]$featurecounts_overlap <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_overlap.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_overlap.counts"),comment.char="#")[,1]),7]
  counts[[i]]$featurecounts_fractionaloverlap <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_fractionaloverlap.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_fractionaloverlap.counts"),comment.char="#")[,1]),7]
  
  counts[[i]]$htseq_union <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_union.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_union.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_intersectionstrict <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionstrict.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionstrict.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_intersectionnonempty <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionnonempty.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionnonempty.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_unionnonuniqueall <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_unionnonuniqueall.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_unionnonuniqueall.counts"),header = F)[,1]),2]
  
  counts[[i]]$kallisto_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".kallisto_default.counts/abundance.tsv"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".kallisto_default.counts/abundance.tsv"))[,1]),4]
  
  counts[[i]]$salmon_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_default.counts/quant.sf"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_default.counts/quant.sf"))[,1]),5]
  counts[[i]]$salmon_optimized <- read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_optimized.counts/quant.sf"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_optimized.counts/quant.sf"))[,1]),5]
}

```

### Parse operon information for E. coli K-12 substr MG1655 from OperonDB

```{R}
operon_page <- readLines("http://operondb.cbcb.umd.edu/cgi-bin/operondb/pairs.cgi?genome_id=376",warn=F)
operon_delimit <- c(grep("coords:",operon_page),length(operon_page))

operon_list <- list()
while(length(operon_delimit) > 1){
	operon_page_subset <- operon_page[(operon_delimit[1]+1):operon_delimit[2]]
	operon_page_subset <- gsub(".*term=","",grep("entrez",operon_page_subset,value = T))
	operon_list <- list.append(operon_list,
								unique(sort(gsub(". .*","",operon_page_subset))))
	operon_delimit <- operon_delimit[-1]
}
operon_list <- operon_list[lapply(operon_list,length)>0]
```

### Read and prep GTF reference files

```{R}
gtf <- read.delim(paste0(REFERENCES_DIR,"/GCF_000005845.2_ASM584v2_genomic.gtf"),header = F)
gene_gtf <- gtf[gtf[,3] == "gene",]
operon_gtf <- read.delim(paste0(REFERENCES_DIR,"/GCF_000005845.2_ASM584v2_genomic.operon.gtf"),header = F)

pseudogenes <- gene_gtf[grep("gene_biotype pseudogene",gene_gtf[,9]),9]

gene_gtf[,9] <- gsub("gene_id b","gene-b",gene_gtf[,9])
gene_gtf[,9] <- gsub(";.*","",gene_gtf[,9])
operon_gtf[,9] <- gsub("gene_id b","gene-b",operon_gtf[,9])
operon_gtf[,9] <- gsub(";.*","",operon_gtf[,9])
pseudogenes <- gsub("gene_id b","gene-b",pseudogenes)
pseudogenes <- gsub(";.*","",pseudogenes)

operon_gtf[grep("operon",operon_gtf[,9]),9] <- gsub("gene_id ","",operon_gtf[grep("operon",operon_gtf[,9]),9])
names(operon_list) <- operon_gtf[grep("operon",operon_gtf[,9]),9]
```

### Conduct hiearchical clustering analysis on counts for simulated datasets

```{R,fig.height = 4, fig.width = 6}
pvclust <- pvclust(log2(counts[[1]]+1), method.hclust="average",
                   method.dist="correlation",
                   nboot=100,
                   iseed = 5)

structure <- hang.dendrogram(as.dendrogram(pvclust$hclust))
structure <- capture.output(str(structure))
structure <- structure[grepl("leaf", structure)]
structure  <- as.numeric(as.character(substr(structure,
             regexpr("h=", structure )+ 3,
             regexpr("  )", structure ))))
# colorkey <- read.table(colorkey.path)
# color <- substr(colnames(df), regexpr("Bm", colnames(df)), regexpr("_.$", colnames(df)) - 1)
# for(i in 1:nrow(colorkey)){
#   color <- gsub(colorkey[i,1], colorkey[i,2], color)
# }
# shape <- c()
# for(i in 1:length(colnames(tpm.de))){
#   if(grepl("agss", colnames(tpm.de)[i])){
#     shape[i] <- 16
#   }else{
#     shape[i] <- 17
#   }
# }
# dendroshape <- shape[result$hclust$order]
dendrocolor <- gsub("_.*","",colnames(counts[[1]]))[pvclust$hclust$order]

dendro.data <- dendro_data(pvclust$hclust)
dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
for(i in 1:nrow(dendro.data)){
  dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
}
dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]

bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro.data$y[duplicated(dendro.data$y)]),
                                            ncol = 2))
for(i in 1:length(dendro.data$y[duplicated(dendro.data$y)])){
  dendro.data.subset <- dendro.data[which(dendro.data$y == dendro.data$y[duplicated(dendro.data$y)][i]),]
  bootstrap.positions[i,1] <- unique(dendro.data.subset$x)
  bootstrap.positions[i,2] <- unique(dendro.data.subset$y)
}

points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                 structure))
pvclust$hclust$labels <- formalize_names(pvclust$hclust$labels)

root.dendro <- ggdendrogram(hang.dendrogram(as.dendrogram(pvclust$hclust)), theme_dendro = T)+
  geom_point(aes(x=points.df[,1], y = points.df[,2], color = dendrocolor), size = 3)+
  guides(color = F)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
for(i in 1:length(pvclust$edges$bp)){
  text <- round(pvclust$edges$bp[i] * 100,0)
  root.dendro <- root.dendro + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.2, y=bootstrap.positions[i,2] + 0.012, size = 2)
}

pdf(paste0(OUTPUT_DIR,"/counts_dendro.pdf"),
    height=4,
    width=6)
print(root.dendro)
dev.off()

png(paste0(OUTPUT_DIR,"/counts_dendro.png"),
    height=4,
    width=6,
    units = "in",res=300)
print(root.dendro)
dev.off()

print(root.dendro)
```

### Identify simulated differentially expressed genes from list of differentially expressed transcripts

```{R}
expected_de <- read.delim(EXPECTED_DE_PATH)
expected_de <- expected_de[expected_de[,3] == T,]
expected_de_up_transcripts <- expected_de[expected_de[,2] == 2,1]
expected_de_down_transcripts <- expected_de[expected_de[,2] == 0.5,1]

length(expected_de_up_transcripts)
length(expected_de_down_transcripts)
```

```{R, eval = F}
[1] 256
[1] 300
```

```{R}
expected_de_up_genes <- expected_de_up_transcripts
expected_de_down_genes <- expected_de_down_transcripts

for(i in grep("operon",expected_de_up_genes)){
  operon <- gsub(";","",expected_de_up_genes[i])
  expected_de_up_genes <- c(as.character(expected_de_up_genes),operon_list[[operon]])
}

for(i in grep("operon",expected_de_down_genes)){
  operon <- gsub(";","",expected_de_down_genes[i])
  expected_de_down_genes <- c(as.character(expected_de_down_genes),operon_list[[operon]])
}

expected_de_up_genes <- unique(paste0("gene-",gsub(";","",expected_de_up_genes[grep("operon",expected_de_up_genes,invert=T)])))
expected_de_down_genes <- unique(paste0("gene-",gsub(";","",expected_de_down_genes[grep("operon",expected_de_down_genes,invert=T)])))

length(expected_de_up_genes)
length(expected_de_down_genes)
```

```{R, eval = F}
[1] 585
[1] 678
```

### Conduct differential expression analysis

#### DESeq2
```{R}
FDRcutoff <- 0.05

deseq2.pairwise.degenes <- list()
deseq2.pairwise.excludegenes <- list()
for(i in 1:ncol(counts[[1]])){
  counts.pairwise <- as.data.frame(cbind(counts[[1]][,i],
                                         counts[[2]][,i],
                                         counts[[3]][,i],
                                         counts[[4]][,i]))
  rownames(counts.pairwise) <- rownames(counts[[1]])
  deseq.groups <- as.data.frame(matrix(nrow=4,
                                       ncol=1))
  deseq.groups[,1] <- as.factor(c("a","a","b","b"))

  rownames(deseq.groups) <- colnames(counts.pairwise)
  colnames(deseq.groups) <- "condition"
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts.pairwise),
					colData = deseq.groups,
					design = ~ condition)
  dds <- DESeq(dds, test="LRT", reduced=~1)
  keep <- rowSums(counts(dds)) >= 10
  deseq2.pairwise.excludegenes[[i]] <-  rownames(counts.pairwise)[which(keep == F)]

  res <- results(dds, cooksCutoff=T)
  res <- res[!is.na(res$padj),]
  deseq2.pairwise.degenes[[i]] <- rownames(res)[res$padj < FDRcutoff]
}
names(deseq2.pairwise.degenes) <- colnames(counts[[1]])
names(deseq2.pairwise.excludegenes) <- colnames(counts[[1]])
```

#### edgeR
```{R}
FDRcutoff <- 0.05

edgeR.pairwise.degenes <- list()
edgeR.pairwise.excludegenes <- list()

for(i in 1:ncol(counts[[1]])){
  counts.pairwise <- as.data.frame(cbind(counts[[1]][,i],
                                         counts[[2]][,i],
                                         counts[[3]][,i],
                                         counts[[4]][,i]))
  rownames(counts.pairwise) <- rownames(counts[[1]])
  groups.pairwise <- c(1,1,2,2)
  #groups.pairwise[,2] <- factor(groups.pairwise[,2],levels=unique(groups.pairwise[,2]))
  
  cpm.cutoff <- 5/min(colSums(counts.pairwise)) * 1000000
  
  y <- DGEList(counts = counts.pairwise, group = groups.pairwise)
  y <- calcNormFactors(y)
  keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups.pairwise))
  edgeR.pairwise.excludegenes[[i]] <- rownames(counts.pairwise)[keep == F]

  y <- y[keep, , keep.lib.sizes = F]
  design <- model.matrix(~groups.pairwise)
  y <- estimateDisp(y , design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
  edgeR.pairwise.degenes[[i]] <- rownames(qlf$table)[qlf$table$padj < FDRcutoff]
}
names(edgeR.pairwise.degenes) <- colnames(counts[[1]])
names(edgeR.pairwise.excludegenes) <- colnames(counts[[1]])
```

### Construct table of DE stats

#### DESeq2
```{R}
deseq2_de_stats <- as.data.frame(matrix(nrow=length(deseq2.pairwise.degenes),
                                 ncol = 0))
for(i in 1:length(deseq2.pairwise.degenes)){
  deseq2_de_stats$'Quantification Method'[i] <- names(deseq2.pairwise.degenes)[i]
  deseq2_de_stats$'DE Genes Excluded by Expression Threshold'[i] <- length(intersect(edgeR.pairwise.excludegenes[[i]],unique(expected_de_up_genes,expected_de_down_genes)))
  deseq2_de_stats$'NDE Genes Excluded by Expression Threshold'[i] <- length(deseq2.pairwise.excludegenes[[i]]) - deseq2_de_stats$'DE Genes Excluded by Expression Threshold'[i]
  deseq2_de_stats$'DE Genes Identified'[i] <- length(deseq2.pairwise.degenes[[i]])
  deseq2_de_stats$'Correct DE Genes'[i] <- length(intersect(expected_de_up_genes,deseq2.pairwise.degenes[[i]])) + length(intersect(expected_de_down_genes,deseq2.pairwise.degenes[[i]]))
  deseq2_de_stats$'Correct Upregulated DE Genes'[i] <- length(intersect(expected_de_up_genes,deseq2.pairwise.degenes[[i]]))
  deseq2_de_stats$'Correct Downregulated DE Genes'[i] <- length(intersect(expected_de_down_genes,deseq2.pairwise.degenes[[i]]))
}
deseq2_de_stats$'False Positive DE Genes' <- deseq2_de_stats$'DE Genes Identified' - deseq2_de_stats$'Correct DE Genes'
deseq2_de_stats$'Missed DE Genes' <- length(unique(c(expected_de_up_genes,expected_de_down_genes))) - deseq2_de_stats$'Correct DE Genes'
deseq2_de_stats$'Correct NDE Genes' <- nrow(counts[[1]]) - deseq2_de_stats$'Correct DE Genes' - deseq2_de_stats$'False Positive DE Genes' - deseq2_de_stats$'Missed DE Genes'

# deseq2_de_stats$correct <- deseq2_de_stats$'Correct Upregulated DE Genes' + deseq2_de_stats$'Correct Downregulated DE Genes' + deseq2_de_stats$'Correct NDE Genes'
# deseq2_de_stats$wrong <- deseq2_de_stats$'False Positive DE Genes' + deseq2_de_stats$'Missed DE Genes'

deseq2_de_stats$'Quantification Method' <- formalize_names(deseq2_de_stats$'Quantification Method')
```

#### edgeR
```{R}
edgeR_de_stats <- as.data.frame(matrix(nrow=length(edgeR.pairwise.degenes),
                                 ncol = 0))
for(i in 1:length(edgeR.pairwise.degenes)){
  edgeR_de_stats$'Quantification Method'[i] <- names(edgeR.pairwise.degenes)[i]
  edgeR_de_stats$'DE Genes Excluded by Expression Threshold'[i] <- length(intersect(edgeR.pairwise.excludegenes[[i]],unique(expected_de_up_genes,expected_de_down_genes)))
  edgeR_de_stats$'NDE Genes Excluded by Expression Threshold'[i] <- length(edgeR.pairwise.excludegenes[[i]]) - edgeR_de_stats$'DE Genes Excluded by Expression Threshold'[i]
  edgeR_de_stats$'DE Genes Identified'[i] <- length(edgeR.pairwise.degenes[[i]])
  edgeR_de_stats$'Correct DE Genes'[i] <- length(intersect(expected_de_up_genes,edgeR.pairwise.degenes[[i]])) + length(intersect(expected_de_down_genes,edgeR.pairwise.degenes[[i]]))
  edgeR_de_stats$'Correct Upregulated DE Genes'[i] <- length(intersect(expected_de_up_genes,edgeR.pairwise.degenes[[i]]))
  edgeR_de_stats$'Correct Downregulated DE Genes'[i] <- length(intersect(expected_de_down_genes,edgeR.pairwise.degenes[[i]]))
}
edgeR_de_stats$'False Positive DE Genes' <- edgeR_de_stats$'DE Genes Identified' - edgeR_de_stats$'Correct DE Genes'
edgeR_de_stats$'Missed DE Genes' <- length(unique(c(expected_de_up_genes,expected_de_down_genes))) - edgeR_de_stats$'Correct DE Genes'
edgeR_de_stats$'Correct NDE Genes' <- nrow(counts[[1]]) - edgeR_de_stats$'Correct DE Genes' - edgeR_de_stats$'False Positive DE Genes' - edgeR_de_stats$'Missed DE Genes'

edgeR_de_stats$'Quantification Method' <- formalize_names(edgeR_de_stats$'Quantification Method')
```

### Plot stacked bar charts from DE stats

#### Plot legend
```{R,fig.height=5,fig.width=8}
deseq2_de_stats.plot_df <- melt(deseq2_de_stats)
plot_categories <- c("Correct DE Genes",
                     "Correct NDE Genes",
                     "Missed DE Genes",
                     "False Positive DE Genes",
                     "DE Genes Excluded by Expression Threshold",
                     "NDE Genes Excluded by Expression Threshold")
deseq2_de_stats.plot_df <- deseq2_de_stats.plot_df[deseq2_de_stats.plot_df[,2] %in% plot_categories,]
deseq2_de_stats.plot_df[,1] <- factor(deseq2_de_stats.plot_df[,1],levels=formalize_names(colnames(counts[[1]]))[pvclust$hclust$order])
deseq2_de_stats.plot_df[,2] <- factor(deseq2_de_stats.plot_df[,2],levels=plot_categories)

deseq2_main.plot <- ggplot(mapping=aes(x=deseq2_de_stats.plot_df $'Quantification Method', 
                                       y=deseq2_de_stats.plot_df[,3], 
                                       fill=deseq2_de_stats.plot_df[,2],
                                       label=deseq2_de_stats.plot_df[,3]))+
  geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="genes",fill="Gene Categories:")+
  scale_y_continuous(expand=c(0.03,0.03))+
  scale_fill_manual(values=rev(c("tan","indianred2","#9B2915","#50A2A7","#688E26","#FAA613")))+
  coord_flip()+
  theme_bw()+
  theme(legend.position="bottom")

legend <- g_legend(deseq2_main.plot)

pdf(paste0(OUTPUT_DIR,"/gene_legend.pdf"),
    height=2,
    width=7)
grid.arrange(legend)
dev.off()

png(paste0(OUTPUT_DIR,"/gene_legend.png"),
    height=2,
    width=7,
	units = "in",res=300)
grid.arrange(legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/gene_legend.png)

#### DESeq2

```{R,fig.height=5,fig.width=8}
deseq2_de_stats.plot_df <- melt(deseq2_de_stats)
plot_categories <- c("Correct DE Genes",
                    "Correct NDE Genes",
                    "Missed DE Genes",
                    "False Positive DE Genes")
deseq2_de_stats.plot_df <- deseq2_de_stats.plot_df[deseq2_de_stats.plot_df[,2] %in% plot_categories,]
deseq2_de_stats.plot_df[,1] <- factor(deseq2_de_stats.plot_df[,1],levels=formalize_names(colnames(counts[[1]]))[pvclust$hclust$order])
deseq2_de_stats.plot_df[,2] <- factor(deseq2_de_stats.plot_df[,2],levels=rev(plot_categories))

deseq2_main.plot <- ggplot(mapping=aes(x=deseq2_de_stats.plot_df $'Quantification Method', 
                                       y=deseq2_de_stats.plot_df[,3], 
                                       fill=deseq2_de_stats.plot_df[,2],
                                       label=deseq2_de_stats.plot_df[,3]))+
  geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="genes",fill="Gene Categories:")+
  scale_y_continuous(expand=c(0.03,0.03))+
  scale_fill_manual(values=c("#9B2915","#50A2A7","#688E26","#FAA613"))+
  guides(fill=F)+
  coord_flip()+
  theme_bw()+
  theme(legend.position="bottom")

deseq2_de_stats.sideplot_df <- melt(deseq2_de_stats)
plot_categories <- c("DE Genes Excluded by Expression Threshold",
                    "NDE Genes Excluded by Expression Threshold")
deseq2_de_stats.sideplot_df <- deseq2_de_stats.sideplot_df[deseq2_de_stats.sideplot_df[,2] %in% plot_categories,]
deseq2_de_stats.sideplot_df[,1] <- factor(deseq2_de_stats.sideplot_df[,1],levels=formalize_names(colnames(counts[[1]]))[pvclust$hclust$order])

deseq2_side.plot <- ggplot(mapping=aes(x=deseq2_de_stats.sideplot_df $'Quantification Method', 
                                             y=deseq2_de_stats.sideplot_df[,3], 
                                             fill=deseq2_de_stats.sideplot_df[,2],
                                             label=deseq2_de_stats.sideplot_df[,3]))+
  geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="genes",fill="Gene Categories:")+
  scale_y_continuous(expand=c(0.05,0.05))+
  scale_fill_manual(values=c("indianred2","tan"))+
  coord_flip()+
  guides(fill=F)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        legend.position="bottom")
```

#### edgeR
```{R,fig.height=5,fig.width=8}
edgeR_de_stats.plot_df <- melt(edgeR_de_stats)
plot_categories <- c("Correct DE Genes",
                    "Correct NDE Genes",
                    "Missed DE Genes",
                    "False Positive DE Genes")
edgeR_de_stats.plot_df <- edgeR_de_stats.plot_df[edgeR_de_stats.plot_df[,2] %in% plot_categories,]
edgeR_de_stats.plot_df[,1] <- factor(edgeR_de_stats.plot_df[,1],levels=formalize_names(colnames(counts[[1]]))[pvclust$hclust$order])
edgeR_de_stats.plot_df[,2] <- factor(edgeR_de_stats.plot_df[,2],levels=rev(plot_categories))

edgeR_main.plot <- ggplot(mapping=aes(x=edgeR_de_stats.plot_df $'Quantification Method', 
                                       y=edgeR_de_stats.plot_df[,3], 
                                       fill=edgeR_de_stats.plot_df[,2],
                                       label=edgeR_de_stats.plot_df[,3]))+
  geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="genes",fill="Gene Categories:")+
  scale_y_continuous(expand=c(0.03,0.03))+
  scale_fill_manual(values=c("#9B2915","#50A2A7","#688E26","#FAA613"))+
  guides(fill=F)+
  coord_flip()+
  theme_bw()+
  theme(legend.position="bottom")

edgeR_de_stats.sideplot_df <- melt(edgeR_de_stats)
plot_categories <- c("DE Genes Excluded by Expression Threshold",
                    "NDE Genes Excluded by Expression Threshold")
edgeR_de_stats.sideplot_df <- edgeR_de_stats.sideplot_df[edgeR_de_stats.sideplot_df[,2] %in% plot_categories,]
edgeR_de_stats.sideplot_df[,1] <- factor(edgeR_de_stats.sideplot_df[,1],levels=formalize_names(colnames(counts[[1]]))[pvclust$hclust$order])

edgeR_side.plot <- ggplot(mapping=aes(x=edgeR_de_stats.sideplot_df $'Quantification Method', 
                                             y=edgeR_de_stats.sideplot_df[,3], 
                                             fill=edgeR_de_stats.sideplot_df[,2],
                                             label=edgeR_de_stats.sideplot_df[,3]))+
  geom_bar(stat="identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="genes",fill="Gene Categories:")+
  scale_y_continuous(expand=c(0.05,0.05))+
  scale_fill_manual(values=c("indianred2","tan"))+
  coord_flip()+
  guides(fill=F)+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        legend.position="bottom")
```

#### Combine DESeq2 and edgeR plots
```{R,fig.height=9,fig.width=8}
pdf(paste0(OUTPUT_DIR,"/sim_de_plot.pdf"),
    height=9,
    width=8)
plot_grid(deseq2_main.plot,deseq2_side.plot,
          edgeR_main.plot,edgeR_side.plot,rel_widths=c(7,1))
dev.off()

png(paste0(OUTPUT_DIR,"/sim_de_plot.png"),
    height=9,
    width=8,
	units = "in",res=300)
plot_grid(deseq2_main.plot,deseq2_side.plot,
          edgeR_main.plot,edgeR_side.plot,rel_widths=c(7,1))
dev.off()

plot_grid(deseq2_main.plot,deseq2_side.plot,
          edgeR_main.plot,edgeR_side.plot,rel_widths=c(7,1))

```

![image](/images/sim_de_plot.png)

### Compare counts for each tool to expected counts from simulation

For each gene, the counts obtained from each tool were divided by the expected number of counts as simulated. For operonic genes, the expected number of counts was obtained by dividing the expected operonic gene counts by the proportion of the operon that the gene covers.

values close to 1 indicate a similar number of counts were obtained compared to the simulated/expected values. Values > 1 indicate a tool is over-counting a gene while values < 1 indicate a tool is undercounting a gene.

```{R}
load(file = EXPECTED_COUNTS_PATH)
rownames(counts_matrix) <- gsub("^b","gene-b",rownames(counts_matrix))
rownames(counts_matrix) <- gsub(";","",rownames(counts_matrix))
counts_matrix <- counts_matrix[!(rownames(counts_matrix) %in% pseudogenes),]

actual_v_expected <- counts[[1]]

for(i in 1:nrow(counts_matrix)){
  if(grepl("operon",rownames(counts_matrix)[i],fixed = T)){
    operon_genes <- unlist(operon_list[grep(rownames(counts_matrix)[i],names(operon_list))])
    operon_length <- operon_gtf[operon_gtf[,9] == rownames(counts_matrix)[i],5] - operon_gtf[operon_gtf[,9] == rownames(counts_matrix)[i],4] + 1
    for(j in 1:length(operon_genes)){
      gene_length <- gene_gtf[gene_gtf[,9] == paste0("gene-",operon_genes[j]),5] - gene_gtf[gene_gtf[,9] == paste0("gene-",operon_genes[j]),4] + 1
      
      actual_v_expected[rownames(actual_v_expected) == paste0("gene-",operon_genes[j]),] <- actual_v_expected[rownames(actual_v_expected) == paste0("gene-",operon_genes[j]),]/(counts_matrix[i,1]*gene_length/operon_length)
    }
    
  }else{
    actual_v_expected[rownames(actual_v_expected) == rownames(counts_matrix)[i],] <- actual_v_expected[which(rownames(actual_v_expected) == rownames(counts_matrix)[i]),]/counts_matrix[i,1]
    
  }
}

print(colMeans(actual_v_expected))
```

```{R, eval = F}
                   fadu_default                       fadu_10em                       fadu_nomm 
                      0.9683484                       0.9673215                       0.9337055 
                express_default               express_optimized           featurecounts_default 
                      0.9385846                       0.9837095                       0.9573941 
          featurecounts_overlap featurecounts_fractionaloverlap                     htseq_union 
                      1.2598515                       1.0586884                       0.8841110 
       htseq_intersectionstrict      htseq_intersectionnonempty         htseq_unionnonuniqueall 
                      0.7334817                       0.8976557                       1.2939391 
               kallisto_default                  salmon_default                salmon_optimized 
                      0.9561832                       0.9581836                       0.9580046 
```

```{R}
log2 <- log2(actual_v_expected)
for(i in 1:ncol(log2)){
  print(paste(colnames(log2)[[i]],
              round(quantile(log2[,i])[2],2),
              round(quantile(log2[,i])[3],2),
              round(quantile(log2[,i])[4],2),
              sep = " | "))
  # print(colnames(log2)[[i]])
  # print(quantile(log2[,i]))
}
```

```{R, eval = F}
[1] "fadu_default | -0.17 | 0 | 0.08"
[1] "fadu_10em | -0.17 | 0 | 0.08"
[1] "fadu_nomm | -0.19 | 0 | 0.08"
[1] "express_default | -0.67 | -0.14 | 0.19"
[1] "express_optimized | -0.31 | -0.04 | 0.15"
[1] "featurecounts_default | -0.19 | 0 | 0.11"
[1] "featurecounts_overlap | 0 | 0.17 | 0.42"
[1] "featurecounts_fractionaloverlap | -0.07 | 0.02 | 0.19"
[1] "htseq_union | -0.47 | -0.07 | 0.05"
[1] "htseq_intersectionstrict | -0.87 | -0.32 | -0.01"
[1] "htseq_intersectionnonempty | -0.44 | -0.06 | 0.05"
[1] "htseq_unionnonuniqueall | 0 | 0.19 | 0.45"
[1] "kallisto_default | -0.28 | -0.01 | 0.06"
[1] "salmon_default | -0.25 | -0.02 | 0.03"
[1] "salmon_optimized | -0.25 | -0.02 | 0.03"
```

#### Plot the distributions of the actual versus expected values for each tool

```{R,fig.height = 6, fig.width = 8}
operon_ave.plot <- melt(actual_v_expected)
operon_ave.plot[,1] <- formalize_names(operon_ave.plot[,1])
operon_ave.plot[,2] <- log2(operon_ave.plot[,2])

operon_ave.plot[,1] <- factor(operon_ave.plot[,1], levels=formalize_names(colnames(actual_v_expected))[pvclust$hclust$order])

operon_ave.plot[,3] <- gsub("\n.*","",operon_ave.plot[,1])

actual_v_expected.plot <- ggplot(mapping=aes(x=operon_ave.plot[,1],y=operon_ave.plot[,2],fill=operon_ave.plot[,3]))+
  geom_violinhalf(scale = "width")+
  geom_text(data = means, aes(label = operon_ave.plot[,2], y = operon_ave.plot[,2] + 0.08))
  geom_boxplot(width=0.2)+
  labs(x="quantification method", y="log2 actual v. expected count ratio")+
  guides(fill = F)+
  coord_flip(ylim = c(-10,10))+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/actual_v_expected_sim_plot.pdf"),
    height=6,
    width=8)
print(actual_v_expected.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/actual_v_expected_sim_plot.png"),
    height=6,
    width=8,
    units = "in",res=300)
print(actual_v_expected.plot)
dev.off()

print(actual_v_expected.plot)
```

![image](/images/actual_v_expected_sim_plot.png)

```{R,fig.height = 6, fig.width = 4}
actual_v_expected_zoom.plot <- ggplot(mapping=aes(x=operon_ave.plot[,1],y=operon_ave.plot[,2],fill=operon_ave.plot[,3]))+
  geom_violinhalf(scale="width")+
  geom_boxplot(width=0.2)+
  labs(x="quantification method", y="log2 actual v. expected count ratio")+
  guides(fill = F)+
  #coord_cartesian()+
  coord_flip(ylim = c(-1,1))+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/actual_v_expected_sim_zoom_plot.pdf"),
    height=6,
    width=4)
print(actual_v_expected_zoom.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/actual_v_expected_sim_zoom_plot.png"),
    height=6,
    width=4,
    units = "in",res=300)
print(actual_v_expected_zoom.plot)
dev.off()

print(actual_v_expected_zoom.plot)
```

![image](/images/actual_v_expected_sim_zoom_plot.png)


# Assess RNA quantification tool performance using actual RNA-Seq data

## Download and construct reference files

### Download E. chaffeensis, E. coli, and wBm reference files

```{bash}
## E. chaffeensis
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/145/GCF_000013145.1_ASM1314v1/GCF_000013145.1_ASM1314v1_genomic.fna.gz
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/145/GCF_000013145.1_ASM1314v1/GCF_000013145.1_ASM1314v1_genomic.gff.gz
gunzip -f "$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.fna.gz
gunzip -f "$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gff.gz

## E. coli
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
gunzip -f "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip -f "$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff.gz

## wBm
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.fna.gz
wget -P "$REFERENCES_DIR" ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.gff.gz
gunzip -f "$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.fna.gz
gunzip -f "$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.gff.gz

## B. malayi
wget -P "$REFERENCES_DIR" ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WS275.genomic.fa.gz
wget -P "$REFERENCES_DIR" ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WS275.annotations.gff3.gz
gunzip -f "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.fa.gz
gunzip -f "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.annotations.gff3.gz
mv "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.fa "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.fna
```

### Create CDS nucleotide fasta files

#### Input Sets
```{bash}
## E. chaffeensis
FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gff

## E. coli
FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff

## wBm
FNA="$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.gff

## B. malayi
FNA="$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.fna
GFF3="$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.annotations.gff3
```

#### Commands
```{bash}
"$SCRIPTS_DIR"/createnuccdsfasta.sh -n "$FNA" -g "$GFF3" -f gene -i ID > "$(echo "$FNA" | sed "s/[.]fna$/.gene.fna/g")" &
```

### Create combined reference files for B. malayi + wBm

```{bash}
cat "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.fna "$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.fna > "$REFERENCES_DIR"/bmalayi_wbm.fna
cat "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.genomic.gene.fna "$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.gene.fna > "$REFERENCES_DIR"/bmalayi_wbm.gene.fna
cat "$REFERENCES_DIR"/b_malayi.PRJNA10729.WS275.annotations.gff3 "$REFERENCES_DIR"/GCF_000008385.1_ASM838v1_genomic.gff > "$REFERENCES_DIR"/bmalayi_wbm.gff
```

### Remove extra wBm entry in combined B. malayi + wBm reference
```{bash}
grep ">" "$REFERENCES_DIR"/bmalayi_wbm.fna | grep -v "wBm_Wolbachia" | sed "s/>//g" | cut -d" "-f1 > "$REFERENCES_DIR"/contigs.list
"$SAMTOOLS_BIN_DIR"/samtools faidx "$REFERENCES_DIR"/bmalayi_wbm.fna -r "$REFERENCES_DIR"/contigs.list > "$REFERENCES_DIR"/temp.fna
mv "$REFERENCES_DIR"/temp.fna "$REFERENCES_DIR"/bmalayi_wbm.fna
rm "$REFERENCES_DIR"/contigs.list
```

## Download FASTQ data from SRA

##### Input Sets
```{bash}
## E. chaffeensis
SRR=SRR1188323
OUTPUT_DIR="$WORKING_DIR"/echaffeensis

## E. coli
SRR=SRR2601722
OUTPUT_DIR="$WORKING_DIR"/ecoli

## wBm
SRR=SRR5192555
OUTPUT_DIR="$WORKING_DIR"/wbm
```

##### Commands)
```{bash}
qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$OUTPUT_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR" --gzip -O "$OUTPUT_DIR"
```

## Quantify FASTQs using alignment-based quantification tools

### Create HISAT2 indexes and align simulated reads

##### Input Sets
```{bash}
THREADS=16

## E. chaffeensis
SRR=SRR1188323
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/echaffeensis/SRR1188323_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/echaffeensis/SRR1188323_2.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/echaffeensis
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.fna

## E. coli
SRR=SRR2601722
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/ecoli/SRR2601722_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/ecoli/SRR2601722_2.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/ecoli
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna

## wBm
SRR=SRR5192555
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/wbm/SRR5192555_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/wbm/SRR5192555_2.fastq.gz
OUTPUT_DIR="$WORKING_DIR"/wbm
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wbm.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/bmalayi_wbm.fna
```

##### Commands
```{bash}
"$HISAT2_BIN_DIR"/hisat2-build --large-index "$NUC_GENOME_FNA" "$NUC_GENOME_FNA"
"$HISAT2_BIN_DIR"/hisat2-build --large-index "$NUC_GENE_FNA" "$NUC_GENE_FNA"

echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENOME_FNA" -1 "$FASTQ1" -2 "$FASTQ2" --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".genome.bam -" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=10G -N hisat2 -wd "$OUTPUT_DIR"

echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENE_FNA" -1 "$FASTQ1" -2 "$FASTQ2" --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".transcript.bam -" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=10G -N hisat2 -wd "$OUTPUT_DIR"
```

### Sort and index BAM files

##### Input Sets
```{bash}
THREADS=4

## E. chaffeensis
SRR=SRR1188323
BAM_DIR="$WORKING_DIR"/echaffeensis

## E. coli
SRR=SRR2601722
BAM_DIR="$WORKING_DIR"/ecoli

## wBm
SRR=SRR5192555
BAM_DIR="$WORKING_DIR"/wbm
```

##### Commands
```{bash}
for SAMPLE in $(find "$BAM_DIR" -name "*[.]bam" | sed "s/.*\\///g" | sort | sed "s/[.].*//g" | uniq | grep -v pseudoalignments)
do
  "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SRR".genome.bam -@ "$THREADS" -o "$BAM_DIR"/"$SRR".genome.sortedbyposition.bam
  #"$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SRR".transcript.bam -@ "$THREADS" -n -o "$BAM_DIR"/"$SRR".transcript.sortedbyname.bam
  "$SAMTOOLS_BIN_DIR"/samtools index "$BAM_DIR"/"$SRR".genome.sortedbyposition.bam -@ "$THREADS"
done
```

### Quantify simulated FASTQ files using alignment-based tools (eXpress, FADU, featureCounts, HTSeq)

##### Input Sets
```{bash}
FEAT_TYPE=gene
ATTR_ID=ID
THREADS=16

## E. chaffeensis
BAM_DIR="$WORKING_DIR"/echaffeensis
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gff
STRANDED=no

## E. coli
BAM_DIR="$WORKING_DIR"/ecoli
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gff
STRANDED=no

## wBm
BAM_DIR="$WORKING_DIR"/wbm
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wbm.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/bmalayi_wbm.fna
GFF3="$REFERENCES_DIR"/bmalayi_wbm.gff
STRANDED=reverse
```

##### Commands
```{bash, eval = F}
for SAMPLE in $(find "$BAM_DIR" -name "*[.]bam" | sed "s/.*\\///g" | sort | sed "s/[.].*//g" | uniq | grep -v pseudoalignments)
do
if [ "$STRANDED" == reverse ]
then

	echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

elif [ "$STRANDED" == yes ]
then

	echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

else

	echo -e ""$EXPRESS_BIN_DIR"/express -o "$WORKING_DIR"/quant/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/
	echo -e ""$EXPRESS_BIN_DIR"/express -B10 --no-bias-correct -o "$WORKING_DIR"/quant/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N express -wd "$WORKING_DIR"/quant/

	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_default.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_10em.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/
	echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".fadu_nomm.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N fadu -wd "$WORKING_DIR"/quant/

	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_union.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionstrict.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_intersectionnonempty.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/quant/"$SAMPLE".htseq_unionnonuniqueall.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N htseq -wd "$WORKING_DIR"/quant/

	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/
	echo -e ""$SUBREAD_BIN_DIR"/featureCounts -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/quant/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N featurecounts -wd "$WORKING_DIR"/quant/

fi
done
```

## Quantify FASTQs using alignment-free quantification tools

### Create reference indexes

##### Input Sets
```{bash, eval = F}
## E. chaffeensis
NUC_GENE_FNA="$REFERENCES_DIR"/canine_echaf.gene.fna

## E. coli
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna

## wBm
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wbm.gene.fna
```

##### Commands
```{bash, eval = F}
"$KALLISTO_BIN_DIR"/kallisto index -i "$NUC_GENE_FNA".kallisto.index --make-unique "$NUC_GENE_FNA"
"$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$NUC_GENE_FNA" -i "$NUC_GENE_FNA".salmon.index 
```

### Quantify simulated FASTQ files using alignment-free tools (kallisto, Salmon)

##### Inputs
```{bash, eval = F}
THREADS=16

## E. chaffeensis
SRR=SRR1188323
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/echaffeensis/SRR1188323_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/echaffeensis/SRR1188323_2.fastq.gz
NUC_GENE_FNA="$REFERENCES_DIR"/canine_echaf.gene.fna
STRANDED=no

## E. coli
SRR=SRR2601722
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/ecoli/SRR2601722_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/ecoli/SRR2601722_2.fastq.gz
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
STRANDED=no

## wBm
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wbm.gene.fna
STRANDED=reverse
SRR=SRR5192555
FASTQ1=/local/projects-t3/EBMAL/mchung_dir/fadu/wbm/SRR5192555_1.fastq.gz
FASTQ2=/local/projects-t3/EBMAL/mchung_dir/fadu/wbm/SRR5192555_2.fastq.gz
```

##### Commands
```{bash, eval = F}
echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$FASTQ1" -2 "$FASTQ2" -p "$THREADS" -o "$WORKING_DIR"/quant/"$SRR".salmon_default.counts --validateMappings" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/

echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$FASTQ1" -2 "$FASTQ2" -p "$THREADS" -o "$WORKING_DIR"/quant/"$SRR".salmon_optimized.counts --validateMappings  --allowDovetail" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/

if [ "$STRANDED" == reverse ]
then
	echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --rf-stranded -o $WORKING_DIR"/quant/"$SRR".kallisto_default.counts "$FASTQ1" "$FASTQ2" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/
elif [ "$STRANDED" == yes ]
then
	echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --fr-stranded -o $WORKING_DIR"/quant/"$SRR".kallisto_default.counts"$FASTQ1" "$FASTQ2" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/

else
	echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" -o $WORKING_DIR"/quant/"$SRR".kallisto_default.counts "$FASTQ1" "$FASTQ2" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N kallisto -wd "$WORKING_DIR"/quant/
fi
```

## Prepare wBm files for more specific analyses

### Split wBm genome BAM file by strandedness

##### Inputs
```{bash, eval = F}
THREADS=4
BAM="$WORKING_DIR"/wbm/SRR5192555.genome.sortedbyposition.bam
STRANDED=reverse
```

##### Commands
```{bash, eval = F}
"$SCRIPTS_BIN_DIR"/split_bam_by_strand.sh -i "$BAM" -s "$STRANDED" -t "$(dirname "$BAM")" -o "$(dirname "$BAM")" -@ "$THREADS"
```

### Index stranded wBm BAM files

##### Inputs
```{bash, eval = F}
BAM_F="$WORKING_DIR"/wbm/SRR5192555.genome.sortedbyposition.f.bam
BAM_R="$WORKING_DIR"/wbm/SRR5192555.genome.sortedbyposition.r.bam
```

##### Commands
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools index "$BAM_F"
"$SAMTOOLS_BIN_DIR"/samtools index "$BAM_R"
```

### Calculate depth for stranded wBm BAM files

##### Inputs
```{bash, eval = F}
CONTIG=NC_006833.1
BAM_F="$WORKING_DIR"/wbm/SRR5192555.genome.sortedbyposition.f.bam
BAM_R="$WORKING_DIR"/wbm/SRR5192555.genome.sortedbyposition.r.bam
```

##### Commands
```{bash, eval = F}
"$SAMTOOLS_BIN_DIR"/samtools depth -aa -d 1000000 -r "$CONTIG" "$BAM_F" > "$BAM_F".depth
"$SAMTOOLS_BIN_DIR"/samtools depth -aa -d 1000000 -r "$CONTIG" "$BAM_R" > "$BAM_R".depth
```

## Determine whether FADU better quantifies compared to other tools in actual data sets

### Set R inputs

```{R}
QUANT_DIR <- "Z:/EBMAL/mchung_dir/fadu/quant"
OUTPUT_DIR <- "C:/Users/MChung.SOM/Documents/plots"

GFF_PATH <- "Z:/EBMAL/mchung_dir/fadu/references/GCF_000008385.1_ASM838v1_genomic.gff"
DEPTH_F_PATH <- "Z:/EBMAL/mchung_dir/fadu/wbm/SRR5192555.genome.sortedbyposition.f.bam.depth"
DEPTH_R_PATH <- "Z:/EBMAL/mchung_dir/fadu/wbm/SRR5192555.genome.sortedbyposition.r.bam.depth"
```

### Load R packages and view sessionInfo

```{R}
library(cowplot)
library(dendextend)
library(DESeq2)
library(ggdendro)
library(gggenes)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(pvclust)
library(reshape2)
library(rlist)
library(see)
library(viridis)

sessionInfo()
```

```{R, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.5.1               viridisLite_0.3.0           see_0.3.0                   rlist_0.4.6.1              
 [5] reshape2_1.4.3              pvclust_2.0-0               gridExtra_2.3               ggrepel_0.8.0              
 [9] gggenes_0.3.2               ggdendro_0.1-20             DESeq2_1.22.2               SummarizedExperiment_1.12.0
[13] DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.55.0          Biobase_2.42.0             
[17] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         IRanges_2.16.0              S4Vectors_0.20.1           
[21] BiocGenerics_0.28.0         dendextend_1.9.0            cowplot_0.9.4               ggplot2_3.2.1              

loaded via a namespace (and not attached):
 [1] bitops_1.0-6           bit64_0.9-7            insight_0.7.1          RColorBrewer_1.1-2     prabclus_2.2-7        
 [6] tools_3.5.1            backports_1.1.3        R6_2.4.0               rpart_4.1-13           Hmisc_4.2-0           
[11] DBI_1.0.0              lazyeval_0.2.2         colorspace_1.4-1       trimcluster_0.1-2.1    nnet_7.3-12           
[16] withr_2.1.2            tidyselect_0.2.5       bit_1.1-14             compiler_3.5.1         htmlTable_1.13.1      
[21] bayestestR_0.4.0       diptest_0.75-7         scales_1.0.0           checkmate_1.9.1        DEoptimR_1.0-8        
[26] mvtnorm_1.0-10         robustbase_0.93-4      genefilter_1.64.0      ggridges_0.5.1         stringr_1.4.0         
[31] digest_0.6.18          foreign_0.8-70         XVector_0.22.0         base64enc_0.1-3        pkgconfig_2.0.2       
[36] htmltools_0.3.6        htmlwidgets_1.3        rlang_0.4.0            rstudioapi_0.10        RSQLite_2.1.1         
[41] mclust_5.4.3           acepack_1.4.1          dplyr_0.8.3            RCurl_1.95-4.12        magrittr_1.5          
[46] modeltools_0.2-22      GenomeInfoDbData_1.2.0 Formula_1.2-3          parameters_0.3.0       Matrix_1.2-14         
[51] Rcpp_1.0.1             munsell_0.5.0          ggfittext_0.6.0        stringi_1.4.3          whisker_0.3-2         
[56] yaml_2.2.0             MASS_7.3-50            zlibbioc_1.28.0        plyr_1.8.4             flexmix_2.3-15        
[61] grid_3.5.1             blob_1.1.1             crayon_1.3.4           lattice_0.20-35        splines_3.5.1         
[66] annotate_1.60.1        locfit_1.5-9.1         knitr_1.22             pillar_1.3.1           fpc_2.1-11.1          
[71] effectsize_0.0.1       geneplotter_1.60.0     XML_3.98-1.19          glue_1.3.1             latticeExtra_0.6-28   
[76] data.table_1.12.2      gtable_0.3.0           purrr_0.3.2            kernlab_0.9-27         assertthat_0.2.1      
[81] xfun_0.6               xtable_1.8-3           class_7.3-14           survival_2.42-3        tibble_2.1.1          
[86] AnnotationDbi_1.44.0   memoise_1.1.0          cluster_2.0.7-   
```

### Load R functions
```{R}
formalize_names <- function(vector){
	vector <- gsub("express_default","eXpress",vector)
	vector <- gsub("express_optimized","eXpress\n-B10 --no-bias-correct",vector)

	vector <- gsub("fadu_default","FADU",vector)
	vector <- gsub("fadu_10em","FADU\n--em_iterations 10",vector)
	vector <- gsub("fadu_nomm","FADU\n--remove_multimapped",vector)

	vector <- gsub("htseq_unionnonuniqueall","HTSeq\n-m union --nonunique all",vector)
	vector <- gsub("htseq_union","HTSeq\n-m union",vector)
	vector <- gsub("htseq_intersectionstrict","HTSeq\n-m intersection-strict",vector)
	vector <- gsub("htseq_intersectionnonempty","HTSeq\n-m intersection-nonempty",vector)
	

	vector <- gsub("featurecounts_default","featureCounts",vector)
	vector <- gsub("featurecounts_overlap","featureCounts\n-O",vector)
	vector <- gsub("featurecounts_fractionaloverlap","featureCounts\n-O --fraction",vector)

	vector <- gsub("salmon_default","Salmon\n--validateMappings",vector)
	vector <- gsub("salmon_optimized","Salmon\n--validateMappings --allowDovetail",vector)

	vector <- gsub("kallisto_default","kallisto",vector)
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

lm_eqn <- function(vector1,vector2){
  df <-as.data.frame(cbind(vector1,vector2))
  m <- lm(df[,2] ~ df[,1], df)
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
```

### Establish list of samples

```{R}
samples <- c("SRR1188323","SRR2601722","SRR5192555")
```

### Store counts data from each tool in a list consisting of dataframes

```{R}
quant_tools <- c("fadu_default",
                 "fadu_10em",
                 "fadu_nomm",
                 "express_default",
                 "express_optimized",
                 "featurecounts_default",
                 "featurecounts_overlap",
                 "featurecounts_fractionaloverlap",
                 "htseq_union",
                 "htseq_intersectionstrict",
                 "htseq_intersectionnonempty",
                 "htseq_unionnonuniqueall",
                 "kallisto_default",
                 "salmon_default",
                 "salmon_optimized")

counts <- list()
for(i in 1:length(samples)){
  counts[[i]] <- as.data.frame(matrix(nrow = nrow(read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))),
                                 ncol = length(quant_tools)))
  rownames(counts[[i]]) <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,1]
  colnames(counts[[i]]) <- quant_tools
  
  counts[[i]]$fadu_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_default.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  counts[[i]]$fadu_10em <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_10em.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  counts[[i]]$fadu_nomm <- read.delim(paste0(QUANT_DIR,"/",samples[i],".fadu_nomm.counts/",samples[i],".genome.sortedbyposition.counts.txt"))[,4]
  
  counts[[i]]$express_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".express_default.counts/results.xprs"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".express_default.counts/results.xprs"))[,2]),8]
  counts[[i]]$express_optimized <- read.delim(paste0(QUANT_DIR,"/",samples[i],".express_optimized.counts/results.xprs"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".express_optimized.counts/results.xprs"))[,2]),8]
  
  counts[[i]]$featurecounts_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_default.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_default.counts"),comment.char="#")[,1]),7]
  counts[[i]]$featurecounts_overlap <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_overlap.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_overlap.counts"),comment.char="#")[,1]),7]
  counts[[i]]$featurecounts_fractionaloverlap <- read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_fractionaloverlap.counts"),comment.char="#")[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".featurecounts_fractionaloverlap.counts"),comment.char="#")[,1]),7]
  
  counts[[i]]$htseq_union <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_union.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_union.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_intersectionstrict <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionstrict.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionstrict.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_intersectionnonempty <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionnonempty.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_intersectionnonempty.counts"),header = F)[,1]),2]
  counts[[i]]$htseq_unionnonuniqueall <- read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_unionnonuniqueall.counts"),header = F)[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".htseq_unionnonuniqueall.counts"),header = F)[,1]),2]
  
  counts[[i]]$kallisto_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".kallisto_default.counts/abundance.tsv"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".kallisto_default.counts/abundance.tsv"))[,1]),4]
  
  counts[[i]]$salmon_default <- read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_default.counts/quant.sf"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_default.counts/quant.sf"))[,1]),5]
  counts[[i]]$salmon_optimized <- read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_optimized.counts/quant.sf"))[match(rownames(counts[[i]]),read.delim(paste0(QUANT_DIR,"/",samples[i],".salmon_optimized.counts/quant.sf"))[,1]),5]
}
```

### Remove B. malayi genes from wBm analysis
```{R}
counts[[3]] <- counts[[3]][grep("Gene:",rownames(counts[[3]]),invert = T),]
```

### Plot FADU log2counts v other tools' log2counts

#### Set plot order 
```{R}
colorder <- c("fadu_default","fadu_10em","fadu_nomm",
              "salmon_optimized","salmon_default","kallisto_default","express_optimized","express_default",
              "featurecounts_default","htseq_intersectionnonempty","htseq_union","htseq_intersectionstrict",
              "htseq_unionnonuniqueall","featurecounts_overlap","featurecounts_fractionaloverlap")
for(i in 1:length(counts)){
  counts[[i]] <- counts[[i]][,match(colorder,colnames(counts[[i]]))]
}
```

#### Plot FADU default v other tools
```{R, fig.height=11, fig.width=4}
plot.list <- list()
for(i in 4:ncol(counts[[1]])){
  for(j in 1:length(counts)){
      plot.df <- as.data.frame(cbind(log2(counts[[j]][,1] + 1),
                                   log2(counts[[j]][,i] + 1)))
      plot.list <- list.append(plot.list,
                   ggplot(plot.df,aes_string(x=plot.df[,1],y=plot.df[,2]))+
                      stat_density_2d(aes(fill = ..level..), geom = "polygon")+
                      scale_fill_viridis()+
                      geom_abline(slope = 1, intercept = 0,color = "black", linetype = "dotted")+
                      annotate("text",x = 0, y = 14, label = lm_eqn(plot.df[,1],plot.df[,2]),parse = T,size = 2.5, hjust = 0)+
                      guides(fill = F)+
                      coord_cartesian(xlim=c(0,15),ylim=c(0,15))+
                      theme_bw()+
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank())
    )
  }
}

pdf(paste0(OUTPUT_DIR,"/actual_fadudefault_countsvcounts.pdf"),
    height=11,
    width=4)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/actual_fadudefault_countsvcounts.png"),
    height=11,
    width=4,
    units = "in",res=300)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

plot_grid(plotlist = plot.list,
          ncol=3)
```

![image](/images/actual_fadudefault_countsvcounts.png)

#### Plot FADU 10em v other tools
```{R fig.height=11, fig.width=4}
plot.list <- list()
for(i in 4:ncol(counts[[1]])){
  for(j in 1:length(counts)){
      plot.df <- as.data.frame(cbind(log2(counts[[j]][,2] + 1),
                                   log2(counts[[j]][,i] + 1)))
      plot.list <- list.append(plot.list,
                   ggplot(plot.df,aes_string(x=plot.df[,1],y=plot.df[,2]))+
                      stat_density_2d(aes(fill = ..level..), geom = "polygon")+
                      scale_fill_viridis()+
                      geom_abline(slope = 1, intercept = 0,color = "black", linetype = "dotted")+
                      annotate("text",x = 0, y = 14, label = lm_eqn(plot.df[,1],plot.df[,2]),parse = T,size = 2.5, hjust = 0)+
                      guides(fill = F)+
                      coord_cartesian(xlim=c(0,15),ylim=c(0,15))+
                      theme_bw()+
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank())
    )
  }
}

pdf(paste0(OUTPUT_DIR,"/actual_fadu10em_countsvcounts.pdf"),
    height=11,
    width=4)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/actual_fadu10em_countsvcounts.png"),
    height=11,
    width=4,
    units = "in",res=300)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

plot_grid(plotlist = plot.list,
          ncol=3)
```

![image](/images/actual_fadu10em_countsvcounts.png)

#### Plot FADU no_multimapped v other tools
```{R, fig.height=11, fig.width=4}
plot.list <- list()
for(i in 4:ncol(counts[[1]])){
  for(j in 1:length(counts)){
      plot.df <- as.data.frame(cbind(log2(counts[[j]][,3] + 1),
                                   log2(counts[[j]][,i] + 1)))
      plot.list <- list.append(plot.list,
                   ggplot(plot.df,aes_string(x=plot.df[,1],y=plot.df[,2]))+
                      stat_density_2d(aes(fill = ..level..), geom = "polygon")+
                      scale_fill_viridis()+
                      geom_abline(slope = 1, intercept = 0,color = "black", linetype = "dotted")+
                      annotate("text",x = 0, y = 14, label = lm_eqn(plot.df[,1],plot.df[,2]),parse = T,size = 2.5, hjust = 0)+
                      guides(fill = F)+
                      coord_cartesian(xlim=c(0,15),ylim=c(0,15))+
                      theme_bw()+
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank())
    )
  }
}

pdf(paste0(OUTPUT_DIR,"/actual_fadunomm_countsvcounts.pdf"),
    height=11,
    width=4)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/actual_fadunomm_countsvcounts.png"),
    height=11,
    width=4,
    units = "in",res=300)
plot_grid(plotlist = plot.list,
          ncol=3)
dev.off()

plot_grid(plotlist = plot.list,
          ncol=3)
```

![image](/images/actual_fadunomm_countsvcounts.png)

#### Create figure legend
```{R,fig.height=2, fig.width=4}
legend.plot <- ggplot(plot.df,aes_string(x=plot.df[,1],y=plot.df[,2]))+
                      stat_density_2d(aes(fill = ..level..), geom = "polygon")+
                      scale_fill_viridis()+
                      geom_abline(slope = 1, intercept = 0,color = "black", linetype = "dotted")+
                      annotate("text",x = 0, y = 14, label = lm_eqn(plot.df[,1],plot.df[,2]),parse = T,size = 2.5, hjust = 0)+
                      coord_cartesian(xlim=c(0,15),ylim=c(0,15))+
                      labs(fill="Level")+
                      theme_bw()+
                      theme(axis.title.x = element_blank(),
                            axis.title.y = element_blank())

level.legend <- g_legend(legend.plot)

pdf(paste0(OUTPUT_DIR,"/countsvcounts_key.pdf"),
    height=2,
    width=4)
grid.arrange(level.legend)
dev.off()

png(paste0(OUTPUT_DIR,"/countsvcounts_key.png"),
    height=2,
    width=4,
    units = "in",res=300)
grid.arrange(level.legend)
dev.off()

grid.arrange(level.legend)
```

![image](/images/countsvcounts_key.png")

### Conduct MA analysis

#### Select quantification methods for MA analysis
```{R}
colorder <- c("fadu_default","fadu_10em",
              "salmon_optimized","express_optimized","featurecounts_default",
              "htseq_union","htseq_unionnonuniqueall","featurecounts_fractionaloverlap")
```

#### Identify the methods over-/under-counted relative to default FADU
```{R}
ma.m.df <- as.data.frame(apply(counts[[3]],2,function(x){return(log2((x+1)/(counts[[3]]$fadu_default+1)))}))
ma.a.df <- as.data.frame(apply(counts[[3]],2,function(x){return(0.5*log2(x+1)+log2(counts[[3]]$fadu_default+1))}))

ma.m.subset.df <- ma.m.df[,which(colnames(ma.a.df) %in% colorder)]
ma.m.simplified.subset.df <- ma.m.subset.df[,2:ncol(ma.m.subset.df)]
ma.m.simplified.subset.df[abs(ma.m.simplified.subset.df) <= 2] <- 0
ma.m.simplified.subset.df[ma.m.simplified.subset.df > 2] <- 1
ma.m.simplified.subset.df[ma.m.simplified.subset.df < -2] <- -1

fadu_overcount_methods <- c()
fadu_undercount_methods <- c()
for(i in 1:nrow(ma.m.simplified.subset.df)){
  if(length(which(ma.m.simplified.subset.df[i,] == -1)) > 0){
    fadu_overcount_methods <- c(fadu_overcount_methods,paste(colnames(ma.m.simplified.subset.df)[which(ma.m.simplified.subset.df[i,] == -1)],collapse = " | "))
  }
  if(length(which(ma.m.simplified.subset.df[i,] == 1)) > 0){
    fadu_undercount_methods <- c(fadu_undercount_methods,paste(colnames(ma.m.simplified.subset.df)[which(ma.m.simplified.subset.df[i,] == 1)],collapse = " | "))
  }
}

rev(sort(table(strsplit(paste(fadu_overcount_methods,collapse = " | "), split = " [|] "))))
```

```{R, eval = F}
      express_optimized        salmon_optimized             htseq_union   featurecounts_default htseq_unionnonuniqueall 
                    122                      68                      67                      29                       2
```

```{R}
rev(sort(table(fadu_overcount_methods)))
```

```{R, eval = F}
                                      salmon_optimized | express_optimized 
                                                                        34 
                                                         express_optimized 
                                                                        33 
salmon_optimized | express_optimized | featurecounts_default | htseq_union 
                                                                        24 
                                           express_optimized | htseq_union 
                                                                        18 
                                                               htseq_union 
                                                                        10 
                        salmon_optimized | express_optimized | htseq_union 
                                                                         8 
                   express_optimized | featurecounts_default | htseq_union 
                                                                         3 
                 express_optimized | htseq_union | htseq_unionnonuniqueall 
                                                                         2 
                    salmon_optimized | featurecounts_default | htseq_union 
                                                                         1 
                                                          salmon_optimized 
                                                                         1 
                                       featurecounts_default | htseq_union 
                                                                         1 
```

```{R}
rev(sort(table(strsplit(paste(fadu_undercount_methods,collapse = " | "), split = " [|] "))))
```

```{R, eval = F}
        htseq_unionnonuniqueall featurecounts_fractionaloverlap                salmon_optimized                     htseq_union 
                             28                              13                               8                               5 
          featurecounts_default               express_optimized 
                              5                               2 
```


```{R}
rev(sort(table(fadu_undercount_methods)))
```

```{R, eval = F}
                                                                                           htseq_unionnonuniqueall 
                                                                                                                15 
                                                         htseq_unionnonuniqueall | featurecounts_fractionaloverlap 
                                                                                                                 8 
                                                                                                  salmon_optimized 
                                                                                                                 6 
                   featurecounts_default | htseq_union | htseq_unionnonuniqueall | featurecounts_fractionaloverlap 
                                                                                                                 3 
salmon_optimized | featurecounts_default | htseq_union | htseq_unionnonuniqueall | featurecounts_fractionaloverlap 
                                                                                                                 2 
                                                                                                 express_optimized 
                                                                                                                 2 
```

#### Identify the genes over-/under-counted relative to default FADU
```{R}
fadu_overcount_genes <- c()
fadu_undercount_genes <- c()

fadu_overcount_genes <- rownames(ma.m.simplified.subset.df)[ma.m.simplified.subset.df$salmon_optimized == -1 &
                                                            ma.m.simplified.subset.df$express_optimized == -1 &
                                                            ma.m.simplified.subset.df$featurecounts_default == -1 &
                                                            ma.m.simplified.subset.df$htseq_union == -1]

fadu_undercount_genes <- rownames(ma.m.simplified.subset.df)[ma.m.simplified.subset.df$htseq_unionnonuniqueall == 1 &
                                                             ma.m.simplified.subset.df$featurecounts_fractionaloverlap == 1]

fadu_overcount_genes
```

```{R, eval = F}
 [1] "gene100" "gene120" "gene167" "gene287" "gene322" "gene335" "gene344" "gene36"  "gene43"  "gene453" "gene477" "gene483"
[13] "gene496" "gene500" "gene545" "gene575" "gene604" "gene672" "gene761" "gene833" "gene835" "gene836" "gene936" "gene989"
```

```{R}
fadu_undercount_genes
```

```{R, eval = F}
 [1] "gene154" "gene23"  "gene335" "gene43"  "gene505" "gene615" "gene761" "gene771" "gene776" "gene779" "gene849" "gene892"
[13] "gene94"
```

```{R}
intersect(fadu_overcount_genes,fadu_undercount_genes)
```

```{R, eval = F}
[1] "gene335" "gene43"  "gene761"
```

#### Plot MA plots
```{R,fig.height=5,fig.width=7}
blue_genes <- c("gene833","gene835","gene836")
red_genes <- c("gene776")

blue_labels <- rownames(counts[[3]])
blue_labels[!(blue_labels %in% blue_genes)] <- ""
blue_genes_coords <- which(rownames(counts[[3]]) %in% blue_genes)

red_labels <- rownames(counts[[3]])
red_labels[!(red_labels %in% red_genes)] <- ""
red_genes_coords <- which(rownames(counts[[3]]) %in% red_genes)

maplot.list <- list()
for(i in 3:length(colorder)){
  plot.df <- as.data.frame(cbind(ma.m.df[,which(colnames(ma.m.df) == colorder[i])],
                                 ma.a.df[,which(colnames(ma.a.df) == colorder[i])]))
  maplot.list <- list.append(maplot.list,
                      ggplot(mapping = aes_string(x=plot.df[,2],y=plot.df[,1]))+
                        geom_hline(mapping=aes(yintercept=2),color = "darkorange", linetype = "dashed", size = 0.25)+
                        geom_hline(mapping=aes(yintercept=-2),color = "darkorange", linetype = "dashed", size = 0.25)+
                        geom_point(size = 0.1)+
                        geom_point(mapping = aes_string(x=plot.df[blue_genes_coords,2],y=plot.df[blue_genes_coords,1]),size= 1,color="blue")+
                        geom_point(mapping = aes_string(x=plot.df[red_genes_coords,2],y=plot.df[red_genes_coords,1]),size= 1,color="red")+
                        geom_text_repel(mapping = aes_string(label="red_labels"),size=3,color="red",ylim = c(5,NA))+
                        geom_text_repel(mapping = aes_string(label="blue_labels"),size=3,color="blue",ylim = c(NA,-5))+
                        #ggtitle(colorder[i])+
                        scale_y_continuous(limits=c(-10,10))+
                        labs(x="A",y="M")+
                        guides(fill = F)+
                        theme_bw()+
                        theme(axis.title.x = element_blank(),
                              axis.title.y = element_blank())
  )
                      
  maplot.list <- list.append(maplot.list,
                      ggplot(mapping = aes_string(x=plot.df[,2],y=plot.df[,1]))+
                        geom_hline(mapping=aes(yintercept=2),color = "darkorange", linetype = "dashed", size = 0.25)+
                        geom_hline(mapping=aes(yintercept=-2),color = "darkorange", linetype = "dashed", size = 0.25)+
                        stat_density_2d(aes(fill = ..level..), geom = "polygon")+
                        scale_fill_viridis()+
                        #ggtitle(colorder[i])+
                        scale_y_continuous(limits=c(-10,10))+
                        labs(x="A",y="M")+
                        guides(fill = F)+
                        theme_bw()+
                        theme(axis.title.x = element_blank(),
                              axis.title.y = element_blank())
  )
}

pdf(paste0(OUTPUT_DIR,"/maplots.pdf"),
    height=2,
    width=4)
plot_grid(plotlist = maplot.list,
          ncol=4)
dev.off()

png(paste0(OUTPUT_DIR,"/maplots.png"),
    height=2,
    width=4,
    units = "in",res=300)
plot_grid(plotlist = maplot.list,
          ncol=4)
dev.off()

plot_grid(plotlist = maplot.list,
          ncol=4)
```

![image](/images/maplots.png")

### Examine counts for genes within a specific operon

#### Set plot order 
```{R}
colorder <- c("salmon_optimized","salmon_default","kallisto_default","express_optimized","express_default",
              "featurecounts_default","fadu_nomm","htseq_intersectionnonempty","htseq_union","htseq_intersectionstrict",
              "htseq_unionnonuniqueall","featurecounts_overlap","featurecounts_fractionaloverlap",
              "fadu_10em","fadu_default")
for(i in 1:length(counts)){
  counts[[i]] <- counts[[i]][,match(colorder,colnames(counts[[i]]))]
}
```

#### Set operon genes and start and stop coordinates
```{R}
operon_genes <- c("gene826","gene827","gene828","gene829","gene830",
                  "gene831","gene832","gene833","gene834","gene835",
                  "gene836")
highlight_operon_genes <- c("gene833","gene835","gene836")

gff <- read.delim(GFF_PATH,comment.char="#",header=F)
gff <- gff[gff[,3] == "gene",]
gff[,10] <- gsub(".*;locus_tag=","",gff[,9])
gff[,10] <- gsub(";.*","",gff[,10])
gff[,11] <- gsub(".*;old_locus_tag=","",gff[,9])
gff[,11] <- gsub(";.*","",gff[,11])
gff[,9] <- gsub(";.*","",gff[,9])
gff[,9] <- gsub("ID=","",gff[,9])

gff.operon <- gff[gff[,9] %in% operon_genes,]

operon_start <- min(c(gff.operon[,4],gff.operon[,5]))
operon_stop <- max(c(gff.operon[,4],gff.operon[,5]))

operon_start <- min(c(gff.operon[,4],gff.operon[,5])) - 0.05*(operon_stop-operon_start)
operon_stop <- max(c(gff.operon[,4],gff.operon[,5])) + 0.05*(operon_stop-operon_start)

strand <- unique(gff.operon[,7])
```

#### Read and subset depth over operon region
```{R}
depth_f <- read.delim(DEPTH_F_PATH,header = F)
depth_r <- read.delim(DEPTH_R_PATH,header = F)

if(strand == "+"){
  depth <- depth_f[depth_f[,2] >= operon_start &
                   depth_f[,2] <= operon_stop,]
}else{
  depth <- depth_r[depth_r[,2] >= operon_start &
                   depth_r[,2] <= operon_stop,]
}
```

#### Calculate normalized relative count values
```{R}
counts.operon <- counts[[3]][rownames(counts[[3]]) %in% operon_genes,]
genelength.operon <- gff.operon[match(operon_genes,gff.operon[,9]),5] - gff.operon[match(operon_genes,gff.operon[,9]),4] + 1

relativecounts.operon <- counts.operon/genelength.operon
for(i in 1:ncol(relativecounts.operon)){
  relativecounts.operon[,i] <- relativecounts.operon[,i]/median(relativecounts.operon[,i])
}
log2relativecounts.operon <- log2(relativecounts.operon)
```

#### Plot depth over operon region
```{R,fig.height=2.5,fig.width=8}
gene_color <- rep("darkgrey",nrow(gff.operon))
gene_color[operon_genes %in% highlight_operon_genes] <- "turquoise2"
gene_label <- rep("",nrow(gff.operon))
gene_label[operon_genes %in% highlight_operon_genes] <- operon_genes[operon_genes %in% highlight_operon_genes]

operon_depth.plot <- ggplot()+
  geom_ribbon(mapping=aes(x=depth[,2],ymin=0,ymax=log10(depth[,3] + 1)),fill="orange")+
  geom_gene_arrow(mapping=aes(xmin=gff.operon[,4],xmax=gff.operon[,5],y=0,forward=gff.operon[,7]),fill=gene_color)+
  geom_text_repel(mapping=aes(x=(gff.operon[,5]+gff.operon[,4])/2,y=0,label=gene_label),ylim = c(0.5,NA),size=3)+
  scale_x_continuous(expand=c(0,0))+
  labs(x="wBm genome position",y="log2 read depth")+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/operon_depth_plot.pdf"),
    height=2.5,
    width=8)
print(operon_depth.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/operon_depth_plot.png"),
    height=2.5,
    width=8,
    units = "in",res=300)
print(operon_depth.plot)
dev.off()

print(operon_depth.plot)
```

![image](/images/operon_depth_plot.png")

#### Plot normalized relative count values of operon region
```{R,fig.height=5,fig.width=8}
log2relativecounts.operon.plot.df <- log2relativecounts.operon
log2relativecounts.operon.plot.df$gene <- rownames(log2relativecounts.operon) 
log2relativecounts.operon.plot.df <- melt(log2relativecounts.operon.plot.df)

log2relativecounts.operon.plot.df$variable <- formalize_names(log2relativecounts.operon.plot.df$variable)
log2relativecounts.operon.plot.df$variable <- factor(log2relativecounts.operon.plot.df$variable,levels=rev(formalize_names(colnames(counts[[3]]))))

log2relativecounts.operon.plot <- ggplot(mapping=aes(x=log2relativecounts.operon.plot.df$gene,
                          y=log2relativecounts.operon.plot.df$variable,
                          fill=log2relativecounts.operon.plot.df$value,
                          label=round(log2relativecounts.operon.plot.df$value,1)))+
  geom_raster()+
  geom_text()+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_gradient2(low = "navyblue", mid = "white", high = "firebrick3",limits=c(-3,3))+
  guides(fill = F)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust=1))

pdf(paste0(OUTPUT_DIR,"/operon_log2relativecounts_plot.pdf"),
    height=2.5,
    width=8)
print(log2relativecounts.operon.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/operon_log2relativecounts_plot.png"),
    height=2.5,
    width=8,
    units = "in",res=300)
print(log2relativecounts.operon.plot)
dev.off()

print(log2relativecounts.operon.plot)
```

![image](/images/operon_log2relativecounts_plot.png")

#### Create figure legend
```{R}
legend.plot <- ggplot(mapping=aes(x=log2relativecounts.operon.plot.df$gene,
                                  y=log2relativecounts.operon.plot.df$variable,
                                  fill=log2relativecounts.operon.plot.df$value,
                                  label=round(log2relativecounts.operon.plot.df$value,1)))+
  geom_raster()+
  geom_text()+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  scale_fill_gradient2(low = "navyblue", mid = "white", high = "firebrick3",limits=c(-3,3))+
  labs(fill = "log2 normalized relative counts")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust=1),
        legend.position = "bottom")

log2relativecounts.legend <- g_legend(legend.plot)

pdf(paste0(OUTPUT_DIR,"/log2relativecounts_key.pdf"),
    height=2,
    width=4)
grid.arrange(log2relativecounts.legend)
dev.off()

png(paste0(OUTPUT_DIR,"/log2relativecounts_key.png"),
    height=2,
    width=4,
    units = "in",res=300)
grid.arrange(log2relativecounts.legend)
dev.off()

grid.arrange(log2relativecounts.legend)
```

![image](/images/log2relativecounts_key.png")

#### Combine plots
```{R,fig.height=8,fig.width=8}
pdf(paste0(OUTPUT_DIR,"/operon_full_plot.pdf"),
    height=8,
    width=8)
plot_grid(plotlist = list(operon_depth.plot,
                          log2relativecounts.operon.plot),
          #align = 'v',
          rel_heights = c(1,3),
          labels="AUTO",
          ncol=1)
dev.off()

png(paste0(OUTPUT_DIR,"/operon_full_plot.png"),
    height=8,
    width=8,
    units = "in",res=300)
plot_grid(plotlist = list(operon_depth.plot,
                          log2relativecounts.operon.plot),
          #align = 'v',
          rel_heights = c(1,3),
          labels="AUTO",
          ncol=1)
dev.off()

plot_grid(plotlist = list(operon_depth.plot,
                          log2relativecounts.operon.plot),
          #align = 'v',
          rel_heights = c(1,3),
          labels="AUTO",
          ncol=1)
```

![image](/images/operon_full_plot.png")

### Examine counts for genes under-counted by FADU

#### Set under-counted genes and start and stop coordinates
```{R}
undercount_genes <- c("gene776","gene777","gene778")
highlight_undercount_genes <- c("gene776")

gff.undercount <- gff[gff[,9] %in% undercount_genes,]

undercount_start <- min(c(gff.undercount[,4],gff.undercount[,5]))
undercount_stop <- max(c(gff.undercount[,4],gff.undercount[,5]))

undercount_start <- min(c(gff.undercount[,4],gff.undercount[,5])) - 0.1*(undercount_stop-undercount_start)
undercount_stop <- max(c(gff.undercount[,4],gff.undercount[,5])) + 0.1*(undercount_stop-undercount_start)

strand <- unique(gff.undercount[,7])
```

#### Read and subset depth over under-counted gene region
```{R}
if(strand == "+"){
  depth <- depth_f[depth_f[,2] >= undercount_start &
                   depth_f[,2] <= undercount_stop,]
}else{
  depth <- depth_r[depth_r[,2] >= undercount_start &
                   depth_r[,2] <= undercount_stop,]
}
```

#### Plot depth over under-counted gene region
```{R,fig.height=2.5,fig.width=8}
gene_color <- rep("darkgrey",nrow(gff.undercount))
gene_color[undercount_genes %in% highlight_undercount_genes] <- "red"
gene_label <- rep("",nrow(gff.undercount))
gene_label[undercount_genes %in% highlight_undercount_genes] <- undercount_genes[undercount_genes %in% highlight_undercount_genes]

undercount_depth.plot <- ggplot()+
  geom_ribbon(mapping=aes(x=depth[,2],ymin=0,ymax=log10(depth[,3] + 1)),fill="orange")+
  geom_gene_arrow(mapping=aes(xmin=gff.undercount[,4],xmax=gff.undercount[,5],y=0,forward=gff.undercount[,7]),fill=gene_color)+
  geom_text_repel(mapping=aes(x=(gff.undercount[,5]+gff.undercount[,4])/2,y=0,label=gene_label),ylim = c(0.5,NA),size=3)+
  scale_x_continuous(expand=c(0,0))+
  labs(x="wBm genome position",y="log2 read depth")+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/undercountgenes_depth_plot.pdf"),
    height=2.5,
    width=8)
print(undercount_depth.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/undercountgenes_depth_plot.png"),
    height=2.5,
    width=8,
    units = "in",res=300)
print(undercount_depth.plot)
dev.off()

print(undercount_depth.plot)
```

![image](/images/undercountgenes_depth_plot.png")

#### Plot count values of operon region
```{R,fig.height=5,fig.width=6}
counts.undercount <- counts[[3]][rownames(counts[[3]]) %in% undercount_genes,]
counts.undercount$genes <- rownames(counts.undercount)

counts.undercount.plot.df <- counts.undercount
counts.undercount.plot.df$gene <- rownames(counts.undercount) 
counts.undercount.plot.df <- melt(counts.undercount.plot.df)

counts.undercount.plot.df$variable <- formalize_names(counts.undercount.plot.df$variable)
counts.undercount.plot.df$variable <- factor(counts.undercount.plot.df$variable,levels=rev(formalize_names(colnames(counts[[3]]))))

counts.undercount.plot <- ggplot(mapping=aes(x=counts.undercount.plot.df$gene,
                          y=counts.undercount.plot.df$variable,
                          label=round(counts.undercount.plot.df$value,1)))+
  geom_raster(fill="white")+
  geom_text()+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  guides(fill = F)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45,hjust=1))

pdf(paste0(OUTPUT_DIR,"/undercountgenes_counts_plot.pdf"),
    height=5,
    width=6)
print(counts.undercount.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/undercountgenes_counts_plot.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(counts.undercount.plot)
dev.off()

print(counts.undercount.plot)
```

![image](/images/undercountgenes_counts_plot.png)

# Generate scripts for benchmarking 

## Create benchmarking scripts directory

```{bash, eval = F}
mkdir "$WORKING_DIR"/benchmarking_scripts
mkdir "$WORKING_DIR"/benchmarking_scripts/outputs/
```

## Create benchmarking scripts

##### Inputs
```{bash, eval = F}
FEAT_TYPE=gene
ATTR_ID=ID

## 1 Core
THREADS=1

## 4 Cores
THREADS=4

### Simulation
OUTPUT_PREFIX=simulation
BAM_DIR="$WORKING_DIR"/bam
READS_DIR="$WORKING_DIR"/simulation/reads
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.gff
STRANDED=yes

### E. chaffeensis
OUTPUT_PREFIX=echaffeensis
BAM_DIR="$WORKING_DIR"/echaffeensis
READS_DIR="$WORKING_DIR"/echaffeensis
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000013145.1_ASM1314v1_genomic.fna
GFF3="$REFERENCES_DIR"/canine_echaf.gff
STRANDED=no

### E. coli
OUTPUT_PREFIX=ecoli
BAM_DIR="$WORKING_DIR"/ecoli
READS_DIR="$WORKING_DIR"/ecoli
NUC_GENE_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.fna
GFF3="$REFERENCES_DIR"/GCF_000005845.2_ASM584v2_genomic.gene.gff
STRANDED=no

### wBm
OUTPUT_PREFIX=wbm
BAM_DIR="$WORKING_DIR"/wbm
READS_DIR="$WORKING_DIR"/wbm
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wbm.gene.fna
NUC_GENOME_FNA="$REFERENCES_DIR"/bmalayi_wbm.fna
GFF3="$REFERENCES_DIR"/bmalayi_wbm.gff
STRANDED=reverse
```

##### Commands
```{bash, eval = F}
for SAMPLE in $(find "$BAM_DIR" -name "*[.]bam" | sed "s/.*\\///g" | sort | sed "s/[.].*//g" | uniq)
do
echo -e ""$HISAT2_BIN_DIR"/hisat2-build -p "$THREADS" --large-index "$NUC_GENOME_FNA" "$NUC_GENOME_FNA"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_hisat2_genome_index_pe"$THREADS".sh
echo -e ""$HISAT2_BIN_DIR"/hisat2-build -p "$THREADS" --large-index "$NUC_GENE_FNA" "$NUC_GENE_FNA"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_hisat2_gene_index_pe"$THREADS".sh

echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENOME_FNA" -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".genome.bam -" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_hisat2_genome_align_pe"$THREADS".sh
echo -e ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -k 200 -X 1000 -x "$NUC_GENE_FNA" -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz --no-spliced-alignment --no-discordant | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".transcript.bam -" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_hisat2_gene_align_pe"$THREADS".sh

echo -e ""$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SAMPLE".genome.bam -@ "$THREADS" -o "$WORKING_DIR"/benchmarking_scripts/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_samtools_genome_sort_pe"$THREADS".sh
echo -e ""$SAMTOOLS_BIN_DIR"/samtools sort "$BAM_DIR"/"$SAMPLE".transcript.bam -@ "$THREADS" -n -o "$WORKING_DIR"/benchmarking_scripts/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_samtools_gene_sort_pe"$THREADS".sh
echo -e ""$SAMTOOLS_BIN_DIR"/samtools index "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -@ "$THREADS"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_samtools_genome_index_pe"$THREADS".sh

echo -e ""$KALLISTO_BIN_DIR"/kallisto index -i "$NUC_GENE_FNA".kallisto.index --make-unique "$NUC_GENE_FNA"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_kallisto_index_pe"$THREADS".sh
echo -e ""$SALMON_BIN_DIR"/salmon index -p "$THREADS" --keepDuplicates -t "$NUC_GENE_FNA" -i "$NUC_GENE_FNA".salmon.index" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_salmon_index_pe"$THREADS".sh

echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".salmon_default.counts --validateMappings" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_salmon_default_pe"$THREADS".sh
echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SAMPLE"_1.fastq.gz -2 "$READS_DIR"/"$SAMPLE"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".salmon_optimized.counts --validateMappings  --allowDovetail"  > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_salmon_optimized_pe"$THREADS".sh

if [ "$STRANDED" == reverse ]
then

echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_default_pe"$THREADS".sh
echo -e ""$EXPRESS_BIN_DIR"/express --rf-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_optimized_pe"$THREADS".sh

echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_default.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_default_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_10em.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_10em_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_nomm.counts -s reverse -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_nomm_pe"$THREADS".sh

echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_union.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_union_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionstrict.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionstrict_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionnonempty.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionnonempty_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s reverse "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_unionnonuniqueall.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_unionnonuniqueall_pe"$THREADS".sh

echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_default_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_overlap_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 2 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_fractionaloverlap_pe"$THREADS".sh

echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --rf-stranded -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_kallisto_default_pe"$THREADS".sh
elif [ "$STRANDED" == yes ]
then

echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_default_pe"$THREADS".sh
echo -e ""$EXPRESS_BIN_DIR"/express --fr-stranded -B10 --no-bias-correct -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_optimized_pe"$THREADS".sh

echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_default.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_default_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_10em.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_10em_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_nomm.counts -s yes -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_nomm_pe"$THREADS".sh

echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_union.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_union_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionstrict.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionstrict_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionnonempty.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionnonempty_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s yes "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_unionnonuniqueall.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_unionnonuniqueall_pe"$THREADS".sh

echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_default_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_overlap_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 1 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_fractionaloverlap_pe"$THREADS".sh

echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" --fr-stranded -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_kallisto_default_pe"$THREADS".sh

else

echo -e ""$EXPRESS_BIN_DIR"/express -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_default.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_default_pe"$THREADS".sh
echo -e ""$EXPRESS_BIN_DIR"/express -B10 --no-bias-correct -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".express_optimized.counts "$NUC_GENE_FNA" "$BAM_DIR"/"$SAMPLE".transcript.sortedbyname.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_express_optimized_pe"$THREADS".sh

echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_default.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID"" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_default_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_10em.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --em_iterations 10" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_10em_pe"$THREADS".sh
echo -e "export JULIA_DEPOT_PATH="$JULIA_LIB_DIR"\nexport JULIA_NUM_THREADS="$THREADS"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -b "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam -g "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".fadu_nomm.counts -s no -f "$FEAT_TYPE" -a "$ATTR_ID" --remove_multimapped" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_fadu_nomm_pe"$THREADS".sh

echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_union.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_union_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-strict --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionstrict.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionstrict_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m intersection-nonempty --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_intersectionnonempty.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_intersectionnonempty_pe"$THREADS".sh
echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_DIR":"$LD_LIBRARY_PATH"\n"$PYTHON_BIN_DIR"/python -m HTSeq.scripts.count -m union --nonunique all --order pos -f bam -t "$FEAT_TYPE" -i "$ATTR_ID" -s no "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam "$GFF3" > "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".htseq_unionnonuniqueall.counts" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_htseq_unionnonuniqueall_pe"$THREADS".sh

echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_default.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_default_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_overlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_overlap_pe"$THREADS".sh
echo -e ""$SUBREAD_BIN_DIR"/featureCounts -T "$THREADS" -p -O --fraction -a "$GFF3" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".featurecounts_fractionaloverlap.counts -t "$FEAT_TYPE" -g "$ATTR_ID" -s 0 "$BAM_DIR"/"$SAMPLE".genome.sortedbyposition.bam" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_featurecounts_fractionaloverlap_pe"$THREADS".sh

echo -e ""$KALLISTO_BIN_DIR"/kallisto quant -i "$NUC_GENE_FNA".kallisto.index -t "$THREADS" -o "$WORKING_DIR"/benchmarking_scripts/outputs/"$SAMPLE".kallisto_default.counts "$READS_DIR"/"$SAMPLE"_1.fastq.gz "$READS_DIR"/"$SAMPLE"_2.fastq.gz" > "$WORKING_DIR"/benchmarking_scripts/"$OUTPUT_PREFIX"_"$SAMPLE"_kallisto_default_pe"$THREADS".sh

fi
done
```

## Remove extra scripts that do not need to be benchmarked

##### Commands
```{bash, eval = F}
rm "$WORKING_DIR"/benchmarking_scripts/*pseudo*
rm "$WORKING_DIR"/benchmarking_scripts/*sample_02*
rm "$WORKING_DIR"/benchmarking_scripts/*sample_03*
rm "$WORKING_DIR"/benchmarking_scripts/*sample_04*
rm "$WORKING_DIR"/benchmarking_scripts/*SRR5189260*
```

## Download Shaun's benchmarking stats

Stats can be found in this Google Sheet: https://docs.google.com/spreadsheets/d/1II4kHQKEsll27-uI4bZC8q3qU0a1Mk4iSAeXfgkpzj4/edit#gid=0

## Plot benchmarking data

### Set R inputs
```{R}
BENCHMARKING_INFO.PATH <- "C:/Users/MChung.SOM/Documents/benchmarking_info.txt"

OUTPUT_DIR <- "C:/Users/MChung.SOM/Documents/plots"
```

### Load R packages and view sessionInfo

```{R}
library(cowplot)
library(ggplot2)

sessionInfo()
```

```{R, eval = F} 
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cowplot_0.9.4 ggplot2_3.2.1

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1       withr_2.1.2      assertthat_0.2.1 crayon_1.3.4     dplyr_0.8.3      grid_3.5.1       R6_2.4.0        
 [8] gtable_0.3.0     magrittr_1.5     scales_1.0.0     pillar_1.3.1     rlang_0.4.0      lazyeval_0.2.2   rstudioapi_0.10 
[15] tools_3.5.1      glue_1.3.1       purrr_0.3.2      munsell_0.5.0    xfun_0.6         yaml_2.2.0       compiler_3.5.1  
[22] pkgconfig_2.0.2  colorspace_1.4-1 knitr_1.22       tidyselect_0.2.5 tibble_2.1.1    
```

### Load R functions
```{R}
formalize_names <- function(vector){
	vector <- gsub("express_default","eXpress",vector)
	vector <- gsub("express_optimized","eXpress\n-B10 --no-bias-correct",vector)

	vector <- gsub("fadu_default","FADU",vector)
	vector <- gsub("fadu_10em","FADU\n--em_iterations 10",vector)
	vector <- gsub("fadu_nomm","FADU\n--remove_multimapped",vector)

	vector <- gsub("htseq_unionnonuniqueall","HTSeq\n-m union --nonunique all",vector)
	vector <- gsub("htseq_union","HTSeq\n-m union",vector)
	vector <- gsub("htseq_intersectionstrict","HTSeq\n-m intersection-strict",vector)
	vector <- gsub("htseq_intersectionnonempty","HTSeq\n-m intersection-nonempty",vector)
	
	vector <- gsub("featurecounts_default","featureCounts",vector)
	vector <- gsub("featurecounts_overlap","featureCounts\n-O",vector)
	vector <- gsub("featurecounts_fractionaloverlap","featureCounts\n-O --fraction",vector)

	vector <- gsub("salmon_default","Salmon\n--validateMappings",vector)
	vector <- gsub("salmon_optimized","Salmon\n--validateMappings --allowDovetail",vector)

	vector <- gsub("kallisto_default","kallisto",vector)
}

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 
```

```{R}
benchmarking_info <- read.delim(BENCHMARKING_INFO.PATH)
benchmarking_info <- benchmarking_info[!(is.na(benchmarking_info[,4])),]
benchmarking_info[,2] <- gsub("[.]",":",benchmarking_info[,2])
```

```{R}
samples <- c("echaffeensis","ecoli","wbm")
```

### Convert wallclock time units to seconds
```{R}
benchmarking_info[,6] <- 0

for(i in 1:nrow(benchmarking_info)){
  hours <- as.numeric(as.character(unlist(strsplit(as.character(benchmarking_info[i,2]),split=":"))[1]))
  min <- as.numeric(as.character(unlist(strsplit(as.character(benchmarking_info[i,2]),split=":"))[2]))
  sec <- as.numeric(as.character(unlist(strsplit(as.character(benchmarking_info[i,2]),split=":"))[3]))

  benchmarking_info[i,6] <- hours*360 + min*60 + sec
}

```

### Parse benchmarking data to create timing data frames for each dataset
```{R}
benchmarking_timingplotdata_pe1.list <- list()
benchmarking_timingplotdata_pe4.list <- list()

for(i in 1:length(samples)){
  benchmarking_info.subset <- benchmarking_info[grep(samples[i],benchmarking_info[,1]),]
  genome_alignindex <- mean(benchmarking_info.subset[grep("hisat2_genome_align_pe1",benchmarking_info.subset[,1]),6]) + mean(benchmarking_info.subset[grep("samtools_genome_sort_pe1",benchmarking_info.subset[,1]),6]) + mean(benchmarking_info.subset[grep("samtools_genome_index_pe1",benchmarking_info.subset[,1]),6])
  gene_alignindex <- mean(benchmarking_info.subset[grep("hisat2_gene_align_pe1",benchmarking_info.subset[,1]),6]) + mean(benchmarking_info.subset[grep("samtools_gene_sort_pe1",benchmarking_info.subset[,1]),6]) + mean(benchmarking_info.subset[grep("hisat2_gene_index_pe1",benchmarking_info.subset[,1]),6])
  kallisto_alignindex <- mean(benchmarking_info.subset[grep("kallisto_index_pe1",benchmarking_info.subset[,1]),6])
  salmon_alignindex <- mean(benchmarking_info.subset[grep("salmon_index_pe1",benchmarking_info.subset[,1]),6])
  
  benchmarking_timingplotdata_pe1.list[[i]] <- as.data.frame(rbind(
    c("salmon_optimized","align",salmon_alignindex),
    c("salmon_optimized","quant",mean(benchmarking_info.subset[grep("salmon_optimized_pe1",benchmarking_info.subset[,1]),6])),
    c("salmon_default","align",salmon_alignindex),
    c("salmon_default","quant",mean(benchmarking_info.subset[grep("salmon_default_pe1",benchmarking_info.subset[,1]),6])),
    c("kallisto_default","align",kallisto_alignindex),
    c("kallisto_default","quant",mean(benchmarking_info.subset[grep("kallisto_default_pe1",benchmarking_info.subset[,1]),6])),
    c("express_optimized","align",gene_alignindex),
    c("express_optimized","quant",mean(benchmarking_info.subset[grep("express_optimized_pe1",benchmarking_info.subset[,1]),6])),
    c("express_default","align",gene_alignindex),
    c("express_default","quant",mean(benchmarking_info.subset[grep("express_default_pe1",benchmarking_info.subset[,1]),6])),
    c("featurecounts_default","align",genome_alignindex),
    c("featurecounts_default","quant",mean(benchmarking_info.subset[grep("featurecounts_default_pe1",benchmarking_info.subset[,1]),6])),  
    c("fadu_nomm","align",genome_alignindex),
    c("fadu_nomm","quant",mean(benchmarking_info.subset[grep("fadu_nomm_pe1",benchmarking_info.subset[,1]),6])),  
    c("htseq_intersectionnonempty","align",genome_alignindex),
    c("htseq_intersectionnonempty","quant",mean(benchmarking_info.subset[grep("htseq_intersectionnonempty_pe1",benchmarking_info.subset[,1]),6])),
    c("htseq_union","align",genome_alignindex),
    c("htseq_union","quant",mean(benchmarking_info.subset[grep("htseq_union_pe1",benchmarking_info.subset[,1]),6])),
    c("htseq_intersectionstrict","align",genome_alignindex),
    c("htseq_intersectionstrict","quant",mean(benchmarking_info.subset[grep("htseq_intersectionstrict_pe1",benchmarking_info.subset[,1]),6])),
    c("htseq_unionnonuniqueall","align",genome_alignindex),
    c("htseq_unionnonuniqueall","quant",mean(benchmarking_info.subset[grep("htseq_unionnonuniqueall_pe1",benchmarking_info.subset[,1]),6])),
    c("featurecounts_overlap","align",genome_alignindex),
    c("featurecounts_overlap","quant",mean(benchmarking_info.subset[grep("featurecounts_overlap_pe1",benchmarking_info.subset[,1]),6])),  
    c("featurecounts_fractionaloverlap","align",genome_alignindex),
    c("featurecounts_fractionaloverlap","quant",mean(benchmarking_info.subset[grep("featurecounts_fractionaloverlap_pe1",benchmarking_info.subset[,1]),6])),  
    c("fadu_10em","align",genome_alignindex),
    c("fadu_10em","quant",mean(benchmarking_info.subset[grep("fadu_10em_pe1",benchmarking_info.subset[,1]),6])),  
    c("fadu_default","align",genome_alignindex),
    c("fadu_default","quant",mean(benchmarking_info.subset[grep("fadu_default_pe1",benchmarking_info.subset[,1]),6]))
  ))
  
  genome_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_genome_align_pe4",benchmarking_info.subset[,1]),6], 
                             benchmarking_info.subset[grep("samtools_genome_sort_pe4",benchmarking_info.subset[,1]),6],
                             benchmarking_info.subset[grep("samtools_genome_index_pe4",benchmarking_info.subset[,1]),6]))
  gene_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_gene_align_pe4",benchmarking_info.subset[,1]),6], 
                           benchmarking_info.subset[grep("samtools_gene_sort_pe4",benchmarking_info.subset[,1]),6],
                           benchmarking_info.subset[grep("hisat2_gene_index_pe4",benchmarking_info.subset[,1]),6]))
  kallisto_alignindex <- max(benchmarking_info.subset[grep("kallisto_index_pe4",benchmarking_info.subset[,1]),6])
  salmon_alignindex <- max(benchmarking_info.subset[grep("salmon_index_pe4",benchmarking_info.subset[,1]),6])
  
  benchmarking_timingplotdata_pe4.list[[i]] <- as.data.frame(rbind(
    c("salmon_optimized","align",salmon_alignindex),
    c("salmon_optimized","quant",mean(benchmarking_info.subset[grep("salmon_optimized_pe4",benchmarking_info.subset[,1]),6])),
    c("salmon_default","align",salmon_alignindex),
    c("salmon_default","quant",mean(benchmarking_info.subset[grep("salmon_default_pe4",benchmarking_info.subset[,1]),6])),
    c("kallisto_default","align",kallisto_alignindex),
    c("kallisto_default","quant",mean(benchmarking_info.subset[grep("kallisto_default_pe4",benchmarking_info.subset[,1]),6])),
    c("express_optimized","align",gene_alignindex),
    c("express_optimized","quant",mean(benchmarking_info.subset[grep("express_optimized_pe4",benchmarking_info.subset[,1]),6])),
    c("express_default","align",gene_alignindex),
    c("express_default","quant",mean(benchmarking_info.subset[grep("express_default_pe4",benchmarking_info.subset[,1]),6])),
    c("featurecounts_default","align",genome_alignindex),
    c("featurecounts_default","quant",mean(benchmarking_info.subset[grep("featurecounts_default_pe4",benchmarking_info.subset[,1]),6])),  
    c("fadu_nomm","align",genome_alignindex),
    c("fadu_nomm","quant",mean(benchmarking_info.subset[grep("fadu_nomm_pe4",benchmarking_info.subset[,1]),6])),  
    c("htseq_intersectionnonempty","align",genome_alignindex),
    c("htseq_intersectionnonempty","quant",mean(benchmarking_info.subset[grep("htseq_intersectionnonempty_pe4",benchmarking_info.subset[,1]),6])),
    c("htseq_union","align",genome_alignindex),
    c("htseq_union","quant",mean(benchmarking_info.subset[grep("htseq_union_pe4",benchmarking_info.subset[,1]),6])),
    c("htseq_intersectionstrict","align",genome_alignindex),
    c("htseq_intersectionstrict","quant",mean(benchmarking_info.subset[grep("htseq_intersectionstrict_pe4",benchmarking_info.subset[,1]),6])),
    c("htseq_unionnonuniqueall","align",genome_alignindex),
    c("htseq_unionnonuniqueall","quant",mean(benchmarking_info.subset[grep("htseq_unionnonuniqueall_pe4",benchmarking_info.subset[,1]),6])),
    c("featurecounts_overlap","align",genome_alignindex),
    c("featurecounts_overlap","quant",mean(benchmarking_info.subset[grep("featurecounts_overlap_pe4",benchmarking_info.subset[,1]),6])),  
    c("featurecounts_fractionaloverlap","align",genome_alignindex),
    c("featurecounts_fractionaloverlap","quant",mean(benchmarking_info.subset[grep("featurecounts_fractionaloverlap_pe4",benchmarking_info.subset[,1]),6])),  
    c("fadu_10em","align",genome_alignindex),
    c("fadu_10em","quant",mean(benchmarking_info.subset[grep("fadu_10em_pe4",benchmarking_info.subset[,1]),6])),  
    c("fadu_default","align",genome_alignindex),
    c("fadu_default","quant",mean(benchmarking_info.subset[grep("fadu_default_pe4",benchmarking_info.subset[,1]),6]))
  ))
}
```

### Parse benchmarking data to create VMem data frames for each dataset
```{R}
benchmarking_vmemplotdata_pe1.list <- list()
benchmarking_vmemplotdata_pe4.list <- list()

for(i in 1:length(samples)){
  benchmarking_info.subset <- benchmarking_info[grep(samples[i],benchmarking_info[,1]),]
  genome_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_genome_align_pe1",benchmarking_info.subset[,1]),4], 
                             benchmarking_info.subset[grep("samtools_genome_sort_pe1",benchmarking_info.subset[,1]),4],
                             benchmarking_info.subset[grep("samtools_genome_index_pe1",benchmarking_info.subset[,1]),4]))
  gene_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_gene_align_pe1",benchmarking_info.subset[,1]),4], 
                           benchmarking_info.subset[grep("samtools_gene_sort_pe1",benchmarking_info.subset[,1]),4],
                           benchmarking_info.subset[grep("hisat2_gene_index_pe1",benchmarking_info.subset[,1]),4]))
  kallisto_alignindex <- max(benchmarking_info.subset[grep("kallisto_index_pe1",benchmarking_info.subset[,1]),4])
  salmon_alignindex <- max(benchmarking_info.subset[grep("salmon_index_pe1",benchmarking_info.subset[,1]),4])
  
  benchmarking_vmemplotdata_pe1.list[[i]] <- as.data.frame(rbind(
    c("salmon_optimized",max(salmon_alignindex,benchmarking_info.subset[grep("salmon_optimized_pe1",benchmarking_info.subset[,1]),4])),
    c("salmon_default",max(salmon_alignindex,benchmarking_info.subset[grep("salmon_default_pe1",benchmarking_info.subset[,1]),4])),
    c("kallisto_default",max(kallisto_alignindex,benchmarking_info.subset[grep("kallisto_default_pe1",benchmarking_info.subset[,1]),4])),
    c("express_optimized",max(gene_alignindex,benchmarking_info.subset[grep("express_optimized_pe1",benchmarking_info.subset[,1]),4])),
    c("express_default",max(gene_alignindex,benchmarking_info.subset[grep("express_default_pe1",benchmarking_info.subset[,1]),4])),
    c("featurecounts_default",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_default_pe1",benchmarking_info.subset[,1]),4])),
    c("fadu_nomm",max(genome_alignindex,benchmarking_info.subset[grep("fadu_nomm_pe1",benchmarking_info.subset[,1]),4])),
    c("htseq_intersectionnonempty",max(genome_alignindex,benchmarking_info.subset[grep("htseq_intersectionnonempty_pe1",benchmarking_info.subset[,1]),4])),
    c("htseq_union",max(genome_alignindex,benchmarking_info.subset[grep("htseq_union_pe1",benchmarking_info.subset[,1]),4])),
    c("htseq_intersectionstrict",max(genome_alignindex,benchmarking_info.subset[grep("htseq_intersectionstrict_pe1",benchmarking_info.subset[,1]),4])),
    c("htseq_unionnonuniqueall",max(genome_alignindex,benchmarking_info.subset[grep("htseq_unionnonuniqueall_pe1",benchmarking_info.subset[,1]),4])),
    c("featurecounts_overlap",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_overlap_pe1",benchmarking_info.subset[,1]),4])),
    c("featurecounts_fractionaloverlap",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_fractionaloverlap_pe1",benchmarking_info.subset[,1]),4])),
    c("fadu_10em",max(genome_alignindex,benchmarking_info.subset[grep("fadu_10em_pe1",benchmarking_info.subset[,1]),4])),
    c("fadu_default",max(genome_alignindex,benchmarking_info.subset[grep("fadu_default_pe1",benchmarking_info.subset[,1]),4]))
  ))
  
  genome_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_genome_align_pe4",benchmarking_info.subset[,1]),4], 
                             benchmarking_info.subset[grep("samtools_genome_sort_pe4",benchmarking_info.subset[,1]),4],
                             benchmarking_info.subset[grep("samtools_genome_index_pe4",benchmarking_info.subset[,1]),4]))
  gene_alignindex <- max(c(benchmarking_info.subset[grep("hisat2_gene_align_pe4",benchmarking_info.subset[,1]),4], 
                           benchmarking_info.subset[grep("samtools_gene_sort_pe4",benchmarking_info.subset[,1]),4],
                           benchmarking_info.subset[grep("hisat2_gene_index_pe4",benchmarking_info.subset[,1]),4]))
  kallisto_alignindex <- max(benchmarking_info.subset[grep("kallisto_index_pe4",benchmarking_info.subset[,1]),4])
  salmon_alignindex <- max(benchmarking_info.subset[grep("salmon_index_pe4",benchmarking_info.subset[,1]),4])
  
  benchmarking_vmemplotdata_pe4.list[[i]] <- as.data.frame(rbind(
    c("salmon_optimized",max(salmon_alignindex,benchmarking_info.subset[grep("salmon_optimized_pe4",benchmarking_info.subset[,1]),4])),
    c("salmon_default",max(salmon_alignindex,benchmarking_info.subset[grep("salmon_default_pe4",benchmarking_info.subset[,1]),4])),
    c("kallisto_default",max(kallisto_alignindex,benchmarking_info.subset[grep("kallisto_default_pe4",benchmarking_info.subset[,1]),4])),
    c("express_optimized",max(gene_alignindex,benchmarking_info.subset[grep("express_optimized_pe4",benchmarking_info.subset[,1]),4])),
    c("express_default",max(gene_alignindex,benchmarking_info.subset[grep("express_default_pe4",benchmarking_info.subset[,1]),4])),
    c("featurecounts_default",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_default_pe4",benchmarking_info.subset[,1]),4])),
    c("fadu_nomm",max(genome_alignindex,benchmarking_info.subset[grep("fadu_nomm_pe4",benchmarking_info.subset[,1]),4])),
    c("htseq_intersectionnonempty",max(genome_alignindex,benchmarking_info.subset[grep("htseq_intersectionnonempty_pe4",benchmarking_info.subset[,1]),4])),
    c("htseq_union",max(genome_alignindex,benchmarking_info.subset[grep("htseq_union_pe4",benchmarking_info.subset[,1]),4])),
    c("htseq_intersectionstrict",max(genome_alignindex,benchmarking_info.subset[grep("htseq_intersectionstrict_pe4",benchmarking_info.subset[,1]),4])),
    c("htseq_unionnonuniqueall",max(genome_alignindex,benchmarking_info.subset[grep("htseq_unionnonuniqueall_pe4",benchmarking_info.subset[,1]),4])),
    c("featurecounts_overlap",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_overlap_pe4",benchmarking_info.subset[,1]),4])),
    c("featurecounts_fractionaloverlap",max(genome_alignindex,benchmarking_info.subset[grep("featurecounts_fractionaloverlap_pe4",benchmarking_info.subset[,1]),4])),
    c("fadu_10em",max(genome_alignindex,benchmarking_info.subset[grep("fadu_10em_pe4",benchmarking_info.subset[,1]),4])),
    c("fadu_default",max(genome_alignindex,benchmarking_info.subset[grep("fadu_default_pe4",benchmarking_info.subset[,1]),4]))
  ))
}
```

### Set plot order
```{R}
colorder <- c("salmon_optimized","salmon_default","kallisto_default","express_optimized","express_default",
              "featurecounts_default","fadu_nomm","htseq_intersectionnonempty","htseq_union","htseq_intersectionstrict",
              "htseq_unionnonuniqueall","featurecounts_overlap","featurecounts_fractionaloverlap",
              "fadu_10em","fadu_default")
```

### Create timing plots
```{R,fig.width=8,fig.height=4}
benchmarking_pe1.timingplot.list <- list()
benchmarking_pe4.timingplot.list <- list()

for(i in 1:length(benchmarking_timingplotdata_pe1.list)){
  benchmarking_timingplotdata_pe1.list[[i]][,1] <- factor(benchmarking_timingplotdata_pe1.list[[i]][,1], levels = rev(colorder))
  
  benchmarking_pe1.timingplot.list[[i]] <- ggplot()+
    geom_bar(mapping=aes_string(x=benchmarking_timingplotdata_pe1.list[[i]][,1],
                                y=as.numeric(as.character(benchmarking_timingplotdata_pe1.list[[i]][,3])),
                                fill=benchmarking_timingplotdata_pe1.list[[i]][,2]),stat="identity")+
    coord_flip()+
    guides(fill = F)+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  
  benchmarking_timingplotdata_pe4.list[[i]][,1] <- factor(benchmarking_timingplotdata_pe4.list[[i]][,1], levels = rev(colorder))
  
  benchmarking_pe4.timingplot.list[[i]] <- ggplot()+
    geom_bar(mapping=aes_string(x=benchmarking_timingplotdata_pe4.list[[i]][,1],
                                y=as.numeric(as.character(benchmarking_timingplotdata_pe4.list[[i]][,3])),
                                fill=benchmarking_timingplotdata_pe4.list[[i]][,2]),stat="identity")+
    coord_flip()+
    guides(fill = F)+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
}

pdf(paste0(OUTPUT_DIR,"/timing_plot.pdf"),
    height=4,
    width=8)
plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/timing_plot.png"),
    height=4,
    width=8,
    units = "in",res=300)
plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]]),
          align = 'v',
          ncol=3)
```

![image](/images/timing_plot.png)

### Create VMem plots
```{R,fig.width=8,fig.height=4}
benchmarking_pe1.vmemplot.list <- list()
benchmarking_pe4.vmemplot.list <- list()

for(i in 1:length(benchmarking_vmemplotdata_pe1.list)){
  benchmarking_vmemplotdata_pe1.list[[i]][,1] <- factor(benchmarking_vmemplotdata_pe1.list[[i]][,1], levels = rev(colorder))
  
  benchmarking_pe1.vmemplot.list[[i]] <- ggplot()+
    geom_bar(mapping=aes_string(x=benchmarking_vmemplotdata_pe1.list[[i]][,1],
                                y=as.numeric(as.character(benchmarking_vmemplotdata_pe1.list[[i]][,2]))),
             fill="#7cae00",
             stat="identity")+
    scale_y_continuous(limits = c(0,5))+
    coord_flip()+
    guides(fill = F)+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  
  
  benchmarking_vmemplotdata_pe4.list[[i]][,1] <- factor(benchmarking_vmemplotdata_pe4.list[[i]][,1], levels = rev(colorder))
  
  benchmarking_pe4.vmemplot.list[[i]] <- ggplot()+
    geom_bar(mapping=aes_string(x=benchmarking_vmemplotdata_pe1.list[[i]][,1],
                                y=as.numeric(as.character(benchmarking_vmemplotdata_pe1.list[[i]][,2]))),
             fill="#7cae00",
             stat="identity")+
    scale_y_continuous(limits = c(0,5))+
    coord_flip()+
    guides(fill = F)+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
}


pdf(paste0(OUTPUT_DIR,"/vmem_plot.pdf"),
    height=4,
    width=8)
plot_grid(plotlist = list(benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]],
                          benchmarking_pe4.vmemplot.list[[1]],
                          benchmarking_pe4.vmemplot.list[[2]],
                          benchmarking_pe4.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/vmem_plot.png"),
    height=4,
    width=8,
    units = "in",res=300)
plot_grid(plotlist = list(benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]],
                          benchmarking_pe4.vmemplot.list[[1]],
                          benchmarking_pe4.vmemplot.list[[2]],
                          benchmarking_pe4.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

plot_grid(plotlist = list(benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]],
                          benchmarking_pe4.vmemplot.list[[1]],
                          benchmarking_pe4.vmemplot.list[[2]],
                          benchmarking_pe4.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
```

![image](/images/vmem_plot.png)

### Combine timing and VMem plots
```{R,fig.width=8,fig.height=8}
pdf(paste0(OUTPUT_DIR,"/timing_vmem_plot.pdf"),
    height=8,
    width=8)
plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]],
                          benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

png(paste0(OUTPUT_DIR,"/timing_vmem_plot.png"),
    height=8,
    width=8,
    units = "in",res=300)
plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]],
                          benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
dev.off()

plot_grid(plotlist = list(benchmarking_pe1.timingplot.list[[1]],
                          benchmarking_pe1.timingplot.list[[2]],
                          benchmarking_pe1.timingplot.list[[3]],
                          benchmarking_pe4.timingplot.list[[1]],
                          benchmarking_pe4.timingplot.list[[2]],
                          benchmarking_pe4.timingplot.list[[3]],
                          benchmarking_pe1.vmemplot.list[[1]],
                          benchmarking_pe1.vmemplot.list[[2]],
                          benchmarking_pe1.vmemplot.list[[3]]),
          align = 'v',
          ncol=3)
```

![image](/images/timing_vmem_plot.png)
