
# Table of Contents
<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Set software and directory paths](#set-software-and-directory-paths)
    - [Software](#software)
    - [Directories](#directories)
    - [Create directories](#create-directories)
    - [Set up reference files](#set-up-reference-files)
        - [Download canine, tick and all E. chaffeensis reference files](#download-canine-tick-and-all-e-chaffeensis-reference-files)
        - [Create combined canine/tick and E. chaffeensis references](#create-combined-caninetick-and-e-chaffeensis-references)
            - [Genomic references](#genomic-references)
            - [CDS references](#cds-references)
- [Identify core Ehrlichia genes between the 8 E. chaffeensis strains and Ehrlichia sp. HF](#identify-core-ehrlichia-genes-between-the-8-e-chaffeensis-strains-and-ehrlichia-sp-hf)
        - [Create CDS FNA files](#create-cds-fna-files)
        - [Create FAA files](#create-faa-files)
        - [Create BLASTP output file for PanOCT](#create-blastp-output-file-for-panoct)
        - [Create gene tags file for PanOCT](#create-gene-tags-file-for-panoct)
        - [Create gene attributes file for PanOCT](#create-gene-attributes-file-for-panoct)
        - [Run PanOCT to find orthologous gene clusters](#run-panoct-to-find-orthologous-gene-clusters)
- [Identify differentially expression genes in canine, tick, and E. chaffeensis strains across the data set](#identify-differentially-expression-genes-in-canine-tick-and-e-chaffeensis-strains-across-the-data-set)
    - [Create SRR mapping file](#create-srr-mapping-file)
    - [Download FASTQs from SRA](#download-fastqs-from-sra)
    - [Find GO terms and InterPro descriptions for canine, tick, and Ehrlichia genes](#find-go-terms-and-interpro-descriptions-for-canine-tick-and-ehrlichia-genes)
        - [Canine/Tick](#caninetick)
            - [Split FASTA files into 1000 sequence chunks](#split-fasta-files-into-1000-sequence-chunks)
            - [Run InterProScan on split FASTA files](#run-interproscan-on-split-fasta-files)
            - [Combine split InterProScan outputs](#combine-split-interproscan-outputs)
        - [Ehrlichia](#ehrlichia)
            - [Run InterProScan](#run-interproscan)
    - [Convert InterProScan outputs to geneinfo files](#convert-interproscan-outputs-to-geneinfo-files)
    - [Create reference indices](#create-reference-indices)
        - [HISAT2](#hisat2)
        - [Salmon](#salmon)
    - [Quantify canine, tick, and E. chaffeensis transcript expression levels](#quantify-canine-tick-and-e-chaffeensis-transcript-expression-levels)
        - [E. chaffeensis](#e-chaffeensis)
            - [Align reads to their respective combined references](#align-reads-to-their-respective-combined-references)
            - [Sort BAM files](#sort-bam-files)
            - [Remove nonsorted BAM files](#remove-nonsorted-bam-files)
            - [Index BAM files](#index-bam-files)
            - [Quantify E. chaffeensis genes from BAM files using FADU](#quantify-e-chaffeensis-genes-from-bam-files-using-fadu)
        - [Canine/Tick](#caninetick-1)
            - [Quantify canine/tick transcripts directly from reads using Salmon](#quantify-caninetick-transcripts-directly-from-reads-using-salmon)
    - [Conduct differential expression analysis](#conduct-differential-expression-analysis)
        - [Canine](#canine)
            - [Set R inputs](#set-r-inputs)
            - [Load R functions](#load-r-functions)
            - [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
            - [Create counts data frame](#create-counts-data-frame)
            - [Create TPM data frame](#create-tpm-data-frame)
            - [Set group levels](#set-group-levels)
            - [Conduct saturation analysis](#conduct-saturation-analysis)
            - [Exclude low count samples](#exclude-low-count-samples)
            - [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally)
            - [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff)
            - [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules)
        - [Tick](#tick)
            - [Set R inputs](#set-r-inputs-1)
            - [Load R functions](#load-r-functions-1)
            - [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
            - [Create counts data frame](#create-counts-data-frame-1)
            - [Create TPM data frame](#create-tpm-data-frame-1)
            - [Set group levels](#set-group-levels-1)
            - [Conduct saturation analysis](#conduct-saturation-analysis-1)
            - [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-1)
            - [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff-1)
            - [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules-1)
        - [Ehrlichia](#ehrlichia-1)
            - [Set R inputs](#set-r-inputs-2)
            - [Load R functions](#load-r-functions-2)
            - [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-2)
            - [Constructs core genome table that consists of only gene clusters consisting of one gene per strain](#constructs-core-genome-table-that-consists-of-only-gene-clusters-consisting-of-one-gene-per-strain)
            - [Create counts data frame](#create-counts-data-frame-2)
            - [Calculate average gene length for all core genes](#calculate-average-gene-length-for-all-core-genes)
            - [Create TPM data frame](#create-tpm-data-frame-2)
            - [Set group levels](#set-group-levels-2)
            - [Conduct saturation analysis](#conduct-saturation-analysis-2)
            - [Exclude low count samples and samples with only 1 replicate](#exclude-low-count-samples-and-samples-with-only-1-replicate)
            - [Identify differentially expressed genes longitudinally](#identify-differentially-expressed-genes-longitudinally-2)
            - [Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff](#conduct-pca-and-hierarchical-clustering-analyses-on-genes-that-passed-the-cpm-cutoff-2)
            - [Divide differentially expressed genes into expression modules](#divide-differentially-expressed-genes-into-expression-modules-2)
            - [Test each module partition for over-represented functional terms](#test-each-module-partition-for-over-represented-functional-terms)

<!-- /MarkdownTOC -->

# Set software and directory paths

For rerunning analyses, all paths in this section must be set by the user.

## Software
```{bash, eval = F}
JULIA_DEPOT_PATH=/home/mattchung/.julia
PYTHON_LIB_PATH=/usr/local/packages/python-3.5/lib

JULIA_BIN_DIR=/usr/local/bin
R_BIN_DIR=/usr/local/packages/r-3.6.0/bin

EMBOSS_BIN_DIR=/usr/local/packages/emboss-6.6.0/bin
FADU_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/FADU_v1.7
HISAT2_BIN_DIR=/usr/local/packages/hisat2-2.1.0
INTERPROSCAN_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/interproscan-5.34-73.0
NCBIBLAST_BIN_DIR=/usr/local/packages/ncbi-blast-2.2.26/bin
PANOCT_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/panoct_v3.23/bin
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.10.0/bin
```

## Directories
```{bash, eval = F}
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/PECHA/
SCRIPTS_DIR=~/scripts/
```

## Create directories
```{bash, eval = F}
mkdir -p "$WORKING_DIR"/bam
mkdir -p "$WORKING_DIR"/fadu
mkdir -p "$WORKING_DIR"/panoct
mkdir -p "$WORKING_DIR"/plots
mkdir -p "$WORKING_DIR"/reads
mkdir -p "$WORKING_DIR"/references
mkdir -p "$WORKING_DIR"/salmon
```

## Set up reference files

### Download canine, tick and all E. chaffeensis reference files

```{bash, eval = F}
## Canine
wget -O "$WORKING_DIR"/references/canine.fna.gz ftp://ftp.ensembl.org/pub/release-99/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.toplevel.fa.gz
wget -O "$WORKING_DIR"/references/canine.cds.fna.gz ftp://ftp.ensembl.org/pub/release-99/fasta/canis_familiaris/cds/Canis_familiaris.CanFam3.1.cds.all.fa.gz

## Tick
wget -O "$WORKING_DIR"/references/tick.fna.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/ixodes_scapularis/dna/Ixodes_scapularis.IscaW1.dna.toplevel.fa.gz
wget -O "$WORKING_DIR"/references/tick.cds.fna.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/ixodes_scapularis/cds/Ixodes_scapularis.IscaW1.cds.all.fa.gz

## E. chaffeensis
wget -O "$WORKING_DIR"/references/Arkansas.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/145/GCF_000013145.1_ASM1314v1/GCF_000013145.1_ASM1314v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/Heartland.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/815/GCF_000632815.1_ASM63281v1/GCF_000632815.1_ASM63281v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/Jax.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/865/GCF_000632865.1_ASM63286v1/GCF_000632865.1_ASM63286v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/Liberty.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/885/GCF_000632885.1_ASM63288v1/GCF_000632885.1_ASM63288v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/Osceola.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/905/GCF_000632905.1_ASM63290v1/GCF_000632905.1_ASM63290v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/StVincent.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/925/GCF_000632925.1_ASM63292v1/GCF_000632925.1_ASM63292v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/Wakulla.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/945/GCF_000632945.1_ASM63294v1/GCF_000632945.1_ASM63294v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/WestPaces.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/965/GCF_000632965.1_ASM63296v1/GCF_000632965.1_ASM63296v1_genomic.fna.gz
wget -O "$WORKING_DIR"/references/HF.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/845/GCF_000632845.1_ASM63284v1/GCF_000632845.1_ASM63284v1_genomic.fna.gz

wget -O "$WORKING_DIR"/references/Arkansas.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/145/GCF_000013145.1_ASM1314v1/GCF_000013145.1_ASM1314v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/Heartland.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/815/GCF_000632815.1_ASM63281v1/GCF_000632815.1_ASM63281v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/Jax.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/865/GCF_000632865.1_ASM63286v1/GCF_000632865.1_ASM63286v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/Liberty.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/885/GCF_000632885.1_ASM63288v1/GCF_000632885.1_ASM63288v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/Osceola.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/905/GCF_000632905.1_ASM63290v1/GCF_000632905.1_ASM63290v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/StVincent.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/925/GCF_000632925.1_ASM63292v1/GCF_000632925.1_ASM63292v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/Wakulla.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/945/GCF_000632945.1_ASM63294v1/GCF_000632945.1_ASM63294v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/WestPaces.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/965/GCF_000632965.1_ASM63296v1/GCF_000632965.1_ASM63296v1_genomic.gff.gz
wget -O "$WORKING_DIR"/references/HF.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/632/845/GCF_000632845.1_ASM63284v1/GCF_000632845.1_ASM63284v1_genomic.gff.gz

gunzip "$WORKING_DIR"/references/*gz
```

### Create combined canine/tick and E. chaffeensis references

#### Genomic references
```{bash, eval = F}
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Arkansas.fna > Arkansas_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Heartland.fna > Heartland_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Jax.fna > Jax_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Liberty.fna > Liberty_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Osceola.fna > Osceola_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/StVincent.fna > StVincent_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/Wakulla.fna > Wakulla_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/WestPaces.fna > WestPaces_canine.combined.fna
cat "$WORKING_DIR"/references/canine.fna "$WORKING_DIR"/references/HF.fna > HF_canine.combined.fna

cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Arkansas.fna > Arkansas_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Heartland.fna > Heartland_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Jax.fna > Jax_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Liberty.fna > Liberty_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Osceola.fna > Osceola_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/StVincent.fna > StVincent_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/Wakulla.fna > Wakulla_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/WestPaces.fna > WestPaces_tick.combined.fna
cat "$WORKING_DIR"/references/tick.fna "$WORKING_DIR"/references/HF.fna > HF_tick.combined.fna
```

#### CDS references
```{bash, eval = F}
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Arkansas.cds.fna > Arkansas_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Heartland.cds.fna > Heartland_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Jax.cds.fna > Jax_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Liberty.cds.fna > Liberty_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Osceola.cds.fna > Osceola_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/StVincent.cds.fna > StVincent_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/Wakulla.cds.fna > Wakulla_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/WestPaces.cds.fna > WestPaces_canine.cds.combined.fna
cat "$WORKING_DIR"/references/canine.cds.fna "$WORKING_DIR"/references/HF.cds.fna > HF_canine.cds.combined.fna

cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Arkansas.cds.fna > Arkansas_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Heartland.cds.fna > Heartland_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Jax.cds.fna > Jax_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Liberty.cds.fna > Liberty_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Osceola.cds.fna > Osceola_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/StVincent.cds.fna > StVincent_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/Wakulla.cds.fna > Wakulla_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/WestPaces.cds.fna > WestPaces_tick.cds.combined.fna
cat "$WORKING_DIR"/references/tick.cds.fna "$WORKING_DIR"/references/HF.cds.fna > HF_tick.cds.combined.fna
```

# Identify core Ehrlichia genes between the 8 E. chaffeensis strains and Ehrlichia sp. HF


### Create CDS FNA files

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
for GFF in "$REFERENCES_DIR"/*gff
do
	FNA="$(echo "$GFF" | sed "s/[.]gff$/.fna/g")"
	"$SAMTOOLS_BIN_DIR"/samtools faidx "$FNA"

	CDS_FNA="$(echo "$GFF" | sed "s/[.]gff$/.cds.fna/g")"
	rm "$CDS_FNA"

	awk -F "\t" '$3 == "gene" {print $0}' "$GFF" | grep "gene_biotype=protein_coding" | while read LINE
	do

	CONTIG="$(echo "$LINE" | awk -F "\t" '{print $1}')"
	START="$(echo "$LINE" | awk -F "\t" '{print $4}')"
	STOP="$(echo "$LINE" | awk -F "\t" '{print $5}')"
	STRAND=$(echo "$LINE" | awk -F "\t" '{print $7}')
	GENE="$(echo "$LINE" | awk -F "\t" '{print $9}' | sed "s/.*ID=//g" | sed "s/;.*//g")"

	if [ "$STRAND" = "+" ]
	then
		"$SAMTOOLS_BIN_DIR"/samtools faidx "$FNA" "$CONTIG":"$START"-"$STOP" | sed -e "s/>.*/>"$GENE"/g" >> "$CDS_FNA"
	else
		"$SAMTOOLS_BIN_DIR"/samtools faidx -i "$FNA" "$CONTIG":"$START"-"$STOP" | sed -e "s/>.*/>"$GENE"/g" >> "$CDS_FNA"
	fi
	done
done
```

### Create FAA files

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
rm "$REFERENCES_DIR"/combined_ehrlichia.faa
for CDS_FNA in "$REFERENCES_DIR"/*cds.fna
do
	"$EMBOSS_BIN_DIR"/transeq -sequence $CDS_FNA -outseq "$(echo "$CDS_FNA" | sed "s/[.]cds.fna/.faa/g")"
done

cat "$WORKING_DIR"/references/*[.]faa > "$REFERENCES_DIR"/combined_ehrlichia.faa
sed -i "s/_1$//g" "$REFERENCES_DIR"/*faa
```

### Create BLASTP output file for PanOCT

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
"$NCBIBLAST_BIN_DIR"/formatdb -i "$WORKING_DIR"/references/combined_ehrlichia.faa -p T
"$NCBIBLAST_BIN_DIR"/blastall -p blastp -i "$WORKING_DIR"/references/combined_ehrlichia.faa -d "$WORKING_DIR"/references/combined_ehrlichia.faa -e 1e-5 -m 8 -o "$WORKING_DIR"/panoct/blastall_out.txt
```

### Create gene tags file for PanOCT

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
cat "$WORKING_DIR"/references/*[.]gff | grep -v "^#" | cut -f1 | sort -n | uniq > "$WORKING_DIR"/panoct/gene_tags.txt
```

### Create gene attributes file for PanOCT

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
rm "$WORKING_DIR"/panoct/gene_atts.txt
awk -F "\t" '$3 == "gene" {print $0}' "$REFERENCES_DIR"/*gff | grep "gene_biotype=protein_coding" | while read LINE
do
	CONTIG="$(echo "$LINE" | awk -F "\t" '{print $1}' )"
	GENE="$(echo "$LINE" | awk -F "\t" '{print $9}' | sed "s/.*ID=//g" | sed "s/;.*//g")"
	START="$(echo "$LINE" | awk -F "\t" '{print $4}' )"
	STOP="$(echo "$LINE" | awk -F "\t" '{print $5}' )"
	PRODUCT="$(grep "$GENE" "$REFERENCES_DIR"/*gff | awk -F "\t" '$3 == "CDS" {print $9}' | sed "s/.*product=//g" | sed "s/;.*//g" | uniq)"

	echo -e ""$CONTIG"\t"$GENE"\t"$START"\t"$STOP"\t"$PRODUCT"\t"$CONTIG""  >> "$WORKING_DIR"/panoct/gene_atts.txt
done
```

### Run PanOCT to find orthologous gene clusters

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
"$PANOCT_BIN_DIR"/panoct.pl -b "$WORKING_DIR"/panoct -t blastall_out.txt -f gene_tags.txt -g gene_atts.txt -Q "$REFERENCES_DIR" -P combined_ehrlichia.faa  -S Y -L 1 -M Y -H Y -V Y -N Y -F 1.33 -G y -c 0,25,50,75,100 -T
```

# Identify differentially expression genes in canine, tick, and E. chaffeensis strains across the data set



## Create SRR mapping file

##### Commands
```{bash, eval = F}
vim "$WORKING_DIR"/pecha_groups.tsv
```

```{bash, eval = F}
SRR1188323	PECHA_Arkansas_mRNA1	#497FCA	16	Arkansas
SRR1188505	PECHA_Arkansas_mRNA2	#497FCA	16	Arkansas
SRR1188516	PECHA_Arkansas_mRNA3	#497FCA	16	Arkansas
SRR1188524	PECHA_Arkansas_totalRNA1	#497FCA	16	Arkansas
SRR1188570	PECHA_Arkansas_totalRNA2	#497FCA	16	Arkansas
SRR1188615	PECHA_Arkansas_totalRNA3	#497FCA	16	Arkansas
SRR1188623	PECHA_DH82_mRNA1	#000000	16	DH82
SRR1188622	PECHA_DH82_mRNA1	#000000	16	DH82
SRR1188624	PECHA_DH82_mRNA2	#000000	16	DH82
SRR1188625	PECHA_DH82_mRNA3	#000000	16	DH82
SRR1188626	PECHA_DH82_totalRNA1	#000000	16	DH82
SRR1188627	PECHA_DH82_totalRNA2	#000000	16	DH82
SRR1188628	PECHA_DH82_totalRNA3	#000000	16	DH82
SRR1188629	PECHA_Heartland_mRNA1	#6c3c83	16	Heartland
SRR1188630	PECHA_Heartland_mRNA2	#6c3c83	16	Heartland
SRR1190411	PECHA_Heartland_mRNA3	#6c3c83	16	Heartland
SRR1190422	PECHA_Heartland_totalRNA1	#6c3c83	16	Heartland
SRR1190421	PECHA_Heartland_totalRNA1	#6c3c83	16	Heartland
SRR1190436	PECHA_Heartland_totalRNA2	#6c3c83	16	Heartland
SRR1190444	PECHA_Heartland_totalRNA3	#6c3c83	16	Heartland
SRR1190443	PECHA_Heartland_totalRNA3	#6c3c83	16	Heartland
SRR1190449	PECHA_HF_mRNA1	#84af6b	16	HF
SRR1190450	PECHA_HF_mRNA2	#84af6b	16	HF
SRR1190456	PECHA_HF_mRNA3	#84af6b	16	HF
SRR1190447	PECHA_HF_totalRNA1	#84af6b	16	HF
SRR1190446	PECHA_HF_totalRNA1	#84af6b	16	HF
SRR1190445	PECHA_HF_totalRNA1	#84af6b	16	HF
SRR1188632	PECHA_HF_totalRNA2	#84af6b	16	HF
SRR1188631	PECHA_HF_totalRNA2	#84af6b	16	HF
SRR1188660	PECHA_HF_totalRNA3	#84af6b	16	HF
SRR1188703	PECHA_ISE6_Arkansas_mRNA1	#497fca	17	Arkansas
SRR1188747	PECHA_ISE6_Arkansas_mRNA4	#497fca	17	Arkansas
SRR1188748	PECHA_ISE6_Arkansas_mRNA5	#497fca	17	Arkansas
SRR1188749	PECHA_ISE6_Arkansas_totalRNA1	#497fca	17	Arkansas
SRR1188750	PECHA_ISE6_Arkansas_totalRNA4	#497fca	17	Arkansas
SRR1188751	PECHA_ISE6_Arkansas_totalRNA5	#497fca	17	Arkansas
SRR1188752	PECHA_ISE6_Heartland_mRNA1	#6c3c83	17	Heartland
SRR1188753	PECHA_ISE6_Heartland_mRNA2	#6c3c83	17	Heartland
SRR1188755	PECHA_ISE6_Heartland_mRNA3	#6c3c83	17	Heartland
SRR1189460	PECHA_ISE6_Heartland_totalRNA1	#6c3c83	17	Heartland
SRR1189488	PECHA_ISE6_Heartland_totalRNA2	#6c3c83	17	Heartland
SRR1189534	PECHA_ISE6_Heartland_totalRNA3	#6c3c83	17	Heartland
SRR1189544	PECHA_ISE6_HF_mRNA1	#84af6b	17	HF
SRR1189545	PECHA_ISE6_HF_mRNA2	#84af6b	17	HF
SRR1189546	PECHA_ISE6_HF_mRNA3	#84af6b	17	HF
SRR1189547	PECHA_ISE6_HF_totalRNA1	#84af6b	17	HF
SRR1189548	PECHA_ISE6_HF_totalRNA2	#84af6b	17	HF
SRR1189549	PECHA_ISE6_HF_totalRNA3	#84af6b	17	HF
SRR1189551	PECHA_ISE6_Jax_mRNA1	#e27660	17	Jacksonville
SRR1189550	PECHA_ISE6_Jax_mRNA1	#e27660	17	Jacksonville
SRR1189553	PECHA_ISE6_Jax_mRNA2	#e27660	17	Jacksonville
SRR1189552	PECHA_ISE6_Jax_mRNA2	#e27660	17	Jacksonville
SRR1189556	PECHA_ISE6_Jax_mRNA3	#e27660	17	Jacksonville
SRR1189555	PECHA_ISE6_Jax_mRNA3	#e27660	17	Jacksonville
SRR1189557	PECHA_ISE6_Jax_totalRNA1	#e27660	17	Jacksonville
SRR1189560	PECHA_ISE6_Jax_totalRNA2	#e27660	17	Jacksonville
SRR1189585	PECHA_ISE6_Jax_totalRNA3	#e27660	17	Jacksonville
SRR1189614	PECHA_ISE6_Liberty_mRNA1	#f0b67f	17	Liberty
SRR1189619	PECHA_ISE6_Liberty_mRNA2	#f0b67f	17	Liberty
SRR1189633	PECHA_ISE6_Liberty_mRNA3	#f0b67f	17	Liberty
SRR1189636	PECHA_ISE6_Liberty_totalRNA1	#f0b67f	17	Liberty
SRR1189644	PECHA_ISE6_Liberty_totalRNA2	#f0b67f	17	Liberty
SRR1189645	PECHA_ISE6_Liberty_totalRNA3	#f0b67f	17	Liberty
SRR1188683	PECHA_ISE6_mRNA1	#000000	17	ISE6
SRR1189646	PECHA_ISE6_mRNA2	#000000	17	ISE6
SRR1189647	PECHA_ISE6_mRNA3	#000000	17	ISE6
SRR1189648	PECHA_ISE6_Osceola_mRNA1	#83c9fc	17	Osceola
SRR1189649	PECHA_ISE6_Osceola_mRNA2	#83c9fc	17	Osceola
SRR1189657	PECHA_ISE6_Osceola_mRNA3	#83c9fc	17	Osceola
SRR1189650	PECHA_ISE6_Osceola_totalRNA1	#83c9fc	17	Osceola
SRR1189651	PECHA_ISE6_Osceola_totalRNA2	#83c9fc	17	Osceola
SRR1189652	PECHA_ISE6_Osceola_totalRNA3	#83c9fc	17	Osceola
SRR1189689	PECHA_ISE6_StVincent_mRNA1	#9966ab	17	St. Vincent
SRR1189688	PECHA_ISE6_StVincent_mRNA1	#9966ab	17	St. Vincent
SRR1189704	PECHA_ISE6_StVincent_mRNA2	#9966ab	17	St. Vincent
SRR1189703	PECHA_ISE6_StVincent_mRNA2	#9966ab	17	St. Vincent
SRR1189711	PECHA_ISE6_StVincent_mRNA3	#9966ab	17	St. Vincent
SRR1189710	PECHA_ISE6_StVincent_mRNA3	#9966ab	17	St. Vincent
SRR1189721	PECHA_ISE6_StVincent_totalRNA1	#9966ab	17	St. Vincent
SRR1189746	PECHA_ISE6_StVincent_totalRNA2	#9966ab	17	St. Vincent
SRR1189747	PECHA_ISE6_StVincent_totalRNA3	#9966ab	17	St. Vincent
SRR1189748	PECHA_ISE6_totalRNA1	#000000	17	ISE6
SRR1189749	PECHA_ISE6_totalRNA2	#000000	17	ISE6
SRR1189750	PECHA_ISE6_totalRNA3	#000000	17	ISE6
SRR1189751	PECHA_ISE6_Wakulla_mRNA1	#BF5073	17	Wakulla
SRR1189752	PECHA_ISE6_Wakulla_mRNA2	#BF5073	17	Wakulla
SRR1189753	PECHA_ISE6_Wakulla_mRNA3	#BF5073	17	Wakulla
SRR1189754	PECHA_ISE6_Wakulla_totalRNA1	#BF5073	17	Wakulla
SRR1189759	PECHA_ISE6_Wakulla_totalRNA2	#BF5073	17	Wakulla
SRR1189779	PECHA_ISE6_Wakulla_totalRNA3	#BF5073	17	Wakulla
SRR1189808	PECHA_ISE6_WestPaces_mRNA1	#c66b9f	17	West Paces
SRR1189807	PECHA_ISE6_WestPaces_mRNA1	#c66b9f	17	West Paces
SRR1189818	PECHA_ISE6_WestPaces_mRNA2	#c66b9f	17	West Paces
SRR1189817	PECHA_ISE6_WestPaces_mRNA2	#c66b9f	17	West Paces
SRR1189835	PECHA_ISE6_WestPaces_mRNA3	#c66b9f	17	West Paces
SRR1189834	PECHA_ISE6_WestPaces_mRNA3	#c66b9f	17	West Paces
SRR1189845	PECHA_ISE6_WestPaces_totalRNA1	#c66b9f	17	West Paces
SRR1189846	PECHA_ISE6_WestPaces_totalRNA2	#c66b9f	17	West Paces
SRR1189847	PECHA_ISE6_WestPaces_totalRNA3	#c66b9f	17	West Paces
SRR1189848	PECHA_Jax_mRNA1	#e27660	16	Jacksonville
SRR1189849	PECHA_Jax_mRNA2	#e27660	16	Jacksonville
SRR1189850	PECHA_Jax_mRNA3	#e27660	16	Jacksonville
SRR1189851	PECHA_Jax_totalRNA1	#e27660	16	Jacksonville
SRR1189853	PECHA_Jax_totalRNA2	#e27660	16	Jacksonville
SRR1189852	PECHA_Jax_totalRNA2	#e27660	16	Jacksonville
SRR1189867	PECHA_Jax_totalRNA3	#e27660	16	Jacksonville
SRR1189889	PECHA_Liberty_mRNA1	#f0b67f	16	Liberty
SRR1189895	PECHA_Liberty_mRNA2	#f0b67f	16	Liberty
SRR1189900	PECHA_Liberty_mRNA3	#f0b67f	16	Liberty
SRR1189912	PECHA_Liberty_totalRNA1	#f0b67f	16	Liberty
SRR1189933	PECHA_Liberty_totalRNA2	#f0b67f	16	Liberty
SRR1189942	PECHA_Liberty_totalRNA3	#f0b67f	16	Liberty
SRR1189951	PECHA_Osceola_mRNA1	#83c9fc	16	Osceola
SRR1189952	PECHA_Osceola_mRNA2	#83c9fc	16	Osceola
SRR1189953	PECHA_Osceola_mRNA3	#83c9fc	16	Osceola
SRR1189955	PECHA_Osceola_totalRNA1	#83c9fc	16	Osceola
SRR1189954	PECHA_Osceola_totalRNA1	#83c9fc	16	Osceola
SRR1189956	PECHA_Osceola_totalRNA2	#83c9fc	16	Osceola
SRR1189963	PECHA_Osceola_totalRNA3	#83c9fc	16	Osceola
SRR1189964	PECHA_StVincent_mRNA1	#9966ab	16	St. Vincent
SRR1189965	PECHA_StVincent_mRNA2	#9966ab	16	St. Vincent
SRR1189966	PECHA_StVincent_mRNA3	#9966ab	16	St. Vincent
SRR1189967	PECHA_StVincent_totalRNA1	#9966ab	16	St. Vincent
SRR1189971	PECHA_StVincent_totalRNA2	#9966ab	16	St. Vincent
SRR1189970	PECHA_StVincent_totalRNA2	#9966ab	16	St. Vincent
SRR1189969	PECHA_StVincent_totalRNA2	#9966ab	16	St. Vincent
SRR1189968	PECHA_StVincent_totalRNA2	#9966ab	16	St. Vincent
SRR1189972	PECHA_StVincent_totalRNA3	#9966ab	16	St. Vincent
SRR1189983	PECHA_Wakulla_mRNA1	#BF5073	16	Wakulla
SRR1189984	PECHA_Wakulla_mRNA2	#BF5073	16	Wakulla
SRR1189985	PECHA_Wakulla_mRNA3	#BF5073	16	Wakulla
SRR1189986	PECHA_Wakulla_totalRNA1	#BF5073	16	Wakulla
SRR1189987	PECHA_Wakulla_totalRNA2	#BF5073	16	Wakulla
SRR1190000	PECHA_Wakulla_totalRNA3	#BF5073	16	Wakulla
SRR1190026	PECHA_WestPaces_mRNA1	#c66b9f	16	West Paces
SRR1190032	PECHA_WestPaces_mRNA2	#c66b9f	16	West Paces
SRR1190040	PECHA_WestPaces_mRNA3	#c66b9f	16	West Paces
SRR1190039	PECHA_WestPaces_mRNA3	#c66b9f	16	West Paces
SRR1190042	PECHA_WestPaces_totalRNA1	#c66b9f	16	West Paces
SRR1190041	PECHA_WestPaces_totalRNA1	#c66b9f	16	West Paces
SRR1190069	PECHA_WestPaces_totalRNA2	#c66b9f	16	West Paces
SRR1190075	PECHA_WestPaces_totalRNA3	#c66b9f	16	West Paces
```

## Download FASTQs from SRA

##### Inputs
```{bash, eval = F}
SRR_MAP="$WORKING_DIR"/pecha_groups.tsv
READS_DIR="$WORKING_DIR"/reads
```

##### Commands
```{bash, eval = F}
cut -f1 "$SRR_MAP" | while read SRR_ID
do
	qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --gzip --split-files "$SRR_ID" -O "$READS_DIR"
done
```

## Find GO terms and InterPro descriptions for canine, tick, and Ehrlichia genes 

### Canine/Tick

#### Split FASTA files into 1000 sequence chunks

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("canine%d.cds.fna",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < "$WORKING_DIR"/references/canine.cds.fna
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("tick%d.cds.fna",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < "$WORKING_DIR"/references/tick.cds.fna
```

#### Run InterProScan on split FASTA files

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
THREADS=16
```

##### Commands
```{bash, eval = F}
ls "$REFERENCES_DIR"/*.fna | grep canine | grep -v combined | grep -v canine.fna | grep -v canine.cds.fna | while read FNA
do
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$INTERPROSCAN_BIN_DIR"/interproscan.sh -i "$FNA" -f tsv -o "$FNA".interproscan.tsv --seqtype n --goterms --iprlookup" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=20G -N interproscan -wd "$REFERENCES_DIR"
done

ls "$REFERENCES_DIR"/*.fna | grep tick | grep -v combined | grep -v tick.fna | grep -v tick.cds.fna | while read FNA
do
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$INTERPROSCAN_BIN_DIR"/interproscan.sh -i "$FNA" -f tsv -o "$FNA".interproscan.tsv --seqtype n --goterms --iprlookup" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=20G -N interproscan -wd "$REFERENCES_DIR"
done
```

#### Combine split InterProScan outputs

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
THREADS=16
```

##### Commands
```{bash, eval = F}
cat "$REFERENCES_DIR"/canine[0-9]*.cds.fna.interproscan.tsv > "$REFERENCES_DIR"/canine.cds.fna.interproscan.tsv
cat "$REFERENCES_DIR"/tick[0-9]*.cds.fna.interproscan.tsv > "$REFERENCES_DIR"/tick.cds.fna.interproscan.tsv

rm "$REFERENCES_DIR"/canine[0-9]*.cds.fna.interproscan.tsv
rm "$REFERENCES_DIR"/tick[0-9]*.cds.fna.interproscan.tsv
```

### Ehrlichia

#### Run InterProScan

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
THREADS=16
```

##### Commands
```{bash, eval = F}
ls "$REFERENCES_DIR"/*cds.fna | grep -v canine | grep -v tick | while read FNA
do
	echo -e "export LD_LIBRARY_PATH="$PYTHON_LIB_PATH":"$LD_LIBRARY_PATH"\n"$INTERPROSCAN_BIN_DIR"/interproscan.sh -i "$FNA" -f tsv -o "$FNA".interproscan.tsv --seqtype n --goterms --iprlookup" | qsub -P jdhotopp-lab -q threaded.q  -pe thread "$THREADS" -l mem_free=20G -N interproscan -wd "$REFERENCES_DIR"
done
```

## Convert InterProScan outputs to geneinfo files

```{bash, eval = F}
ls "$REFERENCES_DIR"/*interproscan.tsv | while read IPRSCAN_OUTPUT
do
	"$R_BIN_DIR"/Rscript "$SCRIPTS_DIR"/interproscan2geneinfo_v2.R "$IPRSCAN_OUTPUT"
done
```

## Create reference indices

### HISAT2

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
ls "$REFERENCES_DIR"/*combined.fna | grep -v cds | while read FNA
do
	qsub -P jdhotopp-lab -l mem_free=20G -N hisat_index -wd "$REFERENCES_DIR" -b y "$HISAT2_BIN_DIR"/hisat2-build --large-index "$FNA" "$FNA"
done

qsub -P jdhotopp-lab -l mem_free=20G -N hisat2_index -wd "$REFERENCES_DIR" -b y "$HISAT2_BIN_DIR"/hisat2-build --large-index "$REFERENCES_DIR"/canine.fna "$REFERENCES_DIR"/canine.fna
qsub -P jdhotopp-lab -l mem_free=20G -N hisat2_index -wd "$REFERENCES_DIR" -b y "$HISAT2_BIN_DIR"/hisat2-build --large-index "$REFERENCES_DIR"/tick.fna "$REFERENCES_DIR"/tick.fna
```

### Salmon

##### Inputs
```{bash, eval = F}
REFERENCES_DIR="$WORKING_DIR"/references/
```

##### Commands
```{bash, eval = F}
ls "$REFERENCES_DIR"/*cds.combined.fna | while read FNA
do
	qsub -P jdhotopp-lab -l mem_free=5G -N salmon_index -wd "$REFERENCES_DIR" -b y "$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$FNA" -i "$FNA".salmon.index
done
qsub -P jdhotopp-lab -l mem_free=5G -N salmon_index -wd "$REFERENCES_DIR" -b y "$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$REFERENCES_DIR"/canine.cds.fna -i "$REFERENCES_DIR"/canine.cds.fna.salmon.index
qsub -P jdhotopp-lab -l mem_free=5G -N salmon_index -wd "$REFERENCES_DIR" -b y "$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$REFERENCES_DIR"/tick.cds.fna -i "$REFERENCES_DIR"/tick.cds.fna.salmon.index
```

## Quantify canine, tick, and E. chaffeensis transcript expression levels

### E. chaffeensis

#### Align reads to their respective combined references

##### Inputs
```{bash, eval = F}
OUTPUT_DIR="$WORKING_DIR"/bam
REFERENCES_DIR="$WORKING_DIR"/references
READS_DIR="$WORKING_DIR"/reads
SRR_MAP="$WORKING_DIR"/pecha_groups.tsv
THREADS=16
```

##### Commands
```{bash, eval = F}
cat "$SRR_MAP" | while read LINE
do
SRR="$(echo "$LINE" | awk -F "\t" '{print $1}')"
ORG1="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $2}')"
ORG2="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $3}' | sed "s/totalRNA.//g")"

if [ "$ORG1"  == "DH82" ]; then
FNA="$REFERENCES_DIR"/canine.fna
elif [ "$ORG1"  == "ISE6" ] && [ "$ORG2"  == "" ]; then
FNA="$REFERENCES_DIR"/tick.fna
elif [ "$ORG1"  == "ISE6" ] && [[ "$ORG2"  =~ "mRNA" ]]; then
FNA="$REFERENCES_DIR"/tick.fna
elif [ "$ORG1"  == "ISE6" ]; then
FNA="$REFERENCES_DIR"/"$ORG2"_tick.combined.fna
else
FNA="$REFERENCES_DIR"/"$ORG1"_canine.combined.fna
fi 

if [ ! -f "$OUTPUT_DIR"/"$SRR".bam ]; then
echo ""$HISAT2_BIN_DIR"/hisat2 -p "$THREADS" -x "$FNA" -1 "$READS_DIR"/"$SRR"_1.fastq.gz -2 "$READS_DIR"/"$SRR"_2.fastq.gz | "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$OUTPUT_DIR"/"$SRR".bam -"  | qsub -q threaded.q  -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=20G -N hisat2 -wd "$OUTPUT_DIR"
fi
done
```

#### Sort BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
THREADS=4
```

##### Commands
```{bash, eval = F}
for BAM in $(find $BAM_DIR -name "*[.]bam" | grep -v "sortedbyposition")
do
if [ ! -f "$OUTPUT_DIR"/"$SRR".sortedbyposition.bam ]; then
qsub -q threaded.q -pe thread "$THREADS" -P jdhotopp-lab -l mem_free=2G -N sort -wd "$BAM_DIR" -b y "$SAMTOOLS_BIN_DIR"/samtools sort "$BAM" -@ "$THREADS" -o "$(echo $BAM | sed "s/[.]bam$/.sortedbyposition.bam"/g)"
fi
done
```

#### Remove nonsorted BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
SRR_MAP="$WORKING_DIR"/pecha_groups.tsv
THREADS=4
```
##### Commands
```{bash, eval = F}
cat "$SRR_MAP" | while read LINE
do
SRR="$(echo "$LINE" | awk -F "\t" '{print $1}')"

if [ -f "$BAM_DIR"/"$SRR".sortedbyposition.bam ] && [ -f "$BAM_DIR"/"$SRR".bam ]; then
rm "$BAM_DIR"/"$SRR".bam
fi
done
```
#### Index BAM files

##### Inputs
```{bash, eval = F}
BAM_DIR="$WORKING_DIR"/bam
```

##### Commands
```{bash, eval = F}
for BAM in $(find $BAM_DIR -name "*[.]sortedbyposition.bam")
do
if [ ! -f "$BAM".bai ]; then
qsub -P jdhotopp-lab -l mem_free=2G -N index -wd "$BAM_DIR" -b y "$SAMTOOLS_BIN_DIR"/samtools index "$BAM"
fi
done

```

#### Quantify E. chaffeensis genes from BAM files using FADU

##### Inputs
```{bash, eval = F}
SRR_MAP="$WORKING_DIR"/pecha_groups.tsv
BAM_DIR="$WORKING_DIR"/bam
FEAT_TYPE="gene"
STRANDEDNESS="no"
ATTR_ID="ID"
OUTPUT_DIR="$WORKING_DIR"/fadu
```

##### Commands
```{bash, eval = F}
cat "$SRR_MAP" | while read LINE
do
SRR="$(echo "$LINE" | awk -F "\t" '{print $1}')"
ORG1="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $2}')"
ORG2="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $3}' | sed "s/totalRNA.//g")"

if [ "$ORG1"  == "DH82" ]; then
GFF=""
elif [ "$ORG1"  == "ISE6" ] && [ "$ORG2"  == "" ]; then
GFF=""
elif [ "$ORG1"  == "ISE6" ] && [[ "$ORG2"  =~ "mRNA" ]]; then
GFF=""
elif [ "$ORG1"  == "ISE6" ]; then
GFF="$REFERENCES_DIR"/"$ORG2".gff
else
GFF="$REFERENCES_DIR"/"$ORG1".gff
fi 

if [ "$GFF"  != "" ]; then
echo -e "export JULIA_DEPOT_PATH="$JULIA_DEPOT_PATH"\n"$JULIA_BIN_DIR"/julia "$FADU_BIN_DIR"/fadu.jl -g "$GFF" -b "$BAM_DIR"/"$SRR".sortedbyposition.bam -o "$OUTPUT_DIR" -s "$STRANDEDNESS" -f "$FEAT_TYPE" -a "$ATTR_ID"" | qsub -P jdhotopp-lab -l mem_free=5G -N fadu -wd "$OUTPUT_DIR"
fi
done
```

### Canine/Tick

#### Quantify canine/tick transcripts directly from reads using Salmon

##### Inputs
```{bash, eval = F}
OUTPUT_DIR="$WORKING_DIR"/salmon
REFERENCES_DIR="$WORKING_DIR"/references
READS_DIR="$WORKING_DIR"/reads
SRR_MAP="$WORKING_DIR"/pecha_groups.tsv
THREADS=4
```

##### Commands
```{bash, eval = F}
cat "$SRR_MAP" | while read LINE
do
	SRR="$(echo "$LINE" | awk -F "\t" '{print $1}')"
	ORG1="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $2}')"
	ORG2="$(echo "$LINE" | awk -F "\t" '{print $2}' | awk -F "_" '{print $3}' | sed "s/totalRNA.//g")"

	if [ "$ORG1"  == "DH82" ]; then
		FNA="$REFERENCES_DIR"/canine.cds.fna
	elif [ "$ORG1"  == "ISE6" ] && [ "$ORG2"  == "" ]; then
		FNA="$REFERENCES_DIR"/tick.cds.fna
	elif [ "$ORG1"  == "ISE6" ] && [[ "$ORG2"  =~ "mRNA" ]]; then
		FNA="$REFERENCES_DIR"/tick.cds.fna
	elif [ "$ORG1"  == "ISE6" ]; then
		FNA="$REFERENCES_DIR"/"$ORG2"_tick.cds.combined.fna
	else
		FNA="$REFERENCES_DIR"/"$ORG1"_canine.cds.combined.fna
	fi 

	if [ ! -f "$OUTPUT_DIR"/"$SRR"/quant.sf ]; then
	echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$FNA".salmon.index --libType A -1 "$READS_DIR"/"$SRR"_1.fastq.gz -2 "$READS_DIR"/"$SRR"_2.fastq.gz -p "$THREADS" -o "$OUTPUT_DIR"/"$SRR" --validateMappings  --allowDovetail" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$OUTPUT_DIR"
	fi
done
```

## Conduct differential expression analysis

### Canine
#### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/PECHA/"
SALMON_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/salmon"
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/PECHA/pecha_groups.tsv"
```

#### Load R functions
```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
                                              ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}

find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
      tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
                                                 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}

functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
                               gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
                                      gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
      freq.all <- functionalterms.list[[i]][j,2]
      freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
                            functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
                            0)
      genes.all <- nrow(geneinfo)
      genes.subset <- nrow(geneinfo.subset)

      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      
      term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

#### Load packages and view sessionInfo
```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.68                  fastcluster_1.1.25          dynamicTreeCut_1.63-1       pvclust_2.0-0              
 [5] gplots_3.0.1.1              ggdendro_0.1-20             FactoMineR_1.42             DESeq2_1.22.2              
 [9] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0         
[13] Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         IRanges_2.16.0             
[17] S4Vectors_0.20.1            BiocGenerics_0.28.0         dendextend_1.12.0           gridExtra_2.3              
[21] ggplot2_3.2.0               vegan_2.5-5                 lattice_0.20-35             permute_0.9-5              
[25] edgeR_3.24.3                limma_3.38.3               

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       htmlTable_1.13.1       XVector_0.22.0         base64enc_0.1-3       
 [5] rstudioapi_0.10        bit64_0.9-7            mvtnorm_1.0-11         AnnotationDbi_1.44.0  
 [9] codetools_0.2-15       splines_3.5.0          leaps_3.0              doParallel_1.0.15     
[13] impute_1.56.0          robustbase_0.93-5      geneplotter_1.60.0     knitr_1.23            
[17] zeallot_0.1.0          Formula_1.2-3          annotate_1.60.1        cluster_2.0.7-1       
[21] GO.db_3.7.0            rrcov_1.4-7            compiler_3.5.0         backports_1.1.4       
[25] assertthat_0.2.1       Matrix_1.2-14          lazyeval_0.2.2         acepack_1.4.1         
[29] htmltools_0.3.6        tools_3.5.0            gtable_0.3.0           glue_1.3.1            
[33] GenomeInfoDbData_1.2.0 dplyr_0.8.3            Rcpp_1.0.2             vctrs_0.2.0           
[37] gdata_2.18.0           preprocessCore_1.44.0  nlme_3.1-137           iterators_1.0.12      
[41] xfun_0.8               stringr_1.4.0          gtools_3.8.1           XML_3.98-1.20         
[45] DEoptimR_1.0-8         zlibbioc_1.28.0        MASS_7.3-51.4          scales_1.0.0          
[49] RColorBrewer_1.1-2     yaml_2.2.0             memoise_1.1.0          rpart_4.1-13          
[53] latticeExtra_0.6-28    stringi_1.4.3          RSQLite_2.1.2          genefilter_1.64.0     
[57] pcaPP_1.9-73           foreach_1.4.7          checkmate_1.9.4        caTools_1.17.1.2      
[61] rlang_0.4.0            pkgconfig_2.0.2        bitops_1.0-6           purrr_0.3.2           
[65] htmlwidgets_1.3        labeling_0.3           bit_1.1-14             tidyselect_0.2.5      
[69] robust_0.4-18.1        magrittr_1.5           R6_2.4.0               Hmisc_4.2-0           
[73] fit.models_0.5-14      DBI_1.0.0              pillar_1.4.2           foreign_0.8-70        
[77] withr_2.1.2            mgcv_1.8-23            survival_2.41-3        scatterplot3d_0.3-41  
[81] RCurl_1.95-4.12        nnet_7.3-12            tibble_2.1.3           crayon_1.3.4          
[85] KernSmooth_2.23-15     viridis_0.5.1          locfit_1.5-9.1         grid_3.5.0            
[89] data.table_1.12.2      blob_1.2.0             digest_0.6.20          flashClust_1.01-2     
[93] xtable_1.8-4           munsell_0.5.0          viridisLite_0.3.0     
```

#### Create counts data frame
```{R}
groups <- read.delim(GROUPS.PATH, header = F)
groups <- groups[intersect(grep("mRNA", groups[,2]),grep("ISE6", groups[,2], invert = T)),]

rownames <- as.character(read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[grep("ENSCAFT", read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[,1]),1])
colnames <- unique(groups[,2])

counts <- as.data.frame(matrix(0,
                               nrow = length(rownames),
                               ncol = length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

for(i in 1:ncol(counts)){
  srr.vector <- groups[groups[,2] == colnames(counts[i]),1]
  for(j in 1:length(srr.vector)){
    counts.subset <- read.delim(paste0(SALMON_OUTPUT.DIR, "/",srr.vector[j],"/quant.sf"))
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),5]
  }
}

write.table(counts,
            paste0(WORKING.DIR,"/canine_counts.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

colSums(counts)
```

```{R, eval = F}
 PECHA_Arkansas_mRNA1  PECHA_Arkansas_mRNA2  PECHA_Arkansas_mRNA3      PECHA_DH82_mRNA1      PECHA_DH82_mRNA2 
          4446675.992           1118448.696           1588787.159            612736.010           1219321.996 
     PECHA_DH82_mRNA3 PECHA_Heartland_mRNA1 PECHA_Heartland_mRNA2 PECHA_Heartland_mRNA3        PECHA_HF_mRNA1 
          2046860.952            620001.016           1074528.431           1429786.216           1748842.155 
       PECHA_HF_mRNA2        PECHA_HF_mRNA3       PECHA_Jax_mRNA1       PECHA_Jax_mRNA2       PECHA_Jax_mRNA3 
          1136318.861           1523891.337           1541799.067              2288.996           1468686.449 
  PECHA_Liberty_mRNA1   PECHA_Liberty_mRNA2   PECHA_Liberty_mRNA3   PECHA_Osceola_mRNA1   PECHA_Osceola_mRNA2 
          1207815.323            984212.022           2720332.698           1011223.949            967659.504 
  PECHA_Osceola_mRNA3 PECHA_StVincent_mRNA1 PECHA_StVincent_mRNA2 PECHA_StVincent_mRNA3   PECHA_Wakulla_mRNA1 
          2091371.252            832188.255             16571.139            995731.667              8867.996 
  PECHA_Wakulla_mRNA2   PECHA_Wakulla_mRNA3 PECHA_WestPaces_mRNA1 PECHA_WestPaces_mRNA2 PECHA_WestPaces_mRNA3 
          1379824.521           1643445.731           1008953.088            973650.126            725770.658 
```

#### Create TPM data frame
```{R}
genelength <- read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))
genelength <- genelength[match(rownames(counts),genelength[,1]),2]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
            paste0(WORKING.DIR,"/canine_tpm.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

dim(tpm)
```

#### Set group levels
```{R}
groups <- unique(groups[2:5])
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
```

#### Conduct saturation analysis
```{R, fig.height=5, fig.width=6}
rarefy.counts <- round(counts,0)
raremax <- round(min(rowSums(t(rarefy.counts))),0)
srare <- rarefy(t(rarefy.counts),raremax) 

rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)

rarefy.df <- as.data.frame(matrix(nrow = 0,
                                  ncol = 5))
rarefy.points.df <- rarefy.df
for(i in 1:length(rarefy.raw.df)){
  steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[i]])))
  detected_genes <- as.numeric(rarefy.raw.df[[i]])
  rarefy.df <- as.data.frame(rbind(rarefy.df,
                                   cbind(as.numeric(steps),as.numeric(detected_genes),as.character(groups[i,1]),as.character(groups[i,2]),groups[i,3])))
  rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
                                          cbind(as.numeric(max(steps)),as.numeric(max(detected_genes)),as.character(groups[i,1]),as.character(groups[i,2],groups[i,3]))))
  
}
rarefy.plot <- ggplot()+
  geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4]),size = 3)+
  guides(colour = F,shape = F)+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x="reads mapping to genes", y="genes detected", color = "Sample")+
  #coord_cartesian(xlim=c(0,100000))+
  theme_bw()


pdf(paste0(WORKING.DIR,"/plots/canine_rarefication_plot.pdf"),
    height=5,
    width=6)
print(rarefy.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_rarefication_plot.png"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

![image](/images/canine_rarefication_plot.png)

#### Exclude low count samples

From the rarefaction analysis, the following samples were excluded (<20k gene counts): "PECHA_Jax_mRNA2", "PECHA_Wakulla_mRNA1", "PECHA_SaintVincent_mRNA2"

```{R}
exclude.samples <- c("PECHA_Jax_mRNA2","PECHA_Wakulla_mRNA1","PECHA_StVincent_mRNA2")
groups <- groups[!(groups[,1] %in% exclude.samples),]

groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))

counts <- counts[,colnames(counts) %in% groups[,1]]
tpm <- tpm[,colnames(tpm) %in% groups[,1]]
```

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = groups[,4])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,4]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,4])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
            paste0(WORKING.DIR,"/canine_counts_edgeR_longitudinal.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```


```{R, eval = F}
[1] "26024 genes excluded with CPM cutoff"
[1] "330 DE genes identified using edgeR longitudinal"
```

#### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

##### Create sample legend for PCA and hierarchical clustering plots

```{R, fig.height = 2, fig.width = 10}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_point(aes(color = groups[,4],shape=groups[,4]), size = 4)+
    scale_shape_manual(values = as.numeric(as.character(groups[,3])))+
    scale_color_manual(values = as.character(unique(groups[,2])))+
    guides(shape = guide_legend(title = "Samples", title.position = "top",nrow=2),
           colour = guide_legend(title = "Samples", title.position = "top",nrow=2))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/canine_pca_hc_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_pca_hc_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/canine_pca_hc_samplekey.png)

##### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

```{R, fig.height=5, fig.width=10}
dendrogram <- as.data.frame(t(scale(t(log2(tpm.keep+1)))))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

structure <- get_dendro_structure(result)
dendro.data <- get_dendro_data(result)
bootstrap.positions <- get_dendro_bootstraps(dendro.data)
  
points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                 structure))
dendrogroups <- groups[,4][result$hclust$order]
dendrocol <- groups[,2][result$hclust$order]
dendroshape <- groups[,3][result$hclust$order]
dendrosize <- colSums(counts.keep)[result$hclust$order]
#dendrosize <- 1

dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
  geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize, shape = dendroshape))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
  guides(colour = guide_legend(ncol = 2), size = F, shape = F)+
  #scale_x_discrete(limits = as.character(dendrolabels))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

for(i in 1:length(result$edges$bp)){
  text <- round(result$edges$bp[i] * 100,0)
  dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.4, y=bootstrap.positions[i,2] + 0.02, size = 2)
}
pdf(paste0(WORKING.DIR,"/plots/canine_tpm_kept_dendrogram.pdf"),
    height=5,
    width=10)
print(dendrogram.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_tpm_kept_dendrogram.png"),
    height=5,
    width=10,
	units = "in",res=300)
print(dendrogram.plot)
dev.off()

print(dendrogram.plot)
```

![image](/images/canine_tpm_kept_dendrogram.png)

##### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

```{R,fig.height=5,fig.width=5}
pca.df <- t(scale(t(log2(tpm.keep + 1))))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.keep) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = dendrocol,size = dendrosize, shape = dendroshape))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
       x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
       y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  guides(color = F,size = F, shape = F)+
  # guides(colour = guide_legend(ncol = 2))+
  scale_color_manual(values = levels(groups[,2]))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/canine_tpm_kept_pca.pdf"),
    height=5,
    width=5)
print(pca.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_tpm_kept_pca.png"),
    height=5,
    width=5,
	units = "in",res=300)
print(pca.plot)
dev.off()

print(pca.plot)
```

![image](/images/canine_tpm_kept_pca.png)

#### Divide differentially expressed genes into expression modules

##### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)

```

![image](/images/canine_tpm_de_wgcna_soft_power_plot.png)

##### Identify expression modules 

```{R, fig.height=5, fig.width=12}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 12, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 12, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/canine_tpm_de_wgcna_merge_eigengenes_plot.png)

##### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
            paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_modules.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/canine_tpm_de_wgcna_eigengene_dendrogram.png)

##### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,2])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)
```

###### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_line(aes(color = groups[,4]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,2])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=3))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/canine_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/canine_hm_samplekey.png)

###### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
                       colours=hmcol,
                       breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/canine_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/canine_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/canine_hm_zscorelog2tpmkey.png)

###### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/canine_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

###### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/canine_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/canine_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)

### Tick

#### Set R inputs
```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/PECHA/"
SALMON_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/salmon"
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/PECHA/pecha_groups.tsv"
```

#### Load R functions
```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
                                              ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}

find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
      tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
                                                 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}

functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
                               gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
                                      gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
      freq.all <- functionalterms.list[[i]][j,2]
      freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
                            functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
                            0)
      genes.all <- nrow(geneinfo)
      genes.subset <- nrow(geneinfo.subset)

      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      
      term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

#### Load packages and view sessionInfo
```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.68                  fastcluster_1.1.25          dynamicTreeCut_1.63-1       pvclust_2.0-0              
 [5] gplots_3.0.1.1              ggdendro_0.1-20             FactoMineR_1.42             DESeq2_1.22.2              
 [9] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0         
[13] Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         IRanges_2.16.0             
[17] S4Vectors_0.20.1            BiocGenerics_0.28.0         dendextend_1.12.0           gridExtra_2.3              
[21] ggplot2_3.2.0               vegan_2.5-5                 lattice_0.20-35             permute_0.9-5              
[25] edgeR_3.24.3                limma_3.38.3               

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       htmlTable_1.13.1       XVector_0.22.0         base64enc_0.1-3       
 [5] rstudioapi_0.10        bit64_0.9-7            mvtnorm_1.0-11         AnnotationDbi_1.44.0  
 [9] codetools_0.2-15       splines_3.5.0          leaps_3.0              doParallel_1.0.15     
[13] impute_1.56.0          robustbase_0.93-5      geneplotter_1.60.0     knitr_1.23            
[17] zeallot_0.1.0          Formula_1.2-3          annotate_1.60.1        cluster_2.0.7-1       
[21] GO.db_3.7.0            rrcov_1.4-7            compiler_3.5.0         backports_1.1.4       
[25] assertthat_0.2.1       Matrix_1.2-14          lazyeval_0.2.2         acepack_1.4.1         
[29] htmltools_0.3.6        tools_3.5.0            gtable_0.3.0           glue_1.3.1            
[33] GenomeInfoDbData_1.2.0 dplyr_0.8.3            Rcpp_1.0.2             vctrs_0.2.0           
[37] gdata_2.18.0           preprocessCore_1.44.0  nlme_3.1-137           iterators_1.0.12      
[41] xfun_0.8               stringr_1.4.0          gtools_3.8.1           XML_3.98-1.20         
[45] DEoptimR_1.0-8         zlibbioc_1.28.0        MASS_7.3-51.4          scales_1.0.0          
[49] RColorBrewer_1.1-2     yaml_2.2.0             memoise_1.1.0          rpart_4.1-13          
[53] latticeExtra_0.6-28    stringi_1.4.3          RSQLite_2.1.2          genefilter_1.64.0     
[57] pcaPP_1.9-73           foreach_1.4.7          checkmate_1.9.4        caTools_1.17.1.2      
[61] rlang_0.4.0            pkgconfig_2.0.2        bitops_1.0-6           purrr_0.3.2           
[65] htmlwidgets_1.3        labeling_0.3           bit_1.1-14             tidyselect_0.2.5      
[69] robust_0.4-18.1        magrittr_1.5           R6_2.4.0               Hmisc_4.2-0           
[73] fit.models_0.5-14      DBI_1.0.0              pillar_1.4.2           foreign_0.8-70        
[77] withr_2.1.2            mgcv_1.8-23            survival_2.41-3        scatterplot3d_0.3-41  
[81] RCurl_1.95-4.12        nnet_7.3-12            tibble_2.1.3           crayon_1.3.4          
[85] KernSmooth_2.23-15     viridis_0.5.1          locfit_1.5-9.1         grid_3.5.0            
[89] data.table_1.12.2      blob_1.2.0             digest_0.6.20          flashClust_1.01-2     
[93] xtable_1.8-4           munsell_0.5.0          viridisLite_0.3.0     
```

#### Create counts data frame
```{R}
groups <- read.delim(GROUPS.PATH, header = F)
groups <- groups[intersect(grep("mRNA", groups[,2]),grep("ISE6", groups[,2])),]

rownames <- as.character(read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[grep("ISCW", read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))[,1]),1])
colnames <- unique(groups[,2])

counts <- as.data.frame(matrix(0,
                               nrow = length(rownames),
                               ncol = length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

for(i in 1:ncol(counts)){
  srr.vector <- groups[groups[,2] == colnames(counts[i]),1]
  for(j in 1:length(srr.vector)){
    counts.subset <- read.delim(paste0(SALMON_OUTPUT.DIR, "/",srr.vector[j],"/quant.sf"))
    counts[,i] <- counts[,i] + counts.subset[match(rownames(counts),counts.subset[,1]),5]
  }
}

write.table(counts,
            paste0(WORKING.DIR,"/tick_counts.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

colSums(counts)
```

```{R, eval = F}
 PECHA_ISE6_Arkansas_mRNA1  PECHA_ISE6_Arkansas_mRNA4  PECHA_ISE6_Arkansas_mRNA5 PECHA_ISE6_Heartland_mRNA1 PECHA_ISE6_Heartland_mRNA2 
                17371256.1                   965294.2                   415261.1                   535804.0                   716451.0 
PECHA_ISE6_Heartland_mRNA3        PECHA_ISE6_HF_mRNA1        PECHA_ISE6_HF_mRNA2        PECHA_ISE6_HF_mRNA3       PECHA_ISE6_Jax_mRNA1 
                  613679.0                   592531.2                   426703.1                  1029629.3                   741714.0 
      PECHA_ISE6_Jax_mRNA2       PECHA_ISE6_Jax_mRNA3   PECHA_ISE6_Liberty_mRNA1   PECHA_ISE6_Liberty_mRNA2   PECHA_ISE6_Liberty_mRNA3 
                 1351246.5                  1105501.0                  1042668.0                  1006309.1                   535050.0 
          PECHA_ISE6_mRNA1           PECHA_ISE6_mRNA2           PECHA_ISE6_mRNA3   PECHA_ISE6_Osceola_mRNA1   PECHA_ISE6_Osceola_mRNA2 
                  363096.0                   594353.0                   533546.0                   618000.0                   619140.0 
  PECHA_ISE6_Osceola_mRNA3 PECHA_ISE6_StVincent_mRNA1 PECHA_ISE6_StVincent_mRNA2 PECHA_ISE6_StVincent_mRNA3   PECHA_ISE6_Wakulla_mRNA1 
                  618716.0                   942616.0                   468259.0                   698992.0                   985561.0 
  PECHA_ISE6_Wakulla_mRNA2   PECHA_ISE6_Wakulla_mRNA3 PECHA_ISE6_WestPaces_mRNA1 PECHA_ISE6_WestPaces_mRNA2 PECHA_ISE6_WestPaces_mRNA3 
                  622454.0                   564877.8                   913495.0                  1779219.9                  2108551.0
```

#### Create TPM data frame
```{R}
genelength <- read.delim(paste0(SALMON_OUTPUT.DIR, "/", groups[1,1],"/quant.sf"))
genelength <- genelength[match(rownames(counts),genelength[,1]),2]

tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
            paste0(WORKING.DIR,"/tick_tpm.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

dim(tpm)
```

#### Set group levels
```{R}
groups <- unique(groups[2:5])
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
```

#### Conduct saturation analysis
```{R, fig.height=5, fig.width=6}
rarefy.counts <- round(counts,0)
raremax <- round(min(rowSums(t(rarefy.counts))),0)
srare <- rarefy(t(rarefy.counts),raremax) 

rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)

rarefy.df <- as.data.frame(matrix(nrow = 0,
                                  ncol = 5))
rarefy.points.df <- rarefy.df
for(i in 1:length(rarefy.raw.df)){
  steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[i]])))
  detected_genes <- as.numeric(rarefy.raw.df[[i]])
  rarefy.df <- as.data.frame(rbind(rarefy.df,
                                   cbind(as.numeric(steps),as.numeric(detected_genes),as.character(groups[i,1]),as.character(groups[i,2]),groups[i,3])))
  rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
                                          cbind(as.numeric(max(steps)),as.numeric(max(detected_genes)),as.character(groups[i,1]),as.character(groups[i,2],groups[i,3]))))
  
}
rarefy.plot <- ggplot()+
  geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4]),size = 3)+
  guides(colour = F,shape = F)+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x="reads mapping to genes", y="genes detected", color = "Sample")+
  #coord_cartesian(xlim=c(0,100000))+
  theme_bw()


pdf(paste0(WORKING.DIR,"/plots/tick_rarefication_plot.pdf"),
    height=5,
    width=6)
print(rarefy.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_rarefication_plot.png"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

![image](/images/tick_rarefication_plot.png)

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = groups[,4])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,4]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,4])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
            paste0(WORKING.DIR,"/tick_counts_edgeR_longitudinal.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```


```{R, eval = F}
[1] "11969 genes excluded with CPM cutoff"
[1] "3862 DE genes identified using edgeR longitudinal"
```

#### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

##### Create sample legend for PCA and hierarchical clustering plots

```{R, fig.height = 2, fig.width = 10}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_point(aes(color = groups[,4],shape=groups[,4]), size = 4)+
    scale_shape_manual(values = as.numeric(as.character(groups[,3])))+
    scale_color_manual(values = as.character(unique(groups[,2])))+
    guides(shape = guide_legend(title = "Samples", title.position = "top",nrow=2),
           colour = guide_legend(title = "Samples", title.position = "top",nrow=2))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/tick_pca_hc_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_pca_hc_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/tick_pca_hc_samplekey.png)

##### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

```{R, fig.height=5, fig.width=10}
dendrogram <- as.data.frame(t(scale(t(log2(tpm.keep+1)))))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

structure <- get_dendro_structure(result)
dendro.data <- get_dendro_data(result)
bootstrap.positions <- get_dendro_bootstraps(dendro.data)
  
points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                 structure))
dendrogroups <- groups[,4][result$hclust$order]
dendrocol <- groups[,2][result$hclust$order]
dendroshape <- groups[,3][result$hclust$order]
dendrosize <- colSums(counts.keep)[result$hclust$order]
#dendrosize <- 1

dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
  geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize, shape = dendroshape))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
  guides(colour = guide_legend(ncol = 2), size = F, shape = F)+
  #scale_x_discrete(limits = as.character(dendrolabels))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

for(i in 1:length(result$edges$bp)){
  text <- round(result$edges$bp[i] * 100,0)
  dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.4, y=bootstrap.positions[i,2] + 0.04, size = 2)
}
pdf(paste0(WORKING.DIR,"/plots/tick_tpm_kept_dendrogram.pdf"),
    height=5,
    width=10)
print(dendrogram.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_tpm_kept_dendrogram.png"),
    height=5,
    width=10,
	units = "in",res=300)
print(dendrogram.plot)
dev.off()

print(dendrogram.plot)
```

![image](/images/tick_tpm_kept_dendrogram.png)
##### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

```{R,fig.height=5,fig.width=5}
pca.df <- t(scale(t(log2(tpm.keep + 1))))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.keep) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = dendrocol,size = dendrosize, shape = dendroshape))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
       x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
       y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  guides(color = F,size = F, shape = F)+
  # guides(colour = guide_legend(ncol = 2))+
  scale_color_manual(values = levels(groups[,2]))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/tick_tpm_kept_pca.pdf"),
    height=5,
    width=5)
print(pca.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_tpm_kept_pca.png"),
    height=5,
    width=5,
	units = "in",res=300)
print(pca.plot)
dev.off()

print(pca.plot)
```

![image](/images/tick_tpm_kept_pca.png)


#### Divide differentially expressed genes into expression modules

##### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)

```

![image](/images/tick_tpm_de_wgcna_soft_power_plot.png)

##### Identify expression modules 

```{R, fig.height=5, fig.width=12}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 12, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 12, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/tick_tpm_de_wgcna_merge_eigengenes_plot.png)

##### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
            paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_modules.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/tick_tpm_de_wgcna_eigengene_dendrogram.png)

##### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,2])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)
```

###### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_line(aes(color = groups[,4]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,2])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=3))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/tick_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/tick_hm_samplekey.png)

###### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
                       colours=hmcol,
                       breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/tick_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/tick_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/tick_hm_zscorelog2tpmkey.png)

###### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/tick_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

###### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/tick_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/tick_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)


### Ehrlichia

#### Set R inputs

```{R}
WORKING.DIR <- "Z:/EBMAL/mchung_dir/PECHA/"
FADU_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/fadu"
INTERPROSCAN_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/references"
PANOCT_OUTPUT.DIR <- "Z:/EBMAL/mchung_dir/PECHA/panoct"
GROUPS.PATH <- "Z:/EBMAL/mchung_dir/PECHA/pecha_groups.tsv"
```

#### Load R functions

```{R}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 

get_dendro_structure <- function(result){
  structure <- hang.dendrogram(as.dendrogram(result$hclust))
  structure <- capture.output(str(structure))
  structure <- structure[grepl("leaf", structure)]
  structure <- as.numeric(as.character(substr(structure, regexpr("h=", structure ) + 3, regexpr("  )", structure))))
  return(structure)
}

get_dendro_data <- function(result){
  dendro.data <- dendro_data(result$hclust)
  dendro.data <- dendro.data$segments[which(dendro.data$segments$y == dendro.data$segments$yend),]
  for(i in 1:nrow(dendro.data)){
    dendro.data$minx[i] <- min(c(dendro.data$x[i], dendro.data$xend[i]))
  }
  dendro.data <- dendro.data[order(as.numeric(as.character(dendro.data$y)), as.numeric(as.character(dendro.data$minx))),]
  return(dendro.data)
}

get_dendro_bootstraps <- function(dendro_data){
  bootstrap.positions <- as.data.frame(matrix(nrow = length(dendro_data$y[duplicated(dendro_data$y)]),
                                              ncol = 2))
  for(i in 1:length(dendro_data$y[duplicated(dendro_data$y)])){
    dendro_data.subset <- dendro_data[which(dendro_data$y == dendro_data$y[duplicated(dendro_data$y)][i]),]
    bootstrap.positions[i,1] <- unique(dendro_data.subset$x)
    bootstrap.positions[i,2] <- unique(dendro_data.subset$y)
  }
  return(bootstrap.positions)
}

find_soft_power <- function(sft){
  df <- as.data.frame(cbind(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]))
  y <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  dy <- diff(y) 
  softpower <- which(abs(dy) < 0.05)[1]
  if(softpower == 1){
    softpower <- which(abs(dy) < 0.05)[2]
  }
  return(softpower)
}

eigengene_invert_id <- function(tpm.de, mergedColors, mergedMEs){
  tpm.de.wgcna <- tpm.de
  tpm.de.wgcna$invert <- T
  tpm.de.wgcna$module <- mergedColors
  for(i in 1:nrow(tpm.de.wgcna)){
    if(cor(t(tpm.de[i,]), mergedMEs[,which(colnames(mergedMEs) == paste0("ME",tpm.de.wgcna$module[i]))], method = "pearson") > 0){
      tpm.de.wgcna$invert[i] <- F
    }
  }
  return(tpm.de.wgcna)
}

wgcna_heatmap_reorder <- function(tpm.de.wgcna){
  clusters <- as.data.frame(table(tpm.de.wgcna$module))
  clusters <- clusters[order(-clusters[,2]),1]
  
  tpm.de.wgcna.reordered <- as.data.frame(matrix(nrow = 0,
                                                 ncol = ncol(tpm.de.wgcna)))
  for(i in 1:length(clusters)){
    tpm.de.wgcna.reordered <- as.data.frame(rbind(tpm.de.wgcna.reordered,
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == F,],
                                                  tpm.de.wgcna[tpm.de.wgcna$module == clusters[i] & tpm.de.wgcna$invert == T,]))
  }
  return(tpm.de.wgcna.reordered)
}

get_heatmap_separators <- function(vector){
  sep <- c()
  for(i in 2:length(unique(vector))){
    sep[length(sep) + 1] <- min(which(vector == unique(vector)[i])) - 1
  }
  return(sep)
}

functionaltermenrichment <- function(genes, geneinfo){
  for(i in 1:ncol(geneinfo)){geneinfo[,i] <- as.character(geneinfo[,i])}
  geneinfo$interpro_description[which(is.na(geneinfo$interpro_description))] <- "No InterPro entry"
  geneinfo$go_biologicalprocess[which(is.na(geneinfo$go_biologicalprocess))] <- "No GO terms for biological process"
  geneinfo$go_cellularcomponent[which(is.na(geneinfo$go_cellularcomponent))] <- "No GO terms for cellular component"
  geneinfo$go_molecularfunction[which(is.na(geneinfo$go_molecularfunction))] <- "No GO terms for molecular function"
  
  functionalterms.list <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo$interpro_description, collapse = "|"),  split = "[|]")))),
                               gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                               gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                               gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  geneinfo.subset <- geneinfo[geneinfo$gene %in% genes,]
  term <- c()
  clusteroccurences <- c()
  genomeoccurences <- c()
  pvalue <- c()
  correctedpvalue <- c()
  oddsratio <- c()

  functionalterms.list.subset <- list(ipr=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$interpro_description, collapse = "|"),  split = "[|]")))),
                                      gobio=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_biologicalprocess, collapse = "|"),  split = "[|]")))),
                                      gocell=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_cellularcomponent, collapse = "|"),  split = "[|]")))),
                                      gomol=as.data.frame(table(unlist(strsplit(paste(geneinfo.subset$go_molecularfunction, collapse = "|"),  split = "[|]")))))
  
  for(i in 1:length(functionalterms.list)){
    for(j in 1:nrow(functionalterms.list[[i]])){
      freq.all <- functionalterms.list[[i]][j,2]
      freq.subset <- ifelse(functionalterms.list[[i]][j,1] %in% functionalterms.list.subset[[i]][,1],
                            functionalterms.list.subset[[i]][functionalterms.list.subset[[i]][,1] == as.character(functionalterms.list[[i]][j,1]),2],
                            0)
      genes.all <- nrow(geneinfo)
      genes.subset <- nrow(geneinfo.subset)

      fisherexact.matrix <- matrix(c(freq.subset, freq.all - freq.subset,
                                     genes.subset - freq.subset, genes.all - genes.subset - freq.all + freq.subset),
                                   nrow = 2,
                                   ncol = 2)
      fisher.test <- fisher.test(fisherexact.matrix)
      
      term[length(term) + 1] <- as.character(functionalterms.list[[i]][j,1])
      clusteroccurences[length(clusteroccurences) + 1] <- as.numeric(as.character(freq.subset))
      genomeoccurences[length(genomeoccurences) + 1] <- as.numeric(as.character(freq.all))
      pvalue[length(pvalue) + 1] <- as.numeric(as.character(fisher.test$p.value))
      correctedpvalue[length(correctedpvalue) + 1] <- p.adjust(as.numeric(as.character(fisher.test$p.value)), method = "fdr", n = nrow(functionalterms.list[[i]]))
      oddsratio[length(oddsratio) + 1] <- as.numeric(as.character(fisher.test$estimate))
    }
  }
  
  terms.df <- as.data.frame(cbind(term,
                                  clusteroccurences,
                                  genomeoccurences,
                                  pvalue,
                                  correctedpvalue,
                                  oddsratio))
  terms.df <- terms.df[order(as.numeric(as.character(terms.df$pvalue))),]
  return(terms.df)
}
```

#### Load packages and view sessionInfo

```{R}
library(dendextend)
library(DESeq2)
library(edgeR)
library(FactoMineR)
library(ggdendro)
library(ggplot2)
library(gplots)
library(gridExtra)
library(pvclust)
library(vegan)
library(WGCNA)

sessionInfo()
```

```{R, eval = F}
R version 3.5.0 (2018-04-23)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] WGCNA_1.68                  fastcluster_1.1.25          dynamicTreeCut_1.63-1       pvclust_2.0-0              
 [5] gplots_3.0.1.1              ggdendro_0.1-20             FactoMineR_1.42             DESeq2_1.22.2              
 [9] SummarizedExperiment_1.12.0 DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0         
[13] Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2         IRanges_2.16.0             
[17] S4Vectors_0.20.1            BiocGenerics_0.28.0         dendextend_1.12.0           gridExtra_2.3              
[21] ggplot2_3.2.0               vegan_2.5-5                 lattice_0.20-35             permute_0.9-5              
[25] edgeR_3.24.3                limma_3.38.3               

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       htmlTable_1.13.1       XVector_0.22.0         base64enc_0.1-3       
 [5] rstudioapi_0.10        bit64_0.9-7            mvtnorm_1.0-11         AnnotationDbi_1.44.0  
 [9] codetools_0.2-15       splines_3.5.0          leaps_3.0              doParallel_1.0.15     
[13] impute_1.56.0          robustbase_0.93-5      geneplotter_1.60.0     knitr_1.23            
[17] zeallot_0.1.0          Formula_1.2-3          annotate_1.60.1        cluster_2.0.7-1       
[21] GO.db_3.7.0            rrcov_1.4-7            compiler_3.5.0         backports_1.1.4       
[25] assertthat_0.2.1       Matrix_1.2-14          lazyeval_0.2.2         acepack_1.4.1         
[29] htmltools_0.3.6        tools_3.5.0            gtable_0.3.0           glue_1.3.1            
[33] GenomeInfoDbData_1.2.0 dplyr_0.8.3            Rcpp_1.0.2             vctrs_0.2.0           
[37] gdata_2.18.0           preprocessCore_1.44.0  nlme_3.1-137           iterators_1.0.12      
[41] xfun_0.8               stringr_1.4.0          gtools_3.8.1           XML_3.98-1.20         
[45] DEoptimR_1.0-8         zlibbioc_1.28.0        MASS_7.3-51.4          scales_1.0.0          
[49] RColorBrewer_1.1-2     yaml_2.2.0             memoise_1.1.0          rpart_4.1-13          
[53] latticeExtra_0.6-28    stringi_1.4.3          RSQLite_2.1.2          genefilter_1.64.0     
[57] pcaPP_1.9-73           foreach_1.4.7          checkmate_1.9.4        caTools_1.17.1.2      
[61] rlang_0.4.0            pkgconfig_2.0.2        bitops_1.0-6           purrr_0.3.2           
[65] htmlwidgets_1.3        labeling_0.3           bit_1.1-14             tidyselect_0.2.5      
[69] robust_0.4-18.1        magrittr_1.5           R6_2.4.0               Hmisc_4.2-0           
[73] fit.models_0.5-14      DBI_1.0.0              pillar_1.4.2           foreign_0.8-70        
[77] withr_2.1.2            mgcv_1.8-23            survival_2.41-3        scatterplot3d_0.3-41  
[81] RCurl_1.95-4.12        nnet_7.3-12            tibble_2.1.3           crayon_1.3.4          
[85] KernSmooth_2.23-15     viridis_0.5.1          locfit_1.5-9.1         grid_3.5.0            
[89] data.table_1.12.2      blob_1.2.0             digest_0.6.20          flashClust_1.01-2     
[93] xtable_1.8-4           munsell_0.5.0          viridisLite_0.3.0     
```

#### Constructs core genome table that consists of only gene clusters consisting of one gene per strain

```{R}
ortho_table <- read.delim(paste0(PANOCT_OUTPUT.DIR,"/matchtable.txt"),header = F,row.names = 1)
ortho_table[ortho_table == "----------"] <- NA
rownames(ortho_table) <- paste0("PANOCT_",rownames(ortho_table))
colnames(ortho_table) <- c("Arkansas","Heartland","HF","Jacksonville","Liberty","Osceola","St. Vincent","Wakulla","West Paces")
ortho_table <- ortho_table[rowSums(!is.na(ortho_table)) == ncol(ortho_table),]

dim(ortho_table)
```

```{R, eval = F}
[1] 795   9
```

#### Create counts data frame

```{R}
groups <- read.delim(GROUPS.PATH, header = F)
groups <- groups[grep("totalRNA", groups[,2]),]
groups <- groups[intersect(grep("DH82_totalRNA",groups[,2],invert = T),grep("ISE6_totalRNA",groups[,2],invert = T)),]

rownames <- rownames(ortho_table)
colnames <- unique(groups[,2])

counts <- as.data.frame(matrix(0,
                               nrow = length(rownames),
                               ncol = length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

for(i in 1:ncol(counts)){
  srr.vector <- groups[groups[,2] == colnames(counts[i]),1]
  for(j in 1:length(srr.vector)){
    counts.subset <- read.delim(paste0(FADU_OUTPUT.DIR, "/",srr.vector[j],".sortedbyposition.counts.txt"))
    counts[,i] <- counts[,i] + counts.subset[match(ortho_table[,colnames(ortho_table) == groups[groups[,1] == srr.vector[j],5]],counts.subset[,1]),4]
  }
}

write.table(counts,
            paste0(WORKING.DIR,"/ehrlichia_counts.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

sort(colSums(counts))
```

```{R, eval = F}
     PECHA_Heartland_totalRNA1 PECHA_ISE6_StVincent_totalRNA1 PECHA_ISE6_StVincent_totalRNA2 
                         10.00                          11.00                          47.07 
PECHA_ISE6_StVincent_totalRNA3   PECHA_ISE6_Liberty_totalRNA3       PECHA_ISE6_Jax_totalRNA3 
                         52.11                          60.32                          93.90 
       PECHA_Liberty_totalRNA1   PECHA_ISE6_Wakulla_totalRNA1   PECHA_ISE6_Wakulla_totalRNA2 
                        113.33                         114.41                         117.64 
  PECHA_ISE6_Liberty_totalRNA2        PECHA_Osceola_totalRNA1            PECHA_Jax_totalRNA2 
                        326.45                         363.98                         375.77 
       PECHA_ISE6_HF_totalRNA2 PECHA_ISE6_WestPaces_totalRNA1       PECHA_ISE6_Jax_totalRNA2 
                        794.11                         798.03                         868.50 
            PECHA_HF_totalRNA3       PECHA_ISE6_Jax_totalRNA1 PECHA_ISE6_WestPaces_totalRNA2 
                        931.12                         935.71                         958.66 
  PECHA_ISE6_Liberty_totalRNA1   PECHA_ISE6_Wakulla_totalRNA3      PECHA_WestPaces_totalRNA1 
                        978.41                        1107.44                        1352.07 
     PECHA_StVincent_totalRNA1 PECHA_ISE6_WestPaces_totalRNA3  PECHA_ISE6_Arkansas_totalRNA1 
                       1520.59                        1580.12                        1587.59 
            PECHA_HF_totalRNA1      PECHA_Heartland_totalRNA2        PECHA_Osceola_totalRNA3 
                       1683.63                        1789.97                        1840.64 
     PECHA_Heartland_totalRNA3             PECHA_HF_totalRNA2        PECHA_Wakulla_totalRNA2 
                       1927.69                        2000.25                        2006.11 
       PECHA_Liberty_totalRNA3  PECHA_ISE6_Arkansas_totalRNA5        PECHA_Wakulla_totalRNA3 
                       2410.99                        2757.66                        3272.84 
     PECHA_StVincent_totalRNA2        PECHA_ISE6_HF_totalRNA3        PECHA_Liberty_totalRNA2 
                       3440.50                        4446.62                        4836.25 
       PECHA_Wakulla_totalRNA1            PECHA_Jax_totalRNA3        PECHA_Osceola_totalRNA2 
                       5078.48                        5428.62                        5583.33 
       PECHA_ISE6_HF_totalRNA1      PECHA_WestPaces_totalRNA3      PECHA_WestPaces_totalRNA2 
                       5678.74                        5994.43                        6079.31 
           PECHA_Jax_totalRNA1       PECHA_Arkansas_totalRNA3  PECHA_ISE6_Arkansas_totalRNA4 
                       6599.09                       16390.22                       22221.51 
      PECHA_Arkansas_totalRNA1      PECHA_StVincent_totalRNA3       PECHA_Arkansas_totalRNA2 
                      24960.98                       27061.62                       27902.70 
```

#### Calculate average gene length for all core genes

```{R}
genelength.ortho_table <- ortho_table
for(i in 1:ncol(ortho_table)){
  counts.subset <- read.delim(paste0(FADU_OUTPUT.DIR, "/",groups[groups[,5] == colnames(ortho_table)[i],1][1],".sortedbyposition.counts.txt"))
  genelength.ortho_table[,i] <- counts.subset[match(ortho_table[,i],counts.subset[,1]),2]
}
genelength <- rowMeans(genelength.ortho_table)
```

#### Create TPM data frame

```{R}
tpm <- counts
for(i in 1:ncol(tpm)){
  tpm[,i] <- tpm[,i]/genelength
  tpm[,i] <- tpm[,i]/(sum(tpm[,i])/1000000)
}

write.table(tpm,
            paste0(WORKING.DIR,"/ehrlichia_tpm.tsv"),
            quote = F,
            col.names = T,
            row.names = T,
            sep = "\t")

dim(tpm)
```

```{R, eval = F}
[1] 795  54
```

#### Set group levels

```{R}
groups <- unique(groups[2:5])
groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))
```

#### Conduct saturation analysis

```{R, fig.height=5, fig.width=6}
counts <- counts[,colSums(counts) != 0]

rarefy.counts <- round(counts,0)
raremax <- round(min(rowSums(t(rarefy.counts))),0)
srare <- rarefy(t(rarefy.counts),raremax) 

rarefy.raw.df <- rarecurve(t(rarefy.counts), step = round(raremax/10,0), sample = raremax)

rarefy.df <- as.data.frame(matrix(nrow = 0,
                                  ncol = 5))
rarefy.points.df <- rarefy.df
for(i in 1:length(rarefy.raw.df)){
  steps <- as.numeric(gsub("N","",names(rarefy.raw.df[[i]])))
  detected_genes <- as.numeric(rarefy.raw.df[[i]])
  rarefy.df <- as.data.frame(rbind(rarefy.df,
                                    cbind(as.numeric(steps),as.numeric(detected_genes),as.character(groups[i,1]),as.character(groups[i,2]),groups[i,3])))
  rarefy.points.df <- as.data.frame(rbind(rarefy.points.df,
                                          cbind(as.numeric(max(steps)),                                                                                                     as.numeric(max(detected_genes)),                                                                                            as.character(groups[i,1]),
                                                as.character(groups[i,2]),
                                                as.character(groups[i,3]))))
  
}
rarefy.plot <- ggplot()+
  geom_line(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  #geom_point(mapping=aes(x=as.numeric(as.character(rarefy.df[,1])), y=as.numeric(as.character(rarefy.df[,2])),group=rarefy.df[,3],color=rarefy.df[,4]))+
  geom_point(mapping=aes(x=as.numeric(as.character(rarefy.points.df[,1])), y=as.numeric(as.character(rarefy.points.df[,2])),group=rarefy.points.df[,3],color=rarefy.points.df[,4],shape=rarefy.points.df[,5]),size = 3)+
  guides(colour = F,shape = F)+
  scale_shape_manual(values = as.numeric(as.character(levels(groups[,3]))))+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x="reads mapping to genes", y="genes detected", color = "Sample")+
  #coord_cartesian(xlim=c(0,100000))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_rarefication_plot.pdf"),
    height=5,
    width=6)
print(rarefy.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_rarefication_plot.png"),
    height=5,
    width=6,
	units = "in",res=300)
print(rarefy.plot)
dev.off()

print(rarefy.plot)
```

![image](/images/ehrlichia_rarefication_plot.png)

#### Exclude low count samples and samples with only 1 replicate

From the rarefaction analysis, these filters were used to exclude samples: (<4k gene counts, approximately 5x the number of core Ehrlichia genes).

```{R}
include.samples <- colnames(counts)[colSums(counts) > 4000]
groups <- groups[groups[,1] %in% include.samples,]

groups[,4] <- as.character(groups[,4])
for(i in 1:nrow(groups)){
  if(grepl("ISE6",groups[i,1])){
    groups[i,4] <- paste0(as.character(groups[i,4]),", Tick")
  }else{
    groups[i,4] <- paste0(as.character(groups[i,4]),", Canine")
  }
}

exclude.samples <- names(table(groups[,4])[table(groups[,4]) == 1])
groups <- groups[!(groups[,4] %in% exclude.samples),]

groups[,1] <- factor(groups[,1], levels = groups[,1])
groups[,2] <- factor(groups[,2], levels = unique(groups[,2]))
groups[,3] <- factor(groups[,3], levels = unique(groups[,3]))
groups[,4] <- factor(groups[,4], levels = unique(groups[,4]))

counts <- counts[,colnames(counts) %in% groups[,1]]
tpm <- tpm[,colnames(tpm) %in% groups[,1]]

dim(tpm)
```

```{R, eval = F}
[1] 795   9
```

#### Identify differentially expressed genes longitudinally

edgeR and DESeq2 are both run with a FDR cutoff of <0.05 and a minimum CPM cutoff of 5 reads in the lowest sequenced sample in the data set.  

```{R}
FDRcutoff <- 0.05
cpm.cutoff <- 5/min(colSums(counts)) * 1000000

y <- DGEList(counts = counts, group = groups[,4])
y <- calcNormFactors(y)
keep <- rowSums(cpm(y) >= cpm.cutoff) >= min(table(groups[,4]))
keep.df <- as.data.frame(table(keep))
print(paste0(keep.df[keep.df[,1] == F,2]," genes excluded with CPM cutoff"))

y <- y[keep, , keep.lib.sizes = F]
design <- model.matrix(~groups[,4])
y <- estimateDisp(y , design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2:ncol(fit))
qlf$table$padj <- p.adjust(qlf$table$PValue, method="BH")
edgeR.longitudinal.degenes <- qlf$table[qlf$table$padj < FDRcutoff,]
print(paste0(nrow(edgeR.longitudinal.degenes)," DE genes identified using edgeR longitudinal"))

write.table(edgeR.longitudinal.degenes,
            paste0(WORKING.DIR,"/ehrlichia_counts_edgeR_longitudinal.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

counts.keep <- counts[keep,]
tpm.keep <- tpm[keep,]
```


```{R, eval = F}
[1] "484 genes excluded with CPM cutoff"
[1] "163 DE genes identified using edgeR longitudinal"
```

#### Conduct PCA and hierarchical clustering analyses on genes that passed the CPM cutoff

##### Create sample legend for PCA and hierarchical clustering plots

```{R, fig.height = 2, fig.width = 10}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_point(aes(color = groups[,4],shape=groups[,3]), size = 4)+
    #scale_shape_manual(values = levels(groups[,3]))+
    scale_color_manual(values = levels(groups[,2]))+
    guides(shape = guide_legend(title = "Samples", title.position = "top",nrow=2),
           colour = guide_legend(title = "Samples", title.position = "top",nrow=2))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_pca_hc_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_pca_hc_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/ehrlichia_pca_hc_samplekey.png)

##### Conduct a hierarchical cluster analysis on the TPM values of all genes that passed CPM cutoff

```{R, fig.height=5, fig.width=10}
dendrogram <- as.data.frame(t(scale(t(log2(tpm.keep+1)))))

result <- pvclust(dendrogram, method.dist="cor", method.hclust="average", nboot=100)

structure <- get_dendro_structure(result)
dendro.data <- get_dendro_data(result)
bootstrap.positions <- get_dendro_bootstraps(dendro.data)
  
points.df <- as.data.frame(cbind(seq(1,length(structure),1),
                                 structure))
dendrogroups <- groups[,4][result$hclust$order]
dendrocol <- groups[,2][result$hclust$order]
dendroshape <- groups[,3][result$hclust$order]
dendrosize <- colSums(counts.keep)[result$hclust$order]
#dendrosize <- 1

dendrogram.plot <- ggdendrogram(hang.dendrogram(as.dendrogram(result$hclust)), theme_dendro = T)+
  geom_point(aes(x=seq(1,length(structure)), y = structure, color = dendrogroups, size = dendrosize, shape = dendroshape))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  scale_color_manual(values = levels(groups[,2]))+
  labs(x = "", y = "", col = "Samples", size = "Reads Mapped\nto Features")+
  guides(colour = guide_legend(ncol = 2), size = F, shape = F)+
  #scale_x_discrete(limits = as.character(dendrolabels))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

for(i in 1:length(result$edges$bp)){
  text <- round(result$edges$bp[i] * 100,0)
  dendrogram.plot <- dendrogram.plot + annotate("text", label = text, x=bootstrap.positions[i,1] + 0.4, y=bootstrap.positions[i,2] + 0.04, size = 2)
}
pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_kept_dendrogram.pdf"),
    height=5,
    width=10)
print(dendrogram.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_kept_dendrogram.png"),
    height=5,
    width=10,
	units = "in",res=300)
print(dendrogram.plot)
dev.off()

print(dendrogram.plot)
```

![image](/images/ehrlichia_tpm_kept_dendrogram.png)
##### Conduct a PCA on the TPM values of all genes that passed CPM cutoff

```{R,fig.height=5,fig.width=5}
pca.df <- t(scale(t(log2(tpm.keep + 1))))
pca.df <- pca.df[rowSums(pca.df == 0) != ncol(pca.df),]
pca <- PCA(as.data.frame(scale(t(pca.df))), graph = FALSE, ncp = ncol(tpm.keep) - 1)

pca.plot <- ggplot()+
  geom_point(aes(x=pca$ind$coord[,1], y=pca$ind$coord[,2], color = dendrocol,size = dendrosize, shape = dendroshape))+
  labs(col = "Samples", size = "Reads Mapped\nto Features", 
       x = paste0("PC1 (", round(pca$eig[1,2],1), "%)"), 
       y = paste0("PC2 (", round(pca$eig[2,2],1), "%)"))+
  guides(color = F,size = F, shape = F)+
  # guides(colour = guide_legend(ncol = 2))+
  scale_color_manual(values = levels(groups[,2]))+
  scale_shape_manual(values= as.numeric(as.character(levels(dendroshape))))+
  theme_bw()

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_kept_pca.pdf"),
    height=5,
    width=5)
print(pca.plot)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_kept_pca.png"),
    height=5,
    width=5,
	units = "in",res=300)
print(pca.plot)
dev.off()

print(pca.plot)
```

![image](/images/ehrlichia_tpm_kept_pca.png)


#### Divide differentially expressed genes into expression modules

##### Find soft power value for WGCNA

```{R}
tpm.de <- tpm[rownames(tpm) %in% rownames(edgeR.longitudinal.degenes),]

wgcna <- as.data.frame(t(tpm.de))
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(wgcna, powerVector = powers, verbose = 5)
softpower <- find_soft_power(sft)

text.color <- rep("black",length(sft$fitIndices[,1]))
text.color[which(sft$fitIndices[,1] == softpower)] <- "red"

scale_independence.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Scale Independence", x = "soft threshold (power)", y = "scale free topology model fit, signed R^2")+
  theme_bw()

mean_connectivity.plot <- ggplot()+
  geom_text(mapping = aes(x = sft$fitIndices[,1], y = sft$fitIndices[,5], label=sft$fitIndices[,1]), color = text.color)+
  labs(title = "Mean Connectivity", x = "soft threshold (power)", y = "mean connectivity")+
  theme_bw()

wgcna_soft_power_plots <- list(scale_independence.plot, mean_connectivity.plot)
lay <- rbind(c(1,2))

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_soft_power_plot.pdf"),
    width = 10, 
    height = 5)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_soft_power_plot.png"),
    width = 10, 
    height = 5,
	units = "in",res=300)
grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)
dev.off()

grid.arrange(grobs = wgcna_soft_power_plots,
             widths = c(5,5),
             heights = c(5),
             layout_matrix = lay)

```

![image](/images/ehrlichia_tpm_de_wgcna_soft_power_plot.png)

##### Identify expression modules 

```{R, fig.height=5, fig.width=12}
adjacency <- adjacency(wgcna, power = softpower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM
geneTree <- hclust(as.dist(dissTOM), method = "average");

minModuleSize <- 1
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(wgcna, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
METree = hclust(as.dist(MEDiss), method = "average")

MEDissThres = 0.25

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_merge_eigengenes_plot.pdf"),
    width = 12, 
    height = 5)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_merge_eigengenes_plot.png"),
    width = 12, 
    height = 5,
    units = "in",res=300)
plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
dev.off()


plot(METree, main = "Clustering of Module Eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
```

![image](/images/ehrlichia_tpm_de_wgcna_merge_eigengenes_plot.png)

##### Merge similar expression modules

```{R}
merge = mergeCloseModules(wgcna, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

tpm.de.wgcna <- eigengene_invert_id(tpm.de, mergedColors, mergedMEs)

write.table(tpm.de.wgcna,
            paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_modules.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_eigengene_dendrogram.pdf"),
    width = 8, 
    height = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_eigengene_dendrogram.png"),
    width = 8, 
    height = 5,
    units = "in",res=300)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("dynamic tree cut", "merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
```

![image](/images/ehrlichia_tpm_de_wgcna_eigengene_dendrogram.png)

##### Plot WGCNA expression modules as a heatmap

```{R}
tpm.de.wgcna <- wgcna_heatmap_reorder(tpm.de.wgcna)

log2tpm.de <- log2(tpm.de.wgcna[,1:(ncol(tpm.de.wgcna) - 2)] + 1)
zscore.log2tpm.de <- as.data.frame(t(scale(t(log2tpm.de))))

hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
rowcol1 <- tpm.de.wgcna$module
rowcol2 <- unlist(lapply(tpm.de.wgcna$invert,function(x){if(x == F){return("grey")}else{return("black")}}))
colcol <- as.character(groups[,2])

rowsep <- get_heatmap_separators(rowcol1)
colsep <- get_heatmap_separators(colcol)
```

###### Create sample legend for WGCNA heatmap

```{R}
legend.plot <- ggplot(mapping=aes(x=groups[,1], y=seq(1,length(groups[,1]),1), group = groups[,4]))+
    geom_line(aes(color = groups[,4]), size = 4)+
    scale_color_manual(values = as.character(unique(groups[,2])))+
    guides(colour = guide_legend(title = "Samples", title.position = "top",nrow=3))+
    theme_bw()+
    theme(legend.position="top",legend.title.align=0.5)

sample.legend <- g_legend(legend.plot)

pdf(paste0(WORKING.DIR,"/plots/ehrlichia_hm_samplekey.pdf"),
    height=2,
    width=10)
grid.arrange(sample.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_hm_samplekey.png"),
    height=2,
    width=10,
    units = "in",res=300)
grid.arrange(sample.legend)
dev.off()

grid.arrange(sample.legend)
```

![image](/images/ehrlichia_hm_samplekey.png)

###### Create z-score log2TPM legend for WGCNA heatmap

```{R, fig,height = 2, fig.width = 7}
hmcol <- colorRampPalette(c("navyblue","white","firebrick3"))(12)
hmcolor.plot <- ggplot() + 
  geom_raster(aes(x=seq(-3,3,0.5), y=seq(-3,3,0.5), fill = seq(-3,3,0.5)))+
  scale_fill_gradientn(name = "z-score log2TPM",
                       colours=hmcol,
                       breaks=c(-3,0,3))+
  theme(legend.position="bottom")+
  guides(fill = guide_colorbar(title.position = "top"))
  

heat.legend <- g_legend(hmcolor.plot)
pdf(paste0(WORKING.DIR,"/plots/ehrlichia_hm_zscorelog2tpmkey.pdf"),
    height=2,
    width=7)
grid.arrange(heat.legend)
dev.off()

png(paste0(WORKING.DIR,"/plots/ehrlichia_hm_zscorelog2tpmkey.png"),
    height=2,
    width=7,
    units = "in",res=300)
grid.arrange(heat.legend)
dev.off()

grid.arrange(heat.legend)
```

![image](/images/ehrlichia_hm_zscorelog2tpmkey.png)

###### Use module assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_zscorelog2tpm_module_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol1,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/ehrlichia_tpm_de_wgcna_zscorelog2tpm_module_heatmap.png)

###### Use inverse assigments as the row color bar

```{R, fig.height=10, fig.width=5}
pdf(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.pdf"),
    width = 5, 
    height = 10)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()
png(paste0(WORKING.DIR,"/plots/ehrlichia_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png"),
    width = 5, 
    height = 10,
    units = "in",res=300)
heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
dev.off()

heatmap.2(as.matrix(zscore.log2tpm.de),
              col=hmcol,
              trace="none",
              labRow=vector(mode = "character", length = nrow(zscore.log2tpm.de)),
              Rowv = F,
              Colv = F,
              RowSideColors=rowcol2,
              ColSideColors=colcol,
              lhei = c(2,8),
              breaks = seq(-3,3,by=.5),
              rowsep = rowsep,
              colsep = colsep,
              dendrogram = "none")
```

![image](/images/ehrlichia_tpm_de_wgcna_zscorelog2tpm_inverse_heatmap.png)

##### Construct geneinfo for Ehrlichia core genome

```{R}
geneinfo <- as.data.frame(matrix(nrow=nrow(ortho_table),
                                 ncol=5))
colnames(geneinfo) <- c("gene","interpro_description","go_biologicalprocess","go_cellularcomponent","go_molecularfunction")
geneinfo[,1] <- rownames(ortho_table)

for(i in 1:ncol(ortho_table)){
  if(colnames(ortho_table)[i] == "Jacksonville"){
    geneinfo.subset <- read.delim(paste0(INTERPROSCAN_OUTPUT.DIR,"/Jax.cds.fna.interproscan.geneinfo.tsv"))
  }else if(colnames(ortho_table)[i] == "St. Vincent"){
    geneinfo.subset <- read.delim(paste0(INTERPROSCAN_OUTPUT.DIR,"/StVincent.cds.fna.interproscan.geneinfo.tsv"))
  }else if(colnames(ortho_table)[i] == "West Paces"){
    geneinfo.subset <- read.delim(paste0(INTERPROSCAN_OUTPUT.DIR,"/WestPaces.cds.fna.interproscan.geneinfo.tsv"))
  }else{
    geneinfo.subset <- read.delim(paste0(INTERPROSCAN_OUTPUT.DIR,"/",colnames(ortho_table)[i],".cds.fna.interproscan.geneinfo.tsv"))
  }
  
  geneinfo[,2] <- paste0(geneinfo[,2],"|",geneinfo.subset[match(ortho_table[,i],geneinfo.subset[,1]),2])
  geneinfo[,3] <- paste0(geneinfo[,3],"|",geneinfo.subset[match(ortho_table[,i],geneinfo.subset[,1]),3])
  geneinfo[,4] <- paste0(geneinfo[,4],"|",geneinfo.subset[match(ortho_table[,i],geneinfo.subset[,1]),4])
  geneinfo[,5] <- paste0(geneinfo[,5],"|",geneinfo.subset[match(ortho_table[,i],geneinfo.subset[,1]),5])
}

geneinfo <- apply(geneinfo,c(1,2),function(x){return(paste(unique(unlist(strsplit(x,split="[|]"))),collapse="|"))})
geneinfo <- gsub("^NA[|]","",geneinfo)
geneinfo[geneinfo == "NA"] <- NA

geneinfo <- as.data.frame(geneinfo)
```

#### Test each module partition for over-represented functional terms

```{R}
terms.colnames <- c("term","clusteroccurences","genomeoccurences","pvalue","correctedpvalue","oddsratio","module","invert")

terms.wgcna <- as.data.frame(matrix(nrow = 0,
                                     ncol = 8))
colnames(terms.wgcna) <- terms.colnames
for(j in 1:length(unique(tpm.de.wgcna$module))){
  terms.wgcna.f <- as.data.frame(cbind(functionaltermenrichment(rownames(tpm.de.wgcna)[tpm.de.wgcna$module == unique(tpm.de.wgcna$module)[j] & tpm.de.wgcna$invert == F],geneinfo),unique(tpm.de.wgcna$module)[j],F))
  terms.wgcna.t <- as.data.frame(cbind(functionaltermenrichment(rownames(tpm.de.wgcna)[tpm.de.wgcna$module == unique(tpm.de.wgcna$module)[j] & tpm.de.wgcna$invert == T],geneinfo),unique(tpm.de.wgcna$module)[j],T))
  
  colnames(terms.wgcna.f) <- terms.colnames
  colnames(terms.wgcna.t) <- terms.colnames
  terms.wgcna <- as.data.frame(rbind(terms.wgcna,
                                     terms.wgcna.f, 
                                     terms.wgcna.t))
}
write.table(terms.wgcna,
            paste0(WORKING.DIR,"/ehrlichia_functionalterms_wgcna.tsv"),
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")
```
