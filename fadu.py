#!/usr/bin/env python

"""
FADU.py - Feature Aggregate Depth Utility
Description - Generate counts of reads that map to non-overlapping portions of genes
Requires: Python 3, Pysam version 0.12.0.1
By: Shaun Adkins (sadkins@som.umaryland.edu)
	Matthew Chung (mattchung@umaryland.edu)

"""
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from argparse import ArgumentParser
from collections import Counter
import multiprocessing as mp
import os
import pysam
import re
import sys
import logging

###########
# Functions (in alphabetical order)
###########

def calc_avg_read_len(bam):
    """ Calculates average read len of all BAM reads """
    name = mp.current_process().name
    logging.debug("{} - Calculating average read length ...".format(name) )
    total_query_len = sum( read.query_length for read in bam.fetch() )
    return round(total_query_len / bam.count())

def calc_depth(depth_dict, bam, strand):
    """ Calculate depth of coverage using 'samtools deptha' """
    name = mp.current_process().name
    logging.debug("{} - Calculating depth of coverage of BAM ...".format(name))
    depth_out = pysam.depth("-aa", "-m", "10000000", bam).splitlines(True)
    depth_out_file = re.sub(r'\.bam', '.depth', bam)
    with open(depth_out_file, 'w') as f:
        for line in depth_out:
            f.write(line)
    store_depth(depth_dict, depth_out, strand)

def calc_readcounts_per_gene(contig_bases, gene_info, depth_dict, out_bam, read_len):
    """ Calculate number of reads that map to the uniq bases of a gene """
    name = mp.current_process().name
    logging.debug("{} - Calculate readcounts per gene for BAM ...".format(name))
    counts_file = re.sub(r'\.bam', '.counts', out_bam)
    with open(counts_file, 'w') as f:
        for gene, fields in gene_info.items():
            (contig, start, stop, strand) = fields
            strand = set_strand(strand)
            # Collect all uniq base positions for this gene, and total the depth for each position
            uniq_coords = [key for key, val in contig_bases[contig][strand].items() if val == gene]
            total_depth = 0
            for coord in uniq_coords:
                assert(coord in depth_dict[contig]), "Coordinate {} for contig {} should have been in the 'samtools depth' output".format(coord, contig)
                total_depth += sum(depth_dict[contig][coord][s] for s in depth_dict[contig][coord] if s == strand)
            readcounts = round(total_depth / read_len, 2)
            f.write("{}\t{}\n".format(gene, readcounts))

def count_uniq_bases_per_gene(contig_bases, gene_info):
    """ For each contig and strand, count the number of bases that belong to a single gene """
    counter = Counter()
    for contig in contig_bases:
        for strand in contig_bases[contig]:
            counter += Counter(val for key, val in contig_bases[contig][strand].items())
    #Address genes that were completely overlapped
    counter.update( {gene: 0 for gene in gene_info if gene not in counter} )
    return counter

def generate_gene_stats(uniq_gene_bases, gene_info, output_dir, gff3_base, stranded_type):
    """ Generate percentage of unique bases per gene """
    logging.debug("Generating overlap stats per gene")
    if stranded_type == "unstranded":
        ext = "unstranded.uniquebp.stats.tsv"
    else:
        ext = "stranded.uniquebp.stats.tsv"
    stats_file = "{}/{}.{}".format(output_dir, gff3_base, ext)
    with open(stats_file, 'w') as f:
        headers = "contig\tstrand\tgene\tlength\tuniq_bases\tpct_uniq\n"
        f.write(headers)
        for key, val in gene_info.items():
            (contig, start, stop, strand) = val
            length = stop - start + 1
            uniq_bases = uniq_gene_bases[key]
            pct_uniq = round(uniq_bases * 100 / length, 2)
            row = ( contig, strand, key, str(length), str(uniq_bases), str(pct_uniq) )
            f.write('\t'.join(row) + "\n")

def get_gene_from_attr(attr_field, ptrn):
    """ Only keep gene ID from attribute section of GFF3 file """
    match = ptrn.search(attr_field)
    assert match.lastindex == 1, "Attribute field found to have more than one match or no matches for ptrn {}".format(ptrn)
    return match.group(1)

def get_start_stop(first, second):
    start = min(first, second)
    stop = max(first, second)
    return (start, stop)

def index_bam(bam):
    name = mp.current_process().name
    logging.debug("{} - Indexing BAM ...".format(name))
    pysam.index(bam)

def merge_bam(bam_list, bam_out):
    """ Merge a list of BAM files """
    name = mp.current_process().name
    logging.debug("{} - Merging BAM files ...".format(name))
    assert len(bam_list), "BAM list has no files!"
    pysam.merge("-f", bam_out, *bam_list)

def parse_gff3(annot_file, is_gff3, stranded_type, feat_type, attr_type):
    """ Parse the GFF3 or GTF annotation file """
    if is_gff3:
        logging.info("Parsing GFF3 file")
        # Attribute field to parse IDs from
        ptrn = re.compile(r';?{0}=(\w+);?'.format(attr_type))
    else:
        logging.info("Parsing GTF file")
        ptrn = re.compile(r';? {0} \"(\w+)\";'.format(attr_type))
    assert(os.stat(annot_file).st_size > 0), "{} was empty".format(annot_file)

    # Dict -> contig_id -> plus/minus -> coord -> gene_id/overlap
    contig_bases = {}
    # Dict -> gene_id -> (contig_id/start/stop/sign)
    gene_info = {}

    stranded = 1
    if stranded_type == "no":
        stranded = 0

    with open(annot_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line == "##FASTA":
                break
            elif line.startswith('#'):
                continue
            entry = re.split(r'\t', line)
            #logging.debug(entry)
            if entry[2] == feat_type:
                contig_id = entry[0]
                gene_id = get_gene_from_attr(entry[8], ptrn)
                assert len(gene_id), "ID for attribute {} was not found for contig {} with feature {} ".format(attr_type, contig_id, feat_type)
                (start, stop) = get_start_stop( int(entry[3]), int(entry[4]) )
                sign = entry[6]
                if not stranded:
                    sign = "+"
                strand = set_strand(sign)
                gene_info[gene_id] = (contig_id, start, stop, sign)

                # Iniitialize strand-specific dicts as contig key is created
                contig_bases.setdefault(contig_id, { 'plus': {}, 'minus': {} })
                update_base_mapping( contig_bases[contig_id], gene_id, strand, start, stop)
    assert len(contig_bases.keys()), "Detected no {} entries in annotation file".format(feat_type)
    return (contig_bases, gene_info)

def process_bam(bam, contig_bases, gene_info, stranded_type, tmp_dir, output_dir):
    bam = bam.rstrip()
    name = mp.current_process().name
    logging.info("{} - Processing BAM file {}".format(name, bam))
    if not os.path.isfile(bam):
        logging.error(bam + " does not seem to exist")
    assert(os.stat(bam).st_size > 0), "{} was empty".format(bam)
    ln_bam = symlink_bam(bam, tmp_dir)
    # Index BAM
    index_bam(ln_bam)
    # Get average read length of BAM
    bam_fh = pysam.AlignmentFile(ln_bam, "rb")
    read_len = calc_avg_read_len(bam_fh)
    bam_fh.close()

    depth_dict = {}
    # Calculate depth of coverage for the BAM file or for split plus/minus BAM files
    if stranded_type == "no":
        calc_depth(depth_dict, ln_bam, "plus")
    else:
        strand_list = ["plus", "minus"]
        (pos_bam, neg_bam) = split_bam_by_strand(ln_bam, stranded_type)
        for idx, split_bam in enumerate([pos_bam, neg_bam]):
            logging.info("{} - Processing {} stranded BAM file".format(name, strand_list[idx]))
            calc_depth(depth_dict, split_bam, strand_list[idx])
    # Calculate readcount stats per gene
    out_bam = output_dir + "/" + os.path.basename(ln_bam)
    calc_readcounts_per_gene(contig_bases, gene_info, depth_dict, out_bam, read_len)
    logging.info("{} - Finished processing BAM file {} !".format(name, bam))

def set_strand(sign):
    if sign == "+":
        return "plus"
    return "minus"

def split_bam_by_strand(bam, strand_type):
    """ Split the input BAM file by plus/minus strands using the bitwise flag """
    name = mp.current_process().name
    logging.debug("{} - Splitting BAM into positive and negative strands".format(name))

    # Forward-stranded assay
    ## Positive Strand:
    ### R1 - forward (96), R2 - revcom (144)
    ## Negative Strand
    ### R1 - revcom (80), R2 - forward (160)
    # Reverse-stranded assay (i.e. Illumina)
    ## Positive Strand
    ### R1 - revcom (80), R2 - forward (160)
    ## Negative Strand
    ### R1 - forward (96), R2 - revcom (144)
    # (Note - Bitflags 1 and 2 are assumed for all of these)

    fwd_flags = {99, 147}
    rev_flags = {83, 163}

    pos_flags = fwd_flags; neg_flags = rev_flags
    if strand_type == "reverse":
        pos_flags = rev_flags; neg_flags = fwd_flags

    logging.debug("{} - Positive strand...".format(name))
    bam_list = []
    for idx, flag in enumerate(pos_flags):
        pos_bam = re.sub(r'\.bam', '.plus{}.bam'.format(idx), bam)
        # Adding 'catch_stdout=False' prevents redirection of output to stdout stream
        pysam.view("-f", str(flag), "-o", pos_bam, bam, catch_stdout=False)
        index_bam(pos_bam)
        bam_list.append(pos_bam)
    pos_bam = re.sub(r'\.bam', '.plus.bam', bam)
    merge_bam(bam_list, pos_bam)

    logging.debug("{} - Negative strand...".format(name))
    bam_list = []
    for idx, flag in enumerate(neg_flags):
        neg_bam = re.sub(r'\.bam', '.minus{}.bam'.format(idx), bam)
        pysam.view("-f", str(flag), "-o", neg_bam, bam, catch_stdout=False)
        index_bam(neg_bam)
        bam_list.append(neg_bam)
    neg_bam = re.sub(r'\.bam', '.minus.bam', bam)
    merge_bam(bam_list, neg_bam)

    return (pos_bam, neg_bam)

def store_depth(depth_dict, depth_out, strand):
    """ Store depth coverage information from 'samtools depth' command """
    name = mp.current_process().name
    logging.debug("{} - Storing 'samtools depth' output".format(name))
    for line in depth_out:
        line = line.rstrip()
        (contig, coord, depth) = re.split(r'\t', line)
        depth_dict.setdefault(contig, {})
        depth_dict[contig].setdefault(coord, {})
        depth_dict[contig][coord][strand] = int(depth)

def symlink_bam(bam, outdir):
    ln_bam = outdir + "/" + os.path.basename(bam)
    if os.path.islink(ln_bam):
        return ln_bam
    os.symlink(bam, ln_bam)
    name = mp.current_process().name
    logging.debug("{} - Symlinking {} to {}".format(name, bam, ln_bam))
    return ln_bam

def update_base_mapping(contig, gene, strand, start, stop):
    """ If coordinate already maps to another gene, it is an overlap """
    for coord in range(start, stop+1):
        str_coord = str(coord)
        if str_coord in contig[strand] and not gene == contig[strand][str_coord]:
            contig[strand][str_coord] = "overlap"
        else:
            contig[strand][str_coord] = gene

########
# Main #
########

def main():
    # Set up options parser and help statement
    description = "FADU - Feature Aggregate Depth Utility\nGenerate counts of reads that map to non-overlapping portions of genes"
    parser = ArgumentParser(description=description)
    bam_group = parser.add_mutually_exclusive_group(required=True)
    bam_group.add_argument("--bam_file", "-b", help="Path to BAM file. Choose between --bam_file or --bam_list.", metavar="/path/to/alignment.bam")
    bam_group.add_argument("--bam_list", "-B", help="List file of BAM-formatted files (no SAM at this time).  Choose between --bam_file or --bam_list.", metavar="/path/to/bam.list")
    gff_group = parser.add_mutually_exclusive_group(required=True)
    gff_group.add_argument("--gff3", "-g", help="Path to GFF3-formatted annotation file. Choose between --gff3 or --gtf.", metavar="/path/to/annotation.gff3")
    gff_group.add_argument("--gtf", "-G", help="Path to GTF-formatted annotation file. Choose between --gff3 or --gtf.", metavar="/path/to/annotation.gtf")
    parser.add_argument("--output_dir", "-o", help="Required. Directory to store output.", metavar="/path/to/output/dir", required=True)
    parser.add_argument("--tmp_dir", "-t", help="Directory to store temporary output.  If not provided, the output directory will serve as the temporary directory also", metavar="/path/to/tmp/dir", required=False)
    parser.add_argument("--stranded", "-s", help="Indicate if BAM reads are from a strand-specific assay", default="yes", choices=['yes', 'no', 'reverse'], required=False)
    parser.add_argument("--feature_type", "-f", help="Which GFF3/GTF feature type (column 3) to obtain readcount statistics for.  Default is 'gene'.  Case-sensitive.", default="gene", required=False)
    parser.add_argument("--attribute_type", "-a", help="Which GFF3/GTF attribute type (column 9) to obtain readcount statistics for.  Default is 'ID'.  Case-sensitive.", default="ID", required=False)
    parser.add_argument("--num_cores", "-n", help="Number of cores to spread processes to when processing BAM list.", default=10, type=check_positive, required=False)
    parser.add_argument("--debug", "-d", help="Set the debug level", default="INFO", metavar="DEBUG/INFO/WARNING/ERROR/CRITICAL")
    args = parser.parse_args()
    check_args(args, parser)

    # Decide if annotation file is gff3 or gtf
    is_gff3 = True
    annot_file = args.gff3
    if args.gtf:
        is_gff3 = False
        annot_file = args.gtf

    # Generate uniq bases per gene in GFF3
    (contig_bases, gene_info) = parse_gff3(annot_file, is_gff3, args.stranded, args.feature_type, args.attribute_type)
    uniq_gene_bases = count_uniq_bases_per_gene(contig_bases, gene_info)
    annot_basename = os.path.basename(annot_file)
    annot_base = os.path.splitext(annot_basename)[0]
    generate_gene_stats(uniq_gene_bases, gene_info, args.output_dir, annot_base, args.stranded)

    # If processing a single bam file...
    if args.bam_file:
        process_bam( args.bam_file, contig_bases, gene_info, args.stranded, args.tmp_dir, args.output_dir )
        return

    # ... otherwise process a list using parallel workers
    logging.info('Creating pool with %d processes\n' % args.num_cores)
    with open(args.bam_list, 'r') as list, mp.Pool(processes=args.num_cores) as p:
        # "apply_async" blocks downstream code from executing after res.get() is called.
        result = [p.apply_async(process_bam, (bam, contig_bases, gene_info, args.stranded, args.tmp_dir, args.output_dir)) for bam in list]
        [res.get() for res in result]

def check_args(args, parser):
    """ Validate the passed arguments """
    log_level = args.debug.upper()
    num_level = getattr(logging, log_level)

    # Verify that our specified log_level has a numerical value associated
    if not isinstance(num_level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    # Create the logger
    logging.basicConfig(level=num_level)

    # Args-specific checks
    if args.gff3 and not os.path.isfile(args.gff3):
        logging.error("GFF3 file does not seem to exist. Please check supplied path.")
    if args.gtf and not os.path.isfile(args.gtf):
        logging.error("GTf file does not seem to exist. Please check supplied path.")
    if args.bam_file and not os.path.isfile(args.bam_file):
        logging.error("BAM file does not seem to exist. Please check supplied path.")
    if args.bam_list and not os.path.isfile(args.bam_list):
        logging.error("BAM list file does not seem to exist. Please check supplied path.")
    if not os.path.isdir(args.output_dir):
        logging.debug("Creating output directory")
        os.mkdir(args.output_dir)
    if not args.tmp_dir:
        logging.debug("--tmp_dir not specified.  Will write temp files to output directory")
        args.tmp_dir = args.output_dir

def check_positive(value):
    """ Verify the integer argument passed in is positive """
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

if __name__ == '__main__':
    mp.freeze_support() # For Windows
    main()
    logging.info("All processes complete!  Exiting.")
    sys.exit(0)
