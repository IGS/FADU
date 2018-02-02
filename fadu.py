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

def adjust_depth(depth_dict, read_pos):
    """ Adjust depth per coordinate position based on inserts and overlaps present in fragment """
    name = mp.current_process().name
    logging.debug("{} - Adjusting depth counts to account for fragments ...".format(name) )
    for query, vals in read_pos.items():
        contig = vals['contig']
        strand = vals['strand']
        (inserts, overlaps) = determine_pair_inserts_overlaps(vals)
        for coords in inserts:
            str_coords = str(coords)
            depth_dict[contig][str_coords][strand] += 1
        for coords in overlaps:
            str_coords = str(coords)
            depth_dict[contig][str_coords][strand] -= 1

def assign_read_to_strand(read, strand_type, pos_fh, neg_fh):
    """ Use the bitwise flags to assign the paired read to the correct strand """

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

    # Forward strand flags
    flag_96 = read.is_read1 and read.mate_is_reverse
    flag_144 = read.is_read2 and read.is_reverse
    # Reverse strand flags
    flag_80 = read.is_read1 and read.is_reverse
    flag_160 = read.is_read2 and read.mate_is_reverse

    pos_flags = {flag_96, flag_144}
    neg_flags = {flag_80, flag_160}
    if strand_type == "reverse":
        pos_flags = {flag_80, flag_160}
        neg_flags = {flag_96, flag_144}

    if any(flags for flags in pos_flags):
        pos_fh.write(read)
        return "plus"
    if any(flags for flags in neg_flags):
        neg_fh.write(read)
        return "minus"

def calc_avg_read_len(bam):
    """ Calculates average read len of all BAM reads """
    name = mp.current_process().name
    logging.debug("{} - Calculating average read length ...".format(name) )
    bam_fh = pysam.AlignmentFile(bam, "rb")
    # Since singletons are kept here, we also keep a paired read if the mate is unmapped
    total_query_len = sum( read.query_length for read in bam_fh.fetch() if not read.is_unmapped)
    num_reads = bam_fh.count()
    bam_fh.close()
    return round(total_query_len / num_reads)

def calc_depth(depth_dict, bam, strand):
    """ Calculate depth of coverage using 'samtools deptha' """
    name = mp.current_process().name
    logging.debug("{} - Calculating depth of coverage of BAM ...".format(name))
    depth_out = pysam.depth("-aa", "-m", "10000000", bam).splitlines(True)
    depth_out_file = re.sub(r'\.bam', '.read.depth', bam)
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
        for gene, val_list in sorted(gene_info.items()):
            # Assuming all coordinates for a given gene feature share the same contig and strand
            (contig, start, stop, strand) = val_list[0]
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

def determine_pair_inserts_overlaps(read_pair):
    """ Keep track of all paired read fragment inserts and overlaps per individual base """

    # Need to make coords 1-based since the depth output is 1-based
    r1start = read_pair['r1start'] + 1
    r2start = read_pair['r2start'] + 1
    r1end = read_pair['r1end'] + 1
    r2end = read_pair['r2end'] + 1

    min_coord = min(r1start, r2start, r1end, r2end)
    max_coord = max(r1start, r2start, r1end, r2end)
    r1_range = set(range(r1start, r1end))
    r2_range = set(range(r2start, r2end))
    combined_range = set(range(min_coord, max_coord))

    # Insert coords will not appear in either read
    inserts = combined_range.difference(*(r1_range, r2_range))
    # Overlapping coords will be common to both reads
    overlaps = r1_range.intersection(r2_range)

    return (inserts, overlaps)

def generate_gene_stats(uniq_gene_bases, gene_info, output_dir, gff3_base, stranded_type):
    """ Generate percentage of unique bases per gene """
    logging.debug("Generating overlap stats per gene")
    stranded = "unstranded" if stranded_type == "no" else "stranded"
    ext = stranded + ".uniquebp.stats.tsv"

    stats_file = "{}/{}.{}".format(output_dir, gff3_base, ext)
    with open(stats_file, 'w') as f:
        headers = "contig\tstrand\tgene\tlength\tuniq_bases\tpct_uniq\n"
        f.write(headers)
        for key, val_list in sorted(gene_info.items()):
            length = 0
            prev_stop = 0
            for elem in val_list:
                overlap = 0
                (contig, start, stop, strand) = elem
                # In some cases, the previous and current coords may have an overlap, possibly due to a frameshift
                if start <= prev_stop:
                    overlap = prev_stop - start + 1
                length += stop - start + 1 - overlap
                prev_stop = stop
            uniq_bases = uniq_gene_bases[key]
            pct_uniq = round(uniq_bases * 100 / length, 2)
            # Assumption is that each feature ID is on the same contig and same strand
            row = ( contig, strand, key, str(length), str(uniq_bases), str(pct_uniq) )
            f.write("\t".join(row) + "\n")

def get_gene_from_attr(attr_field, ptrn):
    """ Only keep gene ID from attribute section of GFF3 file """
    match = ptrn.search(attr_field)
    if not match:
        logging.error("Attribute field '{}' found to have no matches for ptrn {}".format(attr_field, ptrn))
    assert match.lastindex == 2, "Attribute field '{}' found to have more than one match for ptrn {}".format(attr_field, ptrn)
    return match.group(2)

def get_start_stop(first, second):
    start = min(first, second)
    stop = max(first, second)
    return (start, stop)

def index_bam(bam):
    name = mp.current_process().name
    logging.debug("{} - Indexing BAM file {}...".format(name, bam))
    pysam.index(bam)

def parse_bam_for_proper_pairs(bam, read_positions, pp_only, stranded_type, count_by_fragment, **kwargs):
    """ Iterate through the BAM file to get information about the properly paired read alignments """
    bam_fh = pysam.AlignmentFile(bam, "rb")
    ofh = None; pos_ofh = None; neg_ofh = None
    until_eof_flag = False
    strand = "plus"

    # BAM file that will be returned for downstream processing
    working_bam = bam

    if stranded_type != "no":
        pos_bam = re.sub(r'\.bam', '.plus.bam', bam)
        pos_ofh = pysam.AlignmentFile(pos_bam, "wb", template=bam_fh)
        neg_bam = re.sub(r'\.bam', '.minus.bam', bam)
        neg_ofh = pysam.AlignmentFile(neg_bam, "wb", template=bam_fh)

    if pp_only:
        working_bam = re.sub(r'\.bam', '.p_paired.bam', bam)
        ofh = pysam.AlignmentFile(working_bam, "wb", template=bam_fh)

    if kwargs:
        if "until_eof" in kwargs:
            until_eof_flag = True

    # NOTE: if there are lots of reference IDs, may be necessary to change to bam_fh.fetch(until_eof=True)
    for read in bam_fh.fetch(until_eof=until_eof_flag):
        if stranded_type != "no":
            if pp_only and not read.is_proper_pair:
                continue
            strand = assign_read_to_strand(read, stranded_type, pos_ofh, neg_ofh)

        if read.is_proper_pair
            # Store properly paired reads for adjusting depth by fragments later
            if count_by_fragment:
                store_properly_paired_read(read_positions, read, strand)
            # If only keep properly paired reads, write to file
            if pp_only:
                ofh.write(read)

    # Close any open BAM files
    bam_fh.close()
    if pp_only:
        ofh.close()
    if stranded_type != "no":
        pos_ofh.close()
        neg_ofh.close()
    return working_bam

def parse_gff3(annot_file, is_gff3, stranded_type, feat_type, attr_type):
    """ Parse the GFF3 or GTF annotation file """
    if is_gff3:
        logging.info("Parsing GFF3 file")
        # Attribute field to parse IDs from.  Keeping number of matching groups for both GTF and GFF ptrns uniform
        ptrn = re.compile(r'(;)?{0}=([\w|.]+);?'.format(attr_type))
    else:
        logging.info("Parsing GTF file")
        ptrn = re.compile(r'(; )?{0} \"([\w|.]+)\";'.format(attr_type))
    assert(os.stat(annot_file).st_size > 0), "{} was empty".format(annot_file)

    # Dict -> contig_id -> plus/minus -> coord -> gene_id/overlap
    contig_bases = {}
    # Dict -> gene_id -> [ (contig_id/start/stop/sign) ]
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
                    sign = "*"
                strand = set_strand(sign)

                #Store feature ID information to list.  A feature ID can be present in multiple GFF entries
                gene_info.setdefault(gene_id, [])
                gene_info[gene_id].append( (contig_id, start, stop, sign) )

                # Iniitialize strand-specific dicts as contig key is created
                contig_bases.setdefault(contig_id, { 'plus': {}, 'minus': {} })
                update_base_mapping( contig_bases[contig_id], gene_id, strand, start, stop)
    assert len(contig_bases.keys()), "Detected no {} entries in annotation file".format(feat_type)
    return (contig_bases, gene_info)

def process_bam(bam, contig_bases, gene_info, args):
    """ Process a single BAM alignment file for various metrics """
    bam = bam.rstrip()
    name = mp.current_process().name

    # Extract the argument variables that are needed
    stranded_type = args.stranded
    count_by = args.count_by
    pp_only = args.keep_only_properly_paired
    tmp_dir = args.tmp_dir
    output_dir = args.output_dir

    logging.info("{} - Processing BAM file {}".format(name, bam))
    if not os.path.isfile(bam):
        logging.error(bam + " does not seem to exist")
    assert(os.stat(bam).st_size > 0), "{} was empty".format(bam)
    ln_bam = symlink_bam(bam, tmp_dir)
    # Index BAM
    index_bam(ln_bam)

    # Are depth counts fragment-based or read-based?
    count_by_fragment = False
    if count_by.lower() == "fragment":
        count_by_fragment = True

    read_positions = {}
    depth_dict = {}
    read_len = 0

    working_bam = parse_bam_for_proper_pairs(ln_bam, read_positions, pp_only, stranded_type, count_by_fragment)

    # Strandedness determines if the working BAM file needs to be split by strand
    if stranded_type == "no":
        calc_depth(depth_dict, working_bam, "plus")
    else:
        strand_list = ["plus", "minus"]
        pos_bam = re.sub(r'\.bam', '.plus.bam', working_bam)
        neg_bam = re.sub(r'\.bam', '.minus.bam', working_bam)
        for idx, split_bam in enumerate([pos_bam, neg_bam]):
            calc_depth(depth_dict, split_bam, strand_list[idx])

    # Get the average read length of working bam file
    read_len = calc_avg_read_len(working_bam)

    # Adjust depth information for read pair fragments
    if count_by_fragment:
        logging.info("{} - Elected to count by fragments instead of reads".format(name))
        adjust_depth(depth_dict, read_positions)
        write_fragment_depth(depth_dict, working_bam, stranded_type)

    # Calculate readcount stats per gene
    out_bam = output_dir + "/" + os.path.basename(working_bam)
    calc_readcounts_per_gene(contig_bases, gene_info, depth_dict, out_bam, read_len)
    logging.info("{} - Finished processing BAM file {} !".format(name, bam))

def set_strand(sign):
    if sign == "+" or sign == "*":
        return "plus"
    return "minus"

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

def store_properly_paired_read(read_pos, read, strand):
    """ Store alignment information about the current read """
    query_name = read.query_name
    # The contig or chromosome name
    ref_name = read.reference_name

    # NOTE: PySAM coords are 0-based, BAM are 1-based
    # reference start position is always the leftmost coordinate.
    # reference end is one base to the right of the last aligned residue
    if query_name not in read_pos:
        read_pos.setdefault(query_name, {
            'r1start': None,
            'r1end': None,
            'r2start': None,
            'r2end': None,
            'contig': ref_name,
            'strand': strand
        })
        read_pos[query_name]['r1start'] = read.reference_start
        read_pos[query_name]['r1end'] = read.reference_end
    else:
        read_pos[query_name]['r2start'] = read.reference_start
        read_pos[query_name]['r2end'] = read.reference_end


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

def write_fragment_depth(depth_dict, bam, stranded_type):
    """ Write fragment depth to a file, or two files if input BAM was stranded """
    name = mp.current_process().name
    logging.debug("{} - Writing depth file(s) based on fragments ...".format(name))
    if stranded_type == "no":
        depth_out_file = re.sub(r'\.bam', '.fragment.depth', bam)
        with open(depth_out_file, 'w') as ofh:
            for contig, vals in sorted(depth_dict.items()):
                for coord in sorted(vals, key=int):
                    row = [contig, coord, str(vals[coord]['plus'])]
                    ofh.write("\t".join(row) + "\n")
    else:
        strand_list = ["plus", "minus"]
        for strand in strand_list:
            depth_out_file = re.sub(r'\.bam', '.{}.fragment.depth'.format(strand), bam)
            with open(depth_out_file, 'w') as ofh:
                for contig, vals in sorted(depth_dict.items()):
                    for coord in sorted(vals, key=int):
                        row = [contig, coord, str(vals[coord][strand])]
                        ofh.write("\t".join(row) + "\n")

########
# Main #
########

def main():
    # Set up options parser and help statement
    description = "Generate counts of reads that map to non-overlapping portions of genes"
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
    parser.add_argument("--count_by", "-c", help="How to count the reads when performing depth calculations.  Default is 'read'.", default="read", choices=['read', 'fragment'], required=False)
    parser.add_argument("--keep_only_properly_paired", help="Enable flag to remove any reads that are not properly paired from the depth count statistics.", action="store_true", required=False)
    parser.add_argument("--num_cores", "-n", help="Number of cores to spread processes to when processing BAM list.", default=1, type=check_positive, required=False)
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
        process_bam( args.bam_file, contig_bases, gene_info, args )
        return

    # ... otherwise process a list using parallel workers
    logging.info('Creating pool with %d processes\n' % args.num_cores)
    with open(args.bam_list, 'r') as list, mp.Pool(processes=args.num_cores) as p:
        # "apply_async" blocks downstream code from executing after res.get() is called.
        result = [p.apply_async(process_bam, (bam, contig_bases, gene_info, args )) for bam in list]
        [res.get() for res in result]

def check_args(args, parser):
    """ Validate the passed arguments """
    log_level = args.debug.upper()
    num_level = getattr(logging, log_level)

    # Verify that our specified log_level has a numerical value associated
    if not isinstance(num_level, int):
        raise ValueError('Invalid log level: {}'.format(log_level))

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
