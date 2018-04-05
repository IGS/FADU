#!/usr/bin/env python

"""
FADU.py - Feature Aggregate Depth Utility.

Description - Generate counts of reads that map to non-overlapping portions of genes
Requires: Python 3, Pysam version 0.12.0.1
By: Shaun Adkins (sadkins@som.umaryland.edu)
    Matthew Chung (mattchung@umaryland.edu)

"""
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

from collections import Counter
import argparse
import multiprocessing as mp
import os
import re
import sys
import logging
import pysam

###########
# Functions (in alphabetical order)
###########

def adjust_depth(depth_dict, pp_flatfile):
    """Adjust depth per coordinate position based on inserts and overlaps present in fragment."""
    name = mp.current_process().name
    logging.debug("%s - Adjusting depth counts to account for fragments ...", name)
    with open(pp_flatfile, 'r') as f:
        for line in f:
            entry = process_read_pair_line(line)
            contig = entry[0]
            strand = entry[1]
            coords_list = entry[2:6]
            (inserts, overlaps) = determine_pair_inserts_overlaps(coords_list)

            for coord in inserts:
                str_coord = str(coord)
                # Insert coord may originally be depth 0
                depth_dict[contig].setdefault(str_coord, {strand: 0})
                # One strand may have depth but not the other
                depth_dict[contig][str_coord].setdefault(strand, 0)
                depth_dict[contig][str_coord][strand] += 1
            for coord in overlaps:
                str_coord = str(coord)
                depth_dict[contig][str_coord][strand] -= 1

def assign_read_to_strand(read, strand_type, pos_fh, neg_fh):
    """Use the bitwise flags to assign the paired read to the correct strand."""
    # Forward-stranded assay
    ## Positive Strand:
    ### R1 - forward (97), R2 - revcom (145)
    ## Negative Strand
    ### R1 - revcom (81), R2 - forward (161)
    # Reverse-stranded assay (i.e. Illumina)
    ## Positive Strand
    ### R1 - revcom (81), R2 - forward (161)
    ## Negative Strand
    ### R1 - forward (97), R2 - revcom (145)
    # More rare, but flags 65, 129, 113, and 177 follow the same strandedness
    # Singletons will either have flag 0 or 16

    # Forward strand flags
    flag_97 = read.is_read1 and read.mate_is_reverse
    flag_145 = read.is_read2 and read.is_reverse
    flag_65 = read.is_read1 and not (read.is_reverse or read.mate_is_reverse)
    flag_129 = read.is_read2 and not (read.is_reverse or read.mate_is_reverse)
    # Reverse strand flags
    flag_81 = read.is_read1 and read.is_reverse
    flag_161 = read.is_read2 and read.mate_is_reverse
    flag_113 = read.is_read1 and read.is_reverse and read.mate_is_reverse
    flag_177 = read.is_read2 and read.is_reverse and read.mate_is_reverse

    pos_flags = {flag_97, flag_145, flag_65, flag_129}
    neg_flags = {flag_81, flag_161, flag_113, flag_177}
    if strand_type == "reverse":
        pos_flags = {flag_81, flag_161, flag_113, flag_177}
        neg_flags = {flag_97, flag_145, flag_65, flag_129}

    if any(flags for flags in pos_flags):
        pos_fh.write(read)
        return "plus"
    if any(flags for flags in neg_flags):
        neg_fh.write(read)
        return "minus"
    # Read must be a singleton
    if read.is_reverse:
        if strand_type == "reverse":
            neg_fh.write(read)
            return "minus"
        pos_fh.write(read)
        return "plus"
    else:
        if strand_type == "reverse":
            pos_fh.write(read)
            return "plus"
        neg_fh.write(read)
        return "minus"
    logging.warning("Read %s did not get assigned to either strand apparently.", read.query_name)
    return "skip"

def calc_avg_read_len(bam):
    """Calculate average read len of all mapped BAM reads."""
    name = mp.current_process().name
    logging.debug("%s - Calculating average read length ...", name)
    bam_fh = pysam.AlignmentFile(bam, "rb")
    # Unmapped reads will not factor downstream
    total_query_len = sum(read.query_length for read in bam_fh.fetch() if not read.is_unmapped)
    # Skip unmapped reads in counts
    num_reads = bam_fh.count(read_callback='all')
    bam_fh.close()
    avg_read_len = round(total_query_len / num_reads)
    logging.info("%s - Average read length - %d", name, avg_read_len)
    return avg_read_len

def calc_avg_frag_len(bam, pp_flatfile, pp_only):
    """Calculate average fragment length of properly paired reads, incl. read length of reads that are just mapped."""
    name = mp.current_process().name
    logging.debug("%s - Calculating average fragment length ...", name)

    # pp_flatfile stores reads as a pair
    # Multiply both by 2, to give singletons and singly-mapped reads less weight if added to the mix
    num_reads = 0
    total_query_len = 0
    with open(pp_flatfile, 'r') as f:
        for line in f:
            entry = process_read_pair_line(line)
            coords_list = [int(i) for i in entry[2:6]]
            (r1start, r1end, r2start, r2end) = coords_list
            min_coord = min(r1start, r2start, r1end, r2end)
            max_coord = max(r1start, r2start, r1end, r2end)
            frag_len = max_coord - min_coord
            total_query_len += frag_len * 2
            num_reads += 2
    avg_frag_len = round(total_query_len / num_reads)
    logging.info("%s - Average fragment length of just properly paired reads - %d",
                 name, avg_frag_len)

    if not pp_only:
        bam_fh = pysam.AlignmentFile(bam, "rb")
        # Overwrite existing num_reads as this value will include all properly paired reads
        num_reads = bam_fh.count(read_callback='all')
        # Properly paired reads were taken care of above for frag length...
        # get read len of other mapped reads here
        total_query_len += sum(read.query_length
                               for read in bam_fh.fetch() if not (read.is_unmapped or read.is_proper_pair))
        bam_fh.close()
    avg_frag_len = round(total_query_len / num_reads)
    logging.info("%s - Final average fragment length - %d", name, avg_frag_len)
    return avg_frag_len

def calc_depth(depth_dict, bam, strand):
    """Calculate depth of coverage using 'samtools depth'."""
    name = mp.current_process().name
    logging.debug("%s - Calculating depth of coverage of BAM ...", name)
    depth_out = pysam.depth("-aa", "-m", "10000000", bam).splitlines(True)
    depth_outfile = re.sub(r'\.bam', '.read.depth', bam)
    # Write to outfile and reopen that to store depth to dict (saves memory)
    with open(depth_outfile, 'w') as f:
        for line in depth_out:
            f.write(line)
    store_depth(depth_dict, depth_outfile, strand)

def calc_readcounts_per_gene(contig_bases, gene_info, depth_dict, out_bam, read_len):
    """Calculate number of reads that map to the uniq bases of a gene."""
    name = mp.current_process().name
    logging.debug("%s - Calculate readcounts per gene for BAM ...", name)
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
                # Coords without depth are not kept.
                if coord not in depth_dict[contig]:
                    continue
                total_depth += sum(depth_dict[contig][coord][s]
                                   for s in depth_dict[contig][coord] if s == strand)
            readcounts = round(total_depth / read_len, 2)
            f.write("{}\t{}\n".format(gene, readcounts))

def count_uniq_bases_per_gene(contig_bases, gene_info):
    """For each contig and strand, count the number of bases that belong to a single gene."""
    counter = Counter()
    for contig in contig_bases:
        for strand in contig_bases[contig]:
            counter += Counter(val for key, val in contig_bases[contig][strand].items())
    # Address genes that were completely overlapped
    counter.update({gene: 0 for gene in gene_info if gene not in counter})
    return counter

def determine_pair_inserts_overlaps(coords_list):
    """Keep track of all paired read fragment inserts and overlaps per read pair coordinates."""
    coords_list = [int(i) for i in coords_list]
    (r1start, r1end, r2start, r2end) = coords_list
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
    """Generate percentage of unique bases per gene."""
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
                # In some cases, the previous and current coords may have an overlap,
                # possibly due to a frameshift
                if start <= prev_stop:
                    overlap = prev_stop - start + 1
                length += stop - start + 1 - overlap
                prev_stop = stop
            uniq_bases = uniq_gene_bases[key]
            pct_uniq = round(uniq_bases * 100 / length, 2)
            # Assumption is that each feature ID is on the same contig and same strand
            row = (contig, strand, key, str(length), str(uniq_bases), str(pct_uniq))
            f.write("\t".join(row) + "\n")

def get_gene_from_attr(attr_field, ptrn):
    """Only keep gene ID from attribute section of GFF3 file."""
    match = ptrn.search(attr_field)
    if not match:
        raise RuntimeError("Attribute field '%s' found to have no matches"
                          " for ptrn %s", attr_field, ptrn)
    if not match.lastindex == 2:
        raise RuntimeError("Attribute field '{}' found to have"
                           " more than one match for ptrn {}".format(attr_field, ptrn))
    return match.group(2)

def get_start_stop(first, second):
    """Determine start and stop coordinates."""
    start = min(first, second)
    stop = max(first, second)
    return (start, stop)

def index_bam(bam):
    """Index BAM file."""
    name = mp.current_process().name
    logging.debug("%s - Indexing BAM file %s...", name, bam)
    pysam.index(bam)

def parse_bam_for_proper_pairs(
        bam, pp_flatfile, pp_only, stranded_type, count_by_fragment, **kwargs):
    """Iterate through the BAM file to get information about the properly paired read alignments."""
    name = mp.current_process().name
    bam_fh = pysam.AlignmentFile(bam, "rb")
    ofh = None
    pos_ofh = None
    neg_ofh = None
    ppff_ofh = None
    until_eof_flag = False
    strand = "plus"

    # BAM file that will be returned for downstream processing
    working_bam = bam
    if pp_only:
        logging.info("%s - Will remove all reads that are not properly paired ...", name)
        working_bam = re.sub(r'\.bam', '.p_paired.bam', bam)
        ofh = pysam.AlignmentFile(working_bam, "wb", template=bam_fh)
    if stranded_type != "no":
        logging.info("%s - Will split BAM by strand ...", name)
        pos_bam = re.sub(r'\.bam', '.plus.bam', working_bam)
        pos_ofh = pysam.AlignmentFile(pos_bam, "wb", template=bam_fh)
        neg_bam = re.sub(r'\.bam', '.minus.bam', working_bam)
        neg_ofh = pysam.AlignmentFile(neg_bam, "wb", template=bam_fh)
    if count_by_fragment:
        ppff_ofh = open(pp_flatfile, 'w')

    if kwargs:
        if "until_eof" in kwargs:
            until_eof_flag = True

    # Temporary dict to keep track of properly paired reads.
    # Read pair keys will be deleted when mate is found and info is written to file
    read_positions = {}

    # NOTE: if there are lots of reference IDs,
    # may be necessary to change to bam_fh.fetch(until_eof=True)
    for read in bam_fh.fetch(until_eof=until_eof_flag):
        if stranded_type != "no":
            if pp_only and not read.is_proper_pair:
                continue
            strand = assign_read_to_strand(read, stranded_type, pos_ofh, neg_ofh)

        if read.is_proper_pair:
            # Store properly paired reads for adjusting depth by fragments later
            if count_by_fragment:
                store_properly_paired_read(ppff_ofh, read_positions, read, strand)
            # If only keep properly paired reads, write to file
            if pp_only:
                ofh.write(read)

    # Close any open BAM files
    bam_fh.close()
    if count_by_fragment:
        ppff_ofh.close()
    if pp_only:
        ofh.close()
        index_bam(working_bam)
    if stranded_type != "no":
        pos_ofh.close()
        neg_ofh.close()
    return working_bam

def parse_gff3(annot_file, is_gff3, stranded_type, feat_type, attr_type):
    """Parse the GFF3 or GTF annotation file."""
    if is_gff3:
        logging.info("Parsing GFF3 file")
        # Attribute field to parse IDs from.
        # Keeping number of matching groups for both GTF and GFF ptrns uniform
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
            if entry[2] == feat_type:
                contig_id = entry[0]
                gene_id = get_gene_from_attr(entry[8], ptrn)
                if not gene_id:
                    raise AttributeError("ID for attribute {} was not found for contig {} with"
                                         " feature {}".format(attr_type, contig_id, feat_type))
                (start, stop) = get_start_stop(int(entry[3]), int(entry[4]))
                sign = entry[6]
                if not stranded:
                    sign = "*"
                strand = set_strand(sign)

                # Store feature ID information to list.
                # A feature ID can be present in multiple GFF entries
                gene_info.setdefault(gene_id, [])
                gene_info[gene_id].append((contig_id, start, stop, sign))

                # Iniitialize strand-specific dicts as contig key is created
                contig_bases.setdefault(contig_id, {'plus': {}, 'minus': {}})
                update_base_mapping(contig_bases[contig_id], gene_id, strand, start, stop)
    assert contig_bases.keys(), "Detected no {} entries in annotation file".format(feat_type)
    return (contig_bases, gene_info)

def process_bam(bam, contig_bases, gene_info, args):
    """Process a single BAM alignment file for various metrics."""
    bam = bam.rstrip()
    name = mp.current_process().name

    # Extract the argument variables that are needed
    stranded_type = args.stranded
    count_by = args.count_by
    pp_only = args.keep_only_properly_paired
    tmp_dir = args.tmp_dir
    output_dir = args.output_dir

    logging.info("%s - Processing BAM file %s", name, bam)
    if not os.path.isfile(bam):
        raise FileNotFoundError("%s does not seem to exist", bam)
    assert(os.stat(bam).st_size > 0), "{} was empty".format(bam)
    ln_bam = symlink_bam(bam, tmp_dir)
    # Index BAM
    index_bam(ln_bam)

    # Are depth counts fragment-based or read-based?
    count_by_fragment = False
    if count_by.lower() == "fragment":
        count_by_fragment = True

    depth_dict = {}
    read_len = 0

    pp_flatfile = re.sub(r'\.bam', '.properly_paired.tsv', ln_bam)
    working_bam = parse_bam_for_proper_pairs(
        ln_bam, pp_flatfile, pp_only, stranded_type, count_by_fragment)

    # Strandedness determines if the working BAM file needs to be split by strand
    if stranded_type == "no":
        calc_depth(depth_dict, working_bam, "plus")
    else:
        strand_list = ["plus", "minus"]
        pos_bam = re.sub(r'\.bam', '.plus.bam', working_bam)
        neg_bam = re.sub(r'\.bam', '.minus.bam', working_bam)
        for idx, split_bam in enumerate([pos_bam, neg_bam]):
            index_bam(split_bam)
            calc_depth(depth_dict, split_bam, strand_list[idx])

    # Adjust depth information for read pair fragments
    if count_by_fragment:
        logging.info("%s - Elected to count by fragments instead of reads", name)
        adjust_depth(depth_dict, pp_flatfile)
        if stranded_type == "no":
            write_fragment_depth(depth_dict, working_bam, "plus")
        else:
            strand_list = ["plus", "minus"]
            pos_bam = re.sub(r'\.bam', '.plus.bam', working_bam)
            neg_bam = re.sub(r'\.bam', '.minus.bam', working_bam)
            for idx, split_bam in enumerate([pos_bam, neg_bam]):
                write_fragment_depth(depth_dict, split_bam, strand_list[idx])

    if count_by_fragment:
        read_len = calc_avg_frag_len(working_bam, pp_flatfile, pp_only)
    else:
        read_len = calc_avg_read_len(working_bam)

    # Calculate readcount stats per gene
    out_bam = output_dir + "/" + os.path.basename(working_bam)
    calc_readcounts_per_gene(contig_bases, gene_info, depth_dict, out_bam, read_len)
    logging.info("%s - Finished processing BAM file %s !", name, bam)

def process_read_pair_line(line):
    """Process a properly paired read information line."""
    return re.split(r'\t', line.rstrip())

def set_strand(sign):
    """Assign a name to the strand."""
    if sign == "+" or sign == "*":
        return "plus"
    return "minus"

def store_depth(depth_dict, outfile, strand):
    """Store depth coverage information from 'samtools depth' command."""
    name = mp.current_process().name
    logging.debug("%s - Storing 'samtools depth' output", name)
    with open(outfile, 'r') as f:
        for line in f:
            line = line.rstrip()
            (contig, coord, depth) = re.split(r'\t', line)
            depth_dict.setdefault(contig, {})
            int_depth = int(depth)
            # Reducing memory
            if int_depth > 0:
                depth_dict[contig].setdefault(coord, {})
                depth_dict[contig][coord][strand] = int_depth

def store_properly_paired_read(pp_fh, read_pos, read, strand):
    """Store alignment information about the current read."""
    query_name = read.query_name
    # The contig or chromosome name
    contig = read.reference_name

    print ("{}\t{}".format(read.reference_start, read.reference_end))
    # NOTE: PySAM coords are 0-based, BAM are 1-based, so adjust dict to 1-based
    # reference start position is always the leftmost coordinate.
    # reference end is the same (https://www.biostars.org/p/84686/)
    if query_name not in read_pos:
        read_pos.setdefault(query_name, {
            'r1start': read.reference_start + 1,
            'r1end': read.reference_end,
        })
    else:
        r1start = str(read_pos[query_name]['r1start'])
        r1end = str(read_pos[query_name]['r1end'])
        r2start = str(read.reference_start + 1)
        r2end = str(read.reference_end)
        row = (contig, strand, r1start, r1end, r2start, r2end)
        pp_fh.write("\t".join(row) + "\n")
        # Do not need this entry anymore
        read_pos.pop(query_name)


def symlink_bam(bam, outdir):
    """Create symlink for passed in BAM file if it does not exist."""
    ln_bam = outdir + "/" + os.path.basename(bam)
    if os.path.islink(ln_bam):
        return ln_bam
    os.symlink(bam, ln_bam)
    name = mp.current_process().name
    logging.debug("%s - Symlinking %s to %s", name, bam, ln_bam)
    return ln_bam

def update_base_mapping(contig, gene, strand, start, stop):
    """If coordinate already maps to another gene, it is an overlap."""
    for coord in range(start, stop + 1):
        str_coord = str(coord)
        if str_coord in contig[strand] and not gene == contig[strand][str_coord]:
            contig[strand][str_coord] = "overlap"
        else:
            contig[strand][str_coord] = gene

def write_fragment_depth(depth_dict, bam, strand):
    """Write fragment depth to a file, or two files if input BAM was stranded."""
    name = mp.current_process().name
    logging.debug("%s - Writing depth file(s) based on fragments for BAM file %s ...", name, bam)

    bam_fh = pysam.AlignmentFile(bam, "rb")
    if bam_fh.count(read_callback='all') == 0:
        logging.info("%s - BAM file was empty of mapped reads."
                     " Cannot write fragment depth for this file.", name)
        return
    # Get lengths of references associated with BAM file
    assert len(bam_fh.references) == len(bam_fh.lengths), "Different number of reference names and lengths"
    ref_lens = dict(zip(bam_fh.references, bam_fh.lengths))
    bam_fh.close()

    depth_out_file = re.sub(r'\.bam', '.fragment.depth', bam)
    with open(depth_out_file, 'w') as ofh:
        for contig, vals in sorted(depth_dict.items()):
            assert contig in ref_lens, "Contig {} was not in the BAM file list of references.".format(contig)
            for coord in range(1, ref_lens[contig] + 1):
                str_coord = str(coord)
                # Write a 0 if depth was not in the dictionary
                depth = 0
                if str_coord in vals and strand in vals[str_coord]:
                    depth = vals[str_coord][strand]
                row = (contig, str_coord, str(depth))
                ofh.write("\t".join(row) + "\n")

########
# Main #
########

def main():
    """Entryway to program."""
    # Set up options parser and help statement
    description = "Generate counts of reads that map to non-overlapping portions of genes"
    parser = argparse.ArgumentParser(description=description)
    bam_group = parser.add_mutually_exclusive_group(required=True)
    bam_group.add_argument("--bam_file", "-b", metavar="/path/to/alignment.bam",
                           help="Path to BAM file. Choose between --bam_file or --bam_list.")
    bam_group.add_argument("--bam_list", "-B", metavar="/path/to/bam.list",
                           help="List file of BAM-formatted files (no SAM at this time)."
                           " Choose between --bam_file or --bam_list.")
    gff_group = parser.add_mutually_exclusive_group(required=True)
    gff_group.add_argument("--gff3", "-g", metavar="/path/to/annotation.gff3",
                           help="Path to GFF3-formatted annotation file."
                           " Choose between --gff3 or --gtf.")
    gff_group.add_argument("--gtf", "-G", metavar="/path/to/annotation.gtf",
                           help="Path to GTF-formatted annotation file."
                           " Choose between --gff3 or --gtf.")
    parser.add_argument("--output_dir", "-o", metavar="/path/to/output/dir", required=True,
                        help="Required. Directory to store output.")
    parser.add_argument("--tmp_dir", "-t", metavar="/path/to/tmp/dir", required=False,
                        help="Directory to store temporary output. If not provided"
                        " the output directory will serve as the temporary directory also")
    parser.add_argument("--stranded", "-s", default="yes",
                        choices=['yes', 'no', 'reverse'], required=False,
                        help="Indicate if BAM reads are from a strand-specific assay")
    parser.add_argument("--feature_type", "-f", default="gene", required=False,
                        help="Which GFF3/GTF feature type (column 3) to obtain"
                        " readcount statistics for. Default is 'gene'. Case-sensitive.")
    parser.add_argument("--attribute_type", "-a", default="ID", required=False,
                        help="Which GFF3/GTF attribute type (column 9) to obtain"
                        " readcount statistics for. Default is 'ID'. Case-sensitive.")
    parser.add_argument("--count_by", "-c", default="read", choices=['read', 'fragment'],
                        required=False, help="How to count the reads when"
                        " performing depth calculations.  Default is 'read'.")
    parser.add_argument("--keep_only_properly_paired", action="store_true", required=False,
                        help="Enable flag to remove any reads that"
                        " are not properly paired from the depth count statistics.")
    parser.add_argument("--num_cores", "-n", default=1, type=check_positive, required=False,
                        help="Number of cores to spread processes to when processing BAM list.")
    parser.add_argument("--debug", "-d", metavar="DEBUG/INFO/WARNING/ERROR/CRITICAL",
                        default="INFO", help="Set the debug level")
    args = parser.parse_args()
    check_args(args)

    # Decide if annotation file is gff3 or gtf
    is_gff3 = True
    annot_file = args.gff3
    if args.gtf:
        is_gff3 = False
        annot_file = args.gtf

    # Generate uniq bases per gene in GFF3
    (contig_bases, gene_info) = parse_gff3(
        annot_file, is_gff3, args.stranded, args.feature_type, args.attribute_type)
    uniq_gene_bases = count_uniq_bases_per_gene(contig_bases, gene_info)
    annot_basename = os.path.basename(annot_file)
    annot_base = os.path.splitext(annot_basename)[0]
    generate_gene_stats(uniq_gene_bases, gene_info, args.output_dir, annot_base, args.stranded)

    # If processing a single bam file...
    if args.bam_file:
        process_bam(args.bam_file, contig_bases, gene_info, args)
        return

    # ... otherwise process a list using parallel workers
    logging.info('Creating pool with %d processes\n', args.num_cores)
    with open(args.bam_list, 'r') as list, mp.Pool(processes=args.num_cores) as p:
        # "apply_async" blocks downstream code from executing after res.get() is called.
        result = [p.apply_async(process_bam, (bam, contig_bases, gene_info, args)) for bam in list]
        [res.get() for res in result]

def check_args(args):
    """Validate the passed arguments."""
    log_level = args.debug.upper()
    num_level = getattr(logging, log_level)

    # Verify that our specified log_level has a numerical value associated
    if not isinstance(num_level, int):
        raise ValueError('Invalid log level: {}'.format(log_level))

    # Create the logger
    logging.basicConfig(level=num_level)

    # Args-specific checks
    if args.gff3 and not os.path.isfile(args.gff3):
        raise FileNotFoundError("GFF3 file does not seem to exist. Please check supplied path.")
    if args.gtf and not os.path.isfile(args.gtf):
        raise FileNotFoundError("GTf file does not seem to exist. Please check supplied path.")
    if args.bam_file and not os.path.isfile(args.bam_file):
        raise FileNotFoundError("BAM file does not seem to exist. Please check supplied path.")
    if args.bam_list and not os.path.isfile(args.bam_list):
        raise FileNotFoundError("BAM list file does not seem to exist. Please check supplied path.")
    if not os.path.isdir(args.output_dir):
        logging.debug("Creating output directory")
        os.mkdir(args.output_dir)
    if args.tmp_dir:
        if not os.path.isdir(args.tmp_dir):
            logging.debug("Creating temporary directory")
            os.mkdir(args.tmp_dir)
    else:
        logging.debug("--tmp_dir not specified.  Will write temp files to output directory")
        args.tmp_dir = args.output_dir

def check_positive(value):
    """Verify the integer argument passed in is positive."""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue

if __name__ == '__main__':
    mp.freeze_support()  # For Windows
    main()
    logging.info("All processes complete!  Exiting.")
    sys.exit(0)
