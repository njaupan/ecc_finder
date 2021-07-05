#!/usr/bin/env python3
"""
# June 2021
# If using this pipeline please cite : XXXXXXXXXX
#--------------------------------------------------------------------------+    
#                                                   
#	ecc_finder is a tool 
#       to robust and accurate detect extrachromosomal circular DNA (eccDNA)
#                                                        
#--------------------------------------------------------------------------+
#                                                      
#	AUTHOR: panpan ZHANG                            
#	CONTACT: njaupanpan@gmail.com                      
#                                                         
#	LICENSE:                                            
# 	GNU General Public License, Version 3               
#	http://www.gnu.org/licenses/gpl.html  
#                                             
#	VERSION: V.1                    
#                                                                                                       
#--------------------------------------------------------------------------+
"""

import argparse
import sys
import pysam
import re
from collections import Counter

parser = argparse.ArgumentParser(
        description="""Parse a sam/bamfile (from a stream) and output summary stats for each read.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i','--input_sam_or_bam', type=argparse.FileType('r'), default=sys.stdin, nargs='+')
parser.add_argument('-m', '--min_length', type=int, default=None)
parser.add_argument('-a', '--all_alignments', action='store_true',
                    help='Include secondary and supplementary alignments.')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
                    help='Output alignment stats to file instead of stdout.')
parser.add_argument('-s', '--summary', type=argparse.FileType('w'), default=sys.stderr,
                    help='Output summary to file instead of stderr.')

def stats_from_aligned_read(read, references, lengths):
    """Create summary information for an aligned read

    :param read: :class:`pysam.AlignedSegment` object
  
    Get the statistics from sam/bam file, such as Insertion, Deletion, Subsititution and Gap-compressed Identity by the same definition from minimap2 developer Heng Li's blog:
    http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity

     # I get value from https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigar

	M	BAM_CMATCH	0	match
	I	BAM_CINS	1	insertion relative to reference
	D	BAM_CDEL	2	deletion relative to reference
	N	BAM_CREF_SKIP	3	skipped region from the reference
	S	BAM_CSOFT_CLIP	4	soft clip, not aligned but still in sam file
	H	BAM_CHARD_CLIP	5	hard clip, not aligned and not in sam file
	P	BAM_CPAD	6	padding (silent deletion from padded reference)
	=	BAM_CEQUAL	7	I don't know exactly what it is 
	X	BAM_CDIFF	8	I don't know exactly what it is
	B	BAM_CBACK	9	I don't know exactly what it is

    """
    tags = dict(read.tags)
    try:
        tags.get('NM')
    except:
        raise IOError("Read is missing required 'NM' tag. Try running 'samtools fillmd -S - ref.fa'.")


    counts, _ = read.get_cigar_stats()
    match = counts[0]
    ins = counts[1]
    cigar = read.cigartuples 
    counter = Counter(elem[0] for elem in cigar)
    ins2 = counter [1]
    delt = counts[2]
    delt2 =  counter [2]
    # NM is edit distance: NM = INS + DEL + SUB
    tagsNM = tags['NM'] 
    sub = tags['NM'] - ins - delt
    block_length = match +tags['NM']
    blast_iden = 100*float(match)/(block_length)
    gap_compress = ins2 +delt2+sub
    gap_compress_iden = 100*float(match)/(block_length-gap_compress)
    #coverage = 100*float(read.query_alignment_length) / read_length
    direction = '-' if read.is_reverse else '+'

    results = {
	"ref": references[read.reference_id],
        "readID": read.qname,
        #"coverage": coverage,
	"read_len": read.infer_read_length(),
	"length": block_length,
	"blast_iden": blast_iden,
	"gap_compress_iden":gap_compress_iden,
	"qstart": read.query_alignment_start,
        "qend": read.query_alignment_end,
        "direction": direction,	
	"rstart": read.reference_start,
        "rend": read.reference_end,
	"match": match,
        "ins": ins,
	"ins2": ins2,
        "del": delt,
	"del2": delt2,
        "sub": sub,
	"gap_compress":gap_compress,
        #"ref_coverage": 100*float(read.reference_length) / lengths[read.reference_id],
    }

    return results


def main(arguments=None):
    args = parser.parse_args(arguments)

    count = 0
    n_unmapped = 0
    n_short = 0
    headers = ['ref','readID','read_len', 'length','blast_iden', 'gap_compress_iden','qstart', 'qend','direction','rstart', 'rend',  'match','ins', 'del', 'sub','gap_compress']
    args.output.write('\t'.join(headers))
    args.output.write('\n')

    if not isinstance(args.input_sam_or_bam, list):
        args.input_sam_or_bam = ('-',)

    for samfile in (pysam.AlignmentFile(x) for x in args.input_sam_or_bam):
        for read in samfile:
            if read.is_unmapped:
                n_unmapped += 1
                continue
            if not args.all_alignments and (read.is_secondary or read.is_supplementary):
                continue
            results = stats_from_aligned_read(read, samfile.references, samfile.lengths)
            if args.min_length is not None and results['length'] < args.min_length:
                n_short += 1
                continue

            count += 1
            out_row = (str(results[x]) for x in headers)
            args.output.write('\t'.join(out_row))
            args.output.write('\n')
        samfile.close()

    if count == 0:
        raise ValueError('No alignments processed. Check your sam/bam and filtering options.')

    args.summary.write('Mapped/Unmapped/Short: {}/{}/{}\n'.format(count, n_unmapped, n_short))


if __name__ == '__main__':
    main()
