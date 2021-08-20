#!/usr/bin/env python3

"""
# June 2021
# If using this pipeline please cite : XXXXXXXXXX
#--------------------------------------------------------------------------+    
#                                                   
#	ecc_finder is a tool 
#       to detect eccDNA using Illumina and ONT sequencing.  
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
#	VERSION: v1.0.0                  
#                                                                                                       
#--------------------------------------------------------------------------+
"""

import os
import sys
import glob
import argparse
import multiprocessing
import subprocess
import pysam
import numpy as np
import pandas as pd
from eccFinder_lib.utilities import log,run_oae,get_eccFinder_version
from eccFinder_lib.Aligner import Minimap2SAMAligner
from eccFinder_lib.Aligner import bwaAligner
from eccFinder_lib.Spliter import tidehunter


def get_median_read_coverage(file_prefix,output_path, num_threads, overwrite_files):
    log("INFO", "Calculating global read coverage")
    if os.path.isfile(output_path +file_prefix+".bam.stats"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".bam.stats")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".bam.stats")
            stats = pysam.stats("-@", str(num_threads), output_path +file_prefix+".bam")
            with open(output_path +file_prefix+".bam.stats", "w") as f:
                f.write(stats)
    else:
        stats = pysam.stats("-@", str(num_threads), output_path +file_prefix+".bam")
        with open(output_path +file_prefix+".bam.stats", "w") as f:
            f.write(stats)

    # Get the coverage histogram (for 1 to 1k)
    covs = []
    with open(output_path +file_prefix+".bam.stats") as f:
        for line in f:
            if line.startswith("COV"):
                covs.append(int(line.split("\t")[3]))

    # Get the median from the histogram
    covs = np.asarray(covs, dtype=np.int32)

    # Remove the last value, which is a catch-all for coverages > 1k
    covs = covs[:-1]
    mid = sum(covs) // 2
    cs = 0
    for i in range(len(covs)):
        cs += covs[i]
        if cs >= mid:
            return i
    raise ValueError("Unable to calculate read coverage. Check SAM/BAM files and stats file.")


def run_samtools(file_prefix,output_path, num_threads, overwrite_files):
    """ sort, filter and index alignments. """
    if os.path.isfile(output_path +file_prefix+".unit.f.bam"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".unit.f.bam")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".unit.f.bam")
            pysam.sort("-@", str(num_threads),"-n", "-o", output_path +file_prefix+".unit.n.bam", output_path +file_prefix+".unit.sam", catch_stdout=False)
            cmd1 = "samtools view"+ " -@ " + str(num_threads)+ " -h "+ output_path +file_prefix+".unit.n.bam"
            #cmd2 = "awk 'length($10) > 200 || $1 ~ /^@/' " 
            cmd3 = "sed 's|_sub0\\t|_sub0/1\\t|g' |sed 's|_sub1\\t|_sub0/2\\t|g' |sed 's|_sub2\\t|_sub2/1\\t|g'| sed 's|_sub3\\t|_sub2/2\\t|g' > " +output_path +file_prefix+".unit.f.bam"
            fullPipe = "{inS} |{sed} ".format(inS=cmd1, sed=cmd3)
            print(fullPipe)
            os.popen(fullPipe)
    else:
        pysam.sort("-@", str(num_threads),"-n",  "-o", output_path +file_prefix+".unit.n.bam", output_path +file_prefix+".unit.sam", catch_stdout=False)
        cmd1 = "samtools view"+ " -@ " + str(num_threads)+ " -h "+ output_path +file_prefix+".unit.n.bam"
        #cmd2 = "awk 'length($10) > 200 || $1 ~ /^@/' " 
        cmd3 = "sed 's|_sub0\\t|_sub0/1\\t|g' |sed 's|_sub1\\t|_sub0/2\\t|g' |sed 's|_sub2\\t|_sub2/1\\t|g'| sed 's|_sub3\\t|_sub2/2\\t|g' > " +output_path +file_prefix+".unit.f.bam"
        fullPipe = "{inS} |{sed} ".format(inS=cmd1, sed=cmd3)
        print(fullPipe)
        os.popen(fullPipe)

    log("INFO", "Indexing read alignments")
    if os.path.isfile(output_path +file_prefix+".bam.csi"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".bam.csi")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".bam.csi")
            pysam.sort("-@", str(num_threads), "-o",output_path +file_prefix+".bam", output_path +file_prefix+".sam", catch_stdout=False)
            pysam.index("-c", output_path +file_prefix+".bam", catch_stdout=False)
    else:
        pysam.sort("-@", str(num_threads), "-o", output_path +file_prefix+".bam", output_path +file_prefix+".sam", catch_stdout=False)
        pysam.index("-c", output_path +file_prefix+".bam", catch_stdout=False)


def run_sam2paf(file_prefix,output_path, overwrite_files):
    if os.path.isfile(output_path +file_prefix+".paf.bed"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".paf.bed")
        else: 
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".paf.bed")
            cmd1 = "paftools.js sam2paf "+ output_path +file_prefix+".sam"
            cmd2 = "awk -v OFS='\\t' '{print$6,$8,$9,$1,$2,$5}' "  
            cmd3 = "bedtools sort > " +output_path +file_prefix+".paf.bed"
            fullPipe = "{inS} |{awk}|{bed} ".format(inS=cmd1, awk=cmd2, bed=cmd3)
            os.popen(fullPipe)
    else:  
        cmd1 = "paftools.js sam2paf "+ output_path +file_prefix+".sam"
        cmd2 = "awk -v OFS='\\t' '{print$6,$8,$9,$1,$2,$5}' " 
        cmd3 = "bedtools sort > " +output_path +file_prefix+".paf.bed"
        fullPipe = "{inS} |{awk}|{bed} ".format(inS=cmd1, awk=cmd2, bed=cmd3)
        os.popen(fullPipe)

def run_filterBED(file_prefix,output_path, min_read,min_bound,overwrite_files):
    if os.path.isfile(output_path +file_prefix+".paf.bed.tmp2"):
        print("Hi")
    else:  
        bedtools_params= " -wao -f "+ str(min_bound)+" -F "+ str(min_bound)+ " -e "
        bedtools_cmd ="bedtools intersect -a "+ output_path +file_prefix+".site.bed -b " + output_path +file_prefix + ".paf.bed" + bedtools_params+ "> "+output_path +file_prefix+".paf.bed.tmp1"
        print(bedtools_cmd)
        subprocess.call(bedtools_cmd, shell=True) 

        cmd1 = "cat "+ output_path +file_prefix+".paf.bed.tmp1"
        cmd2 = "awk '$10>0'| sort -k7,7 | groupBy -g 1,2,3,7 -c 7,8,10 -o count,distinct,mean |awk '$5>1'" 
        cmd3 = "bedtools sort |groupBy -g 1,2,3 -c 4,5,7 -o count_distinct,sum,median > " +output_path +file_prefix+".paf.bed.tmp2"
        fullPipe = "{inS} |{group}|{bed} ".format(inS=cmd1, group=cmd2, bed=cmd3)
        os.popen(fullPipe)


def run_Genrich(file_prefix,output_path, min_ulen,merge_dist,overwrite_files):
    """ Detecting sites of genomic enrichment. """
    if os.path.isfile(output_path + file_prefix+".site.bed"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".site.bed")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".site.bed")
            GR_params = " -y "
            GR_params += " -l " + str(min_ulen)+" -g " + str(merge_dist)   
            GR_cmd = "Genrich -t "+ output_path +file_prefix+".unit.f.bam" + GR_params+ " -o "+ output_path +file_prefix+".site"
            subprocess.call(GR_cmd, shell=True) 
            cmd1 = "cut -f1-3 " + output_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
            os.popen("{inS} ".format(inS=cmd1))
    else:  
        GR_params = " -y  "
        GR_params += " -l " + str(min_ulen)+" -g " + str(merge_dist)   
        GR_cmd = "Genrich -t "+output_path +file_prefix+".unit.f.bam" + GR_params+ " -o "+ output_path +file_prefix+".site"
        subprocess.call(GR_cmd, shell=True) 
        cmd1 = "cut -f1-3 " + output_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
        os.popen("{inS} ".format(inS=cmd1))

    #for filename in glob.glob(os.path.join(output_path, "*unit*am")):
        #try:
            # Trying to remove a current file
            #os.remove(os.path.join(output_path, filename))
        #except EnvironmentError:
            # You don't have permission to do it
            #pass


def run_TideHunter(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,max_period_size,num_copies,overwrite_files):
    """ Spliting tandem repeats in one long read. """
    if os.path.isfile(output_path +file_prefix+".unit.fa"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".unit.fa")
        else:
            log("INFO", "Overwriting pre-existing file: " +output_path +file_prefix+".unit.fa")
            TH_params = " -c "+ str(num_copies)
            TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P " + str(max_period_size) + " "
            TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -u > "+ output_path +file_prefix+".unit.fa"
            subprocess.call(TH_cmd, shell=True) 
    else:  
        TH_params = " -c "+ str(num_copies)
        TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P " + str(max_period_size) + " "
        TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -u > "+ output_path +file_prefix+".unit.fa"
        print(TH_cmd)
        subprocess.call(TH_cmd, shell=True) 


def main():
    description = "A tool to detect eccDNA loci using ONT sequencing"
    parser = argparse.ArgumentParser(description=description, usage="ecc_finder.py map-ont <reference.idx> <query.fq>")
    parser.add_argument("idx", metavar="<reference.idx>", nargs='?', default="", type=str, help="index file of reference genome")
    parser.add_argument("query", metavar="<query.fq>", nargs='?', default="", type=str, help="query fastq file (uncompressed or bgzipped)")

    map_options = parser.add_argument_group("map options")
    mm2_default = "-ax map-ont" 
    map_options.add_argument('-t', metavar="INT",type=int, default=get_default_thread(),
                            help='number of CPU threads for mapping mode')
    map_options.add_argument("-m", metavar="PATH", type=str, default="minimap2", help="long read executable [minimap2]")
    map_options.add_argument("--mm2-params", metavar="STR", type=str, default=mm2_default, help="minimap2 parameters ['%s']" % mm2_default)
    map_options.add_argument("-l", metavar="INT", type=int, default=200, help="minimum alignment length [200]")
    map_options.add_argument("-g", metavar="INT", type=int, default=100, help="maximum distance between significant site [100]")

    val_options = parser.add_argument_group("validation options")
    val_options.add_argument("-c", metavar="INT", type=int, default=2, help="minimum copy number of tandem repeat in a long read [2]")
    val_options.add_argument("-e", metavar="FLT", type=float, default=0.25, help="maximum allowed divergence rate between two consecutive repeats [0.25]")
    val_options.add_argument("-p", metavar="INT", type=int, default=30, help="minimum period size of tandem repeat (>=2) [30]")
    val_options.add_argument("-P", metavar="INT", type=int, default=100000, help="maximum period size of tandem repeat (<=4294967295) [100K]")
    val_options.add_argument("-v", metavar="INT", type=int, default=10000, help="coverage validation window size [10000]")
    val_options.add_argument("--min-read", metavar="INT", type=int, default=3, help="filter eccDNA loci by unique mapped read number [3]")
    val_options.add_argument("--min-bound", metavar="FLT", type=float, default=0.8, help="filter eccDNA loci by boudary coverage [0.8]")

    out_options = parser.add_argument_group("output options")
    out_options.add_argument("-o", metavar="PATH", type=str, default="eccFinder_output", help="output directory [./eccFinder_output]")
    out_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    out_options.add_argument("-x", type=str, default="ecc.ont", help="add prefix to output [ecc.ont]")

    out_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if not args.idx or not args.query:
        parser.print_help()
        sys.exit("\n** The reference idx and query files are required **")

    log("VERSION", "ecc_finder " + get_eccFinder_version())
    log("CMD", "ecc_finder.py map-ont " + " ".join(sys.argv[1:]))

    idx_file = os.path.abspath(args.idx)
    query_file = os.path.abspath(args.query)

    if not os.path.isfile(idx_file):
        raise FileNotFoundError("Could not find file: %s" % idx_file)

    if not os.path.isfile(query_file):
        raise FileNotFoundError("Could not find file: %s" % query_file)

    num_threads = args.t
    num_copies = args.c
    max_divergence = args.e
    min_period_size = args.p
    max_period_size = args.P
    min_ulen = args.l
    merge_dist = args.g

    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"
    overwrite_files = args.w
    file_prefix = args.x

    # Debugging options
    debug_mode = args.debug
    debug_non_fltrd = output_path + file_prefix + ".debug.unfiltered.sam"
    debug_fltrd = output_path + file_prefix + ".debug.filtered.sam"
    debug_query_info = output_path + file_prefix + ".debug.query.info.txt"

    #Splitting into unit sequences of each tandem repeat for a long read
    TH_params = " -u "
    TH_params += " -t " + str(num_threads) +" -c " + str(num_copies)+" -p " + str(min_period_size)    
    log("INFO", "Splitting into unit sequences of each tandem repeat for a long read")
    run_TideHunter(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,max_period_size,num_copies,overwrite_files)

    # Align the query unit sequencesto the reference.
    mapont_aligner_path = args.m
    mapont_aligner = mapont_aligner_path.split("/")[-1]
    mm2_params = mm2_default
    mm2_params += " -t " + str(num_threads)
    log("INFO", "Mapping the query unit sequences to the reference genome")
    map = Minimap2SAMAligner(idx_file, [output_path +file_prefix+".unit.fa"],mapont_aligner_path, mm2_params,output_path +file_prefix+".unit", in_overwrite=overwrite_files)
    print(map)
    map.run_aligner()
    
    # Align the query raw read to the reference.
    log("INFO", "Mapping the query raw read to the reference genome")
    map_all = Minimap2SAMAligner(idx_file, [query_file],mapont_aligner_path, mm2_params,output_path + file_prefix , in_overwrite=overwrite_files)
    print(map_all)
    map_all.run_aligner()

    # Sort, filter and index the alignments.
    log("INFO", "Sorting, filtering and indexing unit sequence alignments")
    run_samtools(file_prefix,output_path, num_threads, overwrite_files)

    # Detect siginificant sites of genomic enrichment.
    log("INFO", "Detecting sites of genomic enrichment")
    run_Genrich(file_prefix,output_path, min_ulen,merge_dist,overwrite_files)

    # Convert alignment file sam to paf for filtering.
    log("INFO", "Converting alignment file sam to paf for filtering raw read alignments")
    run_sam2paf(file_prefix,output_path, overwrite_files)

    # Boundary coverage thresholds
    min_read= args.min_read
    min_bound = args.min_bound

    # Filtering eccDNA loci by boundary coverage.
    log("INFO", "Filtering eccDNA loci by boundary coverage")
    run_filterBED(file_prefix,output_path, min_read,min_bound,overwrite_files)


def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)



if __name__ == '__main__':
    main()
    
