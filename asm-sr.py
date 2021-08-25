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

from eccFinder_lib.utilities import log,run_oae,get_eccFinder_version
from eccFinder_lib.Aligner import Minimap2SAMAligner
from eccFinder_lib.Aligner import Minimap2Aligner
from eccFinder_lib.Aligner import bwaAligner
from eccFinder_lib.Spliter import tidehunter
from eccFinder_lib.Peaker import genrich


def run_fastp(file_prefix,output_path, query1_file,query2_file,num_threads, overwrite_files):
    """ remove adapter. """
    if os.path.isfile(output_path +file_prefix+".html"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".html")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".html")
            query1_file_out=output_path +file_prefix+".R1.fastq.gz"
            query2_file_out=output_path +file_prefix+".R2.fastq.gz"
            html_file_out=output_path +file_prefix+".html"
            fastp_cmd = "fastp --thread "+ str(num_threads) +" -i " + query1_file+ " -I "+ query2_file +" -o " + query1_file_out+ " -O "+ query2_file_out+" -h "+ html_file_out
            subprocess.call(fastp_cmd, shell=True) 
    else:
        query1_file_out=output_path +file_prefix+".R1.fastq.gz"
        query2_file_out=output_path +file_prefix+".R2.fastq.gz"
        html_file_out=output_path +file_prefix+ ".html"
        fastp_cmd = "fastp --thread "+ str(num_threads) +" -i " + query1_file+ " -I "+ query2_file +" -o " + query1_file_out+ " -O "+ query2_file_out+" -h "+ html_file_out
        subprocess.call(fastp_cmd, shell=True) 

def run_Unicycler(file_prefix,query_file,output_path, num_threads,num_len,overwrite_files):
    """ assembly short read. """
    if os.path.isfile(output_path +"contig.fasta"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +"contig.fasta")
        else:
            log("INFO", "Overwriting pre-existing file: " +output_path +"contig.fasta")
            Spades_params = " -t "+ str(num_threads) +" --min_fasta_length " + str(num_len) 
            Spades_cmd = "unicycler -1 "+ output_path +file_prefix+".R1.fastq.gz"+ " -2 "+ output_path +file_prefix+".R2.fastq.gz" +Spades_params 
            subprocess.call(Spades_cmd, shell=True) 
    else:  
        Spades_params = " -t "+ str(num_threads) +" --min_fasta_length " + str(num_len) 
        Spades_cmd = "unicycler -1 "+ output_path +file_prefix+".R1.fastq.gz"+ " -2 "+ output_path +file_prefix+".R2.fastq.gz" +Spades_params 
        subprocess.call(Spades_cmd, shell=True) 

def main():
    description = "A tool to detect eccDNA loci using ONT sequencing"
    parser = argparse.ArgumentParser(description=description, usage="ecc_finder.py asm-sr <query.fq1> <query.fq2> (option)")
    parser.add_argument("query", metavar="<query.fq>", nargs='?', default="", type=str, help="query fastq file (uncompressed or bgzipped)")

    asm_options = parser.add_argument_group("asm options")   
    asm_options.add_argument('-t', metavar="INT",type=int, default=100,
                            help='number of CPU threads for asmping mode')
    asm_options.add_argument('-l', metavar="INT",type=int, default=get_default_thread(),help='minimum fasta length of assembly[100]')

    out_options = parser.add_argument_group("output options")
    out_options.add_argument("-o", metavar="PATH", type=str, default="eccFinder_output", help="output directory [./eccFinder_output]")
    out_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    out_options.add_argument("-x", type=str, default="ecc.ont", help="add prefix to output [ecc.ont]")
    out_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if not args.query:
        parser.print_help()
        sys.exit("\n** The query files are required **")

    log("VERSION", "ecc_finder " + get_eccFinder_version())
    log("CMD", "ecc_finder.py asm-sr " + " ".join(sys.argv[1:]))

    query_file = os.path.abspath(args.query)
    query1_file = os.path.abspath(args.query1)
    query2_file = os.path.abspath(args.query2)
    
    if not os.path.isfile(query_file):
        raise FileNotFoundError("Could not find file: %s" % query_file)
    
    num_threads = args.t
    num_len = args.l
    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"
    overwrite_files = args.w
    file_prefix = args.x

    # Debugging options
    #debug_mode = args.debug

    #Remove adapter.
    log("INFO", "Remove adapter")
    run_fastp(file_prefix,output_path, query1_file,query2_file,num_threads, overwrite_files)

    #Consensus sequences of each tandem repeat for a long read
    log("INFO", "Run Unicycler assembmy at a wide range of k-mer sizes, evaluating the graph at each one")
    log("INFO", "Which will be the most time-comsuming")
    run_Unicycler(file_prefix,query_file,output_path, num_threads,num_len,overwrite_files)

    log("INFO", "Finished running ecc_finder")
    log("INFO", "Goodbye, have a nice day!")

def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)

if __name__ == '__main__':
    main()
