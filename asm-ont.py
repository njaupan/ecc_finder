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
from eccFinder_lib.Spliter import tidehunter
from eccFinder_lib.Peaker import genrich

def run_TideHunter(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,num_copies,overwrite_files):
    """ Spliting tandem repeats in one long read. """
    if os.path.isfile(output_path +file_prefix+".cons.fa"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".cons.fa")
        else:
            log("INFO", "Overwriting pre-existing file: " +output_path +file_prefix+".cons.fa")
            TH_params = " -c "+ str(num_copies)
            TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
            TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -l |awk '/^>/{print \">Consens\" ++i; next}{print}' > "+ output_path +file_prefix+".cons.fa"       
            subprocess.call(TH_cmd, shell=True) 
    else:  
        TH_params = " -c "+ str(num_copies)
        TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
        TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -l |awk '/^>/{print \">Consens\" ++i; next}{print}' > "+ output_path +file_prefix+".cons.fa"       
        subprocess.call(TH_cmd, shell=True)          

def run_TideHunter2(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,num_copies,five_prime,three_prime,overwrite_files):
    """ Spliting tandem repeats in one long read. """
    if os.path.isfile(output_path +file_prefix+".cons.fa"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".cons.fa")
        else:
            log("INFO", "Overwriting pre-existing file: " +output_path +file_prefix+".cons.fa")
            TH_params = " -c "+ str(num_copies) +" -5 " +five_prime+ " -3 " +three_prime
            TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
            TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -l |awk '/^>/{print \">Consens\" ++i; next}{print}' > "+ output_path +file_prefix+".cons.fa"       
            subprocess.call(TH_cmd, shell=True) 
    else:  
        TH_params = " -c "+ str(num_copies) +" -5 " +five_prime+ " -3 " +three_prime
        TH_params += " -t " + str(num_threads) +" -e " +str(max_divergence) +" -p " + str(min_period_size)+ " -P 1000000 "
        TH_cmd = "TideHunter"+ TH_params+ str(query_file)+ " -l |awk '/^>/{print \">Consens\" ++i; next}{print}' > "+ output_path +file_prefix+".cons.fa"       
        subprocess.call(TH_cmd, shell=True) 

def run_CDHit(file_prefix,output_path, num_threads, identity,memory_limit,seq_length,overwrite_files):
    """ sort, filter and index alignments. """
    if os.path.isfile(output_path +file_prefix+".fasta"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".fasta")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".fasta")
            CDHit_params = "  -c  "+ str(identity) + " -M " +str(memory_limit)+"  -l " +str(seq_length)
            CDHit_params += " -n 10 -d 0 -T " + str(num_threads) 
            CDHit_cmd = "cd-hit-est -i "+ output_path +file_prefix+".cons.fa"+ " -o "+ output_path +file_prefix+".cluster" +CDHit_params
            subprocess.call(CDHit_cmd, shell=True) 

            log("INFO", "Removing singletons")
            cmd1 = "cat "+ output_path +file_prefix+".cluster.clstr |tr '\\n' '\\t' "
            cmd2 = "sed 's/>Cluster/\\n>Cluster/g' |awk '$8 >=1' |tr '\\t' '\\n' >" + output_path +file_prefix+".clster"
            sub = "{inS} |{group}".format(inS=cmd1, group=cmd2)
            ps = subprocess.Popen(sub,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
        
            txt = "cat "+ output_path +file_prefix+".clster |awk '$3 ~/\*/' |cut -f2 -d '>' |awk -F '.' '{print $1}' ) |awk '/^>/{print \">Contig\" ++i; next}{print}' >"
            cmd="seqtk subseq " + output_path +file_prefix+".cluster <(" + txt +output_path +file_prefix+".fasta"
            subprocess.call(cmd, shell=True, executable="/bin/bash") 

    else:
        CDHit_params = "  -c  "+ str(identity) + " -M " +str(memory_limit)+"  -l " +str(seq_length)
        CDHit_params += " -n 10 -d 0 -T " + str(num_threads) 
        CDHit_cmd = "cd-hit-est -i "+ output_path +file_prefix+".cons.fa"+ " -o "+ output_path +file_prefix+".cluster" +CDHit_params
        subprocess.call(CDHit_cmd, shell=True) 

        log("INFO", "Removing singletons")
        cmd1 = "cat "+ output_path +file_prefix+".cluster.clstr |tr '\\n' '\\t' "
        cmd2 = "sed 's/>Cluster/\\n>Cluster/g' |awk '$8 >=1' |tr '\\t' '\\n' >" + output_path +file_prefix+".clster"
        sub = "{inS} |{group}".format(inS=cmd1, group=cmd2)
        ps = subprocess.Popen(sub,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        
        txt = "cat "+ output_path +file_prefix+".clster |awk '$3 ~/\*/' |cut -f2 -d '>' |awk -F '.' '{print $1}' ) |awk '/^>/{print \">Contig\" ++i; next}{print}' >"
        cmd="seqtk subseq " + output_path +file_prefix+".cluster <(" + txt +output_path +file_prefix+".fasta"
        subprocess.call(cmd, shell=True, executable="/bin/bash") 

def main():
    description = "A tool to detect eccDNA loci using ONT sequencing"
    parser = argparse.ArgumentParser(description=description, usage="python ecc_finder.py asm-ont <query.fq> (option)")
    parser.add_argument("query", metavar="<query.fq>", nargs='?', default="", type=str, help="query fastq file (uncompressed or bgzipped)")

    asm_options = parser.add_argument_group("asm options")   
    asm_options.add_argument('-t', metavar="INT",type=int, default=get_default_thread(),
                            help='number of CPU threads for asmping mode')
    asm_options.add_argument("--five-prime",metavar="STR", type=str, help="5' adapter sequence (sense strand) [NULL]")
    asm_options.add_argument("--three-prime",metavar="STR", type=str, help="3' adapter sequence (anti-sense strand) [NULL]")

    val_options = parser.add_argument_group("consensus options")
    val_options.add_argument("-n", metavar="INT", type=int, default=2, help="minimum copy number of tandem repeat in a long read [2]")
    val_options.add_argument("-e", metavar="FLT", type=float, default=0.25, help="maximum allowed divergence rate between two consecutive repeats [0.25]")
    val_options.add_argument("-s", metavar="INT", type=int, default=30, help="minimum period size of tandem repeat (>=2) [30]")
    
    val_options.add_argument("-c", metavar="INT", type=int, default=0.8, help="minimum sequence identity for clustering [0.8]")
    val_options.add_argument("-l", metavar="INT", type=int, default=200, help="minimum length of throw_away_sequences [200]")
    val_options.add_argument("-m", metavar="INT", type=int, default=800, help="memory limit (in MB) for CD-hit clustering program [800]")

    out_options = parser.add_argument_group("output options")
    out_options.add_argument("-o", metavar="PATH", type=str, default="eccFinder_asm_output", help="output directory [./eccFinder_asm_output]")
    out_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    out_options.add_argument("-x", type=str, default="ecc.asm.ont", help="add prefix to output [ecc.asm.ont]")
    out_options.add_argument("--debug", action='store_true', default=False, help=argparse.SUPPRESS)

    args = parser.parse_args()

    if not args.query:
        parser.print_help()
        sys.exit("\n** The query files are required **")

    log("VERSION", "ecc_finder " + get_eccFinder_version())
    log("CMD", "python ecc_finder.py asm-ont " + " ".join(sys.argv[1:]))

    query_file = os.path.abspath(args.query)

    if not os.path.isfile(query_file):
        raise FileNotFoundError("Could not find file: %s" % query_file)
    
    num_threads = args.t
    max_divergence = args.e
    min_period_size = args.s
    num_copies = args.n
    five_prime = args.five_prime
    three_prime = args.three_prime
    seq_length = args.l
    identity = args.c
    memory_limit = args.m

    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"
    overwrite_files = args.w
    file_prefix = args.x

    # Debugging options
    #debug_mode = args.debug

    #consensus sequences of each tandem repeat for a long read
    log("INFO", "Detecting tandem repeat pattern from long reads")
    if five_prime and three_prime:  
        print("Adapter sequences are provided")
        five_prime=os.path.abspath(five_prime)
        three_prime=os.path.abspath(three_prime)
        run_TideHunter2(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,num_copies,five_prime,three_prime,overwrite_files)
    else:
        run_TideHunter(file_prefix,query_file,output_path, num_threads, max_divergence,min_period_size,num_copies,overwrite_files)
        
    # cluster consensus sequences by identity.
    log("INFO", "Clustering by identity")
    run_CDHit(file_prefix,output_path, num_threads,identity,memory_limit,seq_length,overwrite_files)

    log("INFO", "Finished running ecc_finder")
    log("INFO", "Goodbye, have a nice day!")

def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)

if __name__ == '__main__':
    main()
