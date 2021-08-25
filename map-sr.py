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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pybedtools
from pybedtools import BedTool

from eccFinder_lib.utilities import log,run_oae,get_eccFinder_version
from eccFinder_lib.Aligner import Minimap2SAMAligner
from eccFinder_lib.Aligner import Minimap2Aligner
from eccFinder_lib.Aligner import bwaAligner
from eccFinder_lib.Spliter import tidehunter
from eccFinder_lib.Peaker import genrich


def run_fastp(file_prefix,align_path, query1_file,query2_file,num_threads, overwrite_files):
    """ remove adapter. """
    if os.path.isfile(align_path +file_prefix+".html"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + align_path +file_prefix+".html")
        else:
            log("INFO", "Overwriting pre-existing file: " + align_path +file_prefix+".html")
            query1_file_out=align_path +file_prefix+".R1.fastq.gz"
            query2_file_out=align_path +file_prefix+".R2.fastq.gz"
            html_file_out=align_path +file_prefix+".html"
            fastp_cmd = "fastp --thread "+ str(num_threads) +" -i " + query1_file+ " -I "+ query2_file +" -o " + query1_file_out+ " -O "+ query2_file_out+" -h "+ html_file_out
            subprocess.call(fastp_cmd, shell=True) 
    else:
        query1_file_out=align_path +file_prefix+".R1.fastq.gz"
        query2_file_out=align_path +file_prefix+".R2.fastq.gz"
        html_file_out=align_path +file_prefix+ ".html"
        fastp_cmd = "fastp --thread "+ str(num_threads) +" -i " + query1_file+ " -I "+ query2_file +" -o " + query1_file_out+ " -O "+ query2_file_out+" -h "+ html_file_out
        subprocess.call(fastp_cmd, shell=True) 


def run_samtools(file_prefix,align_path, num_threads, overwrite_files):
    """ sort, filter and index alignments. """
    if os.path.isfile(align_path +file_prefix+".bam"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + align_path +file_prefix+".bam")
        else:
            log("INFO", "Overwriting pre-existing file: " + align_path +file_prefix+".bam")
            pysam.sort("-@", str(num_threads),"-n",  "-o", align_path +file_prefix+".bam", align_path +file_prefix+".sam", catch_stdout=False)
    else:     
        pysam.sort("-@", str(num_threads),"-n",  "-o", align_path +file_prefix+".bam", align_path +file_prefix+".sam", catch_stdout=False)

def run_Genrich(file_prefix,align_path,output_path,num_threads, min_peak,max_dist,max_pvalue,overwrite_files):
    """ Detecting sites of genomic enrichment. """
    if os.path.isfile(output_path +file_prefix+".site.bed"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".site.bed")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".site.bed")
            GR_params = " -v "
            GR_params += " -l " + str(min_peak)+" -g " + str(max_dist) +" -p " + str(max_pvalue) 
            GR_cmd = "Genrich -t "+ align_path +file_prefix+".bam" + GR_params+ " -o "+ output_path +file_prefix+".site"
            subprocess.call(GR_cmd, shell=True) 
            cmd1 = "cut -f1-3 " + output_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
            os.popen("{inS} ".format(inS=cmd1))
    else:          
        GR_params = " -v "
        GR_params += " -l " + str(min_peak)+" -g " + str(max_dist) +" -p " + str(max_pvalue)    
        GR_cmd = "Genrich -t "+align_path +file_prefix+".bam" + GR_params+ " -o "+ output_path +file_prefix+".site"
        subprocess.call(GR_cmd, shell=True) 
        cmd1 = "cut -f1-3 " + output_path +file_prefix+".site" + " > " +output_path +file_prefix+".site.bed"
        os.popen("{inS} ".format(inS=cmd1))

def run_bedtoolss(file_prefix,align_path,output_path, num_threads, overwrite_files):
    """ sort, filter and index alignments. """
    if os.path.isfile(output_path +file_prefix+".bam.bed"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".bam.bed")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".bam.bed")
            #pybedtools.BedTool(output_path +file_prefix+".bam").bam_to_bed().saveas(output_path +file_prefix+".bam.bed")
            cmd1 = "bedtools bamtobed -i "+ align_path +file_prefix+".bam" 
            cmd2 = "sed 's|/|\\t|g' |cut -f1-5,7 >" +output_path +file_prefix+".bam.bed"
            beds = "{inS} |{sed} ".format(inS=cmd1, sed=cmd2)
            ps = subprocess.Popen(beds,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
    else:
        #Bamtobed and seperate forward and reverse pair from read ID
        cmd1 = "bedtools bamtobed -i "+ align_path +file_prefix+".bam" 
        cmd2 = "sed 's|/|\\t|g' |cut -f1-5,7 >" +output_path +file_prefix+".bam.bed"
        beds = "{inS} |{sed} ".format(inS=cmd1, sed=cmd2)
        ps = subprocess.Popen(beds,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        output = ps.communicate()[0]

def run_split(file_prefix,output_path, num_threads, overwrite_files):
    if os.path.isfile(output_path +file_prefix+".split.bed"):
        log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".split.bed")
    else:
        df=pybedtools.BedTool(output_path +file_prefix+".bam.bed").to_dataframe()
        #Based on chrom and read ID, sort alignment by base and merge alignment to avoid overlapped reads
        df['chrom'] = df[['chrom','name']].apply(lambda x : '{}__{}'.format(x[0],x[1]), axis=1) 
        x =BedTool.sort(BedTool.from_dataframe(df))
        d =BedTool.merge(x,c = '5,6', o = 'collapse,collapse')
        dM=d.groupby(g = [1], c = [2,3,4,5],o = ['min','max','collapse','collapse'] )
        #Based on chrom and read ID, get strand and pair info and format dataframe
        namess=['chrom','start','end','pair','strand']
        dM=BedTool.to_dataframe(dM,names = namess)
        #Extract split reads, which has 3unique hits, same chromosome and same oriantation for a pair
        splitM =dM[(dM.pair == "1,2,1") & (dM.strand == "+,-,+" ) | 
                (dM.pair =="2,1,2") & (dM.strand == "-,+,-") | 
                (dM.pair =="1,2,1") & (dM.strand == "-,+,-") | 
                (dM.pair =="2,1,2") & (dM.strand == "+,-,+") 
                ]
        splitM[['chrom','read']] = splitM['chrom'].str.split('__',expand=True)
        splitM[['len']] = splitM.end-splitM.start
        splitM =splitM.sort_values(by=['chrom', 'start'])
        splitM.to_csv(output_path +file_prefix+".split.bed", header=False, index = False, sep='\t') 
        
def run_disc(file_prefix,output_path, num_threads, overwrite_files):
    if os.path.isfile(output_path +file_prefix+".disc.bed"):
        log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".disc.bed")
    else:
        raw=output_path +file_prefix+".bam.bed"
        df=pybedtools.BedTool(raw).to_dataframe()
        #Based on chrom and read ID, sort alignment by base and merge alignment to avoid overlapped reads
        df['chrom'] = df[['chrom','name']].apply(lambda x : '{}__{}'.format(x[0],x[1]), axis=1) 
        x =BedTool.sort(BedTool.from_dataframe(df))
        y =BedTool.merge(x,c = '5,6', o = 'collapse,collapse')
        d =BedTool.to_dataframe(y,names = ['chrom','start','end','pair','strand'])
        dM=d[(d.pair != "1,2" ) & (d.pair != "2,1")]
        dM=BedTool.from_dataframe(dM).groupby(g = [1], c = [2,3,4,5],o = ['min','max','collapse','collapse'] )
        #Based on chrom and read ID, get strand and pair info and format dataframe
        namess=['chrom','start','end','pair','strand']
        dM=BedTool.to_dataframe(dM,names = namess)
        pd.options.mode.chained_assignment = None 
        #Extract discordant reads, which has outward facing direction for a pair
        disc=dM[(dM.pair =="1,2") & (dM.strand == "-,+" ) ]
        disc[['chrom','read']] = disc['chrom'].str.split('__',expand=True)
        disc[['len']] = disc.end-disc.start
        disc =disc.sort_values(by=['chrom', 'start'])
        disc.to_csv(output_path +file_prefix+".disc.bed", header=False, index = False, sep='\t') 

def read_files(filelist):
    dfs = []
    for fn in filelist:
        df=pd.read_csv(fn,sep='\t',header=None,names=['Chr', 'start','end','read'])
        dfs.append(df)
    from functools import reduce
    return reduce(lambda  left,right: pd.merge(left,right,on=['Chr', 'start','end'],how='inner'), dfs).fillna(0)

def run_intersect(file_prefix,output_path, num_threads, min_read,min_cov,overwrite_files):
    """ sort, filter and index alignments. """
    if os.path.isfile(output_path +file_prefix+".csv"):
        log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".csv")
    else:
        a=pybedtools.example_bedtool(output_path +file_prefix+".site.bed")
        b=pybedtools.example_bedtool(output_path +file_prefix+".split.bed")
        see=a.intersect(b, wb=True,wa=True)
        s=pybedtools.BedTool(see).to_dataframe()
        s['cov'] = (s.strand-s.score)/(s.end-s.start)*100
        s=s[(s['cov']<=100)]
        s=BedTool.from_dataframe(s).groupby(g = [1,2,3], c = [9],o = ['count_distinct'] ).saveas(output_path +file_prefix+".split.num.bed")

        c=pybedtools.example_bedtool(output_path +file_prefix+".disc.bed")
        see=a.intersect(c, wb=True,wa=True)
        d=pybedtools.BedTool(see).to_dataframe()
        d['cov'] = (d.strand-d.score)/(d.end-d.start)*100
        d=d[ (d['cov']<=100)]
        d=BedTool.from_dataframe(d).groupby(g = [1,2,3], c = [9],o = ['count_distinct'] ).saveas(output_path +file_prefix+".disc.num.bed")
        files = glob.glob(output_path +file_prefix+".split.num.bed")
        files.extend(glob.glob(output_path +file_prefix+".disc.num.bed"))
        df = read_files(files) 
        float_col = df.select_dtypes(include=['float64'])
        for col in float_col.columns.values:
            df[col] = df[col].astype('int64')
        print(df)
        df.to_csv(output_path +file_prefix+".csv", header=False, index = False, sep='\t') 
        
def run_getFasta(output_path, file_prefix ,ref_genome,overwrite_files):
    if os.path.isfile(output_path +file_prefix +".fasta"):
        if not overwrite_files:
            log("INFO", "Retaining pre-existing file: " + output_path +file_prefix+".fasta")
        else:
            log("INFO", "Overwriting pre-existing file: " + output_path +file_prefix+".fasta")
            cmd ="seqtk subseq "+ ref_genome+ " " + output_path +file_prefix +".csv" + "> "+output_path +file_prefix +".fasta"
            subprocess.call(cmd, shell=True) 
    else: 
        cmd ="seqtk subseq "+ ref_genome+ " " + output_path +file_prefix +".csv" + "> "+output_path +file_prefix +".fasta"
        subprocess.call(cmd, shell=True) 

def main():
    description = "A tool to detect eccDNA loci using Illumina read sequencing"
    parser = argparse.ArgumentParser(description=description, usage="ecc_finder.py map-sr <reference.idx> <query.fq1> <query.fq2>")
    parser.add_argument("idx", metavar="<reference.idx>", nargs='?', default="", type=str, help="index file of reference genome")
    parser.add_argument("query1", metavar="<query.fq1>", nargs='?', default="", type=str, help="query fastq forward file (uncompressed or bgzipped)")
    parser.add_argument("query2", metavar="<query.fq2>", nargs='?', default="", type=str, help="query fastq reverse file (uncompressed or bgzipped)")

    map_options = parser.add_argument_group("map options")   
    map_options.add_argument('-t', metavar="INT",type=int, default=get_default_thread(),
                            help='number of CPU threads for mapping mode')
    map_options.add_argument("--aligner", metavar="PATH", type=str, default="bwa", help="short read aligner executable ('bwa', 'minimap2') [bwa]")
    map_options.add_argument("--bwa-params", metavar="STR", type=str, default="mem", help="space delimted bwa parameters ['mem']")
    map_options.add_argument("--minimap2-params", metavar="STR", type=str, default="-ax sr", help="space delimted minimap2 parameters ['-ax sr']")
    map_options.add_argument("-g", metavar="STR", type=str, default="", help="reference genome size larger than 4Gb [yes]")

    peak_options = parser.add_argument_group("peak-calling options")
    peak_options.add_argument("-l", metavar="INT", type=int, default=200, help="minimum length of a peak [200]")
    peak_options.add_argument("-d", metavar="INT", type=int, default=1000, help="maximum distance between signif. sites [1000]")
    peak_options.add_argument("-p", metavar="FLT", type=float, default=0.05, help="maximum p-value [0.05]")

    val_options = parser.add_argument_group("validation options")
    val_options.add_argument("-r", metavar="<reference.fa>", type=str, default="", help="reference genome fasta (uncompressed or bgzipped)")
    val_options.add_argument("--min-read", metavar="INT", type=int, default=3, help="filter locus by unique mapped read number [3]")
    val_options.add_argument("--min-cov", metavar="FLT", type=float, default=1, help="filter locus at regions by raw read coverage (# aligned bases / total bases)")

    out_options = parser.add_argument_group("output options")
    out_options.add_argument("-o", metavar="PATH", type=str, default="eccFinder_output", help="output directory [./eccFinder_output]")
    out_options.add_argument("-w", action='store_true', default=False, help="overwrite intermediate files")
    out_options.add_argument("-x", type=str, default="ecc.ill", help="add prefix to output [ecc.ill]")

    args = parser.parse_args()

    if not args.idx or not args.query1 or not args.query2:
        parser.print_help()
        sys.exit("\n** The reference idx and query files are required **")

    log("VERSION", "ecc_finder " + get_eccFinder_version())
    log("CMD", "ecc_finder.py map-ill " + " ".join(sys.argv[1:]))

    idx_file = os.path.abspath(args.idx)
    query1_file = os.path.abspath(args.query1)
    query2_file = os.path.abspath(args.query2)
    

    if not os.path.isfile(idx_file):
        raise FileNotFoundError("Could not find file: %s" % idx_file)

    if not os.path.isfile(query1_file):
        raise FileNotFoundError("Could not find file: %s" % query1_file)
    if not os.path.isfile(query2_file):
        raise FileNotFoundError("Could not find file: %s" % query2_file)
    
    num_threads = args.t

    min_peak = args.l
    max_dist = args.d
    max_pvalue = args.p

    ref_genome = args.r

    min_read= args.min_read
    min_cov = args.min_cov

    if min_read < 0:
        if min_read != -1:
            raise ValueError("--min-read must be >=3")
    if ref_genome:
        ref_genome=os.path.abspath(ref_genome)

    output_path = args.o
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
    output_path = os.path.abspath(output_path) + "/"

    align_path = os.path.abspath(args.o) + "/align_files/"
    if not os.path.isdir(output_path+"align_files/"):
        os.makedirs(output_path+"align_files/")
    align_path = os.path.abspath(args.o) + "/align_files/"

    overwrite_files = args.w
    file_prefix = args.x

    # Debugging options
    # debug_mode = args.debug
    # debug_non_fltrd = output_path + file_prefix + ".debug.unfiltered.sam"
    # debug_fltrd = output_path + file_prefix + ".debug.filtered.sam"
    # debug_query_info = output_path + file_prefix + ".debug.query.info.txt"

    #Remove adapter.
    log("INFO", "Remove adapter")
    run_fastp(file_prefix,align_path, query1_file,query2_file,num_threads, overwrite_files)

    #Align the query raw read to the reference.
    log("INFO", "Align the query raw read to the reference")
    mapsr_aligner_path = args.aligner
    genome_aligner = mapsr_aligner_path.split("/")[-1]
    if genome_aligner.split("/")[-1] not in {'bwa', 'minimap2'}:
        raise ValueError("Must specify either 'bwa', or 'minimap2' (PATHs allowed) with '--aligner'.")

    bwa_params = args.bwa_params +" -t " + str(num_threads)
    minimap2_params = args.minimap2_params +" -t " + str(num_threads)

    if genome_aligner == "bwa":       
        map_all = bwaAligner(idx_file, [align_path +file_prefix+".R1.fastq.gz",align_path +file_prefix+".R2.fastq.gz"],genome_aligner, bwa_params,align_path + file_prefix , in_overwrite=overwrite_files)
        print(map_all)
        map_all.run_aligner()

    else:
        if args.g=="yes":
            genome_aligner == "minimap2"
            minimap2_params = minimap2_params
            minimap2_params += " --split-prefix=tmp"
            map_all = Minimap2SAMAligner(idx_file, [align_path +file_prefix+".R1.fastq.gz",align_path +file_prefix+".R2.fastq.gz"],genome_aligner, minimap2_params,align_path + file_prefix , in_overwrite=overwrite_files)
            print(map_all)
            map_all.run_aligner()
        else:
            genome_aligner == "minimap2"
            minimap2_params = minimap2_params
            map_all = Minimap2SAMAligner(idx_file, [align_path +file_prefix+".R1.fastq.gz",align_path +file_prefix+".R2.fastq.gz"],genome_aligner, minimap2_params,align_path + file_prefix , in_overwrite=overwrite_files)
            print(map_all)
            map_all.run_aligner()

    # Sort the alignments.
    log("INFO", "Sort the alignments")
    run_samtools(file_prefix,align_path, num_threads, overwrite_files)

    # Check read pair and read pair orientation
    log("INFO", "Check read pair and read pair orientation")
    run_bedtoolss(file_prefix,align_path,output_path, num_threads, overwrite_files)

    # Detect split and discordant reads accordingly
    log("INFO", "Detect split reads")
    run_split(file_prefix,output_path, num_threads, overwrite_files)
    log("INFO", "Detect discordant reads")
    run_disc(file_prefix,output_path, num_threads, overwrite_files)
    # Detect siginificant sites of genomic enrichment.
    log("INFO", "Detecting sites of genomic enrichment")
    run_Genrich(file_prefix,align_path,output_path,num_threads, min_peak,max_dist,max_pvalue,overwrite_files)
    # Filtering eccDNA loci by boundary coverage.
    log("INFO", "Filtering eccDNA locus by number of split and discordant reads")
    run_intersect(file_prefix,output_path, num_threads, min_read,min_cov,overwrite_files)
    log("INFO", "Consider false positives from repetitive loci")
    run_repeat(file_prefix,output_path, num_threads, min_read,min_cov,overwrite_files)
    log("INFO", "Plotting size distribution of detected eccDNA")
    d=pd.read_csv(output_path +file_prefix+".csv", sep='\t',header=None,names=['chr','rstart','rend','num','unit','cov','len'])
    #bins = [100, 200, 400, 600,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]
    x=[k for k in d.len]
    plt.hist(x, color='orange')
    plt.ylabel('Count')
    plt.xlabel('Size distribution')
    plt.show(block=False)
    plt.savefig(output_path +file_prefix+".distribution.png")

    log("INFO", "Finished running ecc_finder")
    log("INFO", "Goodbye, have a nice day!")

def get_default_thread():
    return min(multiprocessing.cpu_count(), 8)

if __name__ == '__main__':
    main()