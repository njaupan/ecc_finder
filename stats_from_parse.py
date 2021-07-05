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
import re
import sys
import argparse 
import pandas as pd
from itertools import groupby
from operator import itemgetter
from pybedtools import BedTool

parser = argparse.ArgumentParser(description="Calculate the statstics of alignment results")
parser.add_argument("-i","--paf",help="paf file name")
parser.add_argument('-m',type=int,help="min align length")
parser.add_argument('-mq', type=int,help="min query length")
args=parser.parse_args()
# Define the search patterns:
pafFileNamePattern = re.compile('(.*)\.paf')

def statsFromPaf(pafFile):
    paf=open(pafFile)
    OutfilePrefix1 =pafFileNamePattern.match(pafFile) # define outfile 
    OutfilePrefix=OutfilePrefix1.group(1)      
    outfile=open( 'paf.stats','w') 
    headers = ['refID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','iden']   
    #outfile.write('\t'.join(headers))
    #outfile.write('\n')
    total=0
    n_primary=0
    for line in paf:         
        parts = line.strip().split("\t")
        total=total+1
        iden =  100.0 *int((parts[9]))/ int((parts[10]))
        #qcov =  100.0 *(int((parts[3]))-int((parts[2])))/ int((parts[1]))
        #rcov =  100.0 *(int((parts[8]))-int((parts[7])))/ int((parts[6]))
        resultss = {
		 "queryID": parts[0],
		 "qlen":  int(parts[1]),
		 "qstart": int(parts[2]),
		 "qend": int(parts[3]),
		 "direction": parts[4],
		 "refID": parts[5],
		 "rlen": int(parts[6]),
		 "rstart": int(parts[7]),
		 "rend": int(parts[8]),
		 "allmatch": int(parts[9]),
		 "iden": iden,

         }
        if args.m is not None and resultss['allmatch'] < args.m:
            n_primary=n_primary+1
            continue
        if args.mq is not None and resultss['qlen'] < args.mq:
            continue
        out_row = (str(resultss[x]) for x in headers)
        outfile.write('\t'.join(out_row))
        outfile.write('\n')
statsFromPaf(str(args.paf))
