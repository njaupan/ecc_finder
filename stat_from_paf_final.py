#!/usr/bin/env python3

import re
import sys
import argparse 
import pandas as pd
from itertools import groupby
from operator import itemgetter
from pybedtools import BedTool

parser = argparse.ArgumentParser(description="Statstics of minimap2 alignment results(.paf files)")
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
    headers = ['refID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','direction','allmatch','blast_iden']   
    #outfile.write('\t'.join(headers))
    #outfile.write('\n')
    total=0
    n_primary=0
    for line in paf:         
        parts = line.strip().split("\t")
        total=total+1
        blast_iden =  100.0 *int((parts[9]))/ int((parts[10]))
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
		 "blast_iden": blast_iden,
		 #"qcov": qcov,
		 #"rcov": rcov,
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
#merge overlapped regions
df=pd.read_csv('paf.stats', skiprows=0,sep='\t',names = ['refID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','direction','allmatch','blast_iden'])
if df.empty:
    print('DataFrame is empty!')
else:
    df['ID']=df['refID']+ '__'+ df['queryID']
    df=df[['ID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','blast_iden'] ] 
    x = BedTool.from_dataframe(df)
    sorted =BedTool.sort(x)
    y =BedTool.merge(sorted, c = '4,5,6,7,8,9',
    o = 'distinct,min,max,distinct,distinct,max',d=1000000)
    df=BedTool.to_dataframe(y,disable_auto_names=True, header=None, 
    names=['ID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','blast_iden'])
    #sum up query coverage
    df['refID']=df['ID'].str.split('__', expand = True)[0]
    df['queryID']=df['ID'].str.split('__', expand = True)[1]
    df=df[['refID','rstart', 'rend','queryID','qstart', 'qend', 'qlen','rlen','blast_iden']]
    df['qcov']=100*(df['qend']-df['qstart'])/df['qlen']
    df['rcov']=100*(df['rend']-df['rstart'])/df['rlen']
    df1=df.groupby(['refID','queryID'])['blast_iden'].mean().reset_index(name='blast_iden_mean').round(2)
    df2=df.groupby(['refID','queryID'])['qcov'].sum().reset_index(name='qcov_sum').round(1)
    df3=df.groupby(['refID','queryID'])['rcov'].sum().reset_index(name='rcov_sum').round(1)
    df=pd.merge(df, df1, on=['refID','queryID'])
    df=pd.merge(df, df2, on=['refID','queryID'])
    df=pd.merge(df, df3, on=['refID','queryID'])  
    df=df.sort_values(by=['rlen'],ascending=False)
    print(df)
    df.to_csv('paf.stats.sum', header=False, index = False, sep='\t')  
    print('make_sum is done')
#split each dictionary based on key
sep = '\t'
with open('paf.stats.sum') as f:
    data = f.read()
split_data = [row.split(sep) for row in data.split('\n')]
gb = groupby(split_data, key=itemgetter(0))
for index, (key, group) in enumerate(gb):
    with open('fresult{}.coords'.format(index), 'w') as f:
        write_data = '\n'.join(sep.join(cell) for cell in group)
        f.write(write_data)
        f.close()
