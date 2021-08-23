#!/usr/bin/env python

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Split the fasta file into individual file with each gene seq")
parser.add_argument('-i', action='store', dest='fasta_file', help='Input fasta file')
result = parser.parse_args()

f_open = open(result.fasta_file, "r")

for record in SeqIO.parse(f_open, "fasta"):
   id = record.id+".fasta"
   seq = record.seq
   id_file = open(id, "w")
   id_file.write(">"+str(record.id)+"\n"+str(seq))
   id_file.close()

f_open.close()

'''
from Bio import SeqIO

original_file = "./original.fasta"
corrected_file = "./corrected.fasta"

with open(original_file) as original, open(corrected_file, 'w') as corrected:
    records = SeqIO.parse(original_file, 'fasta')
    for record in records:
        print record.id            
        if record.id == 'foo':
            record.id = 'bar'
            record.description = 'bar' # <- Add this line
        print record.id 
        SeqIO.write(record, corrected, 'fasta')
'''


