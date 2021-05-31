#!/bin/bash

usage="$(basename "$0") [-h] -r <reference> -i <fastq>

Align Nanopore reads to a reference genome using minimap2.

    -h  show this help text.
    -r  reference, should be a fasta file. (required).
    -i  fastq/fastq.gz input reads (required).
    -P  filter to only primary alignments (i.e. run samtools view -F 2308, 
        defult: do not filter any alignment).
    -t  alignment threads (default: 4).
    -p  output file prefix (default: reads).
    -m  fill MD tag."

PREFIX="reads"
ALIGN_OPTS="-x map-ont"
THREADS=4
FILTER=""
rflag=false
iflag=false
filter_set=0
while getopts ':hr:i:p:Pm:t:' option; do
  case "$option" in
    h  ) echo "$usage" >&2; exit;;
    r  ) rflag=true; REFERENCE=$OPTARG;;
    i  ) iflag=true; INPUT=$OPTARG;;
    P  ) ((filter_set++)); FILTER="-F 2308";;
    m  ) ALIGN_OPTS="${ALIGN_OPTS} --MD";;
    p  ) PREFIX=$OPTARG;;
    t  ) THREADS=$OPTARG;;
    \? ) echo "Invalid option: -${OPTARG}." >&2; exit 1;;
    :  ) echo "Option -$OPTARG requires an argument." >&2; exit 1;;
  esac
done
shift $(($OPTIND - 1))

if ! $iflag || ! $rflag; then
  echo "$usage" >&2;
  echo "-i and -r must be specified." >&2;
  exit 1;
fi

bioawk -c fastx '{ print $name, length($seq) }' ${INPUT} |sort -k2,2nr > ${INPUT}.txt
minimap2 ${ALIGN_OPTS} -t ${THREADS} -a ${REFERENCE} ${INPUT} |
  samtools view -@ ${THREADS} -bS ${FILTER}  -  > ${PREFIX}.mapont.bam
  samtools sort -@ ${THREADS} ${PREFIX}.mapont.bam ${PREFIX}.mapont.sort
rm ${PREFIX}.mapont.bam
samtools index ${PREFIX}.mapont.sort.bam 
bedtools genomecov -bga -split -ibam ${PREFIX}.mapont.sort.bam > ${PREFIX}.mapont.sort.bam.cov.txt
awk '{x+=$4}END {print x/NR}' ${PREFIX}.mapont.sort.bam.cov.txt
