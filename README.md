
# ecc_finder
A robust and accurate detection tool to detect eccDNA using Illumina and ONT sequencing
 
## Table of Contents

- [Introduction](#intro)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [Usage](#Usage)
    - [Long-read-mapping](#Long-read-mapping)
    - [Short-read-mapping](#Short-read-mapping)
    - [Long-read-assembly](#Long-read-assembly)
    - [Short-read-assembly](#Short-read-assembly)
  - [Example](#example)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

ecc_finder is dedicated to 
1) identify eccDNA from Nanopore reads; 
2) identify eccDNA from from Illumina paried end short reads; 
3) provide bona fide locus boundary of eccDNA-producing loci to investigate the origins of eccDNAs.

ecc_finder works well on eccDNA-seq data (either mobilome-seq, Circle-Seq and CIDER-seq) from Arabidopsis, human, and wheat (with genome sizes ranging from 120 Mb to 17 Gb).

Validation 

1) To access ecc_finder accuracy for long read, the confidence score is assigned by tandem repeat pattern from read alignment. Except for satellites, when performing self-alignment, linear reads will not repeat itself while circular reads will be repeated two or more times because it goes through the rolling circle amplification experimentally. Therefore, the circular sequence will have a sub- read alignment in the same direction, and this sub-read alignment will be repeated two or more times on the same boundary. 

2) To access ecc_finder accuracy for short read, the confidence score for the eccDNA locus is assigned to the bona fide locus with an even distribution of split and discordant reads throughout their internal region.

<img width="433" alt="ecc_finder_pipeline" src="https://user-images.githubusercontent.com/8072119/130788739-9d5a5fe8-f971-47c6-926a-8b8e3d1716b4.png">


## Getting Started
## Dependency

- [Minimap2](https://github.com/lh3/minimap2)
- [bwa](https://github.com/lh3/bwa)
- [TideHunter](https://github.com/Xinglab/TideHunter)
- [Genrich](https://github.com/jsh58/Genrich)
- Python 3 

## <a name="install"></a>Installation 

```bash
# install with conda

conda create --name ecc_finder python=3.8
conda install -n ecc_finder -c bioconda -y minimap2 tidehunter samtools genrich bwa seqtk cd-hit

```
## Index
```bash
# before starting, create index file for the reference genome to reduce mapping time.

# Long reads: 
minimap2 -x map-ont -d reference.idx reference.fa -I 4G 

# Short reads:
# You can choose minimap2 short read mapping mode to speed up given a large genome 
# build index accordingly


bwa index reference.bwa.idx reference.fa 

minimap2 -x sr -d reference.mmi.idx reference.fa -I 4G 

```
## <a name="Usage"></a>Usage

```
ecc_finder: Tool for detecting eccDNA loci using Illumina and ONT sequencing.
Version: v1.0.0

usage: ecc_finder.py <command> [options]
    
    Mapping mode:
      map-sr         Call candidate eccDNA loci from paired-end short reads
      map-ont        Call candidate eccDNA loci from Nanopore long reads
    
    Assembly mode:
      asm-sr         Assembly from paired-end short reads
      asm-ont        Assembly from Nanopore long reads
      
    options:
      -v, --version
```

## Docs

Please see the [Wiki](https://github.com/njaupan/ecc_finder/wiki) for detailed documentation.


## <a name="Long-read-mapping"></a>Long-read-mapping

Algorithm and usage in details, please see the [Long-read-mapping](https://github.com/njaupan/ecc_finder/wiki/Long-read-mapping)

```
usage: ecc_finder.py map-ont <reference.idx> <query.fq> -r reference.fasta (option)

```

## <a name="Short-read-mapping"></a>Short-read-mapping

Algorithm and usage in details, please see the [Short-read-mapping](https://github.com/njaupan/ecc_finder/wiki/Short-read-mapping)

Note that, you can choose your favorite short-read aligner (bwa, bowtie2, segemehl or minimap2), and the default is bwa.

```
usage: ecc_finder.py map-sr <reference.idx> <query.fq1> <query.fq2> -r reference.fasta (option)

```
## <a name="Long-read-assembly"></a>Long-read-assembly

Algorithm and usage in details, please see the [Long-read-assembly](https://github.com/njaupan/ecc_finder/wiki/Long-read-assembly)

```
usage: ecc_finder.py asm-ont <query.fq> (option)

```

## <a name="Short-read-assembly"></a>Short-read-assembly

Algorithm and usage in details, please see the [Short-read-assembly](https://github.com/njaupan/ecc_finder/wiki/Short-read-assembly)

```
usage: ecc_finder.py asm-sr <query.fq1> <query.fq2> (option)

```


## Output
All output is in `eccFinder_output`, or whichever directory `-o` specifies.

**ecc_finder.fasta**

>The eccDNA sequence in FASTA format.

**ecc_finder.csv**

>The eccDNA locus in csv format.


|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Reference sequence name                   |
|2  |int   |Reference start on original strand        |
|3  |int   |Reference end on original strand          |
|4  |int   |Circular read number at the locus         |
|5  |int   |Repeat units of all circular reads        |
|6  |int   |Read coverage at the locus                |
|7  |int   |EccDNA sequence length                    |


**ecc_finder.png**

>The Size distribution of detected eccDNA in png format.

<img width="353" alt="Size_distribution" src="https://github.com/njaupan/ecc_finder/blob/main/log/Size_distribution.png">

## <a name="example"></a>Example
#PRJEB46420
https://www.ebi.ac.uk/ena/browser/view/PRJEB46420?show=reads
- Arabidopsis thaliana under heat stress 
- Triticum aestivum 

## Citation
Panpan Z, Haoran P, Christel L, Etienne B, Marie M. Ecc_finder: a robust and accurate tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data 
