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

<img width="453" alt="ecc_finder_pipeline" src="https://user-images.githubusercontent.com/8072119/130788739-9d5a5fe8-f971-47c6-926a-8b8e3d1716b4.png">


## Getting Started
## <a name="install"></a>Quick installation using conda 

```bash

#Step 1. Download the latest ecc_finder:

git clone https://github.com/njaupan/ecc_finder.git

#Step 2. Go to ecc_finder folder

cd ecc_finder

#Step 3. Find the yaml file in the folder and run :

##For Linux64 system:

conda env create -f ecc_finder.yaml

##For MacOS system:

conda env create -f ecc_finder_MAC.yaml

conda activate ecc_finder


```
## Video about how to install and run ecc_finder

You can watch this video:  [ecc_finder_install_Video](https://www.youtube.com/watch?v=nOmYEpaas0w) about the installation and usuage of ecc_finder, which takes 2 minutes to install on a Mac.



## Index
```bash
# before starting, create index file for the reference genome to reduce mapping time.

# Long reads: 
minimap2 -x map-ont -d reference.ont.idx reference.fa -I 4G 

# Short reads:
# You can choose minimap2 short read mapping mode to speed up given a large genome 
# build index accordingly

bwa index reference.fa 

minimap2 -x sr -d reference.sr.idx reference.fa -I 4G 

```
## <a name="Usage"></a>Usage

```
ecc_finder: Tool for detecting eccDNA loci using Illumina and ONT sequencing.
Version: v1.0.0

usage: python ecc_finder.py <command> [options]
    
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

Run Example1: watch the video [Long-read-mapping_Video_example](https://www.youtube.com/watch?v=QUWFUZ19-yw) using the Arabidopsis eccDNA sequencing subsample in the folder test_samples.

```
usage: python ecc_finder.py map-ont <reference.idx> <query.fq> -r reference.fasta (option)

```

## <a name="Short-read-mapping"></a>Short-read-mapping

Algorithm and usage in details, please see the [Short-read-mapping](https://github.com/njaupan/ecc_finder/wiki/Short-read-mapping)

Run Example2: watch the video [Short-read-mapping_Video_example](https://www.youtube.com/watch?v=XnIvlHGYQvw) using the Arabidopsis eccDNA sequencing subsample in the folder test_samples.

Note that because the output from bwa index, it will create 5 files with appendix: reference.fa.bwt, reference.fa.amb, reference.fa.ann, reference.fa.pac and reference.fa.sa

So the <reference.fa> here represents the whole index output for bwa aligner, and "-r reference.fa" represents the reference fasta file

```
usage: 

python ecc_finder.py map-sr <reference.fa> <query.fq1> <query.fq2> -r reference.fa (option)

or to speed up:

python ecc_finder.py map-sr <reference.sr.idx> <query.fq1> <query.fq2> -r reference.fa --aligner minimap2 (option)

```


## <a name="Long-read-assembly"></a>Long-read-assembly

Algorithm and usage in details, please see the [Long-read-assembly](https://github.com/njaupan/ecc_finder/wiki/Long-read-assembly)


Run Example3: watch the video [Long-read-assembly_Video_example](https://www.youtube.com/watch?v=11O_lkUgOwY) using the Arabidopsis eccDNA sequencing subsample in the folder test_samples.

```
usage: python ecc_finder.py asm-ont <query.fq> (option)

```

## <a name="Short-read-assembly"></a>Short-read-assembly

Algorithm and usage in details, please see the [Short-read-assembly](https://github.com/njaupan/ecc_finder/wiki/Short-read-assembly)

Run Example4: watch the video [Short-read-assembly_Video_example](https://www.youtube.com/watch?v=vCyjzNrAYXo) using the Arabidopsis eccDNA sequencing subsample in the folder test_samples.

```
usage: python ecc_finder.py asm-sr <query.fq1> <query.fq2> (option)

```


## Output
All output is in `eccFinder_output` for mapping,  `eccFinder_asm_output` for assembly, or whichever directory `-o` specifies.

**Overview**

|Col|Type  |Description                                                  |
|--:|:----:|:------------------------------------------------------------|
|1  |folder|Folder containing alignment files (align_files)              |
|2  |folder|Folder containing peak calling files (peak_files)            |
|3  |file  |CSV file of the eccDNA loci: ecc.{ont,sr}.csv                |
|4  |file  |FASTA file of the eccDNA sequence: ecc.{ont,sr}.fasta        |
|5  |file  |BED file of genomic enriched sites: ecc.{ont,sr}.site.bed|
|6  |file  |Plot of the size distribution: ecc.{ont,sr}.distribution.png |

**ecc.{ont,sr}.fasta **

>The eccDNA sequence in FASTA format.

**ecc.{ont,sr}.csv **

>The eccDNA locus in csv format.


|Col|Type  |Description                               |
|--:|:----:|:-----------------------------------------|
|1  |string|Reference sequence name                   |
|2  |int   |Reference start on original strand        |
|3  |int   |Reference end on original strand          |
|4  |int   |Circular read number at the locus         |
|5  |int   |Repeat units of all circular reads        |
|6  |int   |EccDNA sequence length                    |


**ecc.{ont,sr}.distribution.png**

>The Size distribution of detected eccDNA in png format.

<img width="353" alt="Size_distribution" src="https://github.com/njaupan/ecc_finder/blob/main/log/Size_distribution.png">

## <a name="example"></a>Example
#PRJEB46420
https://www.ebi.ac.uk/ena/browser/view/PRJEB46420?show=reads
- Arabidopsis thaliana under heat stress 
- Triticum aestivum 

## Citation
Panpan Z, Haoran P, Christel L, Etienne B, Marie M. Ecc_finder: a robust and accurate tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data 
