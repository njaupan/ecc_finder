

# ecc_finder
A robust and accurate detection tool to detect eccDNA using Illumina and ONT sequencing
 
## Table of Contents

- [Introduction](#intro)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [Usage](#Usage)
  - [Example](#example)
- [Limitations](#limit)

## <a name="intro"></a>Introduction

We previously developed the mobilome-seq that selectively sequence extrachromosomal circular DNA (eccDNA) forms of actively transposing transposable elements (TEs) in order to characterize active TEs in any plant or animal tissue (Lanciano et al.  PLoS Genetics, 2017).

ecc_finder is dedicated to 
1) identify eccDNA from Nanopore reads; 
2) identify eccDNA from from Illumina paried end short reads; 
3) provide bona fide locus boundary of eccDNA-producing loci to investigate the origins of eccDNAs.

Validation 

1) To access ecc_finder accuracy for long read, the confidence score is assigned by tandem repeat pattern from read alignment. Except for satellites, when performing self-alignment, linear reads will not repeat itself while circular reads will be repeated two or more times because it goes through the rolling circle amplification experimentally. Therefore, the circular sequence will have a sub- read alignment in the same direction, and this sub-read alignment will be repeated two or more times on the same boundary. 

2) To access ecc_finder accuracy for short read, the confidence score for the eccDNA locus is assigned to the bona fide locus with an even distribution of split and discordant reads throughout their internal region.

<img width="441" alt="ecc_finder_pipeline" src="https://user-images.githubusercontent.com/8072119/124471419-20831e80-dd9d-11eb-89ce-49d5493764d5.png">


## Getting Started
## <a name="install"></a>Installation

- [Minimap2](https://github.com/lh3/minimap2)
- Python 3 
- samtools
- tidehunter

```bash
# install with conda
conda create --name ecc_finder python=3.8
conda install -n ecc_finder -c bioconda -y minimap2 tidehunter samtools genrich bwa seqtk

```
## <a name="Usage"></a>Usage
```
ecc_Finder: Tool for detecting eccDNA loci using Illumina and ONT sequencing.
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
- Long read mapping [map-ont](https://github.com/njaupan/ecc_finder/wiki/map-ont)
- Short read mapping [map-sr](https://github.com/njaupan/ecc_finder/map-sr) 
- Long read assembly [asm-ont](https://github.com/njaupan/ecc_finder/wiki/asm-ont)
- Short read assembly [asm-sr](https://github.com/njaupan/ecc_finder/asm-sr) 
 


## <a name="example"></a>Example
#PRJEB46420
https://www.ebi.ac.uk/ena/browser/view/PRJEB46420?show=reads
- Arabidopsis thaliana under heat stress 
- Triticum aestivum 

## Citation
Panpan Z, Haoran P, Christel L, Etienne B, Marie M. Ecc_finder: a robust and accurate tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data 

