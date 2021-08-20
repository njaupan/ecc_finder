

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

ecc_finder is dedicated to 1) identify eccDNA from Nanopore reads, 2) identify eccDNA from from Illumina paried end short reads, 3) provide clear boundary of eccDNA-producing loci to investigate the origins of eccDNAs.

<img width="641" alt="ecc_finder_pipeline" src="https://user-images.githubusercontent.com/8072119/124471419-20831e80-dd9d-11eb-89ce-49d5493764d5.png">



## Getting Started
## <a name="install"></a>Installation

- [Minimap2](https://github.com/lh3/minimap2)
- Python 3 
- samtools
- tidehunter

```bash
# install with conda
conda create --name ecc_finder python=3.8
conda install -n ecc_finder -c bioconda -y minimap2 tidehunter samtools

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
      -c, --citation  
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

## Docs
## Citation
Lanciano, S., Carpentier, M.C., Llauro, C., Jobet, E., Robakowska-Hyzorek, D., Lasserre, E., Ghesquière, A., Panaud, O. and Mirouze, M., 2017. Sequencing the extrachromosomal circular mobilome reveals retrotransposon activity in plants. PLoS genetics, 13(2), p.e1006630.
