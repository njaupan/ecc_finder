# ecc_finder
A tool to detect eccDNA and chimeric eccDNA using Illumina and ONT sequencing

We previously developed the mobilome-seq that selectively sequence extrachromosomal circular DNA (eccDNA) forms of actively transposing transposable elements (TEs) in order to characterize active TEs in any plant or animal tissue (Lanciano et al.  PLoS Genetics, 2017).
The ecc_finder pipeline is dedicated to the detection of eccDNA and especially chimeric eccDNA in Mobilome-seq by Illumina and ONT sequencing.


<img width="510" alt="ecc_finder_pipeline" src="https://user-images.githubusercontent.com/8072119/124468507-7a81e500-dd99-11eb-8d52-e0d9f3bec733.png">

## Getting Started
## Install required dependencies

- [Minimap2](https://github.com/lh3/minimap2)
- Python 3 
- samtools
- tidehunter

```bash
# install with conda
conda create --name ecc_finder python=3.8
conda install -n ecc_finder -c bioconda -y minimap2 tidehunter samtools

```

## Docs
## Citation
Lanciano, S., Carpentier, M.C., Llauro, C., Jobet, E., Robakowska-Hyzorek, D., Lasserre, E., Ghesquière, A., Panaud, O. and Mirouze, M., 2017. Sequencing the extrachromosomal circular mobilome reveals retrotransposon activity in plants. PLoS genetics, 13(2), p.e1006630.
