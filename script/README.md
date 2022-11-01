# Analysis of pangraph performance

This folder contains all the scripts used to perform the analysis of the performances of pangraph on real and simulated data reported in [our paper](https://github.com/neherlab/pangraph#citing).

The analysis can be reproduced using [Snakemake](https://snakemake.readthedocs.io/en/stable/).

## Dependencies

Running the pipeline requires:

- having the `pangraph` binary available in your path
- having a working [conda](https://docs.conda.io/en/latest/) installation
- having a conda environment with [Snakemake](https://snakemake.readthedocs.io/en/stable/) available.

Moreover the input data needs to be downloaded and prepared as explained in the next section. The analysis and all of the paper's figures can then be reproduced with the command:

```bash
snakemake --use-conda -c4 PF_all
```

N.b.: this will require a high amount of computational time and resources. We also provide an option to submit jobs using SLURM when running the analysis on a cluster:

```bash
snakemake --profile cluster PF_all
```

In this case the name of queues and partitions might have to be adjusted depending on your cluster setup. These are defined in the file `cluster/cluster_config.json`

## Input data preparation

The input data consists of genbank files from RefSeq. A list of accession number for each species can be found in `config/accnums.json`.

For each species, we performed a preliminary analsysis of the data using [PanX](https://github.com/neherlab/pan-genome-analysis). This pipeline groups families of homologous genes. The results of the analysis are stored in a 

PanX stores the results of the analysis in a nested folder structure:

```
data/species/
├── allclusters_final.tsv
├── input_GenBank
│   ├── NC_005042.gbk
│   ├── ...
│   └── NC_009840.gbk
└── vis
    ├── all_gene_alignments.tar.gz
    ├── core_gene_alignments.tar.gz
    ├── coreGenomeTree.json
    ├── geneCluster
    │   ├── GC00000001_11_aa_aln.fa.gz
    │   ├── GC00000001_11_aa_aln_reduced.fa.gz
    │   ├── GC00000001_11_na_aln.fa.gz
    │   ├── GC00000001_11_na_aln_reduced.fa.gz
    │   ├── GC00000001_11.nwk
    │   └── ...
    ├── geneCluster.json
    ├── metaConfiguration.js
    ├── metainfo.tsv
    ├── strainMetainfo.json
    └── strain_tree.nwk
```

For our analysis, we store PanX output for each species in a directory named `panx_data`. This will be used as input data in our analaysis pipeline.

```
panx_data/
├── escherichia_coli
│   ├── allclusters_final.tsv
│   ├── input_GenBank
│   └── vis
├── helicobacter_pylori
│   ├── allclusters_final.tsv
│   ├── input_GenBank
│   └── vis
├── klebsiella_pneumoniae
│   ├── allclusters_final.tsv
│   ├── input_GenBank
│   └── vis
├── mycobacterium_tuberculosis
│   ├── allclusters_final.tsv
│   ├── input_GenBank
│   └── vis
└── prochlorococcus_marinus
    ├── allclusters_final.tsv
    ├── input_GenBank
    └── vis
```