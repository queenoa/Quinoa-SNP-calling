# Quinoa SNP Calling Scripts

Bioinformatics pipeline used for variant calling for quinoa resequencing data, based on the published pipeline by [Abrouk et al., 2020](https://doi.org/10.1038/s41467-020-18329-4), see also https://github.com/IBEXCluster/Wheat-SNPCaller.

## Overview

This repository contains the scripts used for processing raw sequencing reads through variant calling, generating the SNP dataset used in "First application of genomic prediction in quinoa using a statistical and a machine learning approach".

## Pipeline Components

- **pre-processing.sh** - Quality control and read trimming
- **merge_all_using_GATK.sh** - Merging VCF files using GATK
- **merge_all_using_GATK_allsites.sh** - Merging all sites (variant and non-variant)
- **gather_chro_vcf.sh** - Gathering chromosome-level VCF files
- **mergeme.sh** - Helper script for merging operations
- **mergeme_allsites.sh** - Helper script for merging all sites
- **downstream_analysis_generic_v2.0_ALL_allsites.sh** - Downstream variant analysis
- **runme.sh** - Main pipeline execution script

## Requirements

- GATK 4.1.8.0 and 4.0.1.1
- Trimmomatic 0.38
- BWA 0.7.17
- SAMtools 1.8
- Tabix 0.2.6
- Linux/Unix HPC environment with module system

## Usage
```bash
# Basic usage example
bash runme.sh
```

## Citation

If you use this pipeline, please cite:
[Abrouk et al., 2020](https://doi.org/10.1038/s41467-020-18329-4), https://github.com/IBEXCluster/Wheat-SNPCaller




