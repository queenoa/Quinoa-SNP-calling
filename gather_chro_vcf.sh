#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=GatherVCF
#SBATCH --constraint=intel
#SBATCH --nodes=1
##SBATCH --array=1
#SBATCH --mem=115G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --err=GatherVCF-%A_%a.err
#SBATCH --account=c2071

module load gatk/4.1.8.0

gatk GatherVcfs -R /ibex/scratch/projects/c2071/1000quinoa/reference/quinoa_pb_chicago-2-final_PBJELLY_pilon_pseudo_RENAMED_CqUO_Cp_Mt.fasta \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq1B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq1A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq4A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq5B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq5A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq6A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq7A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq2A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq2B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq7B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq4B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq3A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq9A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq3B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq9B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq8B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq6B.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/Cq8A.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/CqUO.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/MK159176.1.SNPs.exclude_filtered.vcf.gz \
-I /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/MK182703.1.SNPs.exclude_filtered.vcf.gz \
-O /ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants/SNPs/CqAllChro_Cp_Mt.SNPs.exclude_filtered.vcf.gz
