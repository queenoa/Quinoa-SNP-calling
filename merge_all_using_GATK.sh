#!/bin/bash
module load gatk/4.1.8.0

export CHR=$1;
export CHUNK=$2;

export REF=/ibex/scratch/projects/c2071/1000quinoa/reference/quinoa_pb_chicago-2-final_PBJELLY_pilon_pseudo_RENAMED_CqUO_Cp_Mt.fasta
export LOCATION=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES
export OUTPUT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_variants
mkdir -p $OUTPUT/SNPs ;
mkdir -p $OUTPUT/INDELs ;
mkdir -p $OUTPUT/LOGs ;

set INPUT_SNP
set INPUT_INDEL
for i in `seq 1 $CHUNK`; 
 do
  INPUT_SNP+="-I $LOCATION/SNPs/hard_filtered_snps.$CHR.$i.vcf "; 
  INPUT_INDEL+="-I $LOCATION/INDELs/hard_filtered_indels.$CHR.$i.vcf "; 
done
echo $INPUT_SNP
echo $INPUT_INDEL
gatk GatherVcfs ${INPUT_SNP} -O $OUTPUT/SNPs/$CHR.SNPs.vcf.gz -R $REF
gatk SelectVariants --variant $OUTPUT/SNPs/$CHR.SNPs.vcf.gz --reference $REF -select-type SNP --exclude-filtered true --output $OUTPUT/SNPs/$CHR.SNPs.exclude_filtered.vcf.gz
gatk GatherVcfs ${INPUT_INDEL} -O $OUTPUT/INDELs/$CHR.INDELs.vcf.gz -R $REF
gatk SelectVariants --variant $OUTPUT/INDELs/$CHR.INDELs.vcf.gz --reference $REF -select-type INDEL --exclude-filtered true --output $OUTPUT/INDELs/$CHR.INDELs.exclude_filtered.vcf.gz
