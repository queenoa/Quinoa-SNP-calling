#!/bin/bash
module load gatk/4.1.8.0

export CHR=$1;
export CHUNK=$2;

export REF=/ibex/scratch/projects/c2071/1000quinoa/reference/quinoa_pb_chicago-2-final_PBJELLY_pilon_pseudo_RENAMED_CqUO_Cp_Mt.fasta
export LOCATION=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES
export OUTPUT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/OUTPUT_ALL_SITES/merge_allsites
mkdir -p $OUTPUT/VCF ;
mkdir -p $OUTPUT/LOGs ;

set INPUT
for i in `seq 1 $CHUNK`; 
 do
  INPUT+="-I $LOCATION/gVCF/$CHR.$i.vcf.gz "; 
done
echo $INPUT
gatk GatherVcfs ${INPUT} -O $OUTPUT/VCF/$CHR.allsites.vcf.gz -R $REF
