#!/bin/bash
#########################################################################################################################################################
# Downstream Analysis, version 2.0															#
#  Last update: May 1, 2020																#
#########################################################################################################################################################
## Software Modules 
module load gatk/4.1.8.0

REF=$1;
INPUT_DIR=$2;
PROJECT=$3;
if [ "$#" -eq 3 ];
then
  if [ ! -f $REF ]; then
    printf " Your Genome reference file does not exist \n" ;
    break;
  fi;
  if [ ! -d $INPUT_DIR ];then
    printf " Your INPUT [ $INPUT_DIR ] directory does not exist \n" ;
    break;
  fi
  if [ -d $PROJECT ];then
   echo "###################### WARNING! "
   echo " Please delete $PROJECT and RERUN the script"
  fi
else
 printf "\033c"
 echo " ***************************************************************************************************************************"
 echo ""
 echo " Run this script with 3 arguments:" 
 echo "      (a) Your Reference file"
 echo "      (b) Your gVCF files directory and "
 echo "      (c) Your output file directory "
 echo "  (absulute path required for all these options) "
 echo ""
 echo ""
 echo "  Example: "
 echo "    ./downstream_analysis_generic_v1.0.sh /ibex/reference/KSL/human_ref/human_g1k_v37.fasta  /ibex/scratch/reyel/1000quinoa/naga/VCF/ /ibex/scratch/reyel/1000quinoa/naga/phase2/downstream_analysis_testing/OUTPUT"
 echo ""
 echo ""
 echo " ***************************************************************************************************************************"
 exit;
fi

## Sample/Data Variables 
export gVCF=${PROJECT}/gVCF;
export logs=${PROJECT}/logs;
export SNPs=${PROJECT}/SNPs;
export INDELs=${PROJECT}/INDELs;
mkdir -p $gVCF;
mkdir -p $SNPs;
mkdir -p $INDELs;
mkdir -p $logs;


## Calculate the Chromosome split 
MaxSize=0;
MaxJobs=2000;
JobSteps=6;

# Find the Max. size of Chromosome 
while IFS=$'\t' read -r -a myREF
do 
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};
 if [ $MaxSize -lt $ChrLen ]; then
  MaxSize=$ChrLen
 fi
done < $REF.fai

# Find the total numbers of Chr
TotalChr=`cat $REF.fai | wc -l `

# Find the best possible Chunk size 
Divide=1;
JobSteps=6;
AssinedJobs=$(( $Divide * $TotalChr * $JobSteps )) ;
echo "$Divide and $AssinedJobs"
while [ ${MaxJobs} -gt $AssinedJobs ];
 do
   Divide=$(( $Divide + 1 ));
   AssinedJobs=$(( $Divide * $TotalChr * $JobSteps )) ;
 done
Divide=$(( Divide - 2 ));
Chunk=$(( $MaxSize / $Divide ));
# echo "total Chr = $TotalChr, Max Chr length = $MaxSize, Chunk number = $Divide, Chunk size = $Chunk "


# Store the all the Chromosome Chunks for reference in a file "submitted_chunks_data.txt"
rm -Rf submitted_chunks_data.txt; 
while IFS=$'\t' read -r -a myREF
do
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};
 Part=$(( $ChrLen / $Chunk )) ;
 tmp=$(( $Part * $Chunk )) ;
 if [ $ChrLen -gt $tmp ]; then
    Part=$(( $Part + 1 ));
    echo "$ChrName split into $Part parts" >> submitted_chunks_data.txt
 else 
    echo "$ChrName split into $Part parts" >> submitted_chunks_data.txt
 fi
done < $REF.fai


###


## Add all the *.g.vcf files for Genotyping 
set INPUT_TMP
for i in `ls -l ${INPUT_DIR}/*.g.vcf.gz | awk '{print $9}'`
do
 INPUT_TMP+="$i -V "; 
done
INPUT=${INPUT_TMP::-4}
 #echo $INPUT
 MEM="115gb"
 CORES=4;
## READ ALL CHROMOSOME and LENGTH from REFERENCE INDEX file 
while IFS=$'\t' read -r -a myREF
do 
 size=1;
 DefSize=$Chunk;
 ChunkSize=$DefSize;
 ChrName=${myREF[0]};
 ChrLen=${myREF[1]};
 echo "******************************************************************************************"
 echo "Chromosome: $ChrName and the size: $ChrLen"
 echo "Chromosome: $ChrName will be split into many chunks for data parallelization"
 echo "Paralle Data split (Chromosome by Chromosome) may takes time ....!"
 echo "Many jobs will be submitted to the scheduler ...!!"
 echo "The following jobs are submitted to the scheduler for the Chromosome: $ChrName"
 echo "******************************************************************************************"

 ## PREPARE CHUNKS for EACH CHROMOSOME 
  for ((Start=1; Start<$ChrLen; Start+=$DefSize))
   do 
     End=$(( $size*$ChunkSize ))
     if [ $End -lt $ChrLen ]
     then
      ## IMPORT Genome into GenomicDB (Chunk by Chunk)
       JOB1_NAME="$ChrName.chunk_$size.GenomicsDBImport"
       JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME} --time=10-00:00:00 --output=$logs/${JOB1_NAME}.%J.out --error=$logs/${JOB1_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB1_CMD="time -p gatk GenomicsDBImport --variant $INPUT --genomicsdb-workspace-path $gVCF/$ChrName.$size --intervals $ChrName:$Start-$End --reader-threads $CORES";
       JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
   echo "$JOB1_NAME with the job id=$JOB1_ID submitted";
      ## CALL GenotypeGVF (Chunk by Chunk)
       JOB2_NAME="$ChrName.chunk_$size.GenotypeGVCFs"
       JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME} --time=10-00:00:00 --output=$logs/${JOB2_NAME}.%J.out --error=$logs/${JOB2_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB2_CMD="time -p gatk GenotypeGVCFs --variant gendb://$gVCF/$ChrName.$size --reference $REF --include-non-variant-sites --output $gVCF/$ChrName.$size.vcf.gz --intervals $ChrName:$Start-$End";
       JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
   echo "$JOB2_NAME with the job id=$JOB2_ID submitted";
      ## Select VARIANT=SNPs (Chunk by Chunk)
       JOB3_NAME="$ChrName.chunk_$size.SNPs.SelectVariants"
       JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME} --time=3-00:00:00 --output=$logs/${JOB3_NAME}.%J.out --error=$logs/${JOB3_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB3_CMD="time -p gatk SelectVariants --variant $gVCF/$ChrName.$size.vcf.gz --reference $REF -select-type SNP --output $SNPs/$ChrName.$size.vcf" ;
       JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
   echo "$JOB3_NAME with the job id=$JOB3_ID submitted";
      ## Call SNPs Filtration (Chunk by Chunk)
       JOB4_NAME="$ChrName.chunk_$size.SNPs.VariantFiltration"
       JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME} --time=3-00:00:00 --output=$logs/${JOB4_NAME}.%J.out --error=$logs/${JOB4_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB4_CMD="time -p gatk VariantFiltration --variant $SNPs/$ChrName.$size.vcf --reference $REF --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0\" --cluster-size 3 --cluster-window-size 10 --filter-name snp_filter --output $SNPs/hard_filtered_snps.$ChrName.$size.vcf" ;
       JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
   echo "$JOB4_NAME with the job id=$JOB4_ID submitted";
      ## Select VARIANT=INDELs (Chunk by Chunk)
       JOB5_NAME="$ChrName.chunk_$size.INDELs.SelectVariants"
       JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME} --time=2-00:00:00 --output=$logs/${JOB5_NAME}.%J.out --error=$logs/${JOB5_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB5_CMD="time -p gatk SelectVariants --variant $gVCF/$ChrName.$size.vcf.gz --reference $REF -select-type INDEL --output $INDELs/$ChrName.$size.vcf" ;
       JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB5_CMD}");
   echo "$JOB5_NAME with the job id=$JOB5_ID submitted";
      ## Call INDELs Filtration (Chunk by Chunk)
       JOB6_NAME="$ChrName.chunk_$size.INDELs.VariantFiltration"
       JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME} --time=2-00:00:00 --output=$logs/${JOB6_NAME}.%J.out --error=$logs/${JOB6_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB6_CMD="time -p gatk VariantFiltration --variant $INDELs/$ChrName.$size.vcf --reference $REF --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\" --filter-name indel_filter  --output $INDELs/hard_filtered_indels.$ChrName.$size.vcf" ; 
       JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
       echo "$JOB6_NAME with the job id=$JOB6_ID submitted";
     ## PREPARE for Next CHUNK
     size=$(( $size + 1 ));
    else   ### Last chunk 
      Start=$(( $End - $ChunkSize +1 ));
      End=$ChrLen;
      ## IMPORT Genome into GenomicDB (Chunk by Chunk)
       JOB1_NAME="$ChrName.chunk_$size.GenomicsDBImport"
       JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME} --time=10-00:00:00 --output=$logs/${JOB1_NAME}.%J.out --error=$logs/${JOB1_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB1_CMD="time -p gatk GenomicsDBImport --variant $INPUT --genomicsdb-workspace-path $gVCF/$ChrName.$size --intervals $ChrName:$Start-$End --reader-threads $CORES";
       JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
       echo "$JOB1_NAME with the job id=$JOB1_ID submitted";
      ## CALL GenotypeGVF (Chunk by Chunk)
       JOB2_NAME="$ChrName.chunk_$size.GenotypeGVCFs"
       JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME} --time=10-00:00:00 --output=$logs/${JOB2_NAME}.%J.out --error=$logs/${JOB2_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB2_CMD="time -p gatk GenotypeGVCFs --variant gendb://$gVCF/$ChrName.$size --reference $REF --include-non-variant-sites --output $gVCF/$ChrName.$size.vcf.gz --intervals $ChrName:$Start-$End";
       JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
       echo "$JOB2_NAME with the job id=$JOB2_ID submitted";
      ## Select VARIANT=SNPs (Chunk by Chunk)
       JOB3_NAME="$ChrName.chunk_$size.SNPs.SelectVariants"
       JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME} --time=3-00:00:00 --output=$logs/${JOB3_NAME}.%J.out --error=$logs/${JOB3_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB3_CMD="time -p gatk SelectVariants --variant $gVCF/$ChrName.$size.vcf.gz --reference $REF -select-type SNP --output $SNPs/$ChrName.$size.vcf" ;
       JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
       echo "$JOB3_NAME with the job id=$JOB3_ID submitted";
      ## Call SNPs Filteration (Chunk by Chunk)
       JOB4_NAME="$ChrName.chunk_$size.SNPs.VariantFiltration"
       JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME} --time=3-00:00:00 --output=$logs/${JOB4_NAME}.%J.out --error=$logs/${JOB4_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB4_CMD="time -p gatk VariantFiltration --variant $SNPs/$ChrName.$size.vcf --reference $REF --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.00 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0\" --cluster-size 3 --cluster-window-size 10 --filter-name snp_filter  --output $SNPs/hard_filtered_snps.$ChrName.$size.vcf" ;
       JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
       echo "$JOB4_NAME with the job id=$JOB4_ID submitted";
      ## Select VARIANT=INDELs (Chunk by Chunk)
       JOB5_NAME="$ChrName.chunk_$size.INDELs.SelectVariants"
       JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME} --time=2-00:00:00 --output=$logs/${JOB5_NAME}.%J.out --error=$logs/${JOB5_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB5_CMD="time -p gatk SelectVariants --variant $gVCF/$ChrName.$size.vcf.gz --reference $REF -select-type INDEL --output $INDELs/$ChrName.$size.vcf" ;
       JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB5_CMD}");
       echo "$JOB5_NAME with the job id=$JOB5_ID submitted";
      ## Call INDELs Filteration (Chunk by Chunk)
       JOB6_NAME="$ChrName.chunk_$size.INDELs.VariantFiltration"
       JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME} --time=2-00:00:00 --output=$logs/${JOB6_NAME}.%J.out --error=$logs/${JOB6_NAME}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
       JOB6_CMD="time -p gatk VariantFiltration --variant $INDELs/$ChrName.$size.vcf --reference $REF --filter-expression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\" --filter-name indel_filter  --output $INDELs/hard_filtered_indels.$ChrName.$size.vcf" ; 
       JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
       echo "$JOB6_NAME with the job id=$JOB6_ID submitted";
    fi ## END OF ALL CHUNKS within the specific CHROMOSOME 
 #  echo "Preparing for the next Chunk:" $ChrName.chunk_$size 
done ## END OF CHROMOSOME (Chr by Chr)
done < $REF.fai
  ## END OF ALL CHROMOSOMEs
