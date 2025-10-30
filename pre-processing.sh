#!/bin/bash
#######################################################################################################################################################
### Trimmomatic + GATK 4 (Picard functionality) + GATK 3.8 HaplotypeCaller  + bgzip ( the *.g.VCF files) & tabix-index based workflow 
#    Naga, Version 2.0 dated 19 May 2019
#    Modification required for any users:
#      1. INPUT (directory location)
#      2. PROJECT (directory location)
#
########################################################################################################################################################
## Software Modules 
module load trimmomatic/0.38 bwa/0.7.17/gnu-6.4.0 samtools/1.8 gatk/4.0.1.1 tabix/0.2.6

## Reference files
export REF=/ibex/scratch/projects/c2071/1000quinoa/reference/quinoa_pb_chicago-2-final_PBJELLY_pilon_pseudo_RENAMED_CqUO_Cp_Mt.fasta

## Sample/Data Variables 
export PROJECT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/reseqlines
export INPUT=/ibex/scratch/projects/c2071/1000quinoa/SNPcallKatsamples/reseqlines/reads
export BAM=${PROJECT}/BAM;
export VCF=${PROJECT}/VCF;
export LOGS=${PROJECT}/LOGS;
export TOTAL_SAMPLE=`ls -lrta $INPUT/*1.fastq.gz | wc -l` 
mkdir -p $BAM;
mkdir -p $VCF;
mkdir -p $LOGS;

echo "========================================"
echo "Total Number of samples: ${TOTAL_SAMPLE}"
echo "========================================"

## For Sample count 
set COUNT=0;

## Workflow steps 
for SAMPLE in `ls $INPUT/*1.fastq.gz`;
do 
  PREFIX=`basename $SAMPLE _1.fastq.gz` ;
  LOCATION=${SAMPLE%/*};
  echo "$PREFIX and $LOCATION" 

 #### Step 1. trimming of reads
    MEM="32gb"
    CORES=4
    JOB1_NAME="Trimming"
    JOB1_TYPE="sbatch --partition=batch --job-name=${JOB1_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB1_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB1_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB1_CMD="time -p java -XX:+UseParallelGC -XX:ParallelGCThreads=${CORES} -jar $TRIMMOMATIC_JAR PE -phred33 $LOCATION/${PREFIX}_1.fastq.gz $LOCATION/${PREFIX}_2.fastq.gz $BAM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.up.1.fast $BAM/$PREFIX.trimmed.P2.fastq $BAM/$PREFIX.up.2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50" ;
    JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
    echo "$PREFIX sample with the job id=$JOB1_ID and Job Name=$JOB1_NAME submitted"
   
 #### Step 2. BWA MEM
    MEM="115gb"
    CORES=16
    JOB2_NAME="bwa-mem"
    JOB2_TYPE="sbatch --partition=batch --job-name=${JOB2_NAME}.${PREFIX} --time=12:00:00 --output=$LOGS/${JOB2_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB2_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB2_CMD="time -p bwa mem -M -k 30 -t $CORES $REF $BAM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.trimmed.P2.fastq | samtools view -@ $CORES -b -S -h -q 30 - | samtools sort - > $BAM/$PREFIX.sorted.bam"
    JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
    echo "$PREFIX sample with the job id=$JOB2_ID and Job Name=$JOB2_NAME submitted"   

 #### Step 3. MarkDuplicates 
    MEM="64gb"
    CORES=1
    JOB3_NAME="MarkDuplicate"
    JOB3_TYPE="sbatch --partition=batch --job-name=${JOB3_NAME}.${PREFIX} --time=4:00:00 --output=$LOGS/${JOB3_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB3_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB3_CMD="time -p gatk MarkDuplicates --INPUT=$BAM/$PREFIX.sorted.bam --METRICS_FILE=$BAM/$PREFIX.metrics.txt --OUTPUT=$BAM/$PREFIX.rmdup.bam"
    JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
    echo "$PREFIX sample with the job id=$JOB3_ID and Job Name=$JOB3_NAME submitted"   

 ##### 4. AddOrReplace
    ## Note: VALIDATION_STRINGENCY=LENIENT missing in GATK 4.0 
    MEM="32gb"
    CORES=1
    JOB4_NAME="AddOrReplace"
    JOB4_TYPE="sbatch --partition=batch --job-name=${JOB4_NAME}.${PREFIX} --time=3:00:00 --output=$LOGS/${JOB4_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB4_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB4_CMD="time -p gatk AddOrReplaceReadGroups --INPUT=$BAM/$PREFIX.rmdup.bam --OUTPUT=$BAM/$PREFIX.rgroup.bam --SORT_ORDER=coordinate --RGSM=$PREFIX --RGPU=none --RGID=1 --RGLB=lib1 --RGPL=Illumina"
    JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
    echo "$PREFIX sample with the job id=$JOB4_ID and Job Name=$JOB4_NAME submitted"   
    
 ##### 5. Samtools Index
    MEM="32gb"
    CORES=1
    JOB5_NAME="Samtool-Index"
    JOB5_TYPE="sbatch --partition=batch --job-name=${JOB5_NAME}.${PREFIX} --time=90:00 --output=$LOGS/${JOB5_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB5_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB5_CMD="time -p samtools index $BAM/$PREFIX.rgroup.bam"
    JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB4_ID} --wrap="${JOB5_CMD}");
    echo "$PREFIX sample with the job id=$JOB5_ID and Job Name=$JOB5_NAME submitted"   


##### 6. HaplotypeCaller
 MEM="115gb"
 CORES=8
 JOB6_NAME="HaplotypeCaller"
 JOB6_TYPE="sbatch --partition=batch --job-name=${JOB6_NAME}.${PREFIX} --time=5-00:00:00 --output=$LOGS/${JOB6_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB6_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
 #JOB6_CMD="time -p java -XX:+UseParallelGC -XX:ParallelGCThreads=${CORES} -jar $GATK -T HaplotypeCaller -nct $CORES -pairHMM VECTOR_LOGLESS_CACHING -R $REF -I $BAM/$PREFIX.rgroup.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 --heterozygosity 0.005 -o $VCF/$PREFIX.snps.indels.g.vcf" ;
 JOB6_CMD="time -p gatk HaplotypeCaller --reference $REF --input $BAM/$PREFIX.rgroup.bam --native-pair-hmm-threads $CORES --emit-ref-confidence GVCF --heterozygosity 0.005 --output $VCF/$PREFIX.snps.indels.g.vcf" ;
 JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
 echo "$PREFIX sample with the job id=$JOB6_ID and Job Name=$JOB6_NAME submitted"   
 
##### 7. Compress the g.VCF file using bgzip
 MEM="32gb"
 CORES=1
 JOB7_NAME="bgzip"
 JOB7_TYPE="sbatch --partition=batch --job-name=${JOB7_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB7_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB7_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
 JOB7_CMD="time -p bgzip -o $VCF/$PREFIX.snps.indels.g.vcf" ;
 JOB7_ID=$(${JOB7_TYPE} --parsable --dependency=afterok:${JOB6_ID} --wrap="${JOB7_CMD}");
 echo "$PREFIX sample with the job id=$JOB7_ID and Job Name=$JOB7_NAME submitted"   

##### 8. Create Tabix-Index for the g.VCF.gz file
 MEM="32gb"
 CORES=1
 JOB8_NAME="tabix"
 JOB8_TYPE="sbatch --partition=batch --job-name=${JOB8_NAME}.${PREFIX} --time=30:00 --output=$LOGS/${JOB8_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB8_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
 JOB8_CMD="time -p tabix -o $VCF/$PREFIX.snps.indels.g.vcf.gz" ;
 JOB8_ID=$(${JOB8_TYPE} --parsable --dependency=afterok:${JOB7_ID} --wrap="${JOB8_CMD}");
 echo "$PREFIX sample with the job id=$JOB8_ID and Job Name=$JOB8_NAME submitted"   
 
 COUNT=$((COUNT + 1))
done
 
echo "========================================"
echo "Total Number of samples submitted: ${COUNT} in 8 Steps"
echo "========================================"
