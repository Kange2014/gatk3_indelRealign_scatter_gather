#!/usr/bin/env bash
#
# gatk3_indelRealign_scatter_gather.sh by Lucius Zheng, 2022/06/09
#
# Scatter-process-gather execution pattern: it will split the
# input into multiple pieces, each of which will be processed in
# parallel, after which they are gathered together in some final
# output.
#

JAVA=/beegfs/work/commercial_test/cupcake/softwares/jdk-1.8.0_211/bin/java
GATK=/beegfs/work/commercial_test/cupcake/softwares/gatk-3.7.0/GenomeAnalysisTK.jar
GATK4=/beegfs/work/commercial_test/cupcake/softwares/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
GATK_BUNDLE=/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19
DATABASE_INDEL=${GATK_BUNDLE}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
DATABASE_SNP=${GATK_BUNDLE}/dbsnp_138.hg19.vcf
PHASE1_INDEL=${GATK_BUNDLE}/1000G_phase1.indels.hg19.sites.vcf
REF_FILE_NEW=/beegfs/work/commercial_test/cupcake/databases/gatk_bundle/2.8/hg19/ucsc.hg19.noconfig.fasta
PANEL_BED=/beegfs/work/commercial_test/cupcake/codes/v2.14.5/database/panel_target_regions/panel15_pro.bed
JAVA_TMP_DIR=/beegfs-ssd/tmp

OUT_DIR="./"
SAMPLE_ID="test"
TIME_LOG=$OUT_DIR/timing.log

NUM_THREADS=8

# test different parallel garbage collection threads

#JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=4
#JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=8"
#JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=8 -XX:+CMSParallelRemarkEnabled"
#JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=4"

# if ssd is available, set java temporary I/O directory to JAVA_TMP_DIR 
#JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREADS} -Djava.io.tmpdir=${JAVA_TMP_DIR}"
JAVA_OPTS="-Xms512M -Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads=${NUM_THREADS}"

# split the panel15_pro target regions to roughly equal-size 4 chunks
#GRP_1='-XL chr1 -XL chr2 -XL chr3 -XL chr4'
#GRP_2='-XL chr5 -XL chr6 -XL chr7 -XL chr8 -XL chr9'
#GRP_3='-XL chr10 -XL chr11 -XL chr12 -XL chr13 -XL chr14 -XL chr15 -XL chr16'
#GRP_4='-XL chr17 -XL chr18 -XL chr19 -XL chr20 -XL chr21 -XL chr22 -XL chrX -XL chrY -XL chrM'

GRP_1='-L chr1 -L chr2 -L chr3 -L chr4'
GRP_2='-L chr5 -L chr6 -L chr7 -L chr8 -L chr9'
GRP_3='-L chr10 -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 -L chr16'
GRP_4='-L chr17 -L chr18 -L chr19 -L chr20 -L chr21 -L chr22 -L chrX -L chrY -L chrM'
#EXCLUDE_GRP='-XL chr1 -XL chr2 -XL chr3 -XL chr4 -XL chr5 -XL chr6 -XL chr7 -XL chr8 -XL chr9 -XL chr10 -XL chr11 -XL chr12 -XL chr13 -XL chr14 -XL chr15 -XL chr16 -XL chr17 -XL chr18 -XL chr19 -XL chr20 -XL chr21 -XL chr22 -XL chrX -XL chrY -XL chrM'

# for samtools use
#GRP_1='chr1 chr2 chr3 chr4'
#GRP_2='chr5 chr6 chr7 chr8 chr9'
#GRP_3='chr10 chr11 chr12 chr13 chr14 chr15 chr16'
#GRP_4='chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM'

# first split bam by chromosomes

START=`date +%s`

for i in {1..4}
do
    GRP=GRP_${i}
    #samtools view -h -b -@ ${NUM_THREADS} ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.bam ${!GRP} > ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.bam
    #samtools index -@ ${NUM_THREADS} ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.bam
    
    # it seems gatk4 PrintReads runs faster: 74 vs 126 seconds for GRP_1 region
    java ${JAVA_OPTS} -jar ${GATK4} PrintReads \
        --input ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.bam \
        --output ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.bam \
        ${!GRP} &
        
done

wait

# RealignerTargetCreator uses -nt/--num_threads to control the number of data 
# threads sent to the processor (acting at the machine level)

for i in {1..4}
do
    ${JAVA} ${JAVA_OPTS} -jar ${GATK} \
        -T RealignerTargetCreator \
        -R ${REF_FILE_NEW} \
        -I ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.bam \
        -o ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.intervals \
        -known ${PHASE1_INDEL} -known ${DATABASE_INDEL} \
        -nt ${NUM_THREADS} \
        -L ${PANEL_BED}
done

END=`date +%s`

echo RealignerTargetCreator ran for $((END - START)) seconds >> $TIME_LOG

START=`date +%s`

# IndelRealigner -> no parallelism, use scatter-gather
# use & and wait for parallel run 

# do not specify -L due to reads discarding

for i in {1..4}
do

    ${JAVA} ${JAVA_OPTS} -jar ${GATK} -T IndelRealigner \
        -R ${REF_FILE_NEW} \
        -I ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.bam \
        -targetIntervals ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.${i}.intervals \
        -o ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.${i}.bam \
        -known ${PHASE1_INDEL} -known ${DATABASE_INDEL} &

done


# ${JAVA} ${JAVA_OPTS} -jar ${GATK} -T IndelRealigner \
        # -R ${REF_FILE_NEW} \
        # -I ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.bam \
        # -targetIntervals ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.intervals \
        # -o ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.5.bam \
        # -known ${PHASE1_INDEL} -known ${DATABASE_INDEL} \
        # ${EXCLUDE_GRP} &

#mv ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.bai ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.bam.bai

wait

samtools merge -f -@ ${NUM_THREADS} ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.bam ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.1.bam ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.2.bam ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.3.bam ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.4.bam

rm ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.[1-9].ba*
rm ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.*.intervals
rm ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.*.ba*

samtools index -@ ${NUM_THREADS} ${OUT_DIR}/${SAMPLE_ID}.sorted.rmdup.realign.bam

END=`date +%s`

echo IndelRealigner ran for $((END - START)) seconds >> $TIME_LOG
