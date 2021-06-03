#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=4G
#$ -l h_vmem=4G
#$ -t 585-599
#$ -l cancergen

#--------
# Input
#---------------------------------------------------------------------------------
pathFile=<TEXT FILE WITH FILE PATHS OF BAMS>
outDir=<WHERE TO SAVE OUTPUT>
refGenome=<DIRECTORY OF HG19 REFERENCE GENOME>
rlib=<PATH TO R LIBRARIES>
#---------------------------------------------------------------------------------

bamPath=$(cat $pathFile | head -n $SGE_TASK_ID | tail -n 1)
sampleName=$(basename $bamPath | sed s/.bam$//)

if [ ! -f $bamPath ]; then
  echo "$bamPath doesn't exist"
  exit 0
fi

if [ -f $outDir/$sampleName/bt2_bam/${sampleName}_bt2.bam ]; then
  echo "$bamPath has already been processed..."
  exit 0
fi

#-----------------------------------------------------------
# Step 1: Extract reads mapping to the mitochondria
#         (not including PCR duplicates)
#         I require that both mates map to the mitochondria
#---------------------------------------------------------------------------------------
mkdir -p $outDir/$sampleName/orig_mito_bam
mkdir -p $outDir/$sampleName/mito_fastq

samtools view -H $bamPath > $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito.sam
samtools view $bamPath 'chrM' | awk -F '\t' '$7=="="{print}' \
  >> $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito.sam
samtools view -bS $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito.sam > \
  $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito.bam
samtools sort -n $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito.bam > \
  $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito_sorted.bam
samtools bam2fq -1 $outDir/$sampleName/mito_fastq/${sampleName}_R1.fastq \
  -2 $outDir/$sampleName/mito_fastq/${sampleName}_R2.fastq \
  $outDir/$sampleName/orig_mito_bam/${sampleName}_orig_mito_sorted.bam

gzip $outDir/$sampleName/mito_fastq/${sampleName}_R1.fastq
gzip $outDir/$sampleName/mito_fastq/${sampleName}_R2.fastq
#---------------------------------------------------------------------------------------


#----------------------------------------------
# Step 2: Align extracted reads using Bowtie2
#------------------------------------------------------------------------------------------------------------
mkdir -p $outDir/$sampleName/temp
mkdir -p $outDir/$sampleName/bt2_bam

bowtie2 \
  --very-sensitive --end-to-end -x $refGenome -1 $outDir/$sampleName/mito_fastq/${sampleName}_R1.fastq.gz \
  -2 $outDir/$sampleName/mito_fastq/${sampleName}_R2.fastq.gz | samtools view -bS - | \
  samtools sort -T $outDir/$sampleName/temp/temp -O BAM -o $outDir/$sampleName/bt2_bam/${sampleName}_bt2.bam -

samtools index $outDir/$sampleName/bt2_bam/${sampleName}_bt2.bam
#------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------
# Step 3: Create a GAlignmentPairs object
#         Count number of reads mapping to mitochondria
#         Get mitochondrial fragment lengths
#---------------------------------------------------------------------------------------------------
mkdir -p $outDir/$sampleName/plots
mkdir -p $outDir/$sampleName/stats

R_LIBS_USER=$rlib Rscript ./make_galp.R \
  $outDir/$sampleName/bt2_bam/${sampleName}_bt2.bam \
  $outDir/$sampleName
#---------------------------------------------------------------------------------------------------
