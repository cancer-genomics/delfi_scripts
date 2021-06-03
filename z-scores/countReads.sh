#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=16G
#$ -l h_vmem=16G
#$ -l h_fsize=250G
#$ -l h_rt=24:00:00

#----------------------
# Description of inputs
#--------------------------------------------------------------------------------------------
# dataFile : the path to file where each line is the path to a bam file
# outDir : the path to a directory to write output files to
# binGR : a GRanges object of bins to count the # of reads mapped to
#--------------------------------------------------------------------------------------------

#--------
# USER INPUT
#------------------------------------------------------------------------------------------------
dataFile=bam-paths.txt
outDir=countReads.R
#----------
#
# Provided
#
binGR=bins50kb.rds
#------------------------------------------------------------------------------------------------

inBamPath=$(cat $dataFile | awk -F '\t' '{print $1}' | grep -e '.bam$' | head -n $SGE_TASK_ID | tail -n 1)
inBam=$(basename $inBamPath)
bamDir=$(dirname $inBamPath)
sampleName=${inBam//.bam}

outFile=$outDir/binCounts/${sampleName}.rds

if [ ! -f $outFile ]; then
  Rscript ./countReads.R $inBam $bamDir $binGR $outDir
fi
