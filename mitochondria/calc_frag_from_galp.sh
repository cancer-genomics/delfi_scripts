#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=14G
#$ -l h_vmem=14G
#$ -t 585-599
#$ -l cancergen

#-------
# Input
#---------------------------------------------------------------------------------
bamPaths=<A PLAIN TEXT FILE CONTAINING BAM FILE PATHS>.txt
rlib=<YOUR R LIBRARY>
#---------------------------------------------------------------------------------

bamFile=$(cat $bamPaths | head -n $SGE_TASK_ID | tail -n 1)

R_LIBS_USER=$rlib Rscript ./calc_frag_from_galp.R $bamFile
